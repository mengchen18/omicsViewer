#' Internal selection bus (unified selection control plane)
#'
#' One canonical record per selection space ("feature", "sample"). Every
#' selection source - scatter figure events (lasso/box/click), the volcano
#' corner auto-selection, the correlation / expression / dynamic heatmaps,
#' feature/sample/expression table row clicks, the gene-set list, snapshot
#' restore and agent state application - reports through
#' \code{selection_report()} / \code{selection_apply()}; every consumer -
#' the table row mirrors, the dynamic-heatmap row filters, the result-space
#' selections, the snapshot state - reads \code{selection_watch()}.
#'
#' The bus mirrors the widget-store philosophy (\code{auxi_widgetStore.R}):
#' the core is an environment of plain values plus per-key reactiveVal
#' invalidation signals, and all protocol logic (report dedupe, mirror
#' policy, transactional apply) operates on plain snapshots, so it is fully
#' unit-testable without a browser.
#'
#' Record fields (per key):
#' \itemize{
#'   \item \code{ids} - character vector; THE canonical semantic selection.
#'     Empty means "no selection" (tables un-filter, result space sees an
#'     empty selection).
#'   \item \code{clicked} - character vector; the last single-click report
#'     (figure clicks). Kept for the scatter snapshot round-trip
#'     (\code{selection_clicked}); \code{ids} is always the effective
#'     selection (selected if non-empty, else clicked).
#'   \item \code{origin} - last writer origin: "figure", "corner",
#'     "clear", "table", "heatmap", "cor_heatmap", "dyn_heatmap",
#'     "gslist", "restore", "system".
#'   \item \code{anchor} - scatter axis signature (see
#'     \code{.scatter_axis_signature}) for figure-space origins. Display
#'     gating only: selection emphasis must not be carried into a
#'     different figure. Never affects propagation.
#'   \item \code{mirror} - TRUE (no filter) or character ids: the table-row
#'     mirror. Table-origin reports leave it unchanged: a DataTable row
#'     click must not re-assert itself through the DT proxy highlight
#'     (the proxy cannot be distinguished from a user click server-side).
#'   \item \code{epoch} - integer; bumps on every applied change.
#' }
#'
#' Adoption rule (uniform, replaces the per-source \code{.xxx_last} guards
#' that previously lived in the data-space module): a source reports only
#' when its INTERACTION REPORT VALUE changes. Module returns recompute for
#' many non-interaction reasons (snapshot attributes, store pushes, DT
#' redraws); \code{selection_report()} skips a report identical to the
#' last one FROM THE SAME ORIGIN, so an echo can never clobber a selection
#' another source has just made, while a genuine value change from any
#' source always wins.
#'
#' @keywords internal
#' @name selectionBusHelpers
NULL

############################################################################
### [1] store construction and ports
############################################################################

.selection_store_origins <- c(
  "figure", "corner", "clear", "table", "heatmap", "cor_heatmap",
  "dyn_heatmap", "gslist", "restore", "system"
)

#' Create a selection bus
#'
#' @param keys Character keys of the selection spaces (default
#'   \code{c("feature", "sample")}).
#' @return An environment holding one record per key (plain value plus a
#'   reactiveVal invalidation signal), per-(key, origin) last-report slots
#'   and a reactive global epoch.
#' @keywords internal
#' @rdname selectionBusHelpers
selection_store_new <- function(keys = c("feature", "sample")) {
  store <- new.env(parent = emptyenv())
  store$keys <- keys
  store$records <- list()   # key -> plain record snapshot
  store$signals <- list()   # key -> reactiveVal (invalidation only)
  store$reports <- list()   # "key\rorigin" -> last report value
  store$global_epoch <- 0L
  store$epoch_rv <- shiny::reactiveVal(0L)
  for (k in keys) {
    store$records[[k]] <- .selection_record_empty()
    store$signals[[k]] <- shiny::reactiveVal(0L)
  }
  store
}

.selection_record_empty <- function() {
  list(ids = character(0), clicked = character(0), origin = "init",
       anchor = NULL, mirror = TRUE, epoch = 0L)
}

#' Bind a selection bus key into a port
#'
#' A port fixes the key so modules report and watch their own space without
#' repeating it (the \code{\link{widget_store_child}} pattern applied to
#' the selection bus).
#'
#' @param store Store from \code{\link{selection_store_new}}.
#' @param key One of the store's keys.
#' @return A list(apply =, report =, watch =, read =, key =, store =).
#' @keywords internal
#' @rdname selectionBusHelpers
selection_port <- function(store, key) {
  stopifnot(!is.null(store), key %in% store$keys)
  list(
    store = store,
    key = key,
    report = function(origin, report, ids, clicked = character(0),
                      anchor = NULL, mirror = NULL)
      selection_report(store, key, origin = origin, report = report,
                       ids = ids, clicked = clicked, anchor = anchor,
                       mirror = mirror),
    apply = function(ids, clicked = character(0), origin = "system",
                     anchor = NULL, mirror = NULL)
      selection_apply(store, key, ids = ids, clicked = clicked,
                      origin = origin, anchor = anchor, mirror = mirror),
    watch = function() selection_watch(store, key),
    read = function() selection_read(store, key)
  )
}

############################################################################
### [2] the transactional protocol
############################################################################

.selection_core <- function(ids, clicked, origin, anchor, mirror, previous) {
  # NULL fields mean "leave unchanged" (anchor/mirror), so the core is
  # resolved against the previous record before the identity comparison.
  list(
    ids = if (is.null(ids)) previous$ids else as.character(ids),
    clicked = if (is.null(clicked)) previous$clicked else as.character(clicked),
    origin = origin,
    anchor = if (is.null(anchor)) previous$anchor else anchor,
    mirror = if (is.null(mirror)) previous$mirror else mirror
  )
}

.selection_commit <- function(store, key, core) {
  previous <- store$records[[key]]
  if (identical(core, previous[names(core)]))
    return(invisible(FALSE))
  record <- c(core, list(epoch = previous$epoch + 1L))
  store$records[[key]] <- record
  store$signals[[key]](record$epoch)  # invalidation signal only
  store$global_epoch <- (store$global_epoch %||% 0L) + 1L
  store$epoch_rv(store$global_epoch)
  invisible(TRUE)
}

#' Report an interaction-driven selection change
#'
#' Reports are deduped per (key, origin) on the report VALUE: an unchanged
#' report is an echo (a module return that recomputed without a user
#' interaction) and is skipped, so it can never revert a newer selection
#' made by another source. A changed report is applied transactionally
#' (identical cores are no-ops).
#'
#' @param store Store from \code{\link{selection_store_new}}.
#' @param key Selection space key.
#' @param origin Origin slot (see \code{\link{selection_store_new}}).
#' @param report The interaction report value compared for dedupe (any
#'   object; compared with \code{identical()}).
#' @param ids Effective selection (selected if non-empty, else clicked).
#' @param clicked Last click report (optional).
#' @param anchor Axis signature for figure-space origins.
#' @param mirror Table-row mirror; NULL leaves it unchanged.
#' @return TRUE when the report was applied.
#' @keywords internal
#' @rdname selectionBusHelpers
selection_report <- function(store, key, origin, report, ids,
                             clicked = character(0), anchor = NULL,
                             mirror = NULL) {
  stopifnot(origin %in% .selection_store_origins)
  slot <- paste(key, origin, sep = "\r")
  if (identical(report, store$reports[[slot]]))
    return(invisible(FALSE))
  store$reports[[slot]] <- report
  selection_apply(store, key, ids = ids, clicked = clicked,
                  origin = origin, anchor = anchor, mirror = mirror)
}

#' Apply a selection change directly (restore / system writers)
#'
#' Bypasses the report dedupe (the writer states the change is real);
#' still dedupes on the resolved record core.
#'
#' @return TRUE when the write changed the record.
#' @keywords internal
#' @rdname selectionBusHelpers
selection_apply <- function(store, key, ids, clicked = character(0),
                            origin = "system", anchor = NULL,
                            mirror = NULL) {
  stopifnot(origin %in% .selection_store_origins, key %in% store$keys)
  core <- .selection_core(ids, clicked, origin, anchor, mirror,
                          store$records[[key]])
  .selection_commit(store, key, core)
}

############################################################################
### [3] reactive reads and snapshot state
############################################################################

#' Reactive view of one selection record
#'
#' @return A reactive that re-derives whenever the key's record changes.
#' @keywords internal
#' @rdname selectionBusHelpers
selection_watch <- function(store, key) {
  shiny::reactive({
    store$signals[[key]]()
    store$records[[key]]
  })
}

#' Plain (non-reactive) read of one selection record
#' @keywords internal
#' @rdname selectionBusHelpers
selection_read <- function(store, key) {
  store$records[[key]]
}

#' Snapshot state of every key (for tests and diagnostics)
#' @keywords internal
#' @rdname selectionBusHelpers
selection_state <- function(store) {
  stats::setNames(lapply(store$keys, function(k) store$records[[k]]),
                  store$keys)
}
