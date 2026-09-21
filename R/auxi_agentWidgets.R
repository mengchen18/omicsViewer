#' Internal generic widget-tier helpers (S3)
#'
#' Thin, registry-driven logic behind the generic assistant tools
#' \code{list_widgets}, \code{get_widget}, and \code{set_widgets}
#' (AGENT_ACCURACY_PLAN.md section 6.3, tier 2). Everything here is a pure
#' function of the canonical widget store, so the whole tier is unit-testable
#' without ellmer or a browser; the ellmer tools in
#' \code{\link{ai_assistant_module}} are one-line wrappers around these.
#'
#' Governing principle: agent-controllability mirrors user-controllability
#' exactly. The list/describe views surface only agent-writable (user
#' editable) bindings, and applies are refused per key for anything else.
#'
#' Callers that touch a live Shiny app must invoke these inside
#' \code{\link[shiny]{isolate}} (or another reactive context): validation
#' reads choices providers, which may themselves read module reactives.
#'
#' @keywords internal
#' @name agentWidgetHelpers
NULL

#' Canonical ids of every agent-visible (user-editable) widget
#'
#' Cheaper than \code{\link{store_registry_view}} when allowed values are
#' not needed: never invokes choices providers, so it is safe in any
#' context.
#'
#' @param store Store (or child view).
#' @return Character vector of canonical ids.
#' @keywords internal
#' @rdname agentWidgetHelpers
agent_widget_ids <- function(store) {
  root <- if (is.null(store$parent)) store else store$parent
  ids <- names(root$bindings)
  ids[vapply(ids, function(k) isTRUE(root$bindings[[k]]$agent_writable),
             logical(1))]
}

#' Registry listing for the generic widget tier
#'
#' Wraps \code{\link{store_registry_view}} and adds each widget's current
#' stored value, so one call is enough for the model to discover ids, kinds,
#' allowed values, and present state.
#'
#' @param store Store (or child view).
#' @param section Optional canonical id prefix (component match: an id is
#'   kept when it equals \code{section} or starts with
#'   \code{"section."}, e.g. \code{"dataspace"} or
#'   \code{"dataspace.expr_heatmap"}).
#' @return \code{widget_count} plus a \code{widgets} list of records
#'   (id, kind, label, help, depends_on, allowed_values, current_value).
#' @keywords internal
#' @rdname agentWidgetHelpers
agent_widget_list <- function(store, section = NULL) {
  if (!is.null(section)) {
    section <- trimws(as.character(section)[1])
    if (!nzchar(section))
      section <- NULL
  }
  view <- store_registry_view(store, prefix = section)
  values <- store_read(store)
  out <- lapply(view, function(r)
    c(r, list(current_value = values[[r$id]])))
  if (length(out))
    names(out) <- vapply(out, function(r) r$id, character(1))
  list(
    section = section,
    widget_count = length(out),
    widgets = out
  )
}

#' Describe one agent-visible widget with its current value
#'
#' @param store Store (or child view).
#' @param id Canonical widget id.
#' @return The registry record plus \code{current_value}.
#' @keywords internal
#' @rdname agentWidgetHelpers
agent_widget_describe <- function(store, id) {
  id <- .agent_trim_scalar(id)
  if (!nzchar(id))
    stop("A widget id is required.")
  record <- store_describe(store, id)
  if (is.null(record)) {
    stop("Unknown or not user-editable widget id: ", id, ".",
         .agent_suggest_text(id, agent_widget_ids(store)))
  }
  c(record, list(current_value = store_read(store, record$id)[[1]]))
}

#' Normalize a generic-tier patch to a named list
#'
#' The ellmer tool declares the patch as a JSON-object string (the tool
#' schema cannot express dynamic keys), so the primary input is JSON text.
#' Direct R callers may pass a named list, an unnamed list of
#' \code{list(id =, value =)} records, or a data.frame with \code{id} and
#' \code{value} columns (the ellmer tibble coercion of array-of-object
#' arguments); all shapes normalize to one named list.
#'
#' @param patch JSON-object string or already-parsed patch structure.
#' @return Named list of canonical id -> proposed value, possibly empty.
#' @keywords internal
#' @rdname agentWidgetHelpers
.agent_widget_normalize_patch <- function(patch) {
  if (is.null(patch))
    return(list())
  if (is.character(patch) && length(patch) == 1L) {
    txt <- trimws(patch)
    if (!nzchar(txt) || txt %in% c("null", "NULL", "{}", "[]"))
      return(list())
    parsed <- tryCatch(
      jsonlite::fromJSON(txt, simplifyVector = FALSE),
      error = function(e)
        stop("The widget patch is not valid JSON: ", conditionMessage(e))
    )
    patch <- parsed
  }
  if (is.data.frame(patch)) {
    if (!all(c("id", "value") %in% names(patch)))
      stop("A patch table requires 'id' and 'value' columns.")
    ids <- as.character(patch$id)
    patch <- as.list(patch$value)
    names(patch) <- ids
  } else if (is.list(patch) && is.null(names(patch))) {
    # unnamed list of {id, value} records
    if (length(patch) && all(vapply(patch, function(e)
      is.list(e) && !is.null(e$id) && !is.null(e$value), logical(1)))) {
      ids <- vapply(patch, function(e) as.character(e$id)[1], character(1))
      patch <- lapply(patch, function(e) e$value)
      names(patch) <- ids
    }
  }
  if (is.null(patch) || !is.list(patch))
    stop("The patch must be a JSON object mapping widget ids to values.")
  if (!length(patch))
    return(list())
  nms <- names(patch)
  if (is.null(nms) || any(is.na(nms)) || any(!nzchar(nms)))
    stop("Every patch entry must be named by a widget id.")
  patch
}

#' Apply a generic-tier widget patch
#'
#' The tier-2 counterpart of the curated apply callbacks: validates ids
#' against the agent-visible registry, then delegates to
#' \code{\link{store_apply}} with \code{origin = "agent"} and per-key
#' resilience (WP2 self-correction: invalid values are rejected with
#' suggestions and the valid remainder still applies).
#'
#' @param store Store (or child view).
#' @param patch JSON-object string or named list (see
#'   \code{\link{.agent_widget_normalize_patch}}).
#' @return Receipt: \code{applied} (ids, dependency-ordered),
#'   \code{applied_values}, \code{unchanged} (no-op ids), and
#'   \code{rejected} (list of \code{list(id, reason)} for unknown,
#'   non-editable, or invalid-value keys; reasons carry closest-match
#'   suggestions).
#' @keywords internal
#' @rdname agentWidgetHelpers
agent_widget_apply <- function(store, patch) {
  patch <- .agent_widget_normalize_patch(patch)
  if (!length(patch))
    stop("The widget patch contains no entries.")

  root <- if (is.null(store$parent)) store else store$parent
  writable <- agent_widget_ids(store)

  known <- list()
  rejected <- list()
  for (id in names(patch)) {
    if (id %in% writable) {
      known[[id]] <- patch[[id]]
      next
    }
    reason <- if (!is.null(root$bindings[[id]]))
      paste0("Widget is not user-editable and cannot be set: ", id, ".")
    else
      paste0("Unknown widget id: ", id, ".",
             .agent_suggest_text(id, writable))
    rejected[[length(rejected) + 1L]] <- list(id = id, reason = reason)
  }

  receipt <- if (length(known))
    store_apply(store, known, origin = "agent", strict = FALSE)
  else
    list(applied = character(), skipped = character(), diff = list())

  if (!is.null(receipt$rejected))
    rejected <- c(receipt$rejected, rejected)

  list(
    applied = receipt$applied,
    applied_values = receipt$diff,
    unchanged = receipt$skipped,
    rejected = rejected
  )
}
