#' Internal universal widget-control-plane store (S1)
#'
#' Canonical widget store implementing the universal control plane designed
#' in AGENT_ACCURACY_PLAN.md section 6. One file, organized by numbered
#' comment sections:
#'
#' \enumerate{
#'   \item store construction and widget registration
#'   \item value validation (per kind, boundary sentinels)
#'   \item the transactional apply protocol (diff, ordering, epochs)
#'   \item acknowledgement and UI-to-store synchronisation
#'   \item snapshot serialisation and restore
#'   \item agent-facing registry views (discovery)
#' }
#'
#' The S1 core is deliberately free of Shiny observers: values live in
#' reactiveVals (usable outside a session), while all protocol logic
#' (validation, diffing, ordering, acks) operates on plain snapshots and is
#' therefore fully unit-testable without a browser. Session glue (binding
#' real inputs, invoking setters, ack subscriptions) lands with S2.
#'
#' Governing principle: agent-controllability mirrors user-controllability
#' exactly. Only user-editable widgets are agent-visible; everything else is
#' at most internal snapshot state and never enters model-facing context.
#'
#' @keywords internal
#' @name widgetStoreHelpers
NULL

############################################################################
### [1] store construction and widget registration
############################################################################

.widget_store_id_pattern <- "^[a-zA-Z][a-zA-Z0-9_]*(\\.[a-zA-Z][a-zA-Z0-9_]*)*$"

.widget_store_kinds <- c(
  "string", "numeric", "integer", "boolean", "enum",
  "select", "select_cascaded", "tabset", "navbar", "checkbox", "slider"
)

#' Create a canonical widget store
#'
#' @return An environment-based store holding bindings, values (each key
#'   keeps a plain value for protocol reads plus a \code{reactiveVal} used
#'   only as an invalidation signal for Shiny consumers - shiny >= 1.14
#'   forbids reading reactiveVals outside reactive contexts), per-key and
#'   global epochs, a reactive epoch counter, pending (unacknowledged)
#'   values, and a bounded override log. Child stores (namespace-prefixed
#'   views) are created with \code{\link{widget_store_child}}.
#' @keywords internal
#' @rdname widgetStoreHelpers
widget_store_new <- function() {
  store <- new.env(parent = emptyenv())
  store$bindings <- list()      # canonical id -> binding record
  store$values <- list()        # canonical id -> reactiveVal holder
  store$epochs <- list()        # canonical id -> integer epoch
  store$pending <- list()       # canonical id -> list(value, epoch)
  store$origins <- list()       # canonical id -> last write origin
  store$global_epoch <- 0L
  store$epoch_rv <- shiny::reactiveVal(0L)  # reactive transaction counter
  store$override_log <- list()  # user-overridden in-flight agent writes
  store
}

#' Create a namespace-prefixed child view of a store
#'
#' Mirrors \code{NS()}: a child prepends \code{prefix} to every id used
#' through it, so modules register and apply under their own namespace
#' without hand-writing full canonical ids.
#'
#' @param store Store from \code{\link{widget_store_new}}.
#' @param prefix Namespace such as \code{"dataspace.feature_space"}.
#' @return A child view (its own environment delegating to the parent).
#' @keywords internal
#' @rdname widgetStoreHelpers
widget_store_child <- function(store, prefix) {
  stopifnot(
    !is.null(store),
    grepl(.widget_store_id_pattern, prefix)
  )
  child <- new.env(parent = emptyenv())
  child$parent <- store
  child$prefix <- prefix
  child
}

#' Describe one widget binding for the store
#'
#' @param id Canonical dotted id (unique per store).
#' @param kind Widget kind; selects the value validator.
#' @param label Short human/model-facing label.
#' @param help One-sentence help text shared by tooltips and the registry.
#' @param depends_on Canonical ids this widget's choices depend on (for
#'   cascaded selects); used to order transactional writes.
#' @param choices_provider Optional function of the effective state
#'   returning the currently allowed values for this widget. The function
#'   receives the full named list of current values (store state overlaid
#'   with any pending patch entries during transactional validation), which
#'   is what makes cascaded selects (subset depends on analysis) validate
#'   correctly inside a single patch.
#' @param values Optional static allowed values for \code{enum}/\code{slider}.
#' @param min,max Optional numeric bounds for \code{numeric}/\code{integer}/
#'   \code{slider}.
#' @param user_editable TRUE (default) when the user can change this widget
#'   in the UI. Governs agent visibility: only user-editable widgets are
#'   agent-writable and agent-discoverable.
#' @param internal TRUE for registered snapshot state that is not a widget
#'   at all (never agent-visible, restorable only by origin
#'   \code{"restore"}).
#' @param setter Optional setter function (S2 session glue); validated but
#'   not invoked by the S1 core.
#' @param getter Optional getter function (S2 session glue).
#' @return A validated binding record (named list).
#' @keywords internal
#' @rdname widgetStoreHelpers
widget_binding <- function(id,
                           kind = c("string", "numeric", "integer", "boolean",
                                    "enum", "select", "select_cascaded",
                                    "tabset", "navbar", "checkbox", "slider"),
                           label = "",
                           help = "",
                           depends_on = character(),
                           choices_provider = NULL,
                           values = NULL,
                           min = NULL, max = NULL,
                           user_editable = TRUE,
                           internal = FALSE,
                           setter = NULL, getter = NULL) {
  kind <- match.arg(kind)
  if (!is.character(id) || length(id) != 1L || !nzchar(id) ||
      !grepl(.widget_store_id_pattern, id))
    stop("Widget id must be a non-empty dotted identifier: ", id)
  label <- trimws(as.character(label)[1])
  help <- trimws(as.character(help)[1])
  if (nchar(help) > 400L)
    stop("Widget help text must be at most 400 characters.")
  if (!is.character(depends_on))
    stop("depends_on must be a character vector of canonical ids.")
  if (any(!nzchar(depends_on)))
    stop("depends_on entries must be non-empty ids.")
  if (!is.logical(user_editable) || length(user_editable) != 1L ||
      is.na(user_editable))
    stop("user_editable must be TRUE or FALSE.")
  if (!is.logical(internal) || length(internal) != 1L || is.na(internal))
    stop("internal must be TRUE or FALSE.")
  # internal state is by definition not user-editable; treat internal=TRUE
  # as forcing user_editable=FALSE instead of erroring (ergonomics).
  user_editable <- isTRUE(user_editable) && !isTRUE(internal)
  for (fn in c(setter, getter)) {
    if (!is.null(fn) && !is.function(fn))
      stop("setter/getter must be functions or NULL.")
  }
  if (!is.null(choices_provider) && !is.function(choices_provider))
    stop("choices_provider must be a function or NULL.")
  if (kind %in% c("enum", "slider") && is.null(values) && kind == "enum")
    stop("enum widgets require static values.")

  list(
    id = id, kind = kind, label = label, help = help,
    depends_on = depends_on,
    choices_provider = choices_provider,
    values = values, min = min, max = max,
    user_editable = isTRUE(user_editable),
    internal = isTRUE(internal),
    setter = setter, getter = getter,
    agent_writable = isTRUE(user_editable) && !isTRUE(internal)
  )
}

.widget_store_key <- function(store, id) {
  if (!is.null(store$parent))
    paste(store$prefix, id, sep = ".")
  else
    id
}

#' Register one or more widget bindings in a store
#'
#' Rejects duplicate ids, unknown dependencies, and dependency cycles so
#' transactional ordering is always well defined.
#'
#' @param store Store (or child view).
#' @param ... Binding records from \code{\link{widget_binding}}.
#' @return The store, invisibly.
#' @keywords internal
#' @rdname widgetStoreHelpers
store_register <- function(store, ...) {
  bindings <- list(...)
  if (!length(bindings))
    return(invisible(store))
  # child views prefix ids, but storage always lives on the root store
  if (!is.null(store$parent)) {
    prefixer <- function(b) {
      b$id <- .widget_store_key(store, b$id)
      if (length(b$depends_on)) {
        # depends_on entries are module-local references: prefix them too,
        # unless already fully qualified for this namespace
        b$depends_on <- vapply(b$depends_on, function(d) {
          if (startsWith(d, paste0(store$prefix, "."))) d
          else .widget_store_key(store, d)
        }, character(1), USE.NAMES = FALSE)
      }
      b
    }
    bindings <- lapply(bindings, prefixer)
    store <- store$parent
  }
  for (b in bindings) {
    if (!is.list(b) || is.null(b$id) || is.null(b$kind))
      stop("Only widget_binding() records can be registered.")
    id <- b$id
    if (id %in% names(store$bindings))
      stop("Widget already registered: ", id)
    store$bindings[[id]] <- b
    store$epochs[[id]] <- 0L
    store$values[[id]] <- list(val = NULL, rv = shiny::reactiveVal(NULL))
  }
  # resolve dependencies after all registrations in this call
  known <- names(store$bindings)
  for (id in known) {
    deps <- store$bindings[[id]]$depends_on
    missing_dep <- setdiff(deps, known)
    if (length(missing_dep))
      stop("Widget ", id, " depends on unregistered id(s): ",
           paste(missing_dep, collapse = ", "))
  }
  # dependency cycle detection (Kahn)
  if (.widget_store_cycle(store))
    stop("Widget dependency graph contains a cycle.")
  invisible(store)
}

.widget_store_cycle <- function(store) {
  ids <- names(store$bindings)
  indeg <- vapply(ids, function(id) {
    length(intersect(store$bindings[[id]]$depends_on, ids))
  }, integer(1))
  ready <- ids[indeg == 0L]
  seen <- character()
  while (length(ready)) {
    cur <- ready[1L]
    ready <- ready[-1L]
    seen <- c(seen, cur)
    for (id in ids) {
      if (cur %in% store$bindings[[id]]$depends_on) {
        indeg[[id]] <- indeg[[id]] - 1L
        if (indeg[[id]] == 0L)
          ready <- c(ready, id)
      }
    }
  }
  length(seen) != length(ids)
}

############################################################################
### [2] value validation (per kind, boundary sentinels)
############################################################################

# Some providers serialise omitted optional strings as literal sentinels
# ("null", "{}", "[]"). Normalise at the store boundary so every widget
# kind benefits instead of each tool patching its own arguments.
.widget_store_sentinel <- function(value) {
  if (is.character(value) && length(value) == 1L &&
      value %in% c("null", "{}", "[]"))
    NULL
  else
    value
}

.widget_store_allowed_values <- function(binding, effective) {
  if (is.function(binding$choices_provider)) {
    tryCatch(binding$choices_provider(effective), error = function(e) NULL)
  } else if (!is.null(binding$values)) {
    binding$values
  } else {
    NULL
  }
}

.widget_store_validate_value <- function(binding, value, effective) {
  value <- .widget_store_sentinel(value)
  if (is.null(value))
    return(list(value = NULL))
  kind <- binding$kind

  if (kind %in% c("string", "select", "select_cascaded", "tabset", "navbar")) {
    if (!is.character(value) || length(value) != 1L || !nzchar(value))
      return(list(error = paste(binding$id, "requires a single non-empty string.")))
    if (kind %in% c("select", "select_cascaded", "tabset", "navbar")) {
      allowed <- .widget_store_allowed_values(binding, effective)
      if (!is.null(allowed) && !value %in% allowed) {
        hint <- .agent_suggest_text(value, allowed)
        return(list(error = paste0(
          "Unknown value for ", binding$id, ": ", value, ".", hint,
          " Allowed: ", paste(head(allowed, 10), collapse = ", "))))
      }
    }
    return(list(value = value))
  }

  if (kind %in% c("numeric", "integer", "slider")) {
    num <- suppressWarnings(as.numeric(value)[1])
    if (is.na(num))
      return(list(error = paste(binding$id, "requires a number.")))
    if (kind == "integer" && num != round(num))
      return(list(error = paste(binding$id, "requires an integer.")))
    if (!is.null(binding$min) && num < binding$min)
      return(list(error = paste0(binding$id, " must be >= ", binding$min, ".")))
    if (!is.null(binding$max) && num > binding$max)
      return(list(error = paste0(binding$id, " must be <= ", binding$max, ".")))
    return(list(value = if (kind == "integer") as.integer(num) else num))
  }

  if (kind %in% c("boolean", "checkbox")) {
    if (is.logical(value) && length(value) == 1L && !is.na(value))
      return(list(value = value))
    txt <- tolower(trimws(as.character(value)[1]))
    if (txt %in% c("true", "false"))
      return(list(value = txt == "true"))
    return(list(error = paste(binding$id, "requires true or false.")))
  }

  if (kind == "enum") {
    value <- as.character(value)[1]
    if (!value %in% binding$values)
      return(list(error = paste0(
        "Unknown value for ", binding$id, ": ", value, ".",
        .agent_suggest_text(value, binding$values),
        " Allowed: ", paste(head(binding$values, 10), collapse = ", "))))
    return(list(value = value))
  }

  list(error = paste("Unsupported widget kind:", kind))
}

############################################################################
### [3] the transactional apply protocol
############################################################################

#' Compute the dependency-ordered write plan for a set of touched ids
#' @keywords internal
#' @rdname widgetStoreHelpers
.store_write_plan <- function(store, touched) {
  ids <- names(store$bindings)
  plan <- character()
  remaining <- intersect(touched, ids)
  # Kahn over the full graph restricted to touched nodes, honouring
  # dependencies on untouched nodes implicitly (they are already correct).
  indeg <- vapply(remaining, function(id) {
    length(intersect(store$bindings[[id]]$depends_on, remaining))
  }, integer(1))
  while (length(remaining)) {
    ready <- remaining[indeg == 0L]
    if (!length(ready))
      stop("Internal error: dependency cycle escaped registration checks.")
    plan <- c(plan, ready)
    remaining <- setdiff(remaining, ready)
    indeg <- indeg[remaining]
    if (length(remaining))
      indeg <- vapply(remaining, function(id) {
        length(intersect(store$bindings[[id]]$depends_on, remaining))
      }, integer(1))
  }
  plan
}

#' Read current store values as a plain named list
#'
#' @param store Store (or child view).
#' @param ids Optional ids to read; default all.
#' @return Named list of current desired values.
#' @keywords internal
#' @rdname widgetStoreHelpers
store_read <- function(store, ids = NULL) {
  if (!is.null(store$parent)) {
    child <- store
    store <- store$parent
    want <- vapply(ids %||% character(), function(k) .widget_store_key(child, k),
                   character(1), USE.NAMES = FALSE)
    if (!length(want))
      want <- names(store$values)
  } else {
    want <- if (is.null(ids)) names(store$values) else ids
  }
  out <- lapply(want, function(id) store$values[[id]]$val)
  names(out) <- want
  out
}

#' Apply a patch transactionally
#'
#' The single entry point for every external write (agent tools, snapshot
#' restore, future automation). Validates the whole patch first, computes
#' the diff against the *store* (never raw widget inputs), then writes
#' touched keys in dependency order, bumping epochs and recording pending
#' (unacknowledged) values. Untouched keys are never written, so applies
#' have no side effects beyond what was requested.
#'
#' @param store Store (or child view).
#' @param patch Named list of canonical id -> proposed value; \code{NULL}
#'   values are treated as omitted sentinels and dropped.
#' @param origin One of \code{"agent"}, \code{"restore"}, \code{"system"}.
#'   Agent-origin writes require agent-writable (user-editable) widgets;
#'   restore-origin writes may also set internal state.
#' @param strict TRUE (default): any invalid value aborts the whole
#'   transaction. FALSE: invalid values are rejected per key (recorded in
#'   \code{receipt$rejected}) and the valid remainder still applies - used
#'   by restores, where one module's transient junk (e.g. an unset
#'   \code{--select--} placeholder) must not veto unrelated keys.
#' @return Invisible receipt: \code{applied} (ordered ids), \code{diff},
#'   \code{epochs}, \code{skipped} (no-op keys), \code{global_epoch},
#'   and (when not strict) \code{rejected}.
#' @keywords internal
#' @rdname widgetStoreHelpers
store_apply <- function(store, patch,
                        origin = c("agent", "restore", "system"),
                        strict = TRUE) {
  origin <- match.arg(origin)
  if (is.null(store$parent)) {
    ids <- names(patch)
  } else {
    parent <- store$parent
    ids <- vapply(names(patch), function(k) .widget_store_key(store, k),
                  character(1), USE.NAMES = FALSE)
    patch <- stats::setNames(patch, ids)
    store <- parent
  }
  if (!is.list(patch))
    stop("Patch must be a named list.")

  # drop explicit NULLs (omitted sentinels)
  entries <- list()
  for (id in ids) {
    value <- .widget_store_sentinel(patch[[id]])
    if (is.null(value)) next
    binding <- store$bindings[[id]]
    if (is.null(binding))
      stop("Unknown widget id: ", id, ".",
           .agent_suggest_text(id, names(store$bindings)))
    if (origin == "agent" && !binding$agent_writable)
      stop("Widget is not user-editable and cannot be set by the agent: ", id)
    entries[[id]] <- value
  }

  # Validate in dependency order against an *effective* view of the state
  # (current values overlaid with earlier patch entries). This is what lets
  # one patch set a cascaded group (analysis -> subset -> variable) whose
  # members are only jointly valid: the subset is checked against the
  # analysis the SAME patch is about to install.
  current <- store_read(store, names(store$bindings))
  plan_order <- .store_write_plan(store, names(entries))
  effective <- current
  rejected <- list()
  for (id in plan_order) {
    checked <- .widget_store_validate_value(store$bindings[[id]], entries[[id]],
                                            effective)
    if (!is.null(checked$error)) {
      if (strict)
        stop(checked$error)
      # keep the effective view unchanged so downstream keys validate
      # against the pre-existing state, and record the rejection
      rejected[[length(rejected) + 1L]] <- list(id = id, reason = checked$error)
      entries[[id]] <- NULL
      next
    }
    entries[[id]] <- checked$value
    effective[[id]] <- checked$value
  }

  # diff against the store: only genuinely changed keys are written
  diff_keys <- names(entries)[!vapply(names(entries), function(id) {
    identical(entries[[id]], current[[id]])
  }, logical(1))]

  plan <- .store_write_plan(store, diff_keys)
  for (id in plan) {
    store$values[[id]]$val <- entries[[id]]
    store$values[[id]]$rv(entries[[id]])
    store$epochs[[id]] <- (store$epochs[[id]] %||% 0L) + 1L
    store$origins[[id]] <- origin
    store$pending[[id]] <- list(value = entries[[id]],
                                epoch = store$epochs[[id]])
  }
  if (length(plan)) {
    store$global_epoch <- store$global_epoch + 1L
    store$epoch_rv(store$global_epoch)  # invalidate reactive consumers once per transaction
  }

  receipt <- list(
    applied = plan,
    diff = entries[plan],
    epochs = store$epochs[plan],
    skipped = setdiff(names(entries), plan),
    global_epoch = store$global_epoch
  )
  if (!strict)
    receipt$rejected <- rejected
  invisible(receipt)
}

############################################################################
### [4] acknowledgement and UI-to-store synchronisation
############################################################################

#' Acknowledge a widget value on behalf of the browser
#'
#' Called by S2 session glue when a widget confirms (or reports) its current
#' value. Clears matching pending entries; a mismatch returns a re-assert
#' signal so the bounded retry loop can re-send the store value.
#'
#' @param store Store (or child view).
#' @param id Canonical id.
#' @param widget_value The value the widget actually reports.
#' @return TRUE when a re-assert (re-send) is needed, otherwise FALSE.
#' @keywords internal
#' @rdname widgetStoreHelpers
store_ack <- function(store, id, widget_value) {
  key <- .widget_store_key(store, id)
  if (is.null(store$parent)) {
    pending <- store$pending[[key]]
  } else {
    pending <- store$parent$pending[[key]]
    store <- store$parent
  }
  if (is.null(pending))
    return(FALSE)
  if (identical(.widget_store_sentinel(widget_value), pending$value)) {
    store$pending[[key]] <- NULL
    FALSE
  } else {
    TRUE
  }
}

#' Mirror a manual user edit into the store
#'
#' UI-to-store sync: updates the stored value with origin \code{"user"},
#' without bumping epochs or marking anything pending. The user always
#' wins: an in-flight (unacknowledged) agent write for the same key is
#' cleared and recorded in the override log.
#'
#' @param store Store (or child view).
#' @param id Canonical id.
#' @param value The value the user chose.
#' @return Invisible TRUE when an agent write was overridden.
#' @keywords internal
#' @rdname widgetStoreHelpers
store_sync_from_ui <- function(store, id, value) {
  key <- .widget_store_key(store, id)
  if (!is.null(store$parent)) store <- store$parent
  binding <- store$bindings[[key]]
  if (is.null(binding))
    return(invisible(FALSE))
  checked <- .widget_store_validate_value(binding, value,
                                          store_read(store, names(store$bindings)))
  value <- if (is.null(checked$error)) checked$value else value
  # acknowledgement: a widget confirming an in-flight external write is NOT
  # a user override - clear the pending entry and keep the write's origin
  pending <- store$pending[[key]]
  if (!is.null(pending) && identical(value, pending$value)) {
    store$pending[[key]] <- NULL
    store$values[[key]]$val <- value
    store$values[[key]]$rv(value)
    return(invisible(FALSE))
  }
  overridden <- !is.null(store$pending[[key]])
  if (overridden) {
    store$override_log <- c(
      store$override_log,
      list(list(id = key, pending = store$pending[[key]]$value,
                user = value,
                epoch = store$epochs[[key]]))
    )
    if (length(store$override_log) > 100L)
      store$override_log <- utils::tail(store$override_log, 100L)
    store$pending[[key]] <- NULL
  }
  store$values[[key]]$val <- value
  store$values[[key]]$rv(value)
  store$origins[[key]] <- "user"
  invisible(overridden)
}

#' React to one widget's store value
#'
#' S2 session glue: returns a reactive expression reading the key's
#' reactiveVal, so effects can depend on store writes directly.
#'
#' @param store Store (or child view).
#' @param id Canonical id.
#' @return A reactive expression yielding the current stored value.
#' @keywords internal
#' @rdname widgetStoreHelpers
store_watch <- function(store, id) {
  key <- .widget_store_key(store, id)
  root <- if (is.null(store$parent)) store else store$parent
  # touch the reactiveVal (the invalidation signal), read the plain value
  shiny::reactive({
    root$values[[key]]$rv()
    root$values[[key]]$val
  })
}

#' Reactive transaction counter for the store
#'
#' Bumped once per applying transaction. Consumers that must re-assert after
#' ANY external write (e.g. a cascade group) watch this instead of tracking
#' per-key epochs.
#'
#' @param store Store (or child view).
#' @return A reactive expression yielding the global epoch.
#' @keywords internal
#' @rdname widgetStoreHelpers
store_epoch <- function(store) {
  root <- if (is.null(store$parent)) store else store$parent
  # epoch_rv stores the global epoch as its value
  shiny::reactive(root$epoch_rv())
}

############################################################################
### [5] snapshot serialisation and restore
############################################################################

#' Serialise the full store state
#'
#' @return A plain list with \code{values} and \code{epochs} for every
#'   registered id (including internal state).
#' @keywords internal
#' @rdname widgetStoreHelpers
store_snapshot <- function(store) {
  if (!is.null(store$parent)) store <- store$parent
  list(
    values = store_read(store),
    epochs = as.list(store$epochs)
  )
}

#' Restore a snapshot through the apply protocol
#'
#' Restore-origin patches may set internal (non-user-editable) state, which
#' is exactly how snapshot round-trips reach keys the agent may never
#' touch. Restores are per-key resilient: a snapshot saved against another
#' dataset may contain values that no longer validate; those keys are
#' reported in \code{receipt$rejected} while the valid remainder still
#' applies. Unknown ids in the snapshot are reported, not applied.
#'
#' @param store Store (or child view).
#' @param snapshot List from \code{\link{store_snapshot}}.
#' @return Invisible receipt from \code{\link{store_apply}} plus
#'   \code{unknown_ids} and (per-key) \code{rejected}.
#' @keywords internal
#' @rdname widgetStoreHelpers
store_restore <- function(store, snapshot) {
  if (!is.null(store$parent)) store <- store$parent
  if (is.null(snapshot) || !is.list(snapshot$values))
    stop("Snapshot must come from store_snapshot().")
  known <- names(store$bindings)
  unknown <- setdiff(names(snapshot$values), known)
  patch <- snapshot$values[intersect(names(snapshot$values), known)]
  receipt <- store_apply(store, patch, origin = "restore", strict = FALSE)
  receipt$unknown_ids <- unknown
  invisible(receipt)
}

############################################################################
### [6] agent-facing registry views (discovery)
############################################################################

#' Registry view of every agent-visible widget
#'
#' Implements the context-hygiene principle: only user-editable widgets
#' appear. Internal state and anything the user cannot change never enters
#' model-facing context. The WP9 capability registry is generated from
#' these records.
#'
#' @param store Store (or child view).
#' @param prefix Optional id prefix filter (component match: an id is kept
#'   when it equals \code{prefix} or starts with \code{"prefix."}, so a
#'   \code{dataspace} section matches \code{dataspace.*} but not
#'   \code{dataspace2.*}).
#' @return A list of records: id, kind, label, help, depends_on, and the
#'   current allowed values when cheap to compute.
#' @keywords internal
#' @rdname widgetStoreHelpers
store_registry_view <- function(store, prefix = NULL) {
  if (!is.null(store$parent)) {
    if (!is.null(prefix))
      prefix <- paste(store$prefix, prefix, sep = ".")
    store <- store$parent
  }
  out <- list()
  for (id in names(store$bindings)) {
    b <- store$bindings[[id]]
    if (!b$agent_writable) next
    if (!is.null(prefix) && !identical(id, prefix) &&
        !startsWith(id, paste0(prefix, "."))) next
    allowed <- NULL
    if (is.function(b$choices_provider)) {
      vals <- store_read(store, names(store$bindings))
      allowed <- tryCatch(b$choices_provider(vals), error = function(e) NULL)
    } else if (!is.null(b$values))
      allowed <- b$values
    out[[length(out) + 1L]] <- list(
      id = id, kind = b$kind, label = b$label, help = b$help,
      depends_on = b$depends_on,
      min = b$min, max = b$max,
      allowed_values = if (is.null(allowed)) NULL else head(allowed, 50)
    )
  }
  out
}

#' Describe one agent-visible widget
#'
#' @param store Store (or child view).
#' @param id Canonical id.
#' @return The registry record, or NULL when the id is unknown or not
#'   agent-visible (never leaked to the model).
#' @keywords internal
#' @rdname widgetStoreHelpers
store_describe <- function(store, id) {
  key <- .widget_store_key(store, id)
  view <- Filter(function(r) identical(r$id, key),
                 store_registry_view(store))
  if (length(view)) view[[1]] else NULL
}
