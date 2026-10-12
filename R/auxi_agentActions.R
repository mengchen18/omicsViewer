############################################################################
### WP-C: run_control server pathway (action tier)
###
### One function sits between the run_control tool / test hook and
### store_action_run, exactly as agent_widget_apply sits between
### set_widgets and store_apply. It owns the settled-policy exclusions:
### the snapshot lifecycle and the dataset selector are PERMANENTLY
### outside the assistant surface (FIX_PLAN policy 1 and 2, 2026-10-11)
### and are rejected by NAME before any registry lookup - they are never
### registered as actions, not even as stubs.
############################################################################

# Controls that are excluded from the assistant surface BY POLICY, not by
# omission: run_control must reject them with the policy reason instead of
# a generic unknown-id error, so a model (or a test) cannot mistake them
# for forgotten coverage. Keys are matched after normalisation, so both
# canonical dotted forms and the DOM-ish dashed forms resolve here.
AGENT_EXCLUDED_CONTROLS <- c(
  "snapshot", "snapshot_save", "snapshot_restore", "snapshot_delete",
  "snapshot_name", "snapshot_include_chat", "selectFile", "select_file",
  "dataset", "dataset_selector"
)

# normalise an arbitrary control reference for policy matching:
# strip an "action:" record prefix, drop leading "app-", unify - and .
.agent_normalize_control_id <- function(x) {
  x <- trimws(as.character(x)[1])
  if (is.na(x) || !nzchar(x)) return("")
  x <- sub("^action:", "", x)
  x <- sub("^app-", "", x)
  chartr("-", ".", x)
}

.agent_control_excluded <- function(x) {
  .agent_normalize_control_id(x) %in% AGENT_EXCLUDED_CONTROLS
}

#' Run one assistant-addressable action (run_control pathway)
#'
#' Mirrors \code{\link{agent_widget_apply}}: normalises the id, enforces
#' the settled policy exclusions (snapshot lifecycle, dataset selector)
#' with a self-describing rejection, then runs
#' \code{\link{store_action_run}} with \code{origin = "agent"}. Hard
#' errors are caught into the receipt's \code{error} field with
#' closest-match suggestions so a model can correct and retry.
#'
#' @param store Canonical widget store (or child view).
#' @param action_id Canonical action id (an "action:" record-id prefix,
#'   as printed by find_controls, is accepted and stripped).
#' @param params Optional named list passed to the handler; actions with
#'   \code{requires_confirm = TRUE} must receive
#'   \code{params$confirm = TRUE}.
#' @return A receipt list: \code{ran} (character, empty on failure),
#'   \code{epoch} (run counter after the run), \code{error} (NULL or a
#'   self-describing reason), \code{excluded_by_policy} (TRUE for the
#'   settled-policy rejection).
#' @keywords internal
#' @rdname agentActionHelpers
agent_action_apply <- function(store, action_id, params = NULL) {
  fail <- function(reason, excluded = FALSE)
    list(ran = character(), epoch = NULL, error = reason,
         excluded_by_policy = excluded)

  if (!is.character(action_id) || length(action_id) != 1L ||
      !nzchar(trimws(action_id)))
    return(fail(paste("run_control requires a single non-empty",
                      "action_id; discover ids with find_controls",
                      "(records of kind 'action').")))

  if (.agent_control_excluded(action_id))
    return(fail(paste0(
      "'", action_id, "' is permanently excluded from the assistant ",
      "surface by policy: the snapshot lifecycle and dataset switching ",
      "are user-only operations. Use the UI controls directly."),
      excluded = TRUE))

  tryCatch({
    receipt <- store_action_run(store, action_id, origin = "agent",
                                params = params)
    list(ran = receipt$ran, epoch = receipt$epoch, error = NULL,
         excluded_by_policy = FALSE)
  }, error = function(e)
    fail(conditionMessage(e)))
}
