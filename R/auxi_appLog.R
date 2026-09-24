# App console logging (debug) ------------------------------------------
#
# Tees everything the running omicsViewer app writes to the R console to
# a per-run log file, without changing what the interactive user sees:
#
#  - print()/cat() output  -> output sink opened with split = TRUE
#    (written to the file AND still echoed to the console)
#  - message()/warning()/error conditions -> calling handlers wrapped
#    around runApp() (logged with a timestamp + level; the usual console
#    display is untouched because the handlers never muffle)
#
# Note on swallowed conditions: calling handlers only see conditions that
# propagate up to the withCallingHandlers() frame. Errors caught by a
# deeper tryCatch() (as Shiny does per observer/reactive) are instead
# captured indirectly through Shiny's own console report ("Warning:
# Error in ..."), which the output sink records.
#
# Default location follows the agent-diagnostic convention
# (tempdir()/omicsviewer-llm-logs): a session-scoped, always-writable
# directory that still survives abnormal termination for post-mortem
# debugging. One file per run (timestamp + pid) so parallel sessions on
# shared servers never clobber each other.
#
# Controls:
#   log_file = <path>        write to an explicit file
#   log_file = FALSE         disable console logging entirely
#   OMICSVIEWER_LOG_DIR      default directory override
#   OMICSVIEWER_LOG = off    disable the default-on behaviour

#' omicsViewer app console-log helpers (internal)
#'
#' Resolve the log destination, start/stop the console tee, and write
#' timestamped condition entries. Not exported; used by
#' \code{\link{omicsViewer}}.
#'
#' @name appConsoleLogHelpers
#' @keywords internal
NULL

.applog_setting <- function(name, default = "") {
  v <- Sys.getenv(name, "")
  if (is.na(v) || !nzchar(v)) default else v
}

.applog_disabled_by_env <- function() {
  tolower(.applog_setting("OMICSVIEWER_LOG", "")) %in%
    c("off", "false", "0", "no")
}

#' @rdname appConsoleLogHelpers
#' @return \code{applog_dir} returns the log directory in effect.
applog_dir <- function() {
  d <- .applog_setting("OMICSVIEWER_LOG_DIR")
  if (!nzchar(d)) d <- file.path(tempdir(), "omicsviewer-logs")
  d
}

#' @rdname appConsoleLogHelpers
#' @return \code{applog_new_file} creates the log directory if needed and
#'   returns the path of a fresh per-run log file.
applog_new_file <- function(dir = applog_dir()) {
  if (!dir.exists(dir))
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(
    dir,
    sprintf("app-%s-p%d.log", format(Sys.time(), "%Y%m%d-%H%M%S"), Sys.getpid())
  )
}

#' @rdname appConsoleLogHelpers
#' @param log_file \code{NULL} (default: auto-named file in
#'   \code{\link{applog_dir}}), an explicit file path, or \code{FALSE} to
#'   disable logging.
#' @return \code{applog_begin} opens the log file, installs the output
#'   sink (split console echo) and writes a header. It returns an
#'   environment with \code{$path}, \code{$con} and a \code{$restore()}
#'   closure, or \code{NULL} when logging is disabled.
applog_begin <- function(log_file = NULL) {
  if (identical(log_file, FALSE) || .applog_disabled_by_env())
    return(NULL)
  if (is.null(log_file) || !nzchar(log_file))
    log_file <- applog_new_file()
  log_file <- path.expand(log_file[1])
  dlog <- dirname(log_file)
  if (!dir.exists(dlog))
    dir.create(dlog, recursive = TRUE, showWarnings = FALSE)

  con <- file(log_file, open = "wt", encoding = "UTF-8")

  cat(
    "== omicsViewer console log ==",
    sprintf("started : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("R       : %s", R.version.string),
    sprintf("shiny   : %s", utils::packageVersion("shiny")),
    sprintf("omicsViewer : %s", utils::packageVersion("omicsViewer")),
    sprintf("pid     : %d", Sys.getpid()),
    "",
    sep = "\n", file = con, append = FALSE
  )
  flush(con)

  sink(con, type = "output", split = TRUE)

  env <- new.env(parent = emptyenv())
  env$path <- log_file
  env$con <- con
  # robust restore: unwind every output-sink frame we own and never let
  # a stack imbalance escape (the app already has a known harmless one)
  env$restore <- function() {
    ok <- TRUE
    tryCatch(
      while (sink.number() > 0L) sink(),
      error = function(e) ok <<- FALSE
    )
    # the message stream is deliberately never diverted here: we never
    # sink type = "message", so unwinding one could pop a user's sink
    tryCatch(
      if (isOpen(env$con)) close(env$con),
      error = function(e) ok <<- FALSE
    )
    invisible(ok)
  }
  env
}

#' @rdname appConsoleLogHelpers
#' @param env logger environment returned by \code{applog_begin}
#'   (\code{NULL} is allowed and a no-op).
#' @return \code{applog_end} restores the console and closes the file;
#'   \code{applog_write} appends one timestamped condition entry.
applog_end <- function(env) {
  if (is.null(env))
    return(invisible(FALSE))
  if (!is.null(env$con)) {
    tryCatch({
      cat(
        sprintf("== ended : %s ==",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
        "",
        sep = "\n", file = env$con, append = TRUE
      )
      flush(env$con)
    }, error = function(e) NULL)
  }
  env$restore()
  invisible(TRUE)
}

#' @rdname appConsoleLogHelpers
#' @param level single character condition level (e.g. "MESSAGE")
#' @param condition a condition caught by a calling handler
applog_write <- function(env, level, condition) {
  if (is.null(env))
    return(invisible(FALSE))
  txt <- conditionMessage(condition)
  cl <- conditionCall(condition)
  if (!is.null(cl)) {
    dc <- deparse(cl)
    if (length(dc) > 0L)
      txt <- paste0(txt, sprintf(" [call: %s]",
        paste(utils::head(dc, 3L), collapse = " ")))
  }
  tryCatch({
    cat(
      sprintf("%s %-7s %s",
              format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"), level, txt),
      file = env$con, append = TRUE, sep = "\n"
    )
    flush(env$con)
  }, error = function(e) invisible(FALSE))
  invisible(TRUE)
}
