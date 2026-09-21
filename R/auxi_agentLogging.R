#' Internal opt-in AI-assistant diagnostic logging
#'
#' Logging is disabled unless explicitly enabled by the administrator or the
#' assistant drawer checkbox. Events are written as one JSON object per line.
#' Logs may contain user prompts, assistant responses, tool arguments, and
#' bounded tool results derived from the loaded dataset; they never contain API
#' keys or other credential values.
#'
#' @return Logger objects, sanitized event payloads, or ellmer turn summaries.
#'
#' @keywords internal
#' @name agentLoggingHelpers
NULL

.agent_log_true <- function(x) {
  tolower(x) %in% c("true", "t", "yes", "y", "on", "1")
}

.agent_log_scalar_setting <- function(name, default = "") {
  x <- Sys.getenv(name)
  if (is.na(x) || !nzchar(x) || is.na(x[1]))
    return(default)
  trimws(x[1])
}

#' Read assistant diagnostic-log settings
#'
#' Supported settings are:
#' \itemize{
#'   \item \code{OMICSVIEWER_LLM_LOG}: true/false, yes/no, on/off, or 1/0.
#'   \item \code{OMICSVIEWER_LLM_LOG_DIR}: writable log directory.
#'   \item \code{OMICSVIEWER_LLM_LOG_MAX_BYTES}: 1 MB--100 MB; default 10 MB.
#' }
#'
#' @return A list containing enabled status, log directory, and size limit.
#' @keywords internal
#' @rdname agentLoggingHelpers
agent_logging_config <- function() {
  enabled <- .agent_log_true(.agent_log_scalar_setting("OMICSVIEWER_LLM_LOG", "false"))

  directory <- .agent_log_scalar_setting("OMICSVIEWER_LLM_LOG_DIR")
  if (!nzchar(directory))
    directory <- file.path(tempdir(), "omicsviewer-llm-logs")

  max_bytes <- suppressWarnings(
    as.integer(.agent_log_scalar_setting("OMICSVIEWER_LLM_LOG_MAX_BYTES", "10485760"))
  )
  if (length(max_bytes) != 1L || is.na(max_bytes))
    max_bytes <- 10485760L
  max_bytes <- max(1024^2L, min(100L * 1024^2L, max_bytes))

  list(
    enabled = enabled,
    directory = directory,
    max_bytes = max_bytes
  )
}

.agent_log_session_label <- function(session) {
  token <- "unknown"
  tryCatch(
    {
      token <- session$token
      if (!is.null(token) && nzchar(token))
        token <- substr(token, 1L, 12L)
    },
    error = function(e) NULL
  )
  paste0(
    format(Sys.time(), "%Y%m%d-%H%M%S", tz = "UTC"),
    "-",
    gsub("[^A-Za-z0-9_-]+", "x", token)
  )
}

#' Create a session diagnostic logger
#'
#' The log file is reserved immediately but the directory and file are not
#' created until logging is enabled. This avoids creating diagnostic files for
#' sessions that never opt in.
#'
#' @param session Shiny session.
#' @param config Configuration from \code{\link{agent_logging_config}}.
#'
#' @return An environment-like internal logger object.
#' @keywords internal
#' @rdname agentLoggingHelpers
agent_logger_new <- function(session, config = agent_logging_config()) {
  directory <- config$directory
  file_name <- paste0("omicsviewer-llm-", .agent_log_session_label(session), ".jsonl")
  path <- file.path(directory, file_name)
  if (file.exists(path))
    path <- make.unique(path, sep = "-")

  logger <- new.env(parent = emptyenv())
  logger$session_label <- .agent_log_session_label(session)
  logger$directory <- directory
  logger$path <- path
  logger$max_bytes <- config$max_bytes
  logger$enabled <- FALSE
  logger$sequence <- 0L
  logger$limit_logged <- FALSE
  logger
}

#' Return a logger's reserved file path
#'
#' @param logger Logger created by \code{\link{agent_logger_new}}.
#' @return Current file path, or NA when logging is disabled and no file exists.
#' @keywords internal
#' @rdname agentLoggingHelpers
agent_logger_path <- function(logger) {
  if (logger$enabled || file.exists(logger$path))
    logger$path
  else
    NA_character_
}

.agent_log_redact_text <- function(x) {
  if (is.null(x) || length(x) == 0)
    return(as.character(x))
  x <- as.character(x)
  x <- gsub("(sk|rk)-[A-Za-z0-9_-]{10,}", "[REDACTED]", x, perl = TRUE)
  x <- gsub("Bearer[[:space:]]+[A-Za-z0-9._~+/=-]{10,}", "Bearer [REDACTED]", x, perl = TRUE)
  x
}

.agent_log_safe_value <- function(x, depth = 0L, max_chars = 20000L) {
  if (depth > 8L)
    return(list(truncated = TRUE, reason = "maximum nesting depth"))

  sensitive <- c(
    "api_key", "apikey", "key", "token", "secret", "password", "credential",
    "credentials", "authorization", "auth"
  )
  if (is.function(x) || inherits(x, "environment") || inherits(x, "R6") ||
      inherits(x, "S7_object") || inherits(x, "ggplot")) {
    return(list(
      object_class = paste(class(x), collapse = "/"),
      omitted = TRUE
    ))
  }

  if (is.raw(x))
    return(list(object_class = "raw", bytes = length(x), omitted = TRUE))

  if (is.data.frame(x)) {
    return(list(
      object_class = "data.frame",
      rows = nrow(x),
      columns = as.list(colnames(x))
    ))
  }

  if (is.list(x)) {
    if (length(x) > 200L) {
      return(list(
        count = length(x),
        values = lapply(utils::head(x, 100L), .agent_log_safe_value, depth = depth + 1L),
        truncated = TRUE
      ))
    }
    values <- lapply(seq_along(x), function(i) {
      # names(x)[[i]] is NULL for completely unnamed lists (e.g. the layer
      # array of a figure spec); guard so the redaction check cannot
      # evaluate to NA and abort the whole event.
      nm <- if (is.null(names(x))) "" else names(x)[[i]]
      if (!is.na(nm) && nzchar(nm) && tolower(nm) %in% sensitive)
        return("[REDACTED]")
      .agent_log_safe_value(x[[i]], depth = depth + 1L)
    })
    if (!is.null(names(x)))
      names(values) <- names(x)
    return(values)
  }

  if (length(x) > 100L) {
    return(list(
      count = length(x),
      values = as.list(utils::head(.agent_log_redact_text(x), 100L)),
      truncated = TRUE
    ))
  }

  if (is.character(x)) {
    x <- .agent_log_redact_text(x)
    too_long <- nzchar(x) & nchar(x, type = "chars", allowNA = TRUE) > max_chars
    if (any(too_long, na.rm = TRUE)) {
      x[too_long] <- paste0(substr(x[too_long], 1L, max_chars), " ... [truncated]")
    }
    return(x)
  }

  if (is.numeric(x) || is.logical(x) || is.factor(x)) {
    if (length(x) <= 20L)
      return(as.vector(x))
    return(list(
      count = length(x),
      values = as.list(utils::head(as.vector(x), 20L)),
      truncated = TRUE
    ))
  }

  if (is.null(x))
    return(NULL)

  .agent_log_redact_text(as.character(x))
}

#' Write one assistant diagnostic event
#'
#' Logging failures are intentionally silent so diagnostics can never break the
#' Shiny session or provider request.
#'
#' @param logger Logger from \code{\link{agent_logger_new}}.
#' @param event Event name.
#' @param details JSON-safe event details.
#' @return Invisible TRUE if an event was written; otherwise invisible FALSE.
#' @keywords internal
#' @rdname agentLoggingHelpers
agent_logger_event <- function(logger, event, details = list()) {
  tryCatch(
    {
      if (!isTRUE(logger$enabled))
        return(invisible(FALSE))

      if (!dir.exists(logger$directory))
        dir.create(logger$directory, recursive = TRUE, showWarnings = FALSE)
      if (!dir.exists(logger$directory))
        return(invisible(FALSE))

      if (file.exists(logger$path) && unname(file.size(logger$path)) >= logger$max_bytes) {
        if (!isTRUE(logger$limit_logged)) {
          logger$limit_logged <- TRUE
          logger$sequence <- logger$sequence + 1L
          line <- jsonlite::toJSON(
            list(
              timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
              sequence = logger$sequence,
              session_id = logger$session_label,
              event = "log_limit_reached",
              max_bytes = logger$max_bytes
            ),
            auto_unbox = TRUE,
            null = "null",
            na = "null"
          )
          cat(line, "\n", file = logger$path, append = TRUE, sep = "")
        }
        logger$enabled <- FALSE
        return(invisible(FALSE))
      }

      logger$sequence <- logger$sequence + 1L
      # Sanitization failures must never silently drop an event: unnamed
      # lists (e.g. the layer array of a figure spec) previously crashed
      # .agent_log_safe_value before this tryCatch, killing the caller's
      # hook and losing tool-request payloads. Fall back to a degraded
      # record that still names the event and the sanitization error.
      safe_details <- tryCatch(
        .agent_log_safe_value(details),
        error = function(e)
          list(sanitization_error = conditionMessage(e), details_omitted = TRUE)
      )
      record <- list(
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
        sequence = logger$sequence,
        session_id = logger$session_label,
        event = event,
        details = safe_details
      )
      line <- jsonlite::toJSON(
        record,
        auto_unbox = TRUE,
        null = "null",
        na = "null",
        dataframe = "columns"
      )
      cat(line, "\n", file = logger$path, append = TRUE, sep = "")
      invisible(TRUE)
    },
    error = function(e) {
      # Last-resort diagnostics on stderr so a failed write is visible in
      # the Shiny process log instead of vanishing.
      try(message("omicsViewer agent log: failed to write event ", event,
                  ": ", conditionMessage(e)), silent = TRUE)
      invisible(FALSE)
    }
  )
}

#' Enable or disable a session diagnostic logger
#'
#' @param logger Logger from \code{\link{agent_logger_new}}.
#' @param enabled Single logical value.
#' @return The logger, invisibly.
#' @keywords internal
#' @rdname agentLoggingHelpers
agent_logger_set_enabled <- function(logger, enabled) {
  enabled <- isTRUE(enabled)
  if (identical(isTRUE(logger$enabled), enabled))
    return(invisible(logger))

  if (enabled) {
    logger$enabled <- TRUE
    agent_logger_event(
      logger,
      "logging_enabled",
      list(log_file = logger$path, max_bytes = logger$max_bytes)
    )
  } else {
    agent_logger_event(logger, "logging_disabled")
    logger$enabled <- FALSE
  }
  invisible(logger)
}

.agent_log_condition <- function(x) {
  if (is.null(x))
    return(NULL)
  list(
    class = paste(class(x), collapse = "/"),
    message = .agent_log_safe_value(conditionMessage(x)),
    omitted_fields = TRUE
  )
}

#' Summarize one ellmer turn for diagnostics
#'
#' Provider-specific raw JSON and rich HTML dependencies are omitted. Tool
#' result display payloads are also omitted to avoid writing large base64
#' figures into diagnostic logs.
#'
#' @param turn An ellmer Turn object.
#' @return A JSON-safe turn summary.
#' @keywords internal
#' @rdname agentLoggingHelpers
agent_log_turn <- function(turn) {
  content_summary <- function(x) {
    out <- list(content_class = paste(class(x), collapse = "/"))
    if (inherits(x, "ellmer::ContentToolRequest")) {
      out$id <- .agent_log_safe_value(x@id)
      out$name <- .agent_log_safe_value(x@name)
      out$arguments <- .agent_log_safe_value(x@arguments)
      return(out)
    }
    if (inherits(x, "ellmer::ContentToolResult")) {
      if (!is.null(x@request)) {
        out$tool_call_id <- .agent_log_safe_value(x@request@id)
        out$tool_name <- .agent_log_safe_value(x@request@name)
      }
      out$error <- if (is.null(x@error)) NULL else .agent_log_condition(x@error)
      out$value <- .agent_log_safe_value(x@value)
      out$display_payload_omitted <- TRUE
      return(out)
    }

    out$text <- tryCatch(
      .agent_log_safe_value(ellmer::contents_text(x)),
      error = function(e) NULL
    )
    out
  }

  out <- list(
    role = .agent_log_safe_value(turn@role),
    text = .agent_log_safe_value(turn@text),
    contents = lapply(turn@contents, content_summary)
  )
  if (inherits(turn, "ellmer::AssistantTurn")) {
    out$tokens <- .agent_log_safe_value(turn@tokens)
    out$cost <- .agent_log_safe_value(turn@cost)
    out$duration <- .agent_log_safe_value(turn@duration)
    out$finish_reason <- .agent_log_safe_value(turn@finish_reason)
    out$provider_json_omitted <- TRUE
  }
  out
}

#' Summarize assistant provider configuration without credentials
#'
#' @param config Provider configuration.
#' @return JSON-safe provider/model information.
#' @keywords internal
#' @rdname agentLoggingHelpers
agent_log_provider_config <- function(config) {
  list(
    provider = config$provider,
    model = if (nzchar(config$model)) config$model else "provider-default",
    base_url = if (nzchar(config$base_url)) config$base_url else "provider-default",
    configured = isTRUE(config$configured),
    key_source = config$key_source,
    api_key_omitted = TRUE
  )
}
