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

#' Classify one assistant error message into the failure taxonomy
#'
#' Ordered regex rules (first match wins) over the lower-cased message.
#' The classes mirror the validation surface of the agent tools
#' (plan section 3, WP4): unknown identifiers, unknown annotation
#' columns, unknown tabs, invalid figure specs, other invalid arguments,
#' missing datasets, the session request limit, the WP12 session
#' budget ceilings (token/cost), and provider-side failures (observed
#' on \code{stream_failure} events).
#'
#' @param message Error message (character, possibly length > 1).
#' @return Single class label.
#' @keywords internal
#' @rdname agentLoggingHelpers
agent_log_error_class <- function(message) {
  msg <- tolower(paste(trimws(as.character(message)), collapse = " "))
  if (!nzchar(msg))
    return("unknown")
  if (grepl("request limit", msg, fixed = TRUE))
    return("request_limit")
  if (grepl("token budget|cost budget", msg))
    return("budget_limit")
  if (grepl("no dataset|metadata are loaded|currently available", msg))
    return("no_dataset")
  if (grepl("provider|http 4|http 5|429|401|403|timeout|timed out|connection|api key|rate limit", msg))
    return("provider_failure")
  if (grepl("unknown .*tab|(data|analysis)[- ]space tab:", msg))
    return("unknown_tab")
  if (grepl("unknown .*(annotation|column)", msg))
    return("unknown_column")
  if (grepl("unknown .*(id|quick view|widget|figure)", msg))
    return("unknown_id")
  if (grepl("figure|geom|layer", msg))
    return("invalid_figure_spec")
  "invalid_argument"
}

.agent_log_na_time <- as.POSIXct(NA_character_, tz = "UTC")
.agent_log_origin <- as.POSIXct("1970-01-01", tz = "UTC")

.agent_log_parse_timestamp <- function(x) {
  x <- trimws(as.character(x)[1])
  if (is.na(x) || !nzchar(x))
    return(.agent_log_na_time)
  as.POSIXct(sub("Z$", "", x), format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")
}

.agent_log_summary_simplify <- function(counts, max_entries = 12L) {
  if (!length(counts))
    return(list())
  counts <- counts[order(-as.integer(counts), names(counts))]
  as.list(utils::head(counts, max_entries))
}

#' Summarize one assistant diagnostic JSONL log
#'
#' Developer tool (plan section 3, WP4): parses one session log and
#' returns event counts, session duration, per-tool call counts, the
#' error taxonomy for failed tool results (plus \code{stream_failure}
#' provider errors), the first-attempt success rate, retry-recovery
#' outcome (same tool re-invoked within \code{retry_window} sequence
#' numbers and succeeding), and the most-rejected argument keys.
#'
#' @param path Path to a \code{.jsonl} assistant log.
#' @param retry_window Sequence-number window in which a same-tool retry
#'   counts as a recovery attempt (default 10).
#' @param max_examples Maximum failed-call examples kept in the summary.
#' @return A classed list with the summary fields; print it for a
#'   console-friendly overview.
#' @keywords internal
#' @rdname agentLoggingHelpers
agent_summarize_log <- function(path, retry_window = 10L, max_examples = 8L) {
  if (!file.exists(path))
    stop("Log file not found: ", path)
  retry_window <- suppressWarnings(as.integer(retry_window)[1])
  if (is.na(retry_window) || retry_window < 1L)
    retry_window <- 10L

  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]

  events <- vector("list", length(lines))
  malformed <- 0L
  for (i in seq_along(lines)) {
    events[[i]] <- tryCatch(
      jsonlite::fromJSON(lines[[i]], simplifyVector = FALSE),
      error = function(e) {
        malformed <<- malformed + 1L
        NULL
      }
    )
  }
  events <- Filter(function(x) !is.null(x) && !is.null(x$event), events)

  event_counts <- table(vapply(events, function(x) as.character(x$event)[1], character(1)))
  timestamps <- vapply(events, function(x) as.character(x$timestamp %||% ""), character(1))
  parsed_times <- unlist(lapply(timestamps[timestamps != ""], .agent_log_parse_timestamp))
  session_start <- if (length(parsed_times) && !all(is.na(parsed_times)))
    .agent_log_origin + min(as.numeric(parsed_times), na.rm = TRUE) else .agent_log_na_time
  session_end <- if (length(parsed_times) && !all(is.na(parsed_times)))
    .agent_log_origin + max(as.numeric(parsed_times), na.rm = TRUE) else .agent_log_na_time

  seq_of <- vapply(events, function(x) {
    v <- suppressWarnings(as.integer(x$sequence)[1])
    if (is.na(v)) 0L else v
  }, integer(1))

  # join tool requests with their results by tool_call_id
  pending <- list()
  pairs <- list()
  orphan_results <- 0L
  unresolved_requests <- 0L
  for (i in seq_along(events)) {
    ev <- events[[i]]$event
    d <- events[[i]]$details
    if (identical(ev, "tool_request")) {
      id <- as.character(d$tool_call_id %||% "")
      if (!nzchar(id)) {
        unresolved_requests <- unresolved_requests + 1L
        next
      }
      pending[[id]] <- list(
        index = i, sequence = seq_of[[i]],
        tool_name = as.character(d$tool_name %||% ""),
        arguments = if (is.list(d$arguments)) d$arguments else list()
      )
    } else if (identical(ev, "tool_result")) {
      id <- as.character(d$tool_call_id %||% "")
      req <- pending[[id]]
      if (is.null(req)) {
        orphan_results <- orphan_results + 1L
        next
      }
      pending[[id]] <- NULL
      error <- d$error
      pairs[[length(pairs) + 1L]] <- list(
        index = i, sequence = seq_of[[i]],
        tool_name = req$tool_name %||% as.character(d$tool_name %||% ""),
        tool_call_id = id, arguments = req$arguments,
        error = !is.null(error),
        error_class = if (is.null(error)) NULL else
          agent_log_error_class(unlist(error$message)),
        error_message = if (is.null(error)) NULL else
          .agent_shorten(paste(unlist(error$message), collapse = " "), 160L)
      )
    }
  }
  unresolved_requests <- unresolved_requests + length(pending)

  tool_counts <- table(vapply(pairs, function(x) x$tool_name, character(1)))
  failed <- Filter(function(x) isTRUE(x$error), pairs)
  succeeded <- Filter(function(x) !isTRUE(x$error), pairs)

  taxonomy <- table(vapply(failed, function(x) x$error_class, character(1)))
  # provider-side failures surface on stream_failure events, not tool
  # results; classify them into the same taxonomy so one view covers all
  # failure modes.
  stream_failures <- Filter(function(x) identical(x$event, "stream_failure"), events)
  stream_classes <- vapply(stream_failures, function(x) {
    msg <- x$details$error$message
    if (is.null(msg)) "unknown" else agent_log_error_class(unlist(msg))
  }, character(1))
  if (length(stream_classes))
    taxonomy <- table(c(
      vapply(failed, function(x) x$error_class, character(1)),
      stream_classes
    ))

  # first-attempt success: per pair (a tool_call_id is one attempt)
  per_tool_rate <- vapply(split(
    vapply(pairs, function(x) !isTRUE(x$error), logical(1)),
    vapply(pairs, function(x) x$tool_name, character(1))
  ), mean, numeric(1))

  # retry recovery: after a failed result, a same-tool request whose own
  # result succeeds within the sequence window
  recovered <- 0L
  for (f in failed) {
    hit <- FALSE
    for (p in pairs) {
      if (p$sequence > f$sequence &&
          identical(p$tool_name, f$tool_name) &&
          !isTRUE(p$error) &&
          (p$sequence - f$sequence) <= retry_window) {
        hit <- TRUE
        break
      }
    }
    if (hit) recovered <- recovered + 1L
  }

  rejected_args <- character()
  for (f in failed) {
    keys <- setdiff(names(f$arguments), "_intent")
    rejected_args <- c(rejected_args, keys)
  }

  examples <- lapply(utils::head(failed, max_examples), function(f)
    list(sequence = f$sequence, tool = f$tool_name, class = f$error_class,
         message = f$error_message))

  structure(
    list(
      path = path,
      session_id = if (length(events)) as.character(events[[1]]$session_id %||% "") else "",
      events_total = length(events) + malformed,
      malformed_lines = malformed,
      event_counts = as.list(event_counts[order(-as.integer(event_counts))]),
      session_start = session_start,
      session_end = session_end,
      duration_seconds = if (!is.na(session_start) && !is.na(session_end))
        round(as.numeric(difftime(session_end, session_start, units = "secs")), 3) else NA_real_,
      tool_calls_total = length(pairs),
      tool_calls_per_tool = as.list(tool_counts[order(-as.integer(tool_counts))]),
      unresolved_requests = unresolved_requests,
      orphan_results = orphan_results,
      tool_failures = length(failed),
      stream_failures = length(stream_failures),
      first_attempt_success_rate = if (length(pairs))
        round(length(succeeded) / length(pairs), 4) else NA_real_,
      first_attempt_success_per_tool = as.list(round(per_tool_rate, 4)),
      error_taxonomy = .agent_log_summary_simplify(taxonomy, max_entries = 12L),
      retry_window = retry_window,
      recovered_retries = if (length(failed)) recovered else 0L,
      retry_recovery_rate = if (length(failed)) round(recovered / length(failed), 4) else NA_real_,
      most_rejected_arguments = .agent_log_summary_simplify(table(rejected_args)),
      failure_examples = examples
    ),
    class = c("omicsViewerAgentLogSummary", "list")
  )
}

#' @method print omicsViewerAgentLogSummary
#' @keywords internal
#' @rdname agentLoggingHelpers
#' @export
print.omicsViewerAgentLogSummary <- function(x, ...) {
  cat("omicsViewer agent log summary\n")
  cat("  file:            ", x$path, "\n", sep = "")
  cat("  events:          ", x$events_total,
      if (x$malformed_lines) paste0("(", x$malformed_lines, " malformed)") else "",
      "\n", sep = " ")
  if (!is.na(x$session_start))
    cat("  session:         ", format(x$session_start),
        "->", format(x$session_end),
        "(", x$duration_seconds, "s)\n", sep = " ")
  cat("  tool calls:      ", x$tool_calls_total,
      "| failed:", x$tool_failures,
      "| first-attempt success:",
      if (is.na(x$first_attempt_success_rate)) "NA"
      else sprintf("%.0f%%", 100 * x$first_attempt_success_rate), "\n", sep = " ")
  if (length(x$tool_calls_per_tool)) {
    cat("  calls per tool:  ",
        paste(names(x$tool_calls_per_tool), x$tool_calls_per_tool, sep = "=", collapse = ", "),
        "\n", sep = "")
  }
  if (length(x$error_taxonomy)) {
    cat("  error taxonomy:  ",
        paste(names(x$error_taxonomy), x$error_taxonomy, sep = "=", collapse = ", "),
        "\n", sep = "")
  }
  if (x$tool_failures)
    cat("  retry recovery:  ", x$recovered_retries, "/", x$tool_failures,
        sprintf("(%.0f%%)", 100 * (x$retry_recovery_rate %||% 0)),
        "within", x$retry_window, "events\n", sep = " ")
  if (length(x$most_rejected_arguments)) {
    cat("  rejected args:   ",
        paste(names(x$most_rejected_arguments), x$most_rejected_arguments,
              sep = "=", collapse = ", "),
        "\n", sep = "")
  }
  for (ex in utils::head(x$failure_examples, 5L)) {
    cat("  failure #", ex$sequence, " [", ex$class, "] ", ex$tool, ": ",
        ex$message, "\n", sep = "")
  }
  invisible(x)
}

#' Summarize every assistant diagnostic log in a directory
#'
#' Convenience wrapper over \code{\link{agent_summarize_log}} for the
#' configured (or given) log directory.
#'
#' @param directory Directory containing \code{.jsonl} logs; defaults to
#'   the configured diagnostic directory.
#' @param retry_window See \code{\link{agent_summarize_log}}.
#' @return Named list of summaries keyed by file name.
#' @keywords internal
#' @rdname agentLoggingHelpers
agent_summarize_logs <- function(directory = agent_logging_config()$directory,
                                 retry_window = 10L) {
  if (!dir.exists(directory) || !length(list.files(directory, pattern = "[.]jsonl$")))
    stop("No assistant .jsonl logs found in: ", directory)
  files <- sort(list.files(directory, pattern = "[.]jsonl$", full.names = TRUE))
  stats::setNames(
    lapply(files, agent_summarize_log, retry_window = retry_window),
    basename(files)
  )
}
