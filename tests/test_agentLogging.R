library(omicsViewer)
library(unittest, quietly = TRUE)
library(ellmer)

agent_logging_config <- omicsViewer:::agent_logging_config
agent_logger_new <- omicsViewer:::agent_logger_new
agent_logger_set_enabled <- omicsViewer:::agent_logger_set_enabled
agent_logger_event <- omicsViewer:::agent_logger_event
agent_logger_path <- omicsViewer:::agent_logger_path
agent_log_turn <- omicsViewer:::agent_log_turn

log_dir <- file.path(tempdir(), paste0("agent-log-", Sys.getpid()))
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
old_log <- Sys.getenv("OMICSVIEWER_LLM_LOG")
old_log_dir <- Sys.getenv("OMICSVIEWER_LLM_LOG_DIR")
Sys.setenv(OMICSVIEWER_LLM_LOG = "true", OMICSVIEWER_LLM_LOG_DIR = log_dir)
on.exit(Sys.setenv(
  OMICSVIEWER_LLM_LOG = old_log,
  OMICSVIEWER_LLM_LOG_DIR = old_log_dir
), add = TRUE)

session <- shiny::MockShinySession$new()
config <- agent_logging_config()
logger <- agent_logger_new(session, config)
agent_logger_set_enabled(logger, TRUE)

ok(ut_cmp_identical(config$enabled, TRUE), "environment enables assistant logging")
ok(
  ut_cmp_identical(startsWith(agent_logger_path(logger), log_dir), TRUE),
  "logger reserves a file in the configured directory"
)

agent_logger_event(
  logger,
  "unit_test_event",
  list(
    message = "hello",
    api_key = "secret-value",
    token = "secret-token",
    text = "key sk-ABCDEFGHIJKLMNOPQRSTUVWXYZ123456",
    values = 1:150,
    nested = list(ok = TRUE)
  )
)

lines <- readLines(agent_logger_path(logger), warn = FALSE)
record <- jsonlite::fromJSON(lines[[length(lines)]], simplifyVector = FALSE)
ok(ut_cmp_identical(record$event, "unit_test_event"), "diagnostic events are JSONL records")
ok(ut_cmp_identical(record$details$message, "hello"), "ordinary event text is retained")
ok(
  ut_cmp_identical(record$details$api_key, "[REDACTED]"),
  "credential-like field names are redacted"
)
ok(
  ut_cmp_identical(record$details$token, "[REDACTED]"),
  "token fields are redacted"
)
ok(
  ut_cmp_identical(grepl("sk-ABCDEFGHIJKLMNOPQRSTUVWXYZ", lines[[length(lines)]], fixed = TRUE), FALSE),
  "credential-like strings are redacted"
)
ok(ut_cmp_identical(record$details$values$truncated, TRUE), "long vectors are bounded")
ok(ut_cmp_identical(record$details$nested$ok, TRUE), "nested ordinary values are retained")

request <- ellmer::ContentToolRequest(
  id = "call-1",
  name = "search_annotations",
  arguments = list(space = "feature", query = "gene")
)
result <- ellmer::ContentToolResult(
  value = list(matching_ids = c("f1", "f2")),
  request = request
)
turn_summary <- agent_log_turn(
  ellmer::UserTurn(list(request, result, ellmer::ContentText("done")))
)
ok(
  ut_cmp_identical(turn_summary$role, "user"),
  "ellmer turn role is logged"
)
ok(
  ut_cmp_identical(turn_summary$contents[[1]]$name, "search_annotations"),
  "tool request names are logged"
)
ok(
  ut_cmp_identical(turn_summary$contents[[2]]$value$matching_ids, c("f1", "f2")),
  "bounded tool result values are logged"
)
ok(
  ut_cmp_identical(turn_summary$contents[[2]]$display_payload_omitted, TRUE),
  "rich tool display payloads are omitted"
)
ok(
  ut_cmp_identical(turn_summary$contents[[3]]$text, "done"),
  "ordinary turn text is logged"
)

agent_logger_set_enabled(logger, FALSE)
agent_logger_event(logger, "should_not_be_written")
lines_after_disable <- readLines(agent_logger_path(logger), warn = FALSE)
ok(
  ut_cmp_identical(any(grepl("should_not_be_written", lines_after_disable, fixed = TRUE)), FALSE),
  "disabled loggers write no further events"
)
ok(
  ut_cmp_identical(tail(lines_after_disable, 1) |> (\(x) jsonlite::fromJSON(x[[1]])$event)(), "logging_disabled"),
  "disabling logger is recorded before logging stops"
)

if (requireNamespace("shinychat", quietly = TRUE)) {
  module_log_dir <- file.path(tempdir(), paste0("agent-log-module-", Sys.getpid()))
  dir.create(module_log_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(OMICSVIEWER_LLM_LOG = "true", OMICSVIEWER_LLM_LOG_DIR = module_log_dir)

  fd <- data.frame(score = 1:4, row.names = paste0("F", 1:4))
  pd <- data.frame(group = c("a", "b"), row.names = paste0("S", 1:2))
  mat <- matrix(1:8, nrow = 4, dimnames = list(rownames(fd), rownames(pd)))
  compact_state <- list(
    dataset = list(id = "d", class = "test", dimensions = c(features = 4L, samples = 2L)),
    active_tabs = list(), selection = list(), available_tabs = list(),
    annotations = list(), quick_views = list(), panels = list(),
    figure_grammar = list(), state_policy = "widget-only"
  )

  shiny::testServer(
    omicsViewer:::ai_assistant_module,
    args = list(
      state = shiny::reactive(compact_state),
      state_available = shiny::reactive(TRUE),
      feature_data = shiny::reactive(fd),
      sample_data = shiny::reactive(pd),
      expression_data = shiny::reactive(mat),
      selected_features = shiny::reactive(character()),
      selected_samples = shiny::reactive(character()),
      apply_state = function(x) list(feature_count = 0L, sample_count = 0L),
      apply_scatter_view = function(...) list()
    ),
    expr = {
      module_log_path <- logger$path
      ok(
        ut_cmp_identical(logger$enabled, TRUE),
        "assistant module starts with environment-enabled diagnostics"
      )
      session$setInputs(enable_logging = FALSE)
      session$flushReact()
      session$flushOutput()
      ok(
        ut_cmp_identical(logger$enabled, FALSE),
        "assistant logging can be disabled from the UI"
      )
    }
  )

  module_file <- list.files(module_log_dir, full.names = TRUE)[1]
  module_lines <- readLines(module_file, warn = FALSE)
  module_events <- vapply(
    module_lines,
    function(line) jsonlite::fromJSON(line)$event,
    character(1)
  )
  ok(
    ut_cmp_identical("assistant_session_start" %in% module_events, TRUE),
    "assistant lifecycle events are logged"
  )
  ok(
    ut_cmp_identical("logging_disabled" %in% module_events, TRUE),
    "disabling logging from the UI is recorded"
  )
  ok(
    ut_cmp_identical(any(grepl("[REDACTED]", module_lines, fixed = TRUE)), FALSE),
    "no credential placeholders are needed when no credential is supplied"
  )
}
