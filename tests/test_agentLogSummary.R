library(omicsViewer)
library(unittest, quietly = TRUE)

agent_summarize_log <- omicsViewer:::agent_summarize_log
agent_summarize_logs <- omicsViewer:::agent_summarize_logs
agent_log_error_class <- omicsViewer:::agent_log_error_class

# ---- pure classifier --------------------------------------------------
ok(
  ut_cmp_identical(agent_log_error_class("Unknown feature ID(s): Foo, Bar."), "unknown_id"),
  "classifier: unknown feature IDs"
)
ok(
  ut_cmp_identical(agent_log_error_class("Unknown feature quick view: volacno."), "unknown_id"),
  "classifier: unknown quick view"
)
ok(
  ut_cmp_identical(agent_log_error_class("Unknown feature X-axis annotation: log.fdrr."), "unknown_column"),
  "classifier: unknown axis annotation"
)
ok(
  ut_cmp_identical(agent_log_error_class("Unknown sample annotation column: groupp."), "unknown_column"),
  "classifier: unknown annotation column"
)
ok(
  ut_cmp_identical(agent_log_error_class("Unknown analysis-space tab: null."), "unknown_tab"),
  "classifier: unknown tab (real Tier B message)"
)
ok(
  ut_cmp_identical(agent_log_error_class("Figure layer point requires mapping(s): x, y"), "invalid_figure_spec"),
  "classifier: invalid figure spec"
)
ok(
  ut_cmp_identical(agent_log_error_class("max_results must be an integer from 1 through 50."), "invalid_argument"),
  "classifier: generic invalid argument"
)
ok(
  ut_cmp_identical(agent_log_error_class("No dataset is currently available to the assistant."), "no_dataset"),
  "classifier: missing dataset"
)
ok(
  ut_cmp_identical(
    agent_log_error_class("Assistant request limit reached for this session (40)."),
    "request_limit"
  ),
  "classifier: request limit (checked before provider failure)"
)
ok(
  ut_cmp_identical(agent_log_error_class("HTTP 429 Too Many Requests"), "provider_failure"),
  "classifier: provider failure"
)
ok(
  ut_cmp_identical(agent_log_error_class(""), "unknown"),
  "classifier: empty message"
)

# ---- synthetic JSONL fixture ------------------------------------------
# Covers every taxonomy class plus one successful retry sequence.
log_dir <- file.path(tempdir(), paste0("agent-log-summary-", Sys.getpid()))
if (dir.exists(log_dir)) unlink(log_dir, recursive = TRUE)
dir.create(log_dir, recursive = TRUE)
log_path <- file.path(log_dir, "omicsviewer-llm-test-session.jsonl")

.event <- function(sequence, seconds, event, details = list()) {
  jsonlite::toJSON(
    list(
      timestamp = sprintf("2026-09-24T10:%02d:%06.3fZ", seconds %/% 60, seconds %% 60),
      sequence = sequence,
      session_id = "test-session",
      event = event,
      details = details
    ),
    auto_unbox = TRUE, null = "null", na = "null"
  )
}
.err <- function(message) list(class = "simpleError/error/condition", message = message)

writeLines(c(
  .event(1, 0, "logging_enabled", list(log_file = log_path)),
  .event(2, 0, "assistant_session_start",
         list(dependencies_available = TRUE, request_limit = 40)),
  .event(3, 1, "user_message", list(message = "switch to the heatmap tab")),
  # taxonomy: unknown_tab + a successful retry of the SAME tool later
  .event(4, 2, "tool_request", list(
    request_index = 1, tool_call_id = "call_1", tool_name = "set_omics_viewer_state",
    arguments = list(data_space_tab = "Heatmep", `_intent` = "switch tab"))),
  .event(5, 2, "tool_result", list(
    request_index = 1, tool_call_id = "call_1", tool_name = "set_omics_viewer_state",
    error = .err("Unknown data-space tab: Heatmep. Closest matches: Heatmap."))),
  .event(6, 4, "tool_request", list(
    request_index = 1, tool_call_id = "call_2", tool_name = "set_omics_viewer_state",
    arguments = list(data_space_tab = "Heatmap", `_intent` = "switch tab"))),
  .event(7, 4, "tool_result", list(
    request_index = 1, tool_call_id = "call_2", tool_name = "set_omics_viewer_state")),
  # taxonomy: unknown_column
  .event(8, 6, "tool_request", list(
    request_index = 2, tool_call_id = "call_3", tool_name = "set_scatter_view",
    arguments = list(space = "feature", x_axis = "log.fdrr",
                     y_axis = "ttest|A_vs_B|log.fdr", `_intent` = "volcano"))),
  .event(9, 6, "tool_result", list(
    request_index = 2, tool_call_id = "call_3", tool_name = "set_scatter_view",
    error = .err("Unknown feature X-axis annotation: log.fdrr."))),
  # taxonomy: unknown_id
  .event(10, 8, "tool_request", list(
    request_index = 3, tool_call_id = "call_4", tool_name = "set_omics_viewer_state",
    arguments = list(features = "GeneN", `_intent` = "select gene"))),
  .event(11, 8, "tool_result", list(
    request_index = 3, tool_call_id = "call_4", tool_name = "set_omics_viewer_state",
    error = .err("Unknown feature ID(s): GeneN."))),
  # taxonomy: invalid_figure_spec
  .event(12, 10, "tool_request", list(
    request_index = 4, tool_call_id = "call_5", tool_name = "create_figure",
    arguments = list(spec = list(geom = "point"), `_intent` = "plot"))),
  .event(13, 10, "tool_result", list(
    request_index = 4, tool_call_id = "call_5", tool_name = "create_figure",
    error = .err("Figure layer point requires mapping(s): x, y"))),
  # taxonomy: invalid_argument
  .event(14, 12, "tool_request", list(
    request_index = 5, tool_call_id = "call_6", tool_name = "search_annotations",
    arguments = list(space = "feature", query = "gene", max_results = 999,
                     `_intent` = "find genes"))),
  .event(15, 12, "tool_result", list(
    request_index = 5, tool_call_id = "call_6", tool_name = "search_annotations",
    error = .err("max_results must be an integer from 1 through 50."))),
  # taxonomy: no_dataset
  .event(16, 14, "tool_request", list(
    request_index = 6, tool_call_id = "call_7", tool_name = "get_omics_viewer_state",
    arguments = list(`_intent` = "look at state"))),
  .event(17, 14, "tool_result", list(
    request_index = 6, tool_call_id = "call_7", tool_name = "get_omics_viewer_state",
    error = .err("No dataset is currently available to the assistant."))),
  # taxonomy: request_limit
  .event(18, 16, "tool_request", list(
    request_index = 41, tool_call_id = "call_8", tool_name = "list_widgets",
    arguments = list(`_intent` = "list widgets"))),
  .event(19, 16, "tool_result", list(
    request_index = 41, tool_call_id = "call_8", tool_name = "list_widgets",
    error = .err("Assistant request limit reached for this session (40)."))),
  # taxonomy: provider_failure (stream_failure, not a tool result)
  .event(20, 18, "stream_failure",
         list(error = .err("HTTP 429 Too Many Requests"))),
  # unresolved request (no matching result) and an orphan result
  .event(21, 20, "tool_request", list(
    request_index = 42, tool_call_id = "call_9", tool_name = "get_widget",
    arguments = list(id = "dataspace.expr_heatmap.heatmap_colors",
                     `_intent` = "describe widget"))),
  .event(22, 22, "tool_result", list(
    request_index = 42, tool_call_id = "call_missing", tool_name = "get_widget")),
  .event(23, 24, "assistant_session_end"),
  # one malformed trailing line
  "{not valid json"
), log_path)

summary <- agent_summarize_log(log_path)

ok(
  ut_cmp_identical(summary$events_total, 24L),
  "summary counts all lines including malformed"
)
ok(
  ut_cmp_identical(summary$malformed_lines, 1L),
  "summary counts malformed lines"
)
ok(
  ut_cmp_identical(summary$event_counts$user_message, 1L) &&
    ut_cmp_identical(summary$event_counts$tool_request, 9L) &&
    ut_cmp_identical(summary$event_counts$tool_result, 9L),
  "summary counts events by type"
)
ok(
  ut_cmp_identical(summary$duration_seconds, 24),
  "summary derives the session duration from timestamps"
)
ok(
  ut_cmp_identical(summary$tool_calls_total, 8L) &&
    ut_cmp_identical(summary$tool_calls_per_tool$set_omics_viewer_state, 3L),
  "summary counts tool calls per tool"
)
ok(
  ut_cmp_identical(summary$tool_failures, 7L),
  "summary counts failed tool results"
)
ok(
  ut_cmp_identical(summary$first_attempt_success_rate, 0.125),
  "summary computes the first-attempt success rate"
)
ok(
  ut_cmp_identical(
    summary$error_taxonomy,
    list(
      invalid_argument = 1L, invalid_figure_spec = 1L, no_dataset = 1L,
      provider_failure = 1L, request_limit = 1L, unknown_column = 1L,
      unknown_id = 1L, unknown_tab = 1L
    )
  ),
  "summary classifies every taxonomy class exactly once"
)
ok(
  ut_cmp_identical(summary$recovered_retries, 1L) &&
    ut_cmp_identical(summary$retry_recovery_rate, 0.1429),
  "summary detects the successful same-tool retry"
)
ok(
  ut_cmp_identical(summary$most_rejected_arguments$space, 2L) &&
    is.null(summary$most_rejected_arguments$`_intent`) &&
    ut_cmp_identical(summary$most_rejected_arguments$data_space_tab, 1L),
  "summary tabulates most-rejected argument keys (intent excluded)"
)
ok(
  ut_cmp_identical(summary$unresolved_requests, 1L) &&
    ut_cmp_identical(summary$orphan_results, 1L),
  "summary reports unresolved requests and orphan results"
)
ok(
  ut_cmp_identical(summary$stream_failures, 1L),
  "summary counts stream failures"
)
ok(
  ut_cmp_identical(summary$failure_examples[[1]]$class, "unknown_tab"),
  "summary keeps bounded failure examples with classes"
)

# retry window: a same-tool success beyond the window does not count
narrow <- agent_summarize_log(log_path, retry_window = 1L)
ok(
  ut_cmp_identical(narrow$recovered_retries, 0L),
  "recovery requires the retry within the sequence window"
)

# print method runs and mentions the key facts
printed <- paste(capture.output(print(summary)), collapse = "\n")
ok(
  grepl("first-attempt success", printed, fixed = TRUE) &&
    grepl("error taxonomy", printed, fixed = TRUE) &&
    grepl("unknown_tab=1", printed, fixed = TRUE),
  "print method renders a console-friendly overview"
)

# ---- directory wrapper -------------------------------------------------
writeLines(
  c(.event(1, 0, "logging_enabled"), .event(2, 1, "assistant_session_end")),
  file.path(log_dir, "omicsviewer-llm-second-session.jsonl")
)
all_logs <- agent_summarize_logs(log_dir)
ok(
  ut_cmp_identical(sort(names(all_logs)),
                   sort(c("omicsviewer-llm-test-session.jsonl",
                          "omicsviewer-llm-second-session.jsonl"))),
  "directory wrapper summarizes every log"
)
ok(
  ut_cmp_identical(all_logs$`omicsviewer-llm-second-session.jsonl`$events_total, 2L),
  "directory wrapper returns per-file summaries"
)
ok(
  ut_cmp_error(agent_summarize_logs(tempfile("no-such-dir")), "No assistant .jsonl logs found"),
  "directory wrapper errors clearly without logs"
)
ok(
  ut_cmp_error(agent_summarize_log(tempfile("no-such-file", fileext = ".jsonl")),
               "Log file not found"),
  "missing file errors clearly"
)

# ---- real archived fixtures still parse --------------------------------
evidence <- file.path("tests", "e2e_agent", "artifacts", "decision-evidence-20260923")
if (dir.exists(evidence)) {
  real <- agent_summarize_log(
    file.path(evidence, "omicsviewer-llm-20260921-211247-5a5dbfd57f42.jsonl")
  )
  ok(
    ut_cmp_identical(real$malformed_lines, 0L) &&
      ut_cmp_identical(real$tool_calls_total, 6L),
    "real archived Tier B log parses with the expected call count"
  )
  ok(
    ut_cmp_identical(real$error_taxonomy$unknown_tab, 2L) &&
      ut_cmp_identical(real$most_rejected_arguments$analysis_space_tab, 2L),
    "real archived log classifies the null-tab failures and rejected args"
  )
} else {
  ok(TRUE, "archived fixtures not present; real-log check skipped")
}
