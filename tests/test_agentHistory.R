# WP11 unit suite: conversation-in-snapshot helpers (auxi_agentAssistant.R)
# and the assistant-module API (snapshot_payload / restore_history).

library(shiny)
library(omicsViewer)
library(unittest, quietly = TRUE)

if (!requireNamespace("ellmer", quietly = TRUE) ||
    !requireNamespace("shinychat", quietly = TRUE)) {
  ok(TRUE, "AI history tests skipped because optional packages are unavailable")
  quit(save = "no", status = 0)
}

library(ellmer)

agent_history_redact <- omicsViewer:::agent_history_redact
agent_transcript_records <- omicsViewer:::agent_transcript_records
agent_history_payload <- omicsViewer:::agent_history_payload
agent_history_restore_payload <- omicsViewer:::agent_history_restore_payload

mk_turns <- function() {
  # ellmer places tool RESULTS in UserTurns (the tool loop appends the
  # result payload to the user side); the assistant turn carries the
  # requests. Fixtures mirror the real placement so downstream consumers
  # (stubber, transcript, snapshot) are exercised on realistic shapes.
  list(
    UserTurn(contents = list(ContentText("Please analyze sk-abcdef0123456789abcdef"))),
    AssistantTurn(contents = list(
      ContentText("Calling a tool"),
      ContentToolRequest("call1", "get_omics_viewer_state", list(`_intent` = "x"))
    )),
    UserTurn(contents = list(
      ContentToolResult(
        value = list(dataset = list(id = "demo")),
        request = ellmer::ContentToolRequest(
          "call1", "get_omics_viewer_state", list(`_intent` = "x")),
        extra = list(display = list(html = "BASE64PREVIEWPLACEHOLDER"))
      )
    )),
    AssistantTurn(contents = list(ContentText("Done. Here is the summary.")))
  )
}

## ------------------------------------------------------------ redaction ----
ok(
  identical(
    agent_history_redact("key sk-abcdef0123456789abcdef here"),
    "key [redacted] here"
  ),
  "OpenAI-style keys are redacted"
)
ok(
  grepl("\\[redacted\\]",
        agent_history_redact("Authorization: Bearer abcdefghijklmnopqr")) &&
    !grepl("abcdefghijklmnopqr",
           agent_history_redact("Authorization: Bearer abcdefghijklmnopqr")),
  "bearer tokens are redacted"
)
ok(
  identical(agent_history_redact("plain text"), "plain text"),
  "ordinary text passes through unchanged"
)

## ----------------------------------------------------------- transcript ----
turns <- mk_turns()
records <- agent_transcript_records(turns)
ok(
  identical(vapply(records, function(r) r$role, character(1)),
            c("user", "assistant", "user", "assistant")),
  "transcript records carry user/assistant roles in order"
)
ok(
  grepl("get_omics_viewer_state", records[[2]]$text),
  "tool calls are summarized in the transcript text"
)
ok(
  grepl("\\[redacted\\]", records[[1]]$text),
  "transcript text is redacted"
)

## ------------------------------------------------------------- payload ----
turns <- mk_turns()
# S7 turns serialize at 27-130 KB each; give the full-payload probe a
# generous budget so it asserts STRUCTURE, not the truncation policy
payload <- agent_history_payload(turns, figures = list(
  list(id = "fig_1", spec = list(data_source = "feature"))
), max_bytes = 4L * 1024L * 1024L)
expected_roles <- vapply(turns, function(t) t@role, character(1))
ok(
  identical(payload$version, 1L) &&
    identical(length(payload$turns), length(turns)) &&
    identical(length(payload$transcript), length(turns)) &&
    identical(vapply(payload$transcript, function(r) r$role, character(1)),
              ifelse(expected_roles == "user", "user", "assistant")) &&
    identical(payload$figures[[1]]$id, "fig_1") &&
    identical(payload$truncated, FALSE),
  "payload carries slim turns, transcript, and the figure registry"
)
tool_result_turn <- Filter(
  function(t) any(vapply(t@contents, function(x)
    inherits(x, "ellmer::ContentToolResult"), logical(1))),
  payload$turns)[[1]]
result_content <- Filter(
  function(x) inherits(x, "ellmer::ContentToolResult"),
  tool_result_turn@contents)[[1]]
ok(
  !any(grepl("BASE64PREVIEWPLACEHOLDER",
             paste(utils::capture.output(str(result_content)), collapse = ""))),
  "display-only payloads (base64 previews) are stripped from saved turns"
)
ok(
  identical(result_content@value, list(dataset = list(id = "demo"))),
  "tool result VALUES survive as inert context"
)
# the structured output-seam copy (extra$data) rides the slim turn when present
structured <- mk_turns()
structured[[3]]@contents[[1]] <- ellmer::ContentToolResult(
  value = structure("{\"figure_id\":\"fig_9\"}", class = "json"),
  request = ellmer::ContentToolRequest("call1", "create_figure", list()),
  extra = list(display = list(html = "x"), data = list(figure_id = "fig_9")))
payload2 <- agent_history_payload(structured, max_bytes = 4L * 1024L * 1024L)
coded <- Filter(function(x) inherits(x, "ellmer::ContentToolResult"),
                payload2$turns[[3]]@contents)[[1]]
ok(
  inherits(coded@value, "json") &&
    identical(coded@extra$data, list(figure_id = "fig_9")),
  "json-string tool values keep their structured extra$data through the slim path"
)
ok(
  identical(agent_history_payload(list(), figures = list()), NULL) &&
    identical(agent_history_payload(list(NULL), figures = list()), NULL),
  "no turns means no payload"
)

# byte cap (todo 3.6): NEWEST turns survive, oldest drop - with DISTINCT
# texts so the assertions cannot pass vacuously (the pre-Stage-2 code kept
# turns 1..k, i.e. the OLDEST, and identical-text tests masked it)
big <- lapply(seq_len(8), function(i)
  UserTurn(contents = list(ContentText(paste(
    "analysis turn", i, paste(rep("x", 600), collapse = ""))))))
all_records <- agent_transcript_records(big)
# S7 turns serialize at ~28 KB each (class metadata per instance), so a
# 150 KB budget keeps a ~5-turn trailing suffix of the 222 KB total
capped <- agent_history_payload(big, figures = list(), max_bytes = 150000L)
kept_from <- length(big) - length(capped$transcript) + 1L
ok(
  capped$truncated && length(capped$transcript) > 1L &&
    kept_from > 1L && kept_from < length(big) &&
    identical(capped$transcript[[1]]$text, all_records[[kept_from]]$text) &&
    identical(capped$transcript[[length(capped$transcript)]]$text,
              all_records[[8]]$text),
  "byte cap keeps the NEWEST turns (first kept is a later turn, last is turn 8)"
)
ok(
  identical(capped$turns[[length(capped$turns)]]@contents[[1]]@text,
            big[[8]]@contents[[1]]@text) &&
    !identical(capped$turns[[1]]@contents[[1]]@text,
               big[[1]]@contents[[1]]@text),
  "kept turns are the trailing suffix, not the leading prefix"
)
tiny_cap <- agent_history_payload(big, figures = list(), max_bytes = 100L)
ok(
  length(tiny_cap$turns) == 2L &&
    identical(tiny_cap$transcript[[2]]$text, all_records[[8]]$text),
  "at least two turns survive even under an absurdly small cap"
)
ok(
  length(agent_history_payload(turns, figures = as.list(1:30))$figures) == 20L,
  "figure registry is capped at 20 entries"
)

## ------------------------------------------- 1.8: nested runtime refs ----
# Live streams attach the registered ToolDef to every ContentToolRequest,
# and ContentToolResult@request keeps that nested reference. The slim path
# must strip BOTH levels so the persisted payload never embeds handler
# closures (the module environment, including credentials), while keeping
# the nested request's id + name (the stubber needs them after restore).
secret_env <- new.env(parent = emptyenv())
secret_env$key <- "sk-slim-path-sentinel"
secret_handler <- function(`_intent`) list(ok = TRUE)
environment(secret_handler) <- secret_env
secret_tool <- ellmer::tool(
  secret_handler, name = "state_tool", description = "closure-heavy",
  arguments = list(`_intent` = ellmer::type_string("intent"))
)
nested_req <- ellmer::ContentToolRequest(
  id = "nested-1", name = "state_tool",
  arguments = list(`_intent` = "x"), tool = secret_tool
)
heavy_turns <- list(
  UserTurn(contents = list(ContentText("call the tool"))),
  AssistantTurn(contents = list(nested_req)),
  UserTurn(contents = list(ellmer::ContentToolResult(
    value = list(ok = TRUE), request = nested_req)))
)
heavy_payload <- agent_history_payload(heavy_turns, figures = list())
heavy_bytes <- serialize(heavy_payload$turns, connection = NULL)
ok(
  identical(length(grepRaw("sk-slim-path-sentinel", heavy_bytes, fixed = TRUE)), 0L),
  "slimmed snapshot turns do not embed handler closures (nested request path)"
)
ok(
  identical(
    heavy_payload$turns[[3]]@contents[[1]]@request@name,
    "state_tool"
  ) && identical(heavy_payload$turns[[3]]@contents[[1]]@request@id, "nested-1"),
  "nested request id + name survive slimming for the restore-time stubber"
)
ok(
  identical(heavy_payload$turns[[3]]@contents[[1]]@request@tool, NULL),
  "nested request carries no ToolDef after slimming"
)

## -------------------------------------------------- restore validation ----
ok(
  identical(agent_history_restore_payload(payload)$version, 1L),
  "valid payloads validate"
)
ok(
  ut_cmp_error(agent_history_restore_payload(NULL), "missing or malformed"),
  "NULL payloads are rejected"
)
ok(
  ut_cmp_error(agent_history_restore_payload(list(version = 2L)), "Unsupported"),
  "unknown versions are rejected"
)
bad_turns <- payload
bad_turns$turns <- c(payload$turns, list("not-a-turn"))
ok(
  ut_cmp_error(agent_history_restore_payload(bad_turns), "turns are malformed"),
  "non-Turn entries are rejected"
)
bad_transcript <- payload
bad_transcript$transcript <- payload$transcript[1:2]
ok(
  ut_cmp_error(agent_history_restore_payload(bad_transcript), "transcript does not match"),
  "transcript/turn mismatch is rejected"
)
bad_figures <- payload
bad_figures$figures <- list(list(no_id = TRUE))
ok(
  ut_cmp_error(agent_history_restore_payload(bad_figures), "figure registry is malformed"),
  "figure entries without ids are rejected"
)
oversize <- payload
oversize$turns <- c(payload$turns, lapply(seq_len(4), function(i)
  UserTurn(contents = list(ContentText(paste(rep("x", 100000), collapse = ""))))))
oversize$transcript <- c(payload$transcript, lapply(seq_len(4), function(i)
  list(role = "user", text = paste(rep("x", 100000), collapse = ""))))
ok(
  ut_cmp_error(agent_history_restore_payload(oversize), "too large"),
  "oversized payloads are rejected at the restore boundary"
)

# saveRDS round-trip: turns are plain files, exactly what .ESS does
tf <- tempfile(fileext = ".ESS")
snapshot <- list(label = "x", assistant = payload)
saveRDS(snapshot, tf)
restored <- readRDS(tf)$assistant
ok(
  identical(agent_history_restore_payload(restored)$version, 1L) &&
    inherits(restored$turns[[1]], "ellmer::UserTurn"),
  "payloads survive saveRDS/readRDS (the .ESS serialization)"
)

## --------------------------------------- module API through chat_server ----
# A real shinychat chat_server (no provider call needed): restore installs
# the transcript into the UI and the full turns as model context.
chatMod <- function(input, output, session) {
  client <- ellmer::chat_openai(
    model = "gpt-4o-mini", api_key = "test-key-not-real"
  )
  chat <- shinychat::chat_server("chat", client, history = FALSE)
  chat$clear(messages = NULL, greeting = FALSE)
  .result <<- tryCatch({
    chat$clear(
      messages = lapply(payload$transcript, function(r)
        list(role = r$role, content = r$text)),
      greeting = FALSE, client_history = "set"
    )
    chat$client$set_turns(payload$turns)
    list(
      n_turns = length(chat$client$get_turns()),
      roles = vapply(chat$client$get_turns(), function(t) t@role, character(1)),
      last_text = ellmer::contents_text(
        utils::tail(chat$client$get_turns(), 1)[[1]])
    )
  }, error = function(e) list(err = conditionMessage(e)))
  invisible(chat)
}
result <- NULL
testServer(chatMod, {})
result <- .result
ok(
  is.null(result$err) &&
    identical(result$n_turns, length(payload$turns)) &&
    setequal(result$roles, c("user", "assistant")),
  "chat restore installs the full turns as model context"
)
ok(
  identical(result$last_text, "Done. Here is the summary."),
  "restored context ends with the saved final answer"
)

## ----------------------------------- assistant module API (no provider) ----
# Without a configured provider the chat never initializes: the API must
# degrade gracefully (NULL payload, figures-only restore).
fd <- data.frame(x = 1:3, row.names = c("G1", "G2", "G3"))
captured_payload <- NULL
api_mod <- function(input, output, session) {
  api <<- omicsViewer:::ai_assistant_module(
    "assistant",
    state = function(sections = NULL) NULL,
    state_available = reactive(TRUE),
    feature_data = reactive(fd),
    sample_data = reactive(fd),
    expression_data = reactive(matrix(1:9, 3)),
    selected_features = reactive(rownames(fd)),
    selected_samples = reactive(rownames(fd)),
    apply_state = function(u) list(),
    apply_scatter_view = function(...) list()
  )
}
testServer(api_mod, {})
ok(
  is.list(api) && setequal(names(api), c("snapshot_payload", "restore_history", "has_conversation")),
  "assistant module exposes the WP11 API surface"
)
ok(
  identical(api$snapshot_payload(), NULL) && !api$has_conversation(),
  "unconfigured assistant yields no conversation payload"
)
restored_figures_only <- NULL
testServer(api_mod, {
  restored_figures_only <<- isTRUE(api$restore_history(payload))
})
ok(
  isTRUE(restored_figures_only),
  "figures-only restore succeeds without a chat"
)
