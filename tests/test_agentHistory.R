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
  list(
    UserTurn(contents = list(ContentText("Please analyze sk-abcdef0123456789abcdef"))),
    AssistantTurn(contents = list(
      ContentText("Calling a tool"),
      ContentToolRequest("call1", "get_omics_viewer_state", list(`_intent` = "x")),
      ContentToolResult(
        value = list(dataset = list(id = "demo")),
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
            c("user", "assistant", "assistant")),
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
payload <- agent_history_payload(turns, figures = list(
  list(id = "fig_1", spec = list(data_source = "feature"))
))
ok(
  identical(payload$version, 1L) &&
    identical(length(payload$turns), 3L) &&
    identical(length(payload$transcript), 3L) &&
    identical(payload$figures[[1]]$id, "fig_1") &&
    identical(payload$truncated, FALSE),
  "payload carries slim turns, transcript, and the figure registry"
)
ok(
  !any(grepl("BASE64PREVIEWPLACEHOLDER",
             vapply(payload$turns[[2]]@contents, function(x)
               paste(utils::capture.output(str(x)), collapse = ""), character(1)))),
  "display-only payloads (base64 previews) are stripped from saved turns"
)
ok(
  identical(payload$turns[[2]]@contents[[3]]@value, list(dataset = list(id = "demo"))),
  "tool result VALUES survive as inert context"
)
ok(
  identical(agent_history_payload(list(), figures = list()), NULL) &&
    identical(agent_history_payload(list(NULL), figures = list()), NULL),
  "no turns means no payload"
)

# byte cap: oldest turns drop first, the last two always survive
big <- lapply(seq_len(8), function(i)
  UserTurn(contents = list(ContentText(paste(rep("analysis ", 600), collapse = "")))))
capped <- agent_history_payload(big, figures = list(), max_bytes = 6000L)
ok(
  capped$truncated && length(capped$turns) < 8L &&
    length(capped$turns) >= 2L &&
    identical(capped$transcript[[length(capped$transcript)]]$text,
              agent_transcript_records(big)[[8]]$text),
  "byte cap drops oldest turns and keeps the newest"
)
tiny_cap <- agent_history_payload(big, figures = list(), max_bytes = 100L)
ok(
  length(tiny_cap$turns) == 2L,
  "at least two turns survive even under an absurdly small cap"
)
ok(
  length(agent_history_payload(turns, figures = as.list(1:30))$figures) == 20L,
  "figure registry is capped at 20 entries"
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
  is.null(result$err) && identical(result$n_turns, 3L) &&
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
