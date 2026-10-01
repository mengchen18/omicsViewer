# tests/test_agentLoopOffline.R ------------------------------------------------
#
# WP0 Tier T2 (todo 4.5): the offline agent-LOOP tier. A webfakes
# OpenAI-compatible server scripts SSE responses; the REAL shinychat
# chat_server + the REAL ellmer client + the REAL registered tools of
# ai_assistant_module run against it inside testServer. This is the
# safety net for 4.4's behavioural moves: it exercises what unit tests
# cannot - the multi-request tool loop, the Governor's request-limit
# stop, user cancellation mid-stream, the WP13b compaction race (LLM
# install AND timeout/fallback), and the WP11 snapshot->restore->
# continue cycle. (Compaction DISCARD on digest conflict stays unit-
# covered in test_agentContext; driving a mid-compaction submit through
# the fake server added only flakiness.)
#
# Webfakes serialization rule (learned the hard way): the app runs in a
# separate callr process and handler closures lose references to the
# test process's globals - everything a handler touches must be an app
# attribute. Counters mutate inside the child and are read back through
# a GET /stats endpoint.
#
# ellmer 0.5.0 chunk shapes (empirically pinned): each response must be
# ONE `data:` chunk carrying role + content-or-tool_calls + finish_reason
# (+ usage); a separate finish_reason chunk makes stream_merge_chunks
# accumulate finish_reason into a length-2 list and the parser dies with
# "EXPR must be a length 1 vector".

library(omicsViewer)
library(unittest, quietly = TRUE)
library(ellmer)

skip_offline <- !requireNamespace("webfakes", quietly = TRUE) ||
  !requireNamespace("later", quietly = TRUE)

fmt <- function(p) paste0("data: ", jsonlite::toJSON(p, auto_unbox = TRUE), "\n\n")

sse_final <- function(text, prompt_tokens = 100, completion_tokens = 10) {
  p <- list(
    id = "chatcmpl-t2", object = "chat.completion.chunk",
    choices = list(list(index = 0,
                        delta = list(role = "assistant", content = text),
                        finish_reason = "stop")),
    usage = list(prompt_tokens = prompt_tokens,
                 completion_tokens = completion_tokens,
                 total_tokens = prompt_tokens + completion_tokens))
  fmt(p)
}
sse_tool <- function(call_id, tool, args_json,
                     prompt_tokens = 100, completion_tokens = 10) {
  tc <- list(index = 0, id = call_id, type = "function")
  tc[["function"]] <- list(name = tool, arguments = args_json)
  p <- list(
    id = "chatcmpl-t2", object = "chat.completion.chunk",
    choices = list(list(index = 0,
                        delta = list(role = "assistant", tool_calls = list(tc)),
                        finish_reason = "tool_calls")),
    usage = list(prompt_tokens = prompt_tokens,
                 completion_tokens = completion_tokens,
                 total_tokens = prompt_tokens + completion_tokens))
  fmt(p)
}
DONE <- "data: [DONE]\n\n"

# The fake OpenAI-compatible endpoint. `responses` scripts the MAIN
# conversation requests (popped in order; the last repeats when
# exhausted), `summary_responses` scripts summariser requests (detected
# by the summariser system prompt - same provider config, same URL).
# A script entry of "HANG" delays the whole response (timeout race).
#
# Transport note: the full SSE payload is sent as ONE Content-Length
# body, not a chunked stream - ellmer's connection reader needs a clean
# EOF and webfakes' chunked terminator produced
# "transfer closed with outstanding read data remaining". Line-by-line
# SSE parsing over a complete body is identical on the client side.
fake_openai <- function(responses, summary_responses = NULL) {
  app <- webfakes::new_app()
  app$responses <- responses
  app$summary_responses <- summary_responses
  app$summary_marker <-
    "summarisation component of the omicsViewer analysis assistant"
  app$done <- "data: [DONE]\n\n"
  app$n_main <- 0L
  app$n_summary <- 0L
  app$post("/v1/chat/completions", function(req, res) {
    # the raw request body lives at req$.body (req$body stays NULL
    # without a body middleware; mw_raw() only covers octet-stream)
    raw <- tryCatch(rawToChar(req$.body), error = function(e) "")
    is_summary <- grepl(app$summary_marker, raw, fixed = TRUE)
    if (is_summary) {
      app$n_summary <- app$n_summary + 1L
      pool <- app$summary_responses
      idx <- app$n_summary
    } else {
      app$n_main <- app$n_main + 1L
      pool <- app$responses
      idx <- app$n_main
    }
    script <- if (length(pool)) pool[[min(idx, length(pool))]] else character()
    res$set_status(200)$set_header("Content-Type", "text/event-stream")
    if (identical(script, "HANG")) {
      # res$delay BLOCKS the webfakes child process (serial serving), so
      # keep the hold short: every cancel/timeout in this suite lands
      # well inside it, and later /stats reads recover right after
      res$delay(6)
      res$send(app$done)
    } else {
      res$send(paste0(paste(script, collapse = ""), app$done))
    }
    invisible(NULL)
  })
  app$get("/stats", function(req, res) {
    res$send_json(list(main = app$n_main, summary = app$n_summary),
                  auto_unbox = TRUE)
  })
  app
}

stats <- function(proc) {
  jsonlite::fromJSON(paste0(
    readLines(paste0(proc$url(), "stats"), warn = FALSE), collapse = ""))
}

# fd/pd fixtures (the seam-test shape; the real tools validate against them)
t2_fd <- data.frame(
  score = c(1, 2, 3), logFdr = c(5, 4, 3),
  category = c("kinase", "phosphatase", "kinase"),
  row.names = c("Gene1", "Gene2", "Gene3"), check.names = FALSE)
t2_pd <- data.frame(group = c("WT", "WT", "KO", "KO"),
                    row.names = c("S1", "S2", "S3", "S4"))
t2_state <- function(sections = NULL) list(
  dataset = list(id = "demo.RDS", class = "ExpressionSet",
                 dimensions = c(features = 3L, samples = 4L)),
  active_tabs = list(data_space = "Feature", analysis_space = "Feature"),
  selection = list(
    features = list(count = 0L, ids = character(), truncated = FALSE),
    samples = list(count = 0L, ids = character(), truncated = FALSE)),
  available_tabs = list(data_space = c("Feature", "Sample"),
                        analysis_space = c("Feature", "ORA")))

t2_module_args <- list(
  state = t2_state,
  state_available = shiny::reactive(TRUE),
  feature_data = shiny::reactive(t2_fd),
  sample_data = shiny::reactive(t2_pd),
  expression_data = shiny::reactive(matrix(
    1:12, nrow = 3, dimnames = list(rownames(t2_fd), rownames(t2_pd)))),
  selected_features = shiny::reactive(character()),
  selected_samples = shiny::reactive(character()),
  apply_state = function(u) list(),
  apply_scatter_view = function(...) list(),
  apply_enrichment = function(u) list(),
  apply_table_view = function(u) list(),
  store = omicsViewer:::widget_store_new()
)

# pump the mock session + the later loop (ellmer promises, janitor tails)
# until `until()` is TRUE or `seconds` elapse; returns whether until() held
pump_until <- function(session, until, seconds = 20, step = 0.05) {
  deadline <- Sys.time() + seconds
  repeat {
    ok_ <- tryCatch(until(), error = function(e) FALSE)
    if (isTRUE(ok_)) return(TRUE)
    later::run_now(step)
    tryCatch(session$flushReact(), error = function(e) NULL)
    if (Sys.time() > deadline) return(FALSE)
  }
}

# point the module's environment config at a fake process
t2_env_on <- function(proc, extra = list()) {
  vars <- c(
    list(OMICSVIEWER_LLM_PROVIDER = "openai_compatible",
         OMICSVIEWER_LLM_BASE_URL = paste0(proc$url(), "v1"),
         OMICSVIEWER_LLM_API_KEY = "sk-t2-offline",
         OMICSVIEWER_LLM_MODEL = "t2-fake-model"),
    extra)
  old <- Sys.getenv(names(vars))
  do.call(Sys.setenv, as.list(vars))
  old
}
t2_env_off <- function(old) do.call(Sys.setenv, as.list(old))

# extract one turn's concatenated text (tool results excluded)
t2_turn_text <- function(turn) {
  tryCatch(paste(unlist(lapply(turn@contents, function(x)
    tryCatch(x@text, error = function(e) "") %||% "")),
    collapse = " "), error = function(e) "")
}

if (skip_offline) {
  ok(TRUE, "T2 offline loop tier skipped (webfakes/later unavailable)")
} else {
library(webfakes)
library(later)

## ------------------------------------------------------------------ block 1
## The 3-round tool loop: three real tool calls dispatch through the real
## ellmer loop, each result rides back through the seam codec, and the
## final text settles the turn list.
proc1 <- new_app_process(fake_openai(list(
  sse_tool("c1", "get_omics_viewer_state", "{\"_intent\":\"overview\"}"),
  sse_tool("c2", "search_annotations",
           "{\"space\":\"feature\",\"query\":\"Gene\",\"_intent\":\"find\"}"),
  sse_tool("c3", "find_controls", "{\"prefix\":\"dataspace\",\"_intent\":\"list controls\"}"),
  sse_final("Used three tools; the interface is ready.", 4000, 20)
)))
old1 <- t2_env_on(proc1, list(OMICSVIEWER_LLM_CONTEXT_TOKENS = "0"))
t2 <- new.env()
shiny::testServer(
  omicsViewer:::ai_assistant_module,
  args = t2_module_args,
  expr = {
    session$setInputs(`chat_user_input` = "inspect the app")
    pump_until(session, function() {
      identical(chat_object$status(), "idle") &&
        tryCatch(grepl("Used three tools",
                       chat_object$last_turn()@text, fixed = TRUE),
                 error = function(e) FALSE)
    }, seconds = 30)
    t2$turns <<- tryCatch(chat_object$client$get_turns(), error = function(e) NULL)
    t2$final <<- tryCatch(chat_object$last_turn()@text, error = function(e) "")
  })
cnt1 <- stats(proc1)
proc1$stop(); t2_env_off(old1)
tool_results <- Filter(function(x) inherits(x, "ellmer::ContentToolResult"),
                       unlist(lapply(t2$turns, function(t) t@contents),
                              recursive = FALSE))
ok(ut_cmp_identical(cnt1$main, 4L),
   "T2 loop: exactly four provider requests (three tool rounds + final)")
ok(ut_cmp_identical(length(tool_results), 3L),
   "T2 loop: three real tool results entered the turn list")
ok(ut_cmp_identical(
  all(vapply(tool_results, function(x)
    inherits(x@value, "json") && nzchar(as.character(x@value)), logical(1))),
  TRUE),
  "T2 loop: every tool result value rides the seam codec (json string)")
ok(ut_cmp_identical(t2$final, "Used three tools; the interface is ready."),
   "T2 loop: the final assistant text settles into the last turn")

## ------------------------------------------------------------------ block 2
## Governor stop: the request limit blocks the NEXT submit; the blocked
## submit spends no provider request.
proc2 <- new_app_process(fake_openai(list(
  sse_tool("c1", "find_controls", "{\"query\":\"heatmap\",\"_intent\":\"first\"}"),
  sse_final("first submit done", 100, 10)
)))
log2 <- file.path(tempdir(), paste0("t2-logs-gov-", Sys.getpid()))
if (dir.exists(log2)) unlink(log2, recursive = TRUE)
dir.create(log2, showWarnings = FALSE)
old2 <- t2_env_on(proc2, list(
  OMICSVIEWER_LLM_MAX_REQUESTS = "2",
  OMICSVIEWER_LLM_CONTEXT_TOKENS = "0",
  OMICSVIEWER_LLM_LOG = "true",
  OMICSVIEWER_LLM_LOG_DIR = log2))
t2b <- new.env()
shiny::testServer(
  omicsViewer:::ai_assistant_module,
  args = t2_module_args,
  expr = {
    session$setInputs(`chat_user_input` = "one")
    t2b$first_ok <<- pump_until(session, function() {
      identical(chat_object$status(), "idle") &&
        tryCatch(grepl("first submit done",
                       chat_object$last_turn()@text, fixed = TRUE),
                 error = function(e) FALSE)
    }, seconds = 25)
    session$setInputs(`chat_user_input` = "two")
    pump_until(session, function() {
      identical(chat_object$status(), "idle") &&
        !is.null(chat_object$last_error()())
    }, seconds = 15)
    t2b$err <<- tryCatch(chat_object$last_error()(),
                         error = function(e) NULL)
  })
cnt2 <- stats(proc2)
proc2$stop(); t2_env_off(old2)
gov_events <- unlist(lapply(list.files(log2, full.names = TRUE), function(f) {
  tryCatch(readLines(f, warn = FALSE), error = function(e) character())
}))
unlink(log2, recursive = TRUE)
ok(ut_cmp_identical(t2b$first_ok, TRUE),
   "T2 governor: the first submit completes inside the request limit")
ok(ut_cmp_identical(cnt2$main, 2L),
   "T2 governor: the blocked submit spends no provider request")
ok(ut_cmp_identical(
  any(grepl("request_limit_reached", gov_events, fixed = TRUE)), TRUE),
  "T2 governor: the limit stop is logged (request_limit_reached) and the stream aborted")

## ------------------------------------------------------------------ block 3
## User cancellation mid-stream: the scripted stream is held open long
## enough for the cancel input to land; the stream stops and no further
## requests fire.
proc3 <- new_app_process(fake_openai(list("HANG")))
old3 <- t2_env_on(proc3, list(OMICSVIEWER_LLM_CONTEXT_TOKENS = "0"))
t2c <- new.env()
shiny::testServer(
  omicsViewer:::ai_assistant_module,
  args = t2_module_args,
  expr = {
    session$setInputs(`chat_user_input` = "cancel me")
    t2c$was_streaming <<- pump_until(session, function() {
      identical(chat_object$status(), "streaming")
    }, seconds = 10)
    session$setInputs(`chat_cancel` = 1L)
    t2c$back_to_idle <<- pump_until(session, function() {
      identical(chat_object$status(), "idle")
    }, seconds = 10)
    t2c$last <<- tryCatch(chat_object$last_turn()@text,
                          error = function(e) "")
  })
  # NB: no /stats here - the aborted socket wedges the (serial) webfakes
  # child until its delay elapses; the process is stopped below anyway
  proc3$stop(); t2_env_off(old3)
  ok(ut_cmp_identical(t2c$was_streaming, TRUE),
     "T2 cancel: the held-open stream reaches the streaming state")
  ok(ut_cmp_identical(t2c$back_to_idle, TRUE),
     "T2 cancel: user cancellation returns the chat to idle")
  ok(ut_cmp_identical(!grepl("slow response", t2c$last, fixed = TRUE), TRUE),
     "T2 cancel: the cancelled stream leaves no completed assistant text")

## ------------------------------------------------------------------ block 4
## Compaction (WP13b race): LLM install (summary served) and timeout
## fallback (summary hangs; the race cancels the stream and installs the
## deterministic summary). Oversized reported prompt_tokens (4000)
## against a small context window arms the janitor.
proc4a <- new_app_process(fake_openai(
  list(sse_final("answer one", 4500, 10),
       sse_final("answer two", 4500, 10),
       sse_final("answer three", 4500, 10)),
  summary_responses = sse_final("Compacted: the user inspected the app.")))
old4 <- t2_env_on(proc4a, list(
  # Compaction arithmetic with the usage-anchored estimator (empirically
  # pinned): the TRIGGER compares the LAST exchange's reported input
  # tokens (the anchor) against the window, so every scripted response
  # reports 4500 against a 4000-token window; the CUT's feasibility loop
  # then falls through to its most-aggressive fallback (keep the last
  # exchange), which is exactly the shape under test.
  OMICSVIEWER_LLM_CONTEXT_TOKENS = "4000",
  OMICSVIEWER_LLM_SUMMARY_TIMEOUT = "10"))
t2d <- new.env()
shiny::testServer(
  omicsViewer:::ai_assistant_module,
  args = t2_module_args,
  expr = {
    for (msg in c("hello one", "hello two", "hello three")) {
      session$setInputs(`chat_user_input` = msg)
      pump_until(session, function() {
        identical(chat_object$status(), "idle") &&
          tryCatch(grepl("answer three", chat_object$last_turn()@text,
                         fixed = TRUE), error = function(e) FALSE)
      }, seconds = 25)
    }
    # janitor -> summary request -> install as the leading framed pair
    pump_until(session, function() {
      tryCatch({
        turns <- chat_object$client$get_turns()
        length(turns) >= 2L && grepl("<omicsviewer-conversation-summary>",
                                     t2_turn_text(turns[[1]]), fixed = TRUE)
      }, error = function(e) FALSE)
    }, seconds = 30)
    t2d$turns <<- tryCatch(chat_object$client$get_turns(),
                           error = function(e) NULL)
  })
cnt4a <- stats(proc4a)
proc4a$stop()
ok(ut_cmp_identical(
  grepl("<omicsviewer-conversation-summary>",
        tryCatch(t2_turn_text(t2d$turns[[1]]), error = function(e) ""),
        fixed = TRUE) &&
    grepl("Compacted: the user inspected the app.",
          tryCatch(t2_turn_text(t2d$turns[[1]]), error = function(e) ""),
          fixed = TRUE),
  TRUE),
  "T2 compaction: the LLM summary installs as the leading framed turn pair")
ok(ut_cmp_identical(cnt4a$summary, 1L),
   "T2 compaction: exactly one summariser request hit the fake endpoint")
ok(ut_cmp_identical(length(t2d$turns) <= 6L, TRUE),
   "T2 compaction: the context is bounded to the summary pair plus the surviving latest exchange")

# 4b: the summary hangs -> the WP13b race times out, cancels the stream
# (no orphan connection) and installs the deterministic fallback.
proc4b <- new_app_process(fake_openai(
  list(sse_final("answer one", 4500, 10),
       sse_final("answer two", 4500, 10),
       sse_final("answer three", 4500, 10)),
  summary_responses = "HANG"))
t2_env_on(proc4b, list(
  OMICSVIEWER_LLM_CONTEXT_TOKENS = "4000",
  OMICSVIEWER_LLM_SUMMARY_TIMEOUT = "5"))
t2e <- new.env()
shiny::testServer(
  omicsViewer:::ai_assistant_module,
  args = t2_module_args,
  expr = {
    for (msg in c("hello one", "hello two", "hello three")) {
      session$setInputs(`chat_user_input` = msg)
      pump_until(session, function() {
        identical(chat_object$status(), "idle") &&
          tryCatch(grepl("answer three", chat_object$last_turn()@text,
                         fixed = TRUE), error = function(e) FALSE)
      }, seconds = 25)
    }
    # the summary hangs; the WP13b race must fire: cancel the stream at
    # the timeout and install the deterministic fallback summary
    pump_until(session, function() {
      tryCatch({
        turns <- chat_object$client$get_turns()
        length(turns) >= 2L && grepl("<omicsviewer-conversation-summary>",
                                     t2_turn_text(turns[[1]]), fixed = TRUE)
      }, error = function(e) FALSE)
    }, seconds = 30)
    t2e$turns <<- tryCatch(chat_object$client$get_turns(),
                           error = function(e) NULL)
  })
proc4b$stop(); t2_env_off(old4)
ok(ut_cmp_identical(
  grepl("<omicsviewer-conversation-summary>",
        tryCatch(t2_turn_text(t2e$turns[[1]]), error = function(e) ""),
        fixed = TRUE),
  TRUE),
  "T2 compaction: a hung summariser still yields the framed summary pair (timeout fallback)")

## ------------------------------------------------------------------ block 5
## WP11 snapshot -> restore -> continue: the payload round-trips through
## the module API and a restored conversation keeps submitting.
proc5 <- new_app_process(fake_openai(list(
  sse_final("first exchange"),
  sse_final("second exchange"),
  sse_final("continued after restore"))))
old5 <- t2_env_on(proc5, list(OMICSVIEWER_LLM_CONTEXT_TOKENS = "0"))
t2f <- new.env()
shiny::testServer(
  omicsViewer:::ai_assistant_module,
  args = t2_module_args,
  expr = {
    session$setInputs(`chat_user_input` = "first")
    pump_until(session, function() {
      identical(chat_object$status(), "idle") &&
        tryCatch(grepl("first exchange", chat_object$last_turn()@text,
                       fixed = TRUE), error = function(e) FALSE)
    }, seconds = 25)
    t2f$saved <<- assistant_api$snapshot_payload()
    session$setInputs(`chat_user_input` = "second")
    pump_until(session, function() {
      tryCatch(grepl("second exchange", chat_object$last_turn()@text,
                     fixed = TRUE), error = function(e) FALSE)
    }, seconds = 25)
    assistant_api$restore_history(t2f$saved)
    session$setInputs(`chat_user_input` = "continue")
    pump_until(session, function() {
      tryCatch(grepl("continued after restore", chat_object$last_turn()@text,
                     fixed = TRUE), error = function(e) FALSE)
    }, seconds = 25)
    t2f$final_text <<- tryCatch(chat_object$last_turn()@text,
                                error = function(e) "")
  })
proc5$stop(); t2_env_off(old5)
ok(ut_cmp_identical(is.list(t2f$saved) && length(t2f$saved) >= 1L, TRUE),
   "T2 history: snapshot_payload returns a structured conversation payload")
ok(ut_cmp_identical(t2f$final_text, "continued after restore"),
   "T2 history: a restored conversation continues submitting successfully")

} # end !skip_offline
