library(omicsViewer)
library(unittest, quietly = TRUE)
library(ellmer)

agent_context_policy <- omicsViewer:::agent_context_policy
agent_stub_history <- omicsViewer:::agent_stub_history
agent_estimate_context_tokens <- omicsViewer:::agent_estimate_context_tokens
agent_compaction_cut <- omicsViewer:::agent_compaction_cut
agent_last_human_turn_index <- omicsViewer:::agent_last_human_turn_index
agent_system_prompt_add_block <- omicsViewer:::agent_system_prompt_add_block
agent_system_prompt_strip_block <- omicsViewer:::agent_system_prompt_strip_block
agent_fallback_summary <- omicsViewer:::agent_fallback_summary
agent_compaction_prompt <- omicsViewer:::agent_compaction_prompt
agent_context_archive_merge <- omicsViewer:::agent_context_archive_merge
agent_context_digest <- omicsViewer:::agent_context_digest
agent_history_payload <- omicsViewer:::agent_history_payload
.agent_context_new_archive <- omicsViewer:::`.agent_context_new_archive`

## -- turn constructors (no provider involved) --------------------------

mk_user <- function(text)
  ellmer::UserTurn(contents = list(ellmer::ContentText(text)))

mk_assistant <- function(text = NULL, tokens = c(100L, 20L, 0L),
                         thinking = NULL, request = NULL) {
  contents <- list()
  if (!is.null(thinking))
    contents <- c(contents, list(ellmer::ContentThinking(thinking = thinking)))
  if (!is.null(request))
    contents <- c(contents, list(request))
  if (!is.null(text))
    contents <- c(contents, list(ellmer::ContentText(text)))
  ellmer::AssistantTurn(contents = contents, tokens = tokens)
}

mk_tool_turn <- function(request, value, error = NULL)
  ellmer::UserTurn(contents = list(
    ellmer::ContentToolResult(value = value, request = request, error = error)
  ))

mk_request <- function(id, name, arguments = list(a = 1))
  ellmer::ContentToolRequest(id = id, name = name, arguments = arguments)

big_state <- function(seed = 1) {
  set.seed(seed)
  list(
    dataset = "demo",
    annotations = as.list(paste0("column_", seq_len(400), "_", seed)),
    scatter_view = list(x = "ttest|A_vs_B|log.fdr", y = NULL)
  )
}

## -- policy --------------------------------------------------------------

withr_local_env <- function(values, code) {
  old <- Sys.getenv(names(values), names = TRUE, unset = NA_character_)
  on.exit(
    if (any(is.na(old)))
      Sys.unsetenv(names(values)[is.na(old)])
    else
      do.call(Sys.setenv, as.list(old[!is.na(old)])),
    add = TRUE
  )
  do.call(Sys.setenv, as.list(values))
  force(code)
}

ok(ut_cmp_identical(agent_context_policy()$tokens, 24000L), "default context token limit")
ok(ut_cmp_identical(agent_context_policy()$tool_result_bytes, 16384L), "default tool-result byte cap")
ok(
  withr_local_env(list(OMICSVIEWER_LLM_CONTEXT_TOKENS = "8000"),
    ut_cmp_identical(agent_context_policy()$tokens, 8000L)),
  "context token limit env override"
)
ok(
  withr_local_env(list(OMICSVIEWER_LLM_CONTEXT_TOKENS = "-5"),
    ut_cmp_identical(agent_context_policy()$tokens, 24000L)),
  "invalid context token limit falls back to default"
)
ok(
  withr_local_env(list(OMICSVIEWER_LLM_COMPACT_TO = "1.5"),
    ut_cmp_identical(agent_context_policy()$compact_to, 0.5)),
  "out-of-range compact_to falls back to default"
)

## -- estimator -----------------------------------------------------------

empty_convo <- list(mk_user("hi"), mk_assistant("hello", tokens = c(1000L, 10L, 0L)))
ok(
  agent_estimate_context_tokens(empty_convo) == 1010,
  "usage anchor alone estimates input+output"
)
grown <- c(empty_convo, list(mk_user(paste(rep("x", 400), collapse = ""))))
est <- agent_estimate_context_tokens(grown)
ok(est > 1010 && est < 1010 + 200, "appended content adds to the anchored estimate")
no_anchor <- list(mk_user("hello"), mk_user("again"))
ok(
  agent_estimate_context_tokens(no_anchor) > 1600,
  "no-anchor estimate includes the fixed tool-schema overhead"
)
ok(ut_cmp_identical(agent_estimate_context_tokens(list()), 0), "empty turn list estimates zero")

## -- stubbing: snapshot supersession -------------------------------------

req1 <- mk_request("call_1", "get_omics_viewer_state")
req2 <- mk_request("call_2", "get_omics_viewer_state")
turns <- list(
  mk_user("what is in the dataset?"),
  mk_assistant(request = req1),
  mk_tool_turn(req1, big_state(1)),
  mk_assistant("here is the overview"),
  mk_user("now show me the panels"),
  mk_assistant(request = req2),
  mk_tool_turn(req2, big_state(2)),
  mk_assistant("here are the panels"),
  mk_user("thanks!")
)

archive <- .agent_context_new_archive()
policy <- agent_context_policy()
st <- agent_stub_history(turns, policy, archive)

ok(isTRUE(st$changed), "stale state snapshot is stubbed")
ok(ut_cmp_identical(st$stub_count, 1L), "exactly one superseded snapshot stubbed")
stub_turn <- st$turns[[3]]
stub_value <- stub_turn@contents[[1]]@value
ok(grepl("get_omics_viewer_state", stub_value, fixed = TRUE), "stub names the tool to re-call")
ok(grepl("ctx_1", stub_value, fixed = TRUE), "stub embeds the archive id")
ok(ut_cmp_identical(stub_turn@contents[[1]]@request@id, "call_1"), "stub preserves the tool-call id")
ok(ut_cmp_identical(archive$values$ctx_1, big_state(1)), "archive keeps the original value")
latest <- st$turns[[7]]@contents[[1]]@value
ok(ut_cmp_identical(latest, big_state(2)), "latest snapshot survives untouched")
ok(ut_cmp_identical(st$turns[[9]], turns[[9]]), "human turns are untouched")

## -- stubbing: protected region ------------------------------------------

archive2 <- .agent_context_new_archive()
fresh <- list(
  mk_user("q1"),
  mk_assistant(request = req1),
  mk_tool_turn(req1, big_state(1)),
  mk_assistant("done"),
  mk_user("q2"),
  mk_assistant(request = req2),
  mk_tool_turn(req2, big_state(2))
)
st2 <- agent_stub_history(fresh, policy, archive2)
ok(isTRUE(st2$changed), "snapshot superseded by a newer fetch is stubbed even when the newer fetch is in the protected exchange")
ok(
  grepl("ctx_1", st2$turns[[3]]@contents[[1]]@value, fixed = TRUE),
  "older snapshot outside the boundary is stubbed"
)
ok(ut_cmp_identical(st2$turns[[7]]@contents[[1]]@value, big_state(2)), "the protected (latest) result itself survives")
only_latest <- list(
  mk_user("q1"),
  mk_assistant(request = req1),
  mk_tool_turn(req1, big_state(1)),
  mk_assistant("done"),
  mk_user("q2"),
  mk_assistant("answer")
)
st3 <- agent_stub_history(only_latest, policy, .agent_context_new_archive())
ok(!isTRUE(st3$changed), "the only state snapshot is never stubbed")

## -- stubbing: figure revision supersession ------------------------------

c1 <- mk_request("call_f1", "create_figure")
u1 <- mk_request("call_f2", "update_figure")
fig_turns <- list(
  mk_user("make a volcano"),
  mk_assistant(request = c1),
  mk_tool_turn(c1, list(figure_id = "fig_1", parent_figure_id = NULL, row_count = 10)),
  mk_assistant("created"),
  mk_user("change the cutoffs"),
  mk_assistant(request = u1),
  mk_tool_turn(u1, list(figure_id = "fig_2", parent_figure_id = "fig_1", row_count = 12)),
  mk_assistant("updated"),
  mk_user("nice")
)
stf <- agent_stub_history(fig_turns, policy, .agent_context_new_archive())
ok(isTRUE(stf$changed), "direct figure revision is stubbed")
ok(
  grepl("fig_1", stf$turns[[3]]@contents[[1]]@value, fixed = TRUE),
  "figure stub names the superseded figure"
)
ok(ut_cmp_identical(stf$turns[[7]]@contents[[1]]@value$figure_id, "fig_2"), "latest revision survives")

indep <- list(
  mk_user("two figures"),
  mk_assistant(request = c1),
  mk_tool_turn(c1, list(figure_id = "fig_1", parent_figure_id = NULL)),
  mk_assistant("one"),
  mk_assistant(request = mk_request("call_f3", "create_figure")),
  mk_tool_turn(mk_request("call_f3", "create_figure"),
               list(figure_id = "fig_3", parent_figure_id = NULL)),
  mk_assistant("two"),
  mk_user("thanks")
)
sti <- agent_stub_history(indep, policy, .agent_context_new_archive())
ok(!isTRUE(sti$changed), "independent figures are not superseded")

## -- stubbing: generic byte cap -------------------------------------------

cap_req <- mk_request("call_big", "search_annotations")
huge <- paste(rep("a", 60000), collapse = "")
capped <- list(
  mk_user("search"),
  mk_assistant(request = cap_req),
  mk_tool_turn(cap_req, huge),
  mk_assistant("found"),
  mk_user("thanks")
)
stc <- agent_stub_history(capped, policy, .agent_context_new_archive())
ok(isTRUE(stc$changed), "oversized tool result is truncated")
ok(
  grepl("Oversized search_annotations", stc$turns[[3]]@contents[[1]]@value, fixed = TRUE),
  "oversize stub names the tool and the reason"
)
ok(
  nchar(stc$turns[[3]]@contents[[1]]@value) < 500,
  "oversize stub is small"
)

## -- stubbing: thinking strip ---------------------------------------------

think_req <- mk_request("call_t1", "search_annotations")
tt <- list(
  mk_user("think hard"),
  mk_assistant(thinking = "long internal reasoning", text = "answer",
               request = think_req),
  mk_tool_turn(think_req, "small result"),
  mk_assistant("done"),
  mk_user("next")
)
stt <- agent_stub_history(tt, policy, .agent_context_new_archive())
ok(
  !any(vapply(stt$turns[[2]]@contents, function(x)
    inherits(x, "ellmer::ContentThinking"), logical(1))),
  "old thinking is stripped from assistant turns"
)
ok(
  any(vapply(stt$turns[[2]]@contents, function(x)
    inherits(x, "ellmer::ContentText"), logical(1))),
  "old assistant text survives the thinking strip"
)
recent_think <- c(tt, list(mk_assistant("final", thinking = "fresh reasoning")))
stt2 <- agent_stub_history(recent_think, policy, .agent_context_new_archive())
ok(
  any(vapply(stt2$turns[[length(recent_think)]]@contents, function(x)
    inherits(x, "ellmer::ContentThinking"), logical(1))),
  "thinking inside the latest exchange survives"
)

## -- stubbing: errors are never stubbed -----------------------------------

err_req <- mk_request("call_e1", "get_omics_viewer_state")
err_req2 <- mk_request("call_e2", "get_omics_viewer_state")
et <- list(
  mk_user("q"),
  mk_assistant(request = err_req),
  mk_tool_turn(err_req, big_state(1), error = simpleError("boom")),
  mk_assistant("it failed"),
  mk_user("q2"),
  mk_assistant(request = err_req2),
  mk_tool_turn(err_req2, big_state(2)),
  mk_assistant("ok"),
  mk_user("q3")
)
ste <- agent_stub_history(et, policy, .agent_context_new_archive())
ok(!is.null(ste$turns[[3]]@contents[[1]]@error), "errored tool results are never stubbed")
ok(
  ut_cmp_identical(ste$turns[[3]]@contents[[1]]@value, big_state(1)),
  "errored tool result value survives"
)

## -- compaction cut --------------------------------------------------------

exchanges <- list()
for (i in 1:4) {
  exchanges <- c(exchanges, list(
    mk_user(paste("question", i, paste(rep("q", 400), collapse = ""))),
    mk_assistant(paste("answer", i), tokens = c(2000L * i, 100L, 0L))
  ))
}
ok(ut_cmp_identical(agent_compaction_cut(exchanges, target_tokens = 1e9), 1L),
   "generous target keeps everything from the first human turn")
cut <- agent_compaction_cut(exchanges, target_tokens = 100)
ok(cut >= 5, "tight target cuts at an exchange boundary, keeping the last exchange whole")
starts <- which(vapply(exchanges, function(t)
  inherits(t, "ellmer::UserTurn") &&
    !any(vapply(t@contents, function(x)
      inherits(x, "ellmer::ContentToolResult"), logical(1))), logical(1)))
ok(cut %in% starts, "cut point is always a human-turn index")
single <- list(mk_user("only"), mk_assistant("exchange"))
ok(is.null(agent_compaction_cut(single, target_tokens = 1)), "nothing to cut returns NULL")

## -- prompt block round trip ----------------------------------------------

base_prompt <- "You are the omicsViewer analysis assistant."
blocked <- agent_system_prompt_add_block(base_prompt, "SUMMARY TEXT")
ok(grepl("SUMMARY TEXT", blocked, fixed = TRUE), "block carries the summary")
ok(ut_cmp_identical(agent_system_prompt_strip_block(blocked), base_prompt),
   "strip recovers the original prompt")
ok(ut_cmp_identical(
  agent_system_prompt_strip_block(agent_system_prompt_add_block(blocked, "OTHER")),
  base_prompt), "add/strip is idempotent across repeated compactions")
ok(ut_cmp_identical(agent_system_prompt_strip_block(base_prompt), base_prompt),
   "strip on a blockless prompt is a no-op")

## -- summaries -------------------------------------------------------------

fb <- agent_fallback_summary(exchanges)
ok(nchar(fb) > 0 && nchar(fb) <= 4200, "fallback summary is non-empty and bounded")
ok(grepl("question 1", fb, fixed = TRUE), "fallback summary reflects the conversation")
cp <- agent_compaction_prompt(exchanges)
ok(grepl("question 1", cp, fixed = TRUE) && grepl("300 words", cp, fixed = TRUE),
   "LLM compaction prompt carries transcript and word budget")

## -- archive merge (snapshot fidelity) -------------------------------------

archive3 <- .agent_context_new_archive()
stubbed_convo <- agent_stub_history(turns, policy, archive3)
merged <- agent_context_archive_merge(stubbed_convo$turns, archive3)
ok(ut_cmp_identical(merged[[3]]@contents[[1]]@value, big_state(1)),
   "archive merge restores the original value for persistence")
ok(ut_cmp_identical(merged[[3]]@contents[[1]]@request@id, "call_1"),
   "archive merge preserves the tool-call id")
no_archive <- agent_context_archive_merge(merged, .agent_context_new_archive())
ok(ut_cmp_identical(no_archive[[3]]@contents[[1]]@value, big_state(1)),
   "merge without archive entries is a no-op")

## -- digest (conflict guard) ----------------------------------------------

ok(
  agent_context_digest(turns) != agent_context_digest(stubbed_convo$turns),
  "digest distinguishes stubbed from original turns"
)
ok(ut_cmp_identical(agent_context_digest(turns), agent_context_digest(turns)),
   "digest is stable for identical turn lists")

## -- human-turn boundary helper --------------------------------------------

ok(ut_cmp_identical(agent_last_human_turn_index(turns), 9L), "last human turn index")
ok(ut_cmp_identical(agent_last_human_turn_index(list()), 0L), "empty list has no human turn")

## -- module integration: restore -> janitor stubs -> snapshot fidelity ------
# A real (unconfigured) assistant module: restoring a payload that carries
# two state snapshots must (a) run the janitor (history_stubbed logged with
# logging enabled) and (b) still snapshot full-fidelity values via the
# in-session archive merge.
if (requireNamespace("shiny", quietly = TRUE) &&
    requireNamespace("ellmer", quietly = TRUE) &&
    requireNamespace("shinychat", quietly = TRUE)) {
  library(shiny)
  log_dir <- file.path(tempdir(), "agent-context-it")
  restored <- NULL
  out <- NULL
  withr_local_env(
    list(
      OMICSVIEWER_LLM_LOG = "true",
      OMICSVIEWER_LLM_LOG_DIR = log_dir
    ),
    {
      int_fd <- data.frame(x = 1:3, row.names = c("G1", "G2", "G3"))
      int_api <- NULL
      int_mod <- function(input, output, session) {
        int_api <<- omicsViewer:::ai_assistant_module(
          "assistant",
          state = function(sections = NULL) NULL,
          state_available = reactive(TRUE),
          feature_data = reactive(int_fd),
          sample_data = reactive(int_fd),
          expression_data = reactive(matrix(1:9, 3)),
          selected_features = reactive(rownames(int_fd)),
          selected_samples = reactive(rownames(int_fd)),
          apply_state = function(u) list(),
          apply_scatter_view = function(...) list()
        )
      }
      testServer(int_mod, {
        # Minimal 7-turn fixture that fits the 512 KB restore budget
        # (ellmer 0.5.0 S7 turns serialize at ~27-130 KB EACH, so the
        # default 256 KB payload cap would truncate a longer fixture).
        # Boundary is turn 5 (q2): the first snapshot (turn 3) is stale,
        # the second (turn 7) sits inside the protected exchange.
        small_turns <- list(
          mk_user("what is in the dataset?"),
          mk_assistant(request = req1),
          mk_tool_turn(req1, list(dataset = "demo", seed = 1)),
          mk_assistant("here is the overview"),
          mk_user("now show me the panels"),
          mk_assistant(request = req2),
          mk_tool_turn(req2, list(dataset = "demo", seed = 2))
        )
        int_payload <- list(
          version = 1L,
          saved_at = "2026-09-26T00:00:00Z",
          turns = small_turns,
          transcript = omicsViewer:::agent_transcript_records(small_turns),
          figures = list(),
          truncated = FALSE
        )
        restored <<- tryCatch(
          int_api$restore_history(int_payload),
          error = function(e) conditionMessage(e)
        )
        out <<- tryCatch(int_api$snapshot_payload(), error = function(e) NULL)
      })
      ok(isTRUE(restored) && isTRUE(int_api$has_conversation()),
         "context fixture restores into the assistant module")
      ok(!is.null(out), "snapshot payload exists after restore")
      ok(
        !is.null(out) && ut_cmp_identical(
          out$turns[[3]]@contents[[1]]@value,
          list(dataset = "demo", seed = 1)
        ),
        "snapshot keeps full-fidelity values although the live context was slimmed"
      )
      log_files <- list.files(log_dir, pattern = "[.]jsonl$", full.names = TRUE)
      stubbed_event <- FALSE
      for (f in log_files) {
        lines <- tryCatch(readLines(f, warn = FALSE), error = function(e) character())
        if (any(grepl("\"history_stubbed\"", lines, fixed = TRUE)))
          stubbed_event <- TRUE
      }
      ok(isTRUE(stubbed_event),
         "restore path runs the context janitor (history_stubbed logged)")
    }
  )
}


## -- WP13b: summary timeout policy + digest parts + promise race ------------

agent_context_digest_parts <- omicsViewer:::agent_context_digest_parts
agent_promise_race <- omicsViewer:::agent_promise_race

ok(
  ut_cmp_identical(agent_context_policy()$summary_timeout, 120),
  "default summary timeout is 120 s"
)
ok(
  withr_local_env(list(OMICSVIEWER_LLM_SUMMARY_TIMEOUT = "9999"),
    ut_cmp_identical(agent_context_policy()$summary_timeout, 120)),
  "summary timeout is hard-capped at 120 s"
)
ok(
  withr_local_env(list(OMICSVIEWER_LLM_SUMMARY_TIMEOUT = "30"),
    ut_cmp_identical(agent_context_policy()$summary_timeout, 30)),
  "summary timeout env override inside the cap"
)
ok(
  withr_local_env(list(OMICSVIEWER_LLM_SUMMARY_TIMEOUT = "2"),
    ut_cmp_identical(agent_context_policy()$summary_timeout, 120)),
  "sub-floor summary timeout falls back to the default"
)
ok(
  ut_cmp_identical(agent_context_policy()$compaction_backoff, 300L),
  "default compaction backoff is 300 s"
)
ok(
  withr_local_env(list(OMICSVIEWER_LLM_COMPACTION_BACKOFF = "0"),
    ut_cmp_identical(agent_context_policy()$compaction_backoff, 0L)),
  "compaction backoff can be disabled with 0"
)

parts_convo <- list(mk_user("hi"), mk_assistant("hello"))
p1 <- agent_context_digest_parts(parts_convo)
ok(ut_cmp_identical(p1$count, 2L), "digest parts count")
ok(length(p1$sizes) == 2 && all(is.finite(p1$sizes)), "digest parts sizes")
ok(
  ut_cmp_identical(agent_context_digest(parts_convo), paste(p1$count, sum(p1$sizes))),
  "digest string stays consistent with its parts"
)
p2 <- agent_context_digest_parts(list(parts_convo[[1]], mk_assistant("hello world")))
ok(
  p2$count == p1$count && p2$sizes[2] != p1$sizes[2],
  "digest parts localise a mutated turn"
)
ok(
  ut_cmp_identical(agent_context_digest_parts(list())$count, 0L),
  "digest parts of an empty list"
)

ok(ut_cmp_identical(agent_context_digest(list()), "empty"),
   "empty digest unchanged after the parts refactor")

## -- promise race --------------------------------------------------------

wait_promise <- function(p, seconds = 2) {
  out <- list(value = NULL, error = NULL, done = FALSE)
  promises::then(
    p,
    function(v) { out$value <<- v; out$done <<- TRUE },
    function(e) { out$error <<- e; out$done <<- TRUE }
  )
  # A single later::run_now() returns as soon as the queue momentarily
  # drains, cutting promise chains mid-hop; poll instead (Sys.sleep
  # advances wall time so future-due callbacks mature).
  deadline <- Sys.time() + seconds
  while (!isTRUE(out$done) && Sys.time() < deadline) {
    later::run_now(0.02)
    if (!isTRUE(out$done)) Sys.sleep(0.01)
  }
  out
}

delayed <- function(value, seconds)
  promises::promise(function(resolve, reject)
    later::later(function() resolve(value), seconds))

res <- wait_promise(agent_promise_race(list(
  promises::promise_resolve(list(method = "llm")),
  delayed(list(method = "timeout"), 0.05)
)))
ok(
  res$done && ut_cmp_identical(res$value$method, "llm"),
  "race: already-resolved promise wins immediately"
)

res <- wait_promise(agent_promise_race(list(
  delayed(list(method = "llm"), 0.2),
  delayed(list(method = "timeout"), 0.02)
)))
ok(
  res$done && ut_cmp_identical(res$value$method, "timeout"),
  "race: earlier-settling timer beats the slow LLM promise"
)

res <- wait_promise(agent_promise_race(list(
  promises::promise(function(resolve, reject)
    later::later(function() reject("boom"), 0.01)),
  delayed(list(method = "timeout"), 0.2)
)))
ok(
  res$done && inherits(res$error, "error") &&
    ut_cmp_identical(conditionMessage(res$error), "boom"),
  "race: first rejection propagates (promises wraps string rejections as errors)"
)
