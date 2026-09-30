library(omicsViewer)
library(unittest, quietly = TRUE)
library(ellmer)

## T0 contract tests (todo 1.3 / 4.5): pin the ellmer / shinychat semantics
## the assistant module depends on, so an upstream upgrade fails LOUDLY here
## instead of silently dead-guarding the janitor/compaction observers again.
##
## Background: shinychat's Chat$status() only ever returns "idle" or
## "streaming" (verified against shinychat 0.5.0 docs/source). The module
## guards compare with "streaming"; if shinychat ever renames the busy
## state, the guards go dead exactly like the old "running" comparisons did.

## ---------------------------------------------------------------- status --
if (requireNamespace("shinychat", quietly = TRUE)) {
  shiny::testServer(
    function(input, output, session) {
      client <- ellmer::chat_openai(api_key = "sk-contract-test")
      chat <- shinychat::chat_server("chat", client, history = FALSE)
      status_now <- chat$status()
      ok(
        ut_cmp_identical(status_now %in% c("idle", "streaming"), TRUE),
        "shinychat status() reports idle/streaming vocabulary"
      )
      ok(
        ut_cmp_identical(status_now, "idle"),
        "a freshly created chat reports idle"
      )
      chat_status <- chat$status
    },
    expr = NULL
  )
}

## ------------------------------------------------- tool value projection --
req <- ContentToolRequest(id = "c1", name = "t", arguments = list(a = 1))
res <- ContentToolResult(value = list(ok = TRUE), request = req)
ok(
  ut_cmp_identical(is.null(contents_text(res)), TRUE),
  "contents_text() of a tool result is NULL (tool outputs are invisible to text consumers)"
)

## tool_string() returns @value verbatim: json-class strings (the planned
## 3.2 output codec) must round-trip untouched.
json_value <- jsonlite::toJSON(list(b = 1, a = c("x", "y")), auto_unbox = TRUE)
json_res <- ContentToolResult(value = json_value, request = req)
ok(
  ut_cmp_identical(ellmer:::`tool_string`(json_res), json_value),
  "tool_string() passes json-class strings through unchanged"
)
ok(
  ut_cmp_identical(
    as.character(jsonlite::toJSON(c(a = 1, b = 2), auto_unbox = TRUE)),
    "[1,2]"
  ),
  "jsonlite drops names on atomic vectors (codec must convert to lists first)"
)
ok(
  ut_cmp_identical(
    as.character(jsonlite::toJSON(as.list(c(a = 1, b = 2)), auto_unbox = TRUE)),
    "{\"a\":1,\"b\":2}"
  ),
  "jsonlite preserves names once vectors are converted to lists"
)

## ------------------------------------------------ runtime ToolDef refs ----
# The janitor strips the ToolDef ellmer attaches to tool requests during a
# live stream. A live-stream attachment assertion needs a real provider and
# is covered by the Tier B harness; here we pin the strip contract instead.
secret <- "sk-runtime-ref-sentinel"
secret_env <- new.env(parent = emptyenv())
secret_env$key <- secret
secret_handler <- function(x) list(x = x)
environment(secret_handler) <- secret_env
tool_def <- ellmer::tool(
  secret_handler,
  name = "echo_secret",
  description = "closes over a sentinel secret",
  arguments = list(x = ellmer::type_string("x"))
)
live_req <- ContentToolRequest(
  id = "c2", name = "echo_secret", arguments = list(x = "a"), tool = tool_def
)
turn <- UserTurn(list(live_req))
stripped <- omicsViewer:::agent_strip_runtime_refs(list(turn))
serialized <- serialize(stripped, NULL)
ok(
  ut_cmp_identical(
    identical(length(grepRaw(secret, serialized, fixed = TRUE)), 0L),
    TRUE
  ),
  "agent_strip_runtime_refs removes handler closures from turns"
)
