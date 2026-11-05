# Stage 0 unit suite for the transport-neutral agent channel
# (auxi_agentChannel.R): JSON-Schema emitter, dispatcher, and
# session-token registry. Pure helper logic needs no browser; the
# ToolDef entry is the verified transport-clean surface (direct call =
# sanitize-in / codec-out, exactly what a model drives).

library(unittest, quietly = TRUE)

if (!requireNamespace("pkgload", quietly = TRUE) ||
    !requireNamespace("ellmer", quietly = TRUE)) {
  ok(TRUE, "agent channel tests skipped because optional packages are unavailable")
  quit(save = "no", status = 0)
}

# The package is NOT reinstalled for this suite; load the working tree
# (the scatterDefaults/snapshotRoundTrip convention).
pkg_root <- if (file.exists("DESCRIPTION")) "." else
  if (file.exists(file.path("..", "DESCRIPTION"))) ".." else "."
pkgload::load_all(pkg_root, quiet = TRUE, attach_testthat = FALSE)
library(ellmer)

agent_tool <- omicsViewer:::agent_tool
agent_channel_toolspecs <- omicsViewer:::agent_channel_toolspecs
agent_channel_dispatch <- omicsViewer:::agent_channel_dispatch
agent_channel_registry_bind <- omicsViewer:::agent_channel_registry_bind
agent_channel_registry_unbind <- omicsViewer:::agent_channel_registry_unbind
agent_channel_registry_list <- omicsViewer:::agent_channel_registry_list
.agent_channel_type_spec <- omicsViewer:::`.agent_channel_type_spec`

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

# Collect every key name reachable in a nested list (D15 assertion).
spec_key_names <- function(x) {
  if (is.list(x))
    unique(c(names(x), unlist(lapply(x, spec_key_names),
                               use.names = FALSE)))
  else
    character()
}

# ---- 1. spec emitter: toy ToolDefs over every used ellmer type ------

# Every type constructor the registry actually uses (verified by the
# static sweep in section 1b): type_string, type_integer, type_number,
# type_boolean, type_enum, type_array(string), type_array(enum),
# type_array(object), type_object with properties/required.
chan_probe <- agent_tool(
  function(space, ids = NULL, counts = NULL, flag = NULL, ratio = NULL,
           sections = NULL, filter = NULL, `_intent`) "ok",
  name = "chan_probe",
  description = "channel schema probe",
  arguments = list(
    space = ellmer::type_enum(c("feature", "sample"), "Annotation space."),
    ids = ellmer::type_array(
      ellmer::type_string("Exact ID."),
      "Optional IDs to scope the call.", required = FALSE),
    counts = ellmer::type_array(
      ellmer::type_integer("Count, from 1 through 50."), required = FALSE),
    flag = ellmer::type_boolean("Optional flag.", required = FALSE),
    ratio = ellmer::type_number(
      "Optional ratio, from 0 through 1.", required = FALSE),
    sections = ellmer::type_array(
      ellmer::type_enum(c("a", "b"), "State section."),
      "Optional sections.", required = FALSE),
    filter = ellmer::type_object(
      "Optional nested filter.",
      column = ellmer::type_string("Exact column name."),
      all = ellmer::type_array(
        ellmer::type_object(
          "One filter leaf.",
          column = ellmer::type_string("Exact column name."),
          op = ellmer::type_enum(c(">", "<"), "Filter operator.")),
        "Combinator leaves.", required = FALSE),
      .required = FALSE),
    `_intent` = ellmer::type_string("Short user-facing reason.")
  )
)

specs <- agent_channel_toolspecs(list(probe = chan_probe))
ok(
  ut_cmp_identical(names(specs), "chan_probe"),
  "specs are keyed by the wire tool name, not the list key"
)
spec <- specs$chan_probe
ok(
  ut_cmp_identical(sort(names(spec)), c("description", "inputSchema", "name")) &&
    ut_cmp_identical(spec$name, "chan_probe") &&
    ut_cmp_identical(spec$description, "channel schema probe"),
  "each spec carries name/description/inputSchema"
)

schema <- spec$inputSchema
ok(
  ut_cmp_identical(schema$type, "object") &&
    setequal(names(schema$properties),
          c("counts", "filter", "flag", "ids", "ratio", "sections",
            "space", "_intent")) &&
    setequal(schema$required, c("space", "_intent")),
  "inputSchema is the arguments object: properties + required array"
)
ok(
  ut_cmp_identical(is.null(schema$description), TRUE) &&
    ut_cmp_identical(schema$additionalProperties, FALSE),
  "empty object descriptions are omitted; additionalProperties FALSE lands (ellmer 0.5.0 default)"
)
ok(
  ut_cmp_identical(schema$properties$space$type, "string") &&
    ut_cmp_identical(schema$properties$space$enum, c("feature", "sample")) &&
    ut_cmp_identical(schema$properties$space$description, "Annotation space."),
  "enum values land with their description"
)
ok(
  ut_cmp_identical(schema$properties$ids$type, "array") &&
    ut_cmp_identical(schema$properties$ids$items$type, "string") &&
    ut_cmp_identical(schema$properties$ids$description,
                     "Optional IDs to scope the call."),
  "type_array(type_string) lands as array/items with description"
)
ok(
  ut_cmp_identical(schema$properties$sections$items$type, "string") &&
    ut_cmp_identical(schema$properties$sections$items$enum, c("a", "b")),
  "type_array(type_enum) lands the nested enum on items"
)
ok(
  ut_cmp_identical(schema$properties$counts$items$type, "integer") &&
    ut_cmp_identical(schema$properties$flag$type, "boolean") &&
    ut_cmp_identical(schema$properties$ratio$type, "number"),
  "integer/boolean/number basic types map to their JSON types"
)
nested <- schema$properties$filter
ok(
  ut_cmp_identical(nested$type, "object") &&
    ut_cmp_identical(nested$required, "column") &&
    ut_cmp_identical(nested$additionalProperties, FALSE) &&
    ut_cmp_identical(nested$properties$all$items$type, "object") &&
    ut_cmp_identical(nested$properties$all$items$properties$op$enum, c(">", "<")),
  "nested objects carry their own required/additionalProperties/enum"
)
ok(
  ut_cmp_identical(
    any(spec_key_names(schema) %in%
          c("minimum", "maximum", "exclusiveMinimum", "exclusiveMaximum")),
    FALSE) &&
    grepl("from 0 through 1", schema$properties$ratio$description,
          fixed = TRUE),
  "numeric ranges stay description-only (D15); no min/max schema keys anywhere"
)
# JSON-ready at the edge: toJSON round-trip projects the spec cleanly
edge <- tryCatch(
  jsonlite::toJSON(specs, auto_unbox = TRUE, null = "null"),
  error = function(e) e)
ok(
  !inherits(edge, "error") &&
    identical(
      jsonlite::fromJSON(edge)$chan_probe$inputSchema$type,
      "object"),
  "specs serialize through toJSON(auto_unbox) at the edge without error"
)
ok(
  ut_cmp_identical(agent_channel_toolspecs(list()), list()),
  "an empty tool list emits an empty spec list"
)
ok(
  ut_cmp_error(agent_channel_toolspecs(list(x = 1)),
               "not an ellmer::ToolDef"),
  "non-ToolDef elements fail loudly"
)
ok(
  ut_cmp_error(.agent_channel_type_spec(ellmer::type_ignore()),
               "Unsupported ellmer type"),
  "unknown ellmer types are rejected, never silently degraded"
)

# ---- 1b. static coverage: every registry constructor has a branch ----

# Defensive sweep: read the registry (and the figure grammar it embeds
# via agent_figure_spec_schema), extract every type_* constructor call,
# and require a converter branch for each. A new ellmer type appearing
# in the registry fails HERE first, not at a channel client.
src_paths <- file.path(
  pkg_root,
  c("R/auxi_agentToolRegistry.R", "R/auxi_agentFigureGrammar.R"))
src_paths <- src_paths[file.exists(src_paths)]
registry_src <- paste(
  vapply(src_paths, function(f) paste(readLines(f, warn = FALSE),
                                      collapse = "\n"), ""),
  collapse = "\n")
used <- unique(regmatches(
  registry_src,
  gregexpr("\\btype_[a-z]+\\s*\\(", registry_src))[[1]])
used <- unique(sub("\\s*\\($", "", used))

representatives <- list(
  type_string = function() ellmer::type_string("d", required = FALSE),
  type_integer = function() ellmer::type_integer("d", required = FALSE),
  type_number = function() ellmer::type_number("d", required = FALSE),
  type_boolean = function() ellmer::type_boolean("d", required = FALSE),
  type_enum = function() ellmer::type_enum(c("a", "b"), "d", required = FALSE),
  type_array = function() ellmer::type_array(
    ellmer::type_string("i"), "d", required = FALSE),
  type_object = function() ellmer::type_object(
    "d", x = ellmer::type_string("x"))
)
uncovered <- setdiff(used, names(representatives))
ok(
  ut_cmp_identical(uncovered, character()),
  paste("every type constructor used by the registry has an emitter",
        "branch (used: ", paste(sort(used), collapse = ", "), ")", sep = "")
)
for (constructor in intersect(used, names(representatives))) {
  converted <- tryCatch(
    .agent_channel_type_spec(representatives[[constructor]]()),
    error = function(e) e)
  ok(
    !inherits(converted, "error") && !is.null(converted$type),
    paste("constructor", constructor, "converts to a typed schema fragment")
  )
}

# The REAL figure spec schema (deeply nested object/array/enum shapes
# actually embedded in the registry) converts completely.
figure_spec <- tryCatch(
  .agent_channel_type_spec(
    omicsViewer:::agent_figure_spec_schema(
      required = FALSE, layers_required = FALSE)),
  error = function(e) e)
ok(
  !inherits(figure_spec, "error") &&
    identical(figure_spec$type, "object") &&
    identical(figure_spec$properties$layers$items$type, "object") &&
    identical(figure_spec$properties$layers$items$properties$geom$type,
              "string") &&
    length(figure_spec$properties$layers$items$properties$geom$enum) > 0L,
  "the real nested figure schema converts with geoms landing as enums"
)

# ---- 2. dispatcher round-trips ---------------------------------------

last_call <- NULL
chan_echo <- agent_tool(
  function(space, n = 2L, ids = NULL, `_intent`) {
    last_call <<- list(space = space, n = n, ids = ids,
                       intent = `_intent`)
    ellmer::ContentToolResult(
      value = list(space = space, n = n, ids = ids))
  },
  name = "chan_echo",
  description = "channel dispatch probe",
  arguments = list(
    space = ellmer::type_enum(c("feature", "sample"), "Space."),
    n = ellmer::type_integer("Count, from 1 through 50.", required = FALSE),
    ids = ellmer::type_array(ellmer::type_string("id."), required = FALSE),
    `_intent` = ellmer::type_string("intent.")
  )
)
chan_tools <- list(echo = chan_echo)

# (a) plain args: dispatcher result identical to the direct ToolDef call
d <- agent_channel_dispatch(
  chan_tools, "chan_echo",
  list(space = "feature", n = 5L, ids = list("F1", "F2"),
       `_intent` = "probe"))
via_dispatch <- last_call
last_call <<- NULL
direct <- chan_echo(space = "feature", n = 5L, ids = list("F1", "F2"),
                    `_intent` = "probe")
via_direct <- last_call
ok(
  ut_cmp_identical(d$ok, TRUE) && is.null(d$error),
  "plain dispatch succeeds"
)
ok(
  ut_cmp_identical(via_dispatch, via_direct),
  "the handler sees identical canonical args through dispatch and direct call"
)
ok(
  ut_cmp_identical(d$value, direct),
  "dispatch returns exactly the direct ToolDef result (codec included)"
)
ok(
  ut_cmp_identical(
    agent_channel_dispatch(chan_tools, "echo", list(
      space = "sample", `_intent` = "w"))$ok, TRUE),
  "lookup also accepts the list key"
)

# (b) omitted optional args: handler defaults apply
d_b <- agent_channel_dispatch(
  chan_tools, "chan_echo", list(space = "sample", `_intent` = "w"))
ok(
  d_b$ok &&
    ut_cmp_identical(d_b$value@extra$data$n, 2L) &&
    is.null(d_b$value@extra$data$ids),
  "omitted optionals take the handler defaults"
)

# (c) explicit NULL args are dropped before the call (mirror agent_tool)
d_c <- agent_channel_dispatch(
  chan_tools, "chan_echo",
  list(space = "sample", n = NULL, ids = NULL, `_intent` = "w"))
ok(
  d_c$ok && ut_cmp_identical(d_c$value@extra$data, d_b$value@extra$data),
  "explicit NULL args behave exactly like omitted args"
)

# (d) sentinel strings are treated as omitted (AGENT_SENTINEL_STRINGS)
sentinel_results <- lapply(
  omicsViewer:::AGENT_SENTINEL_STRINGS,
  function(s) agent_channel_dispatch(
    chan_tools, "chan_echo",
    list(space = "sample", n = s, ids = s, `_intent` = "w")))
ok(
  all(vapply(sentinel_results, function(r) r$ok, logical(1))) &&
    all(vapply(sentinel_results,
               function(r) identical(r$value@extra$data$n, 2L) &&
                 is.null(r$value@extra$data$ids),
               logical(1))),
  "every sentinel string arg is treated as omitted (defaults apply)"
)

# (e) handler errors: ok=FALSE, message captured, nothing escapes
chan_boom <- agent_tool(
  function(`_intent`) stop("chan boom: detonated"),
  name = "chan_boom", description = "error probe",
  arguments = list(`_intent` = ellmer::type_string("intent."))
)
d_e <- agent_channel_dispatch(
  list(boom = chan_boom), "chan_boom", list(`_intent` = "x"))
ok(
  ut_cmp_identical(d_e$ok, FALSE) && is.null(d_e$value) &&
    grepl("chan boom: detonated", d_e$error, fixed = TRUE),
  "handler errors return ok=FALSE with the condition message"
)
d_e2 <- agent_channel_dispatch(
  chan_tools, "chan_echo",
  list(space = "feature", bogus = "field", `_intent` = "x"))
ok(
  ut_cmp_identical(d_e2$ok, FALSE) && is.null(d_e2$value) &&
    grepl("bogus", d_e2$error, fixed = TRUE) &&
    grepl("space", d_e2$error, fixed = TRUE) &&
    grepl("n", d_e2$error, fixed = TRUE),
  "unknown argument names are rejected before the call, naming the valid ones"
)
d_e2b <- agent_channel_dispatch(
  chan_tools, "chan_echo", list(spa = "feature", `_intent` = "x"))
ok(
  ut_cmp_identical(d_e2b$ok, FALSE) &&
    grepl("spa", d_e2b$error, fixed = TRUE) &&
    grepl("space", d_e2b$error, fixed = TRUE),
  "partial-match argument names (spa for space) are rejected, not silently matched"
)
ok(
  ut_cmp_identical(
    agent_channel_dispatch(
      chan_tools, "chan_echo",
      list(space = "feature", `_intent` = "exact names"))$ok, TRUE),
  "exact argument names still dispatch (I9 does not over-reject)"
)
d_e3 <- agent_channel_dispatch(chan_tools, "chan_echo", NULL)
ok(
  ut_cmp_identical(d_e3$ok, FALSE) &&
    grepl("missing", d_e3$error),
  "missing required args return ok=FALSE (dispatch with no args)"
)

# (h) string tool errors: an errored ContentToolResult whose @error is
# a plain string (ellmer's S7 union allows NULL/character/condition,
# and ellmer itself produces string errors) must not trip
# conditionMessage() - the dispatcher returns ok=FALSE with the string.
chan_strerr <- agent_tool(
  function(`_intent`) ellmer::ContentToolResult(
    value = NULL, error = "plain string failure"),
  name = "chan_strerr", description = "string error probe",
  arguments = list(`_intent` = ellmer::type_string("intent."))
)
d_h <- agent_channel_dispatch(
  list(strerr = chan_strerr), "chan_strerr", list(`_intent` = "x"))
ok(
  ut_cmp_identical(d_h$ok, FALSE) && is.null(d_h$value) &&
    grepl("plain string failure", d_h$error, fixed = TRUE),
  "string tool errors return ok=FALSE carrying the string; no condition escapes"
)

# (f) unknown tool name: ok=FALSE naming the available tools
d_f <- agent_channel_dispatch(chan_tools, "no_such_tool",
                              list(`_intent` = "x"))
ok(
  ut_cmp_identical(d_f$ok, FALSE) && is.null(d_f$value) &&
    grepl("chan_echo", d_f$error, fixed = TRUE) &&
    grepl("Unknown tool", d_f$error, fixed = TRUE),
  "unknown tool names list the available tools"
)
ok(
  ut_cmp_identical(
    agent_channel_dispatch(chan_tools, "", list())$ok, FALSE) &&
    ut_cmp_identical(
      agent_channel_dispatch(chan_tools, "chan_echo", "not a list")$ok,
      FALSE) &&
    ut_cmp_identical(
      agent_channel_dispatch(chan_tools, "chan_echo", list("unnamed"))$ok,
      FALSE),
  "malformed name/args requests are rejected without conditions escaping"
)

# (g) oversize handler value: codec truncation marker in the result
big_ids <- paste0("gene", 1:800)
chan_big <- agent_tool(
  function(`_intent`) ellmer::ContentToolResult(value = list(ids = big_ids)),
  name = "chan_big", description = "truncation probe",
  arguments = list(`_intent` = ellmer::type_string("intent."))
)
d_g <- withr_local_env(
  list(OMICSVIEWER_LLM_TOOL_OUTPUT_BYTES = "2048"),
  agent_channel_dispatch(list(big = chan_big), "chan_big",
                         list(`_intent` = "w")))
ok(
  d_g$ok &&
    grepl("Tool result truncated", d_g$value@value, fixed = TRUE) &&
    grepl("Narrow the request", d_g$value@value, fixed = TRUE),
  "oversize values carry the codec truncation marker (env knob respected)"
)
ok(
  ut_cmp_identical(length(d_g$value@extra$data$ids), 800L),
  "the structured original survives truncation in extra$data"
)
ok(
  ut_cmp_identical(
    Sys.getenv("OMICSVIEWER_LLM_TOOL_OUTPUT_BYTES", ""), ""),
  "the env knob is unset after the truncation probe"
)

# ---- 3. registry lifecycle -------------------------------------------

reg_env <- omicsViewer:::`.agent_channel_registry`
tokens_used <- c("chan_tok_1", "chan_tok_2", "chan_tok_3")
on.exit(for (t in tokens_used) agent_channel_registry_unbind(t),
       add = TRUE)

ok(
  ut_cmp_identical(agent_channel_registry_list(), character()) ||
    !any(tokens_used %in% agent_channel_registry_list()),
  "registry starts without this suite's tokens"
)
vis <- withVisible(agent_channel_registry_bind("chan_tok_1", list(a = 1)))
ok(
  ut_cmp_identical(vis$value, "chan_tok_1") &&
    ut_cmp_identical(vis$visible, FALSE),
  "bind returns the token invisibly"
)
ok(
  "chan_tok_1" %in% agent_channel_registry_list() &&
    ut_cmp_identical(reg_env$chan_tok_1, list(a = 1)),
  "bind stores the entry keyed by token"
)
agent_channel_registry_bind("chan_tok_1", list(b = 2))
ok(
  ut_cmp_identical(reg_env$chan_tok_1, list(b = 2)) &&
    ut_cmp_identical(
      sum(agent_channel_registry_list() == "chan_tok_1"), 1L),
  "double-binding the same token overwrites (single entry)"
)
agent_channel_registry_unbind("chan_tok_1")
ok(
  !("chan_tok_1" %in% agent_channel_registry_list()) &&
    !exists("chan_tok_1", envir = reg_env, inherits = FALSE),
  "unbind removes the entry"
)
ok(
  ut_cmp_identical(
    tryCatch({agent_channel_registry_unbind("chan_never_bound"); "silent"},
             error = function(e) conditionMessage(e)),
    "silent"),
  "unbinding an absent token is silent"
)
ended <- 0L
agent_channel_registry_bind("chan_tok_2", list(), on_ended = function() {
  ended <<- ended + 1L
})
ok(
  ut_cmp_identical(ended, 1L) && "chan_tok_2" %in%
    agent_channel_registry_list(),
  "on_ended wiring runs exactly once, immediately after storing"
)
# a throwing on_ended fails the bind AND leaves no orphaned entry (F4)
bind_fail <- tryCatch(
  agent_channel_registry_bind(
    "chan_tok_3", list(a = 1),
    on_ended = function() stop("on_ended exploded")),
  error = function(e) e)
ok(
  inherits(bind_fail, "error") &&
    grepl("on_ended exploded", conditionMessage(bind_fail), fixed = TRUE) &&
    !("chan_tok_3" %in% agent_channel_registry_list()) &&
    !exists("chan_tok_3", envir = reg_env, inherits = FALSE),
  "a throwing on_ended fails the bind and removes the just-stored entry"
)
ok(
  ut_cmp_error(agent_channel_registry_bind(1L, list()), "single non-empty string"),
  "non-string tokens are rejected"
)
