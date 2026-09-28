library(omicsViewer)
library(unittest, quietly = TRUE)

# WP15: canonical tool-argument transport (the agent boundary seam).
# Boundary-level tests, not consumer-level: the shape handlers receive
# through the REAL ellmer dispatch is asserted directly, so an ellmer
# upgrade that changes convert = FALSE semantics fails here instead of
# in a user session.

if (!requireNamespace("ellmer", quietly = TRUE) ||
    !requireNamespace("shinychat", quietly = TRUE)) {
  ok(TRUE, "agent tool seam tests skipped because optional packages are unavailable")
  quit(save = "no", status = 0)
}

agent_args_sanitize <- omicsViewer:::agent_args_sanitize
agent_tool <- omicsViewer:::agent_tool

# ---- 1. agent_args_sanitize canonicalization rules ------------------------

ok(
  ut_cmp_identical(agent_args_sanitize("null"), NULL) &&
    ut_cmp_identical(agent_args_sanitize("NULL"), NULL) &&
    ut_cmp_identical(agent_args_sanitize("{}"), NULL) &&
    ut_cmp_identical(agent_args_sanitize("[]"), NULL),
  "top-level sentinel strings become NULL"
)
ok(
  ut_cmp_identical(agent_args_sanitize(list(a = "null", b = "Gene1")),
                   list(a = NULL, b = "Gene1")),
  "sentinels inside objects become NULL entries"
)
ok(
  ut_cmp_identical(agent_args_sanitize(list(a = list(b = "x", c = "[]"))),
                   list(a = list(b = "x", c = NULL))),
  "sentinels inside nested objects are dropped to NULL entries"
)
ok(
  ut_cmp_identical(
    agent_args_sanitize(list(sections = list("annotations", "panels"))),
    list(sections = c("annotations", "panels"))
  ),
  "scalar string arrays become atomic vectors"
)
ok(
  ut_cmp_identical(
    agent_args_sanitize(list(values = list(1, 2.5))),
    list(values = c(1, 2.5))
  ),
  "scalar numeric arrays become atomic numeric vectors"
)
ok(
  ut_cmp_identical(
    agent_args_sanitize(list(features = list("F1"))),
    list(features = "F1")
  ),
  "single-element scalar arrays become length-1 vectors"
)
ok(
  ut_cmp_identical(agent_args_sanitize(list(x = list())), list(x = list())),
  "empty arrays stay explicit empty lists (distinct from omitted)"
)
ok(
  ut_cmp_identical(
    agent_args_sanitize(list(x = list(NULL), y = list(NULL, NULL))),
    list(x = list(), y = list())
  ),
  "all-null arrays (the [null] artifact) stay explicit empty lists"
)
# the min-rule distinction: an explicit empty array reaches the store
# as an empty list (rejected by min-bound multi_selects), while an
# omitted or null-valued key never reaches it at all
ok(
  ut_cmp_identical(
    agent_args_sanitize(list(columns = list(), theme = "null")),
    list(columns = list(), theme = NULL)
  ),
  "empty arrays survive while sentinel strings become named NULLs"
)
ok(
  ut_cmp_identical(
    agent_args_sanitize(list(all = list(list(column = "a"), list(column = "b")))),
    list(all = list(list(column = "a"), list(column = "b")))
  ),
  "arrays of objects stay lists of named lists (the nested combinator shape)"
)
ok(
  ut_cmp_identical(
    agent_args_sanitize(list(mixed = list("a", list(column = "b")))),
    list(mixed = list("a", list(column = "b")))
  ),
  "mixed arrays keep list shape"
)
ok(
  ut_cmp_identical(agent_args_sanitize(list(a = 1L, b = TRUE, c = "x")),
                   list(a = 1L, b = TRUE, c = "x")),
  "atomic scalars pass through untouched"
)
ok(
  ut_cmp_identical(agent_args_sanitize(list(a = c("x", "y"))),
                   list(a = c("x", "y"))),
  "already-atomic vectors pass through untouched"
)
ok(
  ut_cmp_identical(agent_args_sanitize(list(a = NA)), list(a = NA)),
  "NA is data, not a sentinel, and survives"
)

# idempotence: sanitizing a canonical document is the identity
canonical_doc <- list(
  spec = list(
    geom = "point",
    filter = NULL,
    params = list(color = "#bdbdbd", alpha = 0.4),
    all = list(list(column = "a", op = ">", value = 1),
               list(column = "b", op = "<", value = 2))
  ),
  ids = c("F1", "F2")
)
ok(
  ut_cmp_identical(agent_args_sanitize(canonical_doc), canonical_doc),
  "sanitize is idempotent on canonical documents"
)

deep <- list()
for (i in seq_len(70)) deep <- list(deep)
ok(
  ut_cmp_error(agent_args_sanitize(deep), "at most 64 levels"),
  "pathologically deep arguments are rejected, not recursed"
)

# ---- 2. the agent_tool wrapper (direct + ellmer dispatch) -----------------

captured <- NULL
record_tool <- agent_tool(
  function(spec = NULL, filter = NULL, ids = NULL, `_intent`) {
    captured <<- list(spec = spec, filter = filter, ids = ids,
                      intent = `_intent`)
    "recorded"
  },
  name = "seam_record",
  description = "records sanitized arguments",
  arguments = list(
    spec = ellmer::type_string("s", required = FALSE),
    filter = ellmer::type_object(
      "f",
      column = ellmer::type_string("c", required = FALSE),
      all = ellmer::type_array(
        ellmer::type_object("leaf",
                            column = ellmer::type_string("c", required = FALSE)),
        required = FALSE),
      .required = FALSE
    ),
    ids = ellmer::type_array(ellmer::type_string("id"), required = FALSE),
    `_intent` = ellmer::type_string("intent")
  )
)
ok(
  ut_cmp_identical(record_tool@convert, FALSE),
  "agent_tool registers with convert = FALSE"
)

# the wire payload from the live failure class: absent nested object,
# nested combinator array, scalar string array, sentinel strings
wire <- jsonlite::fromJSON(paste0(
  '{"_intent":"x","ids":["F1","F2"],',
  '"spec":"null",',
  '"filter":{"all":[{"column":"a"},{"column":"b"}]}}'
), simplifyVector = FALSE)
request <- ellmer::ContentToolRequest(
  id = "c1", name = "seam_record", arguments = wire)
request@tool <- record_tool
result <- ellmer:::invoke_tool(request)
ok(
  ut_cmp_identical(result@error, NULL) &&
    ut_cmp_identical(result@value, "recorded"),
  "the failing live payload class dispatches without error"
)
ok(
  ut_cmp_identical(captured$ids, c("F1", "F2")) &&
    ut_cmp_identical(captured$spec, NULL) &&
    ut_cmp_identical(captured$filter,
                     list(all = list(list(column = "a"), list(column = "b")))),
  "handlers receive the canonical document through real ellmer dispatch"
)

# direct calls on the ToolDef run through the same wrapper
captured <<- NULL
direct <- record_tool(spec = "s", filter = NULL, ids = list("F1"), `_intent` = "y")
ok(
  ut_cmp_identical(captured$ids, "F1") &&
    ut_cmp_identical(captured$spec, "s") &&
    ut_cmp_identical(captured$filter, NULL),
  "direct ToolDef calls pass through the identical seam"
)

# unknown arguments are still a self-correctable error (ellmer's
# convert = TRUE unused-argument check preserved at the wrapper)
extra <- jsonlite::fromJSON(
  '{"_intent":"x","bogus":"field"}', simplifyVector = FALSE)
request2 <- ellmer::ContentToolRequest(
  id = "c2", name = "seam_record", arguments = extra)
request2@tool <- record_tool
result2 <- ellmer:::invoke_tool(request2)
ok(
  ut_cmp_identical(is.null(result2@error), FALSE) &&
    grepl("unused argument", conditionMessage(result2@error)),
  "unknown tool arguments error for model self-correction"
)

# ---- 3. round-trip property test over the REAL tool schemas --------------

# For every registered tool schema and generated payload variant, what
# the handler receives through real ellmer dispatch must equal the
# canonical payload: the tested shape IS the shipped shape.
seam_sample <- function(type, variant) {
  if (S7::S7_inherits(type, ellmer::TypeArray)) {
    if (identical(variant, "empty")) return(list())
    lapply(seq_len(2L), function(...)
      seam_sample(type@items,
                  if (identical(variant, "empty")) "full" else variant))
  } else if (S7::S7_inherits(type, ellmer::TypeObject)) {
    out <- list()
    for (nm in names(type@properties)) {
      prop <- type@properties[[nm]]
      if (isTRUE(prop@required)) {
        out[[nm]] <- seam_sample(prop, variant)
      } else if (identical(variant, "full")) {
        out[[nm]] <- seam_sample(prop, "full")
      } else if (identical(variant, "nulls")) {
        out[[nm]] <- NULL
      } else if (identical(variant, "sentinel")) {
        out[[nm]] <- if (S7::S7_inherits(prop, ellmer::TypeObject) ||
                         S7::S7_inherits(prop, ellmer::TypeArray)) "{}" else "null"
      }
      # "omitted": leave the property out entirely
    }
    out
  } else if (S7::S7_inherits(type, ellmer::TypeEnum)) {
    type@values[[1]]
  } else {
    switch(type@type,
           string = "seam", number = 1.5, integer = 2L, boolean = TRUE, "seam")
  }
}

fd <- data.frame(
  score = c(1, 2, 3), logFdr = c(5, 4, 3),
  category = c("kinase", "phosphatase", "kinase"),
  row.names = c("Gene1", "Gene2", "Gene3"), check.names = FALSE)
pd <- data.frame(group = c("WT", "WT", "KO", "KO"),
                 row.names = c("S1", "S2", "S3", "S4"))
state_builder <- function(sections = NULL) list(
  dataset = list(id = "demo.RDS", class = "ExpressionSet",
                 dimensions = c(features = 3L, samples = 4L)),
  active_tabs = list(data_space = "Feature", analysis_space = "Feature"),
  selection = list(
    features = list(count = 0L, ids = character(), truncated = FALSE),
    samples = list(count = 0L, ids = character(), truncated = FALSE)),
  available_tabs = list(data_space = c("Feature", "Sample"),
                        analysis_space = c("Feature", "ORA")))

variants <- c("full", "nulls", "omitted", "sentinel", "empty")
roundtrip_failures <- character()

shiny::testServer(
  omicsViewer:::ai_assistant_module,
  args = list(
    state = state_builder,
    state_available = shiny::reactive(TRUE),
    feature_data = shiny::reactive(fd),
    sample_data = shiny::reactive(pd),
    expression_data = shiny::reactive(matrix(
      1:12, nrow = 3, dimnames = list(c("Gene1", "Gene2", "Gene3"), rownames(pd)))),
    selected_features = shiny::reactive(character()),
    selected_samples = shiny::reactive(character()),
    apply_state = function(u) list(),
    apply_scatter_view = function(...) list(),
    apply_enrichment = function(u) list(),
    apply_table_view = function(u) list(),
    store = omicsViewer:::widget_store_new()
  ),
  expr = {
    tools <- chat_object$client$get_tools()
    for (tool_name in names(tools)) {
      schema <- tools[[tool_name]]@arguments
      for (variant in variants) {
        payload <- seam_sample(schema, variant)
        expected <- agent_args_sanitize(payload)
        # absent top-level arguments are dropped at dispatch (ellmer's
        # own convert = TRUE null-dropping), so the capture never sees
        # them as explicit NULL entries
        expected <- expected[!vapply(expected, is.null, logical(1))]
        capture <- agent_tool(
          function(...) ellmer::ContentToolResult(value = list(...)),
          description = "capture",
          name = paste0("capture_", tool_name),
          arguments = schema@properties)
        req <- ellmer::ContentToolRequest(
          id = "rt", name = capture@name, arguments = payload)
        req@tool <- capture
        res <- tryCatch(
          ellmer:::invoke_tool(req),
          error = function(e) e)
        received <- if (inherits(res, "error")) res else res@value
        if (inherits(res, "error")) {
          roundtrip_failures <<- c(roundtrip_failures,
            paste(tool_name, variant, "ERROR:", conditionMessage(res)))
        } else if (!isTRUE(identical(received, expected))) {
          diff <- paste(
            tool_name, variant, "GOT",
            paste(capture.output(str(received)), collapse = " | "),
            "EXPECTED",
            paste(capture.output(str(expected)), collapse = " | "))
          roundtrip_failures <<- c(roundtrip_failures, substr(diff, 1L, 400L))
        }
      }
    }
    all_convert_false <<- all(vapply(
      tools, function(t) isFALSE(t@convert), logical(1)))
    tool_names <<- names(tools)
  })

ok(
  ut_cmp_identical(roundtrip_failures, character()),
  paste("every tool schema round-trips every payload variant through",
        "the real ellmer dispatch unchanged")
)
ok(
  ut_cmp_identical(all_convert_false, TRUE),
  "every registered tool travels on convert = FALSE"
)
ok(
  ut_cmp_identical(length(tool_names) >= 14L, TRUE),
  "the property sweep covers the full registered tool surface"
)

# ---- 4. golden traffic replay: the live failing payloads ------------------

# tests/fixtures/agent/update_figure_seq{35,55}.json are the exact
# payloads from live session 20260927-064108 (demo.RDS, glm via
# OpenRouter): 6/6 update_figure calls rejected with
# "Figure filter requires a column." because ellmer materialized the
# absent layer filter as a 1-row all-NA tibble and the nested all
# combinator as list(<tibble>). Both must render first-try through the
# seam, on the demo dataset, with the exact payloads.

demo <- readRDS(file.path("inst", "extdata", "demo.RDS"))
demo_fd <- Biobase::fData(demo)
demo_pd <- Biobase::pData(demo)
demo_mat <- Biobase::exprs(demo)

replay_errors <- list()
replay_filters <- list()
replay_warnings <- list()

shiny::testServer(
  omicsViewer:::ai_assistant_module,
  args = list(
    state = state_builder,
    state_available = shiny::reactive(TRUE),
    feature_data = shiny::reactive(demo_fd),
    sample_data = shiny::reactive(demo_pd),
    expression_data = shiny::reactive(demo_mat),
    selected_features = shiny::reactive(rownames(demo_fd)[1:5]),
    selected_samples = shiny::reactive(character()),
    apply_state = function(u) list(),
    apply_scatter_view = function(...) list(),
    apply_enrichment = NULL, apply_table_view = NULL, store = NULL
  ),
  expr = {
    tools <- chat_object$client$get_tools()
    seed <- tools$create_figure(
      template = "volcano",
      x = "ttest|RE_vs_LE|mean.diff", y = "ttest|RE_vs_LE|log.fdr",
      `_intent` = "seed the session figure fig_1"
    )
    for (fx in c("update_figure_seq35.json", "update_figure_seq55.json")) {
      payload <- jsonlite::fromJSON(
        file.path("tests", "fixtures", "agent", fx), simplifyVector = FALSE)
      req <- ellmer::ContentToolRequest(
        id = "replay", name = "update_figure", arguments = payload)
      req@tool <- tools$update_figure
      res <- ellmer:::invoke_tool(req)
      # NB: assigning NULL via [[<- would REMOVE the entry (R list
      # semantics), so success is recorded as the string "ok".
      replay_errors[[fx]] <<- if (is.null(res@error)) "ok" else
        conditionMessage(res@error)
      replay_filters[[fx]] <<- if (is.null(res@error)) list(
        layer2 = res@value$spec$layers[[2]]$filter,
        layer3 = res@value$spec$layers[[3]]$filter
      ) else NULL
      replay_warnings[[fx]] <<- if (is.null(res@error))
        (res@value$warnings %||% character()) else character("failed")
    }
  })

for (fx in names(replay_errors)) {
  ok(
    ut_cmp_identical(replay_errors[[fx]], "ok"),
    paste("golden replay:", fx, "renders first-try (live 20260927-064108 rejection)")
  )
  ok(
    ut_cmp_identical(
      length(replay_filters[[fx]]$layer2$all), 2L) &&
      ut_cmp_identical(
        replay_filters[[fx]]$layer2$all[[1]]$column,
        "ttest|RE_vs_LE|mean.diff") &&
      ut_cmp_identical(
        replay_filters[[fx]]$layer3$all[[2]]$column,
        "ttest|RE_vs_LE|log.fdr"),
    paste("golden replay:", fx, "carries the intended quadrant filters")
  )
  ok(
    ut_cmp_identical(
      any(grepl("filter", replay_warnings[[fx]], ignore.case = TRUE,
                fixed = FALSE)),
      FALSE),
    paste("golden replay:", fx, "raises no filter warnings")
  )
}
