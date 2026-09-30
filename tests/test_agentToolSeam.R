library(omicsViewer)
library(unittest, quietly = TRUE)

# WP15: canonical tool-argument transport (the agent boundary seam).
# Boundary-level tests, not consumer-level: the shape handlers receive
# through the REAL ellmer dispatch is asserted directly, so an ellmer
# upgrade that changes convert = FALSE semantics fails here instead of
# in a user session.

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

# ---- 1b. sentinel strings never disappear inside DATA arrays (3.3) ------
# The sentinel rule exists for omitted OBJECT FIELDS. Inside id/value
# arrays the same strings are legitimate data (a feature literally named
# "NULL" was silently dropped pre-Stage-2), so array elements keep them.
ok(
  ut_cmp_identical(
    agent_args_sanitize(list(features = list("TP53", "NULL", "null"))),
    list(features = c("TP53", "NULL", "null"))
  ),
  "sentinel-equal strings inside id arrays survive sanitization"
)
ok(
  ut_cmp_identical(
    agent_args_sanitize(list(values = list("[]", "a"), x = "[]")),
    list(values = c("[]", "a"), x = NULL)
  ),
  "array elements keep sentinels while object fields still drop them"
)
ok(
  ut_cmp_identical(
    agent_args_sanitize(list(layers = list(list(geom = "point", filter = "null",
                                               color = "null")))),
    list(layers = list(list(geom = "point", filter = NULL, color = NULL)))
  ),
  "objects nested inside arrays still apply the sentinel rule per field"
)
ok(
  ut_cmp_identical(
    agent_args_sanitize(list(x = list(list(a = "null"), "null"))),
    list(x = list(list(a = NULL), "null"))
  ),
  "mixed arrays keep sentinel strings as elements, objects drop per field"
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
        received <- if (inherits(res, "error")) res else
          (if (!is.null(res@extra$data)) res@extra$data else res@value)
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

# ---- 4b. output codec at the seam (3.2) ----------------------------------
# Every non-string, non-error tool value is serialized ONCE at the seam:
# names preserved (ellmer's toJSON(auto_unbox) drops named-vector
# names), NULL as null (never the {} sentinel), size capped, and the
# structured original kept in extra$data for internal consumers.
codec_tool <- agent_tool(
  function(payload = NULL, `_intent`) ellmer::ContentToolResult(
    value = payload),
  description = "codec probe",
  name = "seam_codec",
  arguments = list(
    payload = ellmer::type_object(
      "Arbitrary structured payload to echo back.",
      dimensions = ellmer::type_array(
        ellmer::type_number("n"), "numbers", required = FALSE),
      note = ellmer::type_string("s", required = FALSE),
      .required = FALSE
    ),
    `_intent` = ellmer::type_string("intent")
  )
)
req <- ellmer::ContentToolRequest(
  id = "c3", name = "seam_codec",
  arguments = list(`_intent` = "x",
                   payload = list(dimensions = c(features = 2702, samples = 60),
                                  note = "demo")))
req@tool <- codec_tool
coded <- ellmer:::invoke_tool(req)
ok(
  inherits(coded@value, "json") &&
    grepl('"dimensions":{"features":2702', coded@value, fixed = TRUE) &&
    !grepl('"dimensions":[2702', coded@value, fixed = TRUE),
  "named vectors keep their names through the seam (wire: object, not [2702,60])"
)
ok(
  identical(as.character(coded@value),
            as.character(jsonlite::toJSON(
              list(dimensions = list(features = 2702, samples = 60), note = "demo"),
              auto_unbox = TRUE, null = "null"))),
  "the model-facing value is exactly our serialization"
)
ok(
  ut_cmp_identical(coded@extra$data$dimensions, c(features = 2702, samples = 60)),
  "extra$data keeps the structured original for internal consumers"
)
# tool_string (ellmer's wire projection) returns the json value verbatim
ok(
  identical(ellmer:::tool_string(coded), coded@value),
  "tool_string passes json-class values through unchanged"
)

# NULL value serializes as literal json null, never {}
req_null <- ellmer::ContentToolRequest(
  id = "c4", name = "seam_codec",
  arguments = list(`_intent` = "x", payload = NULL))
req_null@tool <- codec_tool
coded_null <- ellmer:::invoke_tool(req_null)
ok(
  identical(as.character(coded_null@value), "null") &&
    is.null(coded_null@extra$data),
  "NULL values become literal json null (no {} object sentinel)"
)

# oversized values hit the truncation marker; extra$data stays intact
big_payload <- list(ids = paste0("gene", 1:4000))
big_tool <- agent_tool(
  function(`_intent`) ellmer::ContentToolResult(value = big_payload),
  description = "oversized probe", name = "seam_big",
  arguments = list(`_intent` = ellmer::type_string("intent"))
)
req_big <- ellmer::ContentToolRequest(
  id = "c5", name = "seam_big", arguments = list(`_intent` = "x"))
req_big@tool <- big_tool
coded_big <- ellmer:::invoke_tool(req_big)
ok(
  grepl("Tool result truncated", coded_big@value, fixed = TRUE) &&
    grepl("Narrow the request", coded_big@value, fixed = TRUE) &&
    nchar(coded_big@value, type = "bytes") < 12000L,
  "oversized results are replaced by a bounded truncation marker"
)
ok(
  ut_cmp_identical(length(coded_big@extra$data$ids), 4000L),
  "the structured original survives truncation in extra$data"
)

# env knob clamps
ok(
  ut_cmp_identical(omicsViewer:::.agent_tool_output_max_bytes(), 8192L) &&
    ut_cmp_identical(
      withr_local_env(list(OMICSVIEWER_LLM_TOOL_OUTPUT_BYTES = "100000"),
        omicsViewer:::.agent_tool_output_max_bytes()), 65536L) &&
    ut_cmp_identical(
      withr_local_env(list(OMICSVIEWER_LLM_TOOL_OUTPUT_BYTES = "10"),
        omicsViewer:::.agent_tool_output_max_bytes()), 1024L),
  "tool-output byte cap default and clamp bounds"
)

# ---- 5. golden traffic replay: the live failing payloads ------------------

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
        layer2 = res@extra$data$spec$layers[[2]]$filter,
        layer3 = res@extra$data$spec$layers[[3]]$filter
      ) else NULL
      replay_warnings[[fx]] <<- if (is.null(res@error))
        (res@extra$data$warnings %||% character()) else "failed"
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

# ---- 5. scalar-argument arrays raise named correctable errors (1.7) ------
# A length>1 value for a scalar tool argument previously crashed the scalar
# helpers with "'length = N' in coercion to 'logical(1)'" -- a cryptic
# provider-facing error. The helpers must name the argument instead.
search_annotations_fn <- omicsViewer:::agent_search_annotations
ok(
  ut_cmp_error(
    search_annotations_fn(space = "feature", query = c("score", "cat"), fd, pd),
    "Expected a single string for query"
  ),
  "array payload for a scalar search query names the argument"
)
ok(
  ut_cmp_error(
    omicsViewer:::agent_summarize_annotation(
      space = "feature", column = c("score", "category"), fd, pd),
    "Expected a single string for column"
  ),
  "array payload for a scalar annotation column names the argument"
)
ok(
  ut_cmp_error(
    omicsViewer:::agent_normalize_scatter_view(
      space = "feature", x_axis = c("a|b|c", "d|e|f"), y_axis = "a|b|c",
      feature_columns = c("a|b|c")),
    "Expected a single string for x_axis"
  ),
  "array payload for a scatter axis names the argument"
)
ok(
  ut_cmp_error(
    omicsViewer:::agent_figure_template_spec(
      template = c("volcano", "scatter"), feature_data = fd, sample_data = pd),
    "Expected a single string for template"
  ),
  "array payload for a figure template names the argument"
)
ok(
  ut_cmp_error(
    omicsViewer:::agent_normalize_figure_spec(
      list(geom = c("point", "bar"), x = "score"), feature_data = fd,
      sample_data = pd),
    "geom"
  ),
  "figure geom arrays are rejected with the field named"
)
ok(
  ut_cmp_identical(
    omicsViewer:::.agent_trim_scalar(c("  x ", "y"), arg = NULL),
    "x"
  ),
  "unnamed trim-scalar callers keep the historic first-element behaviour"
)
