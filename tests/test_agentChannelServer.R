# Stage 2 server suite for the agent channel: the localhost JSON-RPC
# endpoint (R/auxi_agentChannelServer.R) driven at three layers -
#
#   [6]  the PURE handler (agent_channel_handle): initialize /
#        tools/list / tools/call happy paths against a Stage-1-style
#        constructed registry, every error code (parse / batch / id /
#        unknown method / no session / arg + tool errors / rate cap),
#        the I() array-unboxing obligation, the value@value-only wire
#        projection (F3), and token hygiene.
#   [7]  the bind hook (agent_channel_session_attach/detach): inert
#        when OMICSVIEWER_AGENT_CHANNEL is unset (incl. a full
#        ai_assistant_module boot via testServer), first-bind-wins,
#        onEnded unbind + last-unbind server stop.
#   [8]  loopback socket smoke (skipped without the `websocket`
#        package): real httpuv server + real WS client, wrong-token
#        refusal, initialize/tools/list/tools/call round trips, binary
#        frame refusal, token never on the wire.
#
# Run: Rscript tests/test_agentChannelServer.R   (from the repo root)
#
# Env-var discipline: the channel is NEVER enabled outside a
# withr-style local scope (plain Sys.setenv with on.exit restore).
# Registry construction reuses the proven tests/test_agentChannelLive.R
# patterns (decisions I12-I17); duplication is deliberate - the Stage-1
# file is frozen by the Stage-2 edit-scope rule.

library(unittest, quietly = TRUE)

if (!requireNamespace("pkgload", quietly = TRUE) ||
    !requireNamespace("ellmer", quietly = TRUE) ||
    !requireNamespace("shinychat", quietly = TRUE)) {
  ok(TRUE, "agent channel server tests skipped because optional packages are unavailable")
  quit(save = "no", status = 0)
}

pkg_root <- if (file.exists("DESCRIPTION")) "." else
  if (file.exists(file.path("..", "DESCRIPTION"))) ".." else "."
pkgload::load_all(pkg_root, quiet = TRUE, attach_testthat = FALSE)
library(ellmer)

for (nm in c(
  "widget_store_new", "widget_store_child", "widget_binding", "widget_action",
  "store_register", "store_register_action", "store_seed", "store_apply",
  "store_read", "store_snapshot", "agent_widget_apply",
  "agent_tool_registry", "agent_channel_dispatch", "agent_channel_toolspecs",
  "agent_channel_handle", "agent_channel_state_new",
  "agent_channel_server_start", "agent_channel_server_stop",
  "agent_channel_server_info", "agent_channel_session_attach",
  "agent_channel_session_detach", "agent_channel_registry_bind",
  "agent_channel_registry_unbind", "agent_channel_registry_list",
  "agent_normalize_state_update", "agent_normalize_scatter_view",
  "agent_normalize_enrichment_update", "agent_normalize_table_view_update",
  "agent_logger_new", "agent_logger_event", "agent_logger_set_enabled",
  "trisetter", "build_app_state", "agent_compact_state",
  "agent_annotation_catalog", "ai_assistant_module",
  ".agent_channel_wire_spec", ".agent_channel_token_tag",
  ".agent_channel_new_token"))
  assign(nm, get(nm, envir = asNamespace("omicsViewer")))

# ==========================================================================
# [1]-[4] Registry harness (duplicate of the Stage-1 constructor, trimmed
#         of its equivalence machinery - see tests/test_agentChannelLive.R
#         and STATUS decisions I12-I17)
# ==========================================================================

demo_path <- system.file("extdata", "demo.RDS", package = "omicsViewer")
if (!nzchar(demo_path))
  demo_path <- file.path(pkg_root, "inst", "extdata", "demo.RDS")
demo <- readRDS(demo_path)
fd <- Biobase::fData(demo)
pd <- Biobase::pData(demo)
mat <- Biobase::exprs(demo)

set.seed(1L)
fids <- rownames(fd)
gs_ids <- paste0("gs_set", c("A", "B", "C", "D"))
gs <- data.frame(
  featureId = factor(sample(fids, 300L, replace = TRUE), levels = fids),
  gsId = factor(sample(gs_ids, 300L, replace = TRUE), levels = gs_ids),
  weight = 1L)
attr(fd, "GS") <- gs

ftri <- trisetter(meta = fd, expr = mat, combine = "feature")
ptri <- trisetter(meta = pd, expr = mat, combine = "pheno")
quick_views <- list(feature = attr(fd, "quickViews"),
                    sample = attr(pd, "quickViews"))

DS_TABS <- c("Feature", "Feature table", "Sample", "Sample table", "Cor",
             "Heatmap", "Dynamic heatmap", "Expression", "GSList")
ANALYSIS_TABS <- c("Feature", "ORA", "fGSEA", "Geneshot", "Sample")

feature_data <- reactive(fd)
sample_data <- reactive(pd)
expression_data <- reactive(mat)
selected_features <- reactiveVal(character())
selected_samples <- reactiveVal(character())

make_store <- function() {
  store <- widget_store_new()
  sds <- widget_store_child(store, "dataspace")
  store_register(
    sds,
    widget_binding("active_tab", "navbar", label = "Data-space tab",
                   help = "Active tab of the data-space navbar",
                   values = DS_TABS))
  store_seed(sds, list(active_tab = "Feature"))

  tri_providers <- function(tri, prefix) {
    kx1 <- paste(prefix, "x_analysis", sep = ".")
    kx2 <- paste(prefix, "x_subset", sep = ".")
    ky1 <- paste(prefix, "y_analysis", sep = ".")
    ky2 <- paste(prefix, "y_subset", sep = ".")
    list(
      a1 = function(v) unique(tri[, 1]),
      a2 = function(v) unique(tri[tri[, 1] %in% v[[kx1]], 2]),
      a3 = function(v) tri[tri[, 1] %in% v[[kx1]] & tri[, 2] %in% v[[kx2]], 3],
      b1 = function(v) unique(tri[, 1]),
      b2 = function(v) unique(tri[tri[, 1] %in% v[[ky1]], 2]),
      b3 = function(v) tri[tri[, 1] %in% v[[ky1]] & tri[, 2] %in% v[[ky2]], 3])
  }
  bind_scatter <- function(child, providers) {
    store_register(
      child,
      widget_binding("x_analysis", "select", label = "X-axis analysis",
        help = "Analysis category for the X axis",
        choices_provider = providers$a1),
      widget_binding("x_subset", "select", label = "X-axis subset",
        help = "Subset within the X-axis analysis", depends_on = "x_analysis",
        choices_provider = providers$a2),
      widget_binding("x_variable", "select_cascaded", label = "X-axis variable",
        help = "Annotation variable plotted on the X axis",
        depends_on = c("x_analysis", "x_subset"),
        choices_provider = providers$a3),
      widget_binding("y_analysis", "select", label = "Y-axis analysis",
        help = "Analysis category for the Y axis",
        choices_provider = providers$b1),
      widget_binding("y_subset", "select", label = "Y-axis subset",
        help = "Subset within the Y-axis analysis", depends_on = "y_analysis",
        choices_provider = providers$b2),
      widget_binding("y_variable", "select_cascaded", label = "Y-axis variable",
        help = "Annotation variable plotted on the Y axis",
        depends_on = c("y_analysis", "y_subset"),
        choices_provider = providers$b3),
      widget_binding("axis_mode", "enum", label = "Axis mode",
        help = "Quick-view shortcuts vs the custom axis selector",
        values = c("quick", "custom")))
    split3 <- function(nm) strsplit(nm, "|", fixed = TRUE)[[1]]
    dx <- split3("PCA|All|PC1(10.5%)")
    dy <- split3("PCA|All|PC2(7.2%)")
    store_seed(child, list(
      x_analysis = dx[[1]], x_subset = dx[[2]], x_variable = dx[[3]],
      y_analysis = dy[[1]], y_subset = dy[[2]], y_variable = dy[[3]],
      axis_mode = "quick"))
  }
  bind_scatter(widget_store_child(store, "dataspace.feature_space"),
               tri_providers(ftri, "dataspace.feature_space"))
  bind_scatter(widget_store_child(store, "dataspace.sample_space"),
               tri_providers(ptri, "dataspace.sample_space"))

  srs <- widget_store_child(store, "resultspace")
  store_register(
    srs,
    widget_binding("analyst_tab", "navbar", label = "Analysis tab",
                   help = "Active tab of the analysis navbar",
                   choices_provider = function(v) ANALYSIS_TABS))
  store_seed(srs, list(analyst_tab = "Feature"))

  stf <- widget_store_child(store, "dataspace.tab_feature")
  col_choices <- function(v) colnames(fd)
  store_register(
    stf,
    widget_binding("multi_selection", "boolean", label = "Multiple selection",
      help = "Allow selecting more than one row"),
    widget_binding("columns", "multi_select", label = "Shown columns",
      help = "Table columns currently displayed", min = 1L,
      choices_provider = col_choices),
    widget_binding("page", "integer", label = "Table page",
      help = "1-based page currently visible", min = 1L),
    widget_binding("column_filters", "mapping", label = "Column filters",
      help = "Per-column search patterns",
      choices_provider = col_choices))
  store_register_action(
    stf,
    widget_action("clear", label = "Clear table selection",
      help = "Show all rows again",
      handler = function(params) invisible(TRUE),
      input_id = "tab_feature-clear"))
  store_seed(stf, list(multi_selection = FALSE,
                       columns = utils::head(colnames(fd), 5L)))

  shm <- widget_store_child(store, "dataspace.expr_heatmap")
  store_register(
    shm,
    widget_binding("heatmap_colors", "select", label = "Heatmap color panel",
      help = "Diverging palette for the heatmap cells",
      values = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy",
                 "RdYlBu", "RdYlGn")))
  store_seed(shm, list(heatmap_colors = "RdBu"))
  store
}

make_deps <- function() {
  store <- make_store()
  apply_state <- function(update) {
    validated <- agent_normalize_state_update(
      update = update, data_tabs = DS_TABS, analysis_tabs = ANALYSIS_TABS,
      feature_ids = rownames(isolate(feature_data())),
      sample_ids = rownames(isolate(sample_data())))
    tab_patch <- list()
    if (!is.null(validated$data_space_tab))
      tab_patch[["dataspace.active_tab"]] <- validated$data_space_tab
    if (!is.null(validated$analysis_space_tab))
      tab_patch[["resultspace.analyst_tab"]] <- validated$analysis_space_tab
    if (length(tab_patch))
      store_apply(store, tab_patch, origin = "agent", strict = TRUE)
    list(data_space_tab = validated$data_space_tab,
         analysis_space_tab = validated$analysis_space_tab,
         feature_count = length(validated$features %||% character()),
         sample_count = length(validated$samples %||% character()))
  }
  apply_scatter_view <- function(space, quick_view_id = NULL, x_axis = NULL,
                                 y_axis = NULL) {
    view <- agent_normalize_scatter_view(
      space = space, quick_view_id = quick_view_id, x_axis = x_axis,
      y_axis = y_axis, quick_views = quick_views,
      feature_columns = colnames(isolate(feature_data())),
      sample_columns = colnames(isolate(sample_data())))
    prefix <- if (view$space == "feature") "dataspace.feature_space" else
      "dataspace.sample_space"
    split_axis <- function(axis) strsplit(axis, "|", fixed = TRUE)[[1]]
    patch <- stats::setNames(
      c(as.list(split_axis(view$x_axis)), as.list(split_axis(view$y_axis))),
      c(paste(prefix, c("x_analysis", "x_subset", "x_variable"), sep = "."),
        paste(prefix, c("y_analysis", "y_subset", "y_variable"), sep = ".")))
    shiny::isolate(store_apply(store, patch, origin = "agent", strict = TRUE))
    view
  }
  apply_enrichment <- function(update) {
    validated <- agent_normalize_enrichment_update(update,
                                                   isolate(feature_data()))
    receipt <- shiny::isolate(store_apply(store, validated$patch,
                                          origin = "agent", strict = FALSE))
    list(method = validated$method, applied = receipt$applied)
  }
  apply_table_view <- function(update) {
    validated <- agent_normalize_table_view_update(update)
    receipt <- shiny::isolate(store_apply(store, validated$patch,
                                          origin = "agent", strict = FALSE))
    list(table = validated$table, applied = receipt$applied)
  }
  state <- function(sections = NULL) {
    shiny::isolate({
      full <- build_app_state(
        dataset = demo, dataset_id = "demo.RDS",
        data_status = list(eset_active_tab =
          store_read(store, "dataspace.active_tab")[[1]] %||% "Feature"),
        result_status = list(analyst_active_tab =
          store_read(store, "resultspace.analyst_tab")[[1]] %||% "Feature"),
        selected_features = selected_features(),
        selected_samples = selected_samples(),
        label = "AI assistant current state")
      wanted <- sections %||% character()
      agent_compact_state(
        state = full,
        annotations = if ("annotations" %in% wanted)
          agent_annotation_catalog(fd, pd),
        quick_views = quick_views,
        available_tabs = list(data_space = DS_TABS,
                              analysis_space = ANALYSIS_TABS),
        figure_grammar = if ("figure_grammar" %in% wanted)
          agent_figure_grammar(),
        sections = sections,
        store = store)
    })
  }
  list(store = store, apply_state = apply_state,
       apply_scatter_view = apply_scatter_view,
       apply_enrichment = apply_enrichment,
       apply_table_view = apply_table_view, state = state)
}

domain <- shiny:::MockShinySession$new()
logger <- agent_logger_new(domain)
figures <- reactiveVal(list())
figure_counter <- reactiveVal(0L)
figure_directory <- file.path(tempdir(), "agent-channel-server")
dir.create(figure_directory, showWarnings = FALSE, recursive = TRUE)

deps <- make_deps()
registry <- agent_tool_registry(
  state = deps$state,
  feature_data = feature_data,
  sample_data = sample_data,
  expression_data = expression_data,
  selected_features = selected_features,
  selected_samples = selected_samples,
  apply_state = deps$apply_state,
  apply_scatter_view = deps$apply_scatter_view,
  apply_enrichment = deps$apply_enrichment,
  apply_table_view = deps$apply_table_view,
  store = deps$store,
  figures = figures,
  figure_counter = figure_counter,
  figure_directory = figure_directory,
  session_domain = domain,
  logger = logger)

ok(
  ut_cmp_identical(length(registry$tools), 14L),
  "harness constructs the real 14-tool registry"
)

# ==========================================================================
# helpers
# ==========================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

# withr-style scoped env (plain Sys.setenv + on.exit INSIDE the wrapper,
# so vars live exactly for the duration of `expr`)
with_env <- function(vars, expr) {
  old <- Sys.getenv(names(vars), unset = NA_character_)
  var_names <- names(vars)
  do.call(Sys.setenv, as.list(vars))
  on.exit({
    for (i in seq_along(var_names)) {
      if (is.na(old[[i]]))
        Sys.unsetenv(var_names[[i]])
      else
        do.call(Sys.setenv, stats::setNames(list(old[[i]]), var_names[[i]]))
    }
  })
  eval(substitute(expr), envir = parent.frame())
}

stub_session <- function(token) {
  env <- new.env(parent = emptyenv())
  env$token <- token
  env$ended_callback <- NULL
  env$onEnded <- function(cb) env$ended_callback <- cb
  env
}

parse_response <- function(json)
  jsonlite::fromJSON(json, simplifyVector = FALSE)

err_code <- function(json) parse_response(json)$error$code

call_msg <- function(id, name, args = list(`_intent` = "server test")) {
  paste0("{\"jsonrpc\":\"2.0\",\"id\":", jsonlite::toJSON(id, auto_unbox = TRUE),
         ",\"method\":\"tools/call\",\"params\":{\"name\":\"", name,
         "\",\"args\":", jsonlite::toJSON(args, auto_unbox = TRUE,
                                          null = "null"), "}}")
}

req_msg <- function(id, method) {
  paste0("{\"jsonrpc\":\"2.0\",\"id\":", jsonlite::toJSON(id, auto_unbox = TRUE),
         ",\"method\":\"", method, "\"}")
}

find_free_port <- function() {
  repeat {
    candidate <- sample(20000:60000, 1L)
    free <- tryCatch({
      probe <- httpuv::startServer("127.0.0.1", candidate, list())
      httpuv::stopServer(probe)
      TRUE
    }, error = function(e) FALSE)
    if (free)
      return(candidate)
  }
}

# ==========================================================================
# [5] env knobs + token utilities
# ==========================================================================

ok(
  ut_cmp_identical(
    with_env(c(OMICSVIEWER_AGENT_CHANNEL = "1"),
      get(".agent_channel_enabled", envir = asNamespace("omicsViewer"))()), TRUE) &&
    ut_cmp_identical(
      with_env(c(OMICSVIEWER_AGENT_CHANNEL = ""),
        get(".agent_channel_enabled", envir = asNamespace("omicsViewer"))()), FALSE),
  "env gate: truthy values enable, empty disables"
)

ok(
  ut_cmp_identical(
    with_env(c(OMICSVIEWER_AGENT_CHANNEL = "0",
               OMICSVIEWER_AGENT_CHANNEL_PORT = "9123",
               OMICSVIEWER_AGENT_CHANNEL_MAX_CALLS = "0"),
      {
        ns <- asNamespace("omicsViewer")
        c(enabled = get(".agent_channel_enabled", envir = ns)(),
          port = get(".agent_channel_port_default", envir = ns)(),
          max_calls = get(".agent_channel_max_calls", envir = ns)())
      }),
    c(enabled = FALSE, port = 9123L, max_calls = 0L)),
  "env gate recognizes falsy words; port and max-calls knobs parse (0 = unlimited)"
)

tok <- .agent_channel_new_token()
tag <- .agent_channel_token_tag(tok)
ok(
  ut_cmp_identical(nchar(tok) >= 32L, TRUE) &&
    ut_cmp_identical(grepl("^[0-9a-f]+$", tok), TRUE) &&
    ut_cmp_identical(nchar(tag), 16L) &&
    ut_cmp_identical(tag, .agent_channel_token_tag(tok)) &&
    ut_cmp_identical(tag == .agent_channel_token_tag(.agent_channel_new_token()), FALSE),
  "tokens are >=32 hex chars; tags are stable 16-hex and differ per token"
)

# ==========================================================================
# [6] pure handler: happy paths against the real registry
# ==========================================================================

live_state <- agent_channel_state_new(tools = registry$tools, logger = logger,
                                      session_token = domain$token)

## 6.1 initialize -----------------------------------------------------------

resp <- agent_channel_handle(live_state, req_msg(1, "initialize"))
ok(ut_cmp_identical(err_code(resp), NULL), "initialize answers without error")
init <- parse_response(resp)$result
ok(
  ut_cmp_identical(init$protocol, "omicsviewer-agent-channel") &&
    ut_cmp_identical(init$jsonrpc, "2.0") &&
    ut_cmp_identical(init$channel_schema_version, 1L) &&
    ut_cmp_identical(init$tool_tier, "read") &&
    ut_cmp_identical(nzchar(init$package_version), TRUE),
  "initialize result carries protocol, jsonrpc, schema version (D8), read tier, package version"
)

## 6.2 tools/list: tier filter, versions echo, I() arrays -------------------

resp <- agent_channel_handle(live_state, req_msg(2, "tools/list"))
ok(ut_cmp_identical(err_code(resp), NULL), "tools/list answers without error")
listed <- parse_response(resp)$result
ok(
  ut_cmp_identical(
    vapply(listed$tools, function(spec) spec$name, character(1)),
    c("get_omics_viewer_state", "search_annotations", "summarize_annotation",
      "find_controls", "describe_control", "get_figure")),
  "tools/list exposes exactly the six read-tier tools in tier order"
)
ok(
  ut_cmp_identical(listed$channel_schema_version, 1L) &&
    ut_cmp_identical(nzchar(listed$package_version), TRUE),
  "tools/list echoes versions (D8)"
)
req_array <- parse_response(resp)$result$tools[[1]]$inputSchema$required
ok(
  ut_cmp_identical(is.list(req_array), TRUE) &&
    ut_cmp_identical(req_array[[1]], "_intent"),
  "length-1 required array survives the wire as an array (I() obligation)"
)
ok(
  ut_cmp_identical(
    any(vapply(listed$tools, function(spec)
      identical(spec$name, "set_widgets"), logical(1))), FALSE),
  "write-tier tools are absent from tools/list in Stage 2"
)

# direct unit probe of the I() wrapper: length-1 enum and required
wire <- .agent_channel_wire_spec(list(
  type = "object", required = "solo",
  properties = list(x = list(type = "string", enum = "only"))))
ok(
  grepl('"required":\\["solo"\\]', jsonlite::toJSON(wire, auto_unbox = TRUE), perl = TRUE) &&
    grepl('"enum":\\["only"\\]', jsonlite::toJSON(wire, auto_unbox = TRUE), perl = TRUE),
  "wire-spec wrapper keeps length-1 enum/required as JSON arrays (auto_unbox edge)"
)

## 6.3 tools/call happy path + value@value-only projection (F3) -------------

resp <- agent_channel_handle(live_state, call_msg(3, "search_annotations",
  list(space = "feature", query = "Gene.name", max_results = 5L,
       `_intent` = "find the name column")))
ok(ut_cmp_identical(err_code(resp), NULL), "tools/call search_annotations answers without error")
content <- parse_response(resp)$result$content
ok(
  ut_cmp_identical(length(content), 1L) &&
    ut_cmp_identical(content[[1]]$type, "text") &&
    ut_cmp_identical(is.character(content[[1]]$text), TRUE),
  "tools/call result wraps one text content block"
)
direct <- agent_channel_dispatch(registry$tools, "search_annotations",
                                 list(space = "feature", query = "Gene.name",
                                      max_results = 5L,
                                      `_intent` = "find the name column"))
ok(
  ut_cmp_identical(content[[1]]$text, as.character(direct$value@value)),
  "wire text is exactly value@value (byte-identical to the direct dispatch)"
)
ok(
  ut_cmp_identical(grepl("display", resp, fixed = TRUE), FALSE) &&
    ut_cmp_identical(grepl("value_preview", resp, fixed = TRUE), FALSE) &&
    ut_cmp_identical(grepl("extra", resp, fixed = TRUE), FALSE),
  "extra$display / extra$data payload never crosses the wire (F3 projection)"
)
ok(
  ut_cmp_identical(
    is.list(jsonlite::fromJSON(content[[1]]$text, simplifyVector = FALSE)),
    TRUE),
  "the content text itself is the codec's bounded JSON serialization"
)

resp <- agent_channel_handle(live_state, call_msg(4, "get_omics_viewer_state",
  list(`_intent` = "overview")))
ok(
  ut_cmp_identical(err_code(resp), NULL) &&
    grepl("scatter_view", resp, fixed = TRUE),
  "tools/call get_omics_viewer_state returns the compact state overview"
)

resp <- agent_channel_handle(live_state, call_msg(5, "describe_control",
  list(id = "dataspace.feature_space.x_variable", `_intent` = "x axis")))
ok(
  ut_cmp_identical(err_code(resp), NULL) &&
    grepl("select_cascaded", resp, fixed = TRUE),
  "tools/call describe_control returns the live widget record"
)

# arrays as JSON args (features list) survive the JSON-decode + sanitizer
resp <- agent_channel_handle(live_state, call_msg(6, "summarize_annotation",
  list(space = "sample", column = "General|All|Origin", `_intent` = "origin groups")))
ok(
  ut_cmp_identical(err_code(resp), NULL) &&
    grepl('"non_missing":60',
          parse_response(resp)$result$content[[1]]$text, fixed = TRUE),
  "tools/call summarize_annotation summarizes a real demo column (60 observed)"
)

## 6.4 error codes -----------------------------------------------------------

ok(
  ut_cmp_identical(err_code(agent_channel_handle(live_state, "{oops")), -32700L) &&
    ut_cmp_identical(err_code(agent_channel_handle(live_state, "  ")), -32700L),
  "parse errors answer -32700"
)
ok(
  ut_cmp_identical(err_code(agent_channel_handle(live_state, "5")), -32600L) &&
    ut_cmp_identical(
      err_code(agent_channel_handle(live_state, req_msg(7, "initialize") |> (\(m) paste0("[", m, "]"))())),
      -32600L) &&
    ut_cmp_identical(
      err_code(agent_channel_handle(live_state,
        gsub("\"id\":7,", "", req_msg(7, "initialize"), fixed = TRUE))), -32600L),
  "non-object messages, batch arrays, and missing id answer -32600"
)
ok(
  ut_cmp_identical(
    err_code(agent_channel_handle(live_state, req_msg(8, "no/such/method"))),
    -32601L),
  "unknown method answers -32601"
)
ok(
  ut_cmp_identical(
    err_code(agent_channel_handle(live_state,
      call_msg(9, "set_widgets", list(patch = "{}")))), -32602L),
  "write-tier tools are refused on the read tier (-32602)"
)
resp <- agent_channel_handle(live_state, call_msg(10, "search_annotations",
  list(spa = "feature", `_intent` = "typo")))
ok(
  ut_cmp_identical(err_code(resp), -32602L) &&
    grepl("Unknown argument 'spa'", resp, fixed = TRUE),
  "exact-name argument validation surfaces as -32602 with the offending arg (I9)"
)
ok(
  ut_cmp_identical(
    err_code(agent_channel_handle(live_state, call_msg(11, "no_such_tool"))),
    -32602L),
  "unknown tool answers -32602 listing available tools"
)
resp <- agent_channel_handle(live_state,
  paste0("{\"jsonrpc\":\"2.0\",\"id\":12,\"method\":\"tools/call\",",
         "\"params\":{\"name\":\"find_controls\",\"args\":[1,2]}}"))
ok(
  ut_cmp_identical(err_code(resp), -32602L),
  "array-shaped args answer -32602 (params.args must be an object)"
)

# omitted / empty `args` is NOT a params error (S2-review F1): the call
# must reach the DISPATCHER, exactly like an explicit argument-less
# call. NB R's lazy missing-argument semantics: no tool handler ever
# forces `_intent` in its body, so a dispatcher-layer rejection names
# the first REQUIRED-AND-USED formal (space, query, id, ...) - while
# tools whose only required arg is `_intent` outright SUCCEED with an
# argument-less call. Both outcomes are dispatcher behavior; what F1
# pins is that neither ever surfaces the args-object guard error.
resp <- agent_channel_handle(live_state,
  paste0("{\"jsonrpc\":\"2.0\",\"id\":40,\"method\":\"tools/call\",",
         "\"params\":{\"name\":\"search_annotations\"}}"))
ok(
  ut_cmp_identical(err_code(resp), -32602L) &&
    grepl("is missing, with no default", resp, fixed = TRUE) &&
    !grepl("must be a JSON object", resp, fixed = TRUE),
  "tools/call with OMITTED args reaches the dispatcher (required-arg missing error), not the args-object guard (F1)"
)
resp_omitted <- agent_channel_handle(live_state,
  paste0("{\"jsonrpc\":\"2.0\",\"id\":41,\"method\":\"tools/call\",",
         "\"params\":{\"name\":\"find_controls\"}}"))
resp_empty <- agent_channel_handle(live_state,
  paste0("{\"jsonrpc\":\"2.0\",\"id\":42,\"method\":\"tools/call\",",
         "\"params\":{\"name\":\"find_controls\",\"args\":{}}}"))
ok(
  ut_cmp_identical(err_code(resp_omitted), NULL) &&
    ut_cmp_identical(err_code(resp_empty), NULL) &&
    ut_cmp_identical(
      parse_response(resp_omitted)$result$content[[1]]$text,
      parse_response(resp_empty)$result$content[[1]]$text),
  "omitted and explicit {} args dispatch the same argument-less call: find_controls answers both (F1)"
)
ok(
  ut_cmp_identical(
    err_code(agent_channel_handle(live_state,
      paste0("{\"jsonrpc\":\"2.0\",\"id\":43,\"method\":\"tools/call\",",
             "\"params\":{\"name\":\"find_controls\",",
             "\"args\":{\"_intent\":\"list the controls\"}}}"))), NULL),
  "omitted-args semantics stay intact: args carrying only _intent succeeds"
)
resp <- agent_channel_handle(live_state,
  paste0("{\"jsonrpc\":\"2.0\",\"id\":44,\"method\":\"tools/call\",",
         "\"params\":{}}"))
ok(
  ut_cmp_identical(err_code(resp), -32602L) &&
    grepl("`name` must be a single non-empty tool name", resp, fixed = TRUE) &&
    !grepl("needs a params object", resp, fixed = TRUE),
  "empty params object falls through to the name check (symmetric empty-object semantics)"
)

## 6.5 no-session states ------------------------------------------------------

no_session <- agent_channel_state_new()
ok(
  ut_cmp_identical(err_code(agent_channel_handle(no_session, req_msg(13, "initialize"))), NULL),
  "initialize works without a bound session"
)
ok(
  ut_cmp_identical(
    err_code(agent_channel_handle(no_session, req_msg(14, "tools/list"))), -32000L) &&
    ut_cmp_identical(
      err_code(agent_channel_handle(no_session, call_msg(15, "find_controls"))),
      -32000L),
  "tools/list and tools/call answer -32000 with no session bound"
)

## 6.6 rate cap (D6) ----------------------------------------------------------

capped <- agent_channel_state_new(tools = registry$tools, max_calls = 2L)
ok(
  ut_cmp_identical(err_code(agent_channel_handle(capped, call_msg(21, "find_controls"))), NULL) &&
    ut_cmp_identical(err_code(agent_channel_handle(capped, call_msg(22, "find_controls"))), NULL) &&
    ut_cmp_identical(capped$calls, 2L),
  "the per-connection counter tracks executed tools/call only"
)
resp <- agent_channel_handle(capped, call_msg(23, "find_controls"))
ok(
  ut_cmp_identical(err_code(resp), -32000L) &&
    grepl("Rate cap exceeded", resp, fixed = TRUE) &&
    grepl("2 tool calls", resp, fixed = TRUE),
  "exceeding the cap answers -32000 naming the cap"
)
ok(
  ut_cmp_identical(err_code(agent_channel_handle(capped, req_msg(24, "initialize"))), NULL) &&
    ut_cmp_identical(capped$calls, 2L),
  "the connection stays usable after a cap hit (non-tools/call methods keep answering)"
)
unlimited <- agent_channel_state_new(tools = registry$tools, max_calls = 0L)
ok(
  ut_cmp_identical(
    local({
      okays <- vapply(25:27, function(i)
        is.null(err_code(agent_channel_handle(
          unlimited, call_msg(i, "find_controls")))), logical(1))
      all(okays)
    }), TRUE),
  "max_calls = 0 means unlimited"
)

## 6.7 logger hygiene: tags, never tokens -------------------------------------

log_dir <- file.path(tempdir(), "agent-channel-server-logs")
unlink(log_dir, recursive = TRUE, force = TRUE)
log_domain <- shiny:::MockShinySession$new()
log_logger <- agent_logger_new(
  log_domain,
  list(enabled = TRUE, directory = log_dir, max_bytes = 10485760L))
agent_logger_set_enabled(log_logger, TRUE)
logged_state <- agent_channel_state_new(tools = registry$tools,
                                        logger = log_logger,
                                        session_token = log_domain$token)
invisible(agent_channel_handle(logged_state, req_msg(31, "initialize")))
invisible(agent_channel_handle(logged_state, call_msg(32, "find_controls")))
log_lines <- readLines(log_logger$path)
ok(
  ut_cmp_identical(
    any(grepl("channel_request_start", log_lines, fixed = TRUE)) &&
      any(grepl("channel_request_end", log_lines, fixed = TRUE)), TRUE),
  "request start/end events flow through the session logger (D5)"
)
ok(
  ut_cmp_identical(
    any(grepl(log_domain$token, log_lines, fixed = TRUE)), FALSE),
  "no log line carries the full session token (the session_id label uses only a 12-char prefix by the logger's own design)"
)

# attach/detach events through an ENABLED logger: tag present, token absent
with_env(c(OMICSVIEWER_AGENT_CHANNEL = "1",
           OMICSVIEWER_AGENT_CHANNEL_PORT = as.character(find_free_port())), {
  s <- stub_session("log-hygiene-session-token")
  ok(
    ut_cmp_identical(
      agent_channel_session_attach(
        session = s, tools = registry$tools, session_domain = domain,
        logger = log_logger), TRUE),
    "attach binds with an enabled logger (channel_bind + channel_server_start events)"
  )
  agent_channel_session_detach("log-hygiene-session-token")
})
log_lines <- readLines(log_logger$path)
ok(
  ut_cmp_identical(
    any(grepl("channel_bind", log_lines, fixed = TRUE)) &&
      any(grepl("channel_server_start", log_lines, fixed = TRUE)) &&
      any(grepl(.agent_channel_token_tag("log-hygiene-session-token"),
                log_lines, fixed = TRUE)), TRUE) &&
    ut_cmp_identical(
      any(grepl("log-hygiene-session-token", log_lines, fixed = TRUE)), FALSE),
  "bind/start events carry the token tag, never the session token"
)
server_token_line <- grep("channel_server_start", log_lines, fixed = TRUE, value = TRUE)
server_token <- agent_channel_server_info()$token
ok(
  ut_cmp_identical(
    if (length(server_token_line) && !is.null(server_token))
      !grepl(server_token, server_token_line[[1]], fixed = TRUE)
    else TRUE, TRUE),
  "the pairing token appears only via its tag in channel_server_start events"
)

# ==========================================================================
# [7] bind hook: attach/detach lifecycle
# ==========================================================================

# 7.1 inert when the env is unset -------------------------------------------

# ensure a clean slate regardless of earlier sections
for (existing_token in agent_channel_registry_list())
  agent_channel_registry_unbind(existing_token)
agent_channel_server_stop()

ok(
  ut_cmp_identical(
    with_env(c(OMICSVIEWER_AGENT_CHANNEL = ""),
      agent_channel_session_attach(
        session = stub_session("stub-session-token-1"),
        tools = registry$tools, session_domain = domain, logger = logger)), FALSE) &&
    ut_cmp_identical(length(agent_channel_registry_list()), 0L) &&
    ut_cmp_identical(agent_channel_server_info()$running, FALSE),
  "attach is inert with the env unset: nothing binds, no server starts"
)

# 7.2 attach, first-bind-wins, detach ----------------------------------------

free_port <- find_free_port()
bound_token <- NULL
ok(
  ut_cmp_identical(
    with_env(c(OMICSVIEWER_AGENT_CHANNEL = "1",
               OMICSVIEWER_AGENT_CHANNEL_PORT = as.character(free_port)),
      agent_channel_session_attach(
        session = stub_session("stub-session-token-A"),
        tools = registry$tools, session_domain = domain, logger = logger)),
    TRUE),
  "attach binds the first session when the channel is enabled"
)
bound_token <- if (length(agent_channel_registry_list()))
  agent_channel_registry_list()[[1]] else NULL
info <- agent_channel_server_info()
ok(
  ut_cmp_identical(bound_token, "stub-session-token-A") &&
    ut_cmp_identical(info$running, TRUE) &&
    ut_cmp_identical(info$port, free_port) &&
    ut_cmp_identical(nchar(info$token) >= 32L, TRUE) &&
    ut_cmp_identical(info$tier, "read"),
  "the server started on the requested port with a fresh pairing token"
)
ok(
  ut_cmp_identical(
    with_env(c(OMICSVIEWER_AGENT_CHANNEL = "1"),
      agent_channel_session_attach(
        session = stub_session("stub-session-token-B"),
        tools = registry$tools, session_domain = domain, logger = logger)),
    FALSE) &&
    ut_cmp_identical(agent_channel_registry_list(), "stub-session-token-A"),
  "second attach warns and does NOT rebind (first-bind-wins)"
)
ok(
  ut_cmp_identical(.agent_channel_token_tag("stub-session-token-A") !=
                     .agent_channel_token_tag("stub-session-token-B"), TRUE),
  "tag sanity for the refusal warning path"
)
# fire the wired onEnded callback: unbind + last-unbind server stop
agent_channel_session_detach(bound_token)
ok(
  ut_cmp_identical(length(agent_channel_registry_list()), 0L) &&
    ut_cmp_identical(agent_channel_server_info()$running, FALSE),
  "onEnded unbind evicts the session and stops the server with the last unbind"
)

# 7.3 full module boot is inert with the env unset ---------------------------

module_boot_bound <- NULL
shiny::testServer(
  function(input, output, session) {
    ai_assistant_module(
      id = "ai", state = deps$state,
      state_available = reactive(TRUE),
      feature_data = feature_data, sample_data = sample_data,
      expression_data = expression_data,
      selected_features = selected_features,
      selected_samples = selected_samples,
      apply_state = deps$apply_state,
      apply_scatter_view = deps$apply_scatter_view,
      apply_enrichment = deps$apply_enrichment,
      apply_table_view = deps$apply_table_view,
      store = deps$store)
  },
  expr = {
    # force a flush so moduleServer bodies (registry + attach hook) run
    session$flushReact()
    module_boot_bound <<- agent_channel_registry_list()
  })
ok(
  ut_cmp_identical(length(module_boot_bound), 0L) &&
    ut_cmp_identical(agent_channel_server_info()$running, FALSE),
  "ai_assistant_module boots with the channel unset: nothing binds (inert hook)"
)

# 7.4 attach reports server-start failure (F2): port pre-occupied -----

startfail_dir <- file.path(tempdir(), "agent-channel-startfail-logs")
unlink(startfail_dir, recursive = TRUE, force = TRUE)
startfail_domain <- shiny:::MockShinySession$new()
startfail_logger <- agent_logger_new(
  startfail_domain,
  list(enabled = TRUE, directory = startfail_dir, max_bytes = 10485760L))
agent_logger_set_enabled(startfail_logger, TRUE)

attach_with_busy_port <- function(logger) {
  probe_port <- find_free_port()
  probe <- httpuv::startServer("127.0.0.1", probe_port, list())
  on.exit(try(httpuv::stopServer(probe), silent = TRUE), add = TRUE)
  with_env(c(OMICSVIEWER_AGENT_CHANNEL = "1",
             OMICSVIEWER_AGENT_CHANNEL_PORT = as.character(probe_port)),
    agent_channel_session_attach(
      session = stub_session("stub-session-token-F2"),
      tools = registry$tools, session_domain = domain, logger = logger))
}
attach_start_failure <- attach_with_busy_port(startfail_logger)
startfail_lines <- readLines(startfail_logger$path)
ok(
  ut_cmp_identical(attach_start_failure, FALSE) &&
    ut_cmp_identical(agent_channel_registry_list(), "stub-session-token-F2") &&
    ut_cmp_identical(agent_channel_server_info()$running, FALSE),
  "attach returns FALSE when the port is taken; the session stays bound, no server (F2)"
)
ok(
  ut_cmp_identical(
    any(grepl("channel_server_start_failed", startfail_lines, fixed = TRUE)), TRUE) &&
    ut_cmp_identical(
      any(grepl("channel_bind", startfail_lines, fixed = TRUE)), FALSE),
  "a failed start logs channel_server_start_failed and NO channel_bind success event (F2)"
)
ok(
  ut_cmp_identical(
    any(grepl(.agent_channel_token_tag("stub-session-token-F2"),
              startfail_lines, fixed = TRUE)), TRUE) &&
    ut_cmp_identical(
      any(grepl("stub-session-token-F2", startfail_lines, fixed = TRUE)), FALSE),
  "the failure event carries the 16-hex tag, never the session token"
)
agent_channel_session_detach("stub-session-token-F2")

# ==========================================================================
# [8] loopback socket smoke (real httpuv server + real WS client)
# ==========================================================================

if (!requireNamespace("websocket", quietly = TRUE)) {
  ok(TRUE, "socket smoke skipped: the websocket package is not installed")
} else {
  with_env(c(OMICSVIEWER_AGENT_CHANNEL = "1",
             OMICSVIEWER_AGENT_CHANNEL_PORT = as.character(smoke_port <- find_free_port())), {
  for (existing_token in agent_channel_registry_list())
    agent_channel_registry_unbind(existing_token)
  agent_channel_server_stop()

  agent_channel_registry_bind(
    "smoke-session-token",
    list(tools = registry$tools, session_domain = domain, logger = logger))
  ok(
    ut_cmp_identical(agent_channel_server_start(), TRUE),
    "server starts on a random free port with a bound session"
  )
  token <- agent_channel_server_info()$token

  wait_for <- function(probe, timeout = 10) {
    deadline <- Sys.time() + timeout
    while (!isTRUE(probe()) && Sys.time() < deadline)
      later::run_now(0.05)
    isTRUE(probe())
  }

  make_client <- function(url) {
    cli <- new.env(parent = emptyenv())
    cli$open <- FALSE
    cli$closed <- FALSE
    cli$messages <- character()
    ws <- websocket::WebSocket$new(url)
    ws$onOpen(function(ev) cli$open <- TRUE)
    ws$onClose(function(ev) cli$closed <- TRUE)
    ws$onMessage(function(ev) cli$messages <- c(cli$messages, ev$data))
    cli$ws <- ws
    cli
  }
  rpc <- function(cli, msg) {
    n <- length(cli$messages)
    cli$ws$send(msg)
    ok(wait_for(function() length(cli$messages) > n),
       paste0("socket round trip answered: ", substr(msg, 1, 40), "..."))
    cli$messages[length(cli$messages)]
  }

  cli <- make_client(paste0("ws://127.0.0.1:", smoke_port, "/", token))
  ok(wait_for(function() cli$open), "client opens the token-paired connection")

  resp <- rpc(cli, req_msg(1, "initialize"))
  ok(
    ut_cmp_identical(err_code(resp), NULL) &&
      ut_cmp_identical(parse_response(resp)$result$protocol,
                       "omicsviewer-agent-channel"),
    "socket initialize answers the protocol result"
  )
  resp <- rpc(cli, req_msg(2, "tools/list"))
  ok(
    ut_cmp_identical(err_code(resp), NULL) &&
      ut_cmp_identical(length(parse_response(resp)$result$tools), 6L),
    "socket tools/list returns the six read-tier tools"
  )
  ok(
    ut_cmp_identical(
      is.list(parse_response(resp)$result$tools[[1]]$inputSchema$required), TRUE),
    "socket tools/list keeps length-1 required arrays as arrays"
  )
  resp <- rpc(cli, call_msg(3, "search_annotations",
                            list(space = "feature", query = "Origin")))
  ok(
    ut_cmp_identical(err_code(resp), NULL) &&
      ut_cmp_identical(
        parse_response(resp)$result$content[[1]]$type, "text"),
    "socket tools/call search_annotations returns a text content result"
  )

  # binary frame -> -32700 (connection stays)
  cli$ws$send(charToRaw("not json"))
  ok(
    wait_for(function() length(cli$messages) >= 4L) &&
      ut_cmp_identical(err_code(cli$messages[[4]]), -32700L),
    "binary frames answer -32700 (text-frame protocol)"
  )
  resp <- rpc(cli, req_msg(5, "initialize"))
  ok(ut_cmp_identical(err_code(resp), NULL),
     "the connection still answers after a binary frame")

  # token hygiene: the full token never crosses the wire in either direction
  ok(
    ut_cmp_identical(
      any(grepl(token, cli$messages, fixed = TRUE)), FALSE),
    "no response ever contains the pairing token"
  )

  # wrong-token connection: no messages, connection terminated (the client
  # onOpen may race the server-side close, so the open flag is not asserted)
  bad <- make_client(paste0("ws://127.0.0.1:", smoke_port, "/wrongtoken"))
  ok(
    wait_for(function() bad$closed || bad$ws$readyState() %in% c(2L, 3L), 5) &&
      ut_cmp_identical(length(bad$messages), 0L),
    "wrong-token connection is refused: no messages, connection closed"
  )

  # cleanup
  try(cli$ws$close(), silent = TRUE)
  agent_channel_registry_unbind("smoke-session-token")
  agent_channel_server_stop()
  ok(
    ut_cmp_identical(agent_channel_server_info()$running, FALSE),
    "socket cleanup stopped the server"
  )
  })
}

# final hygiene: nothing bound, no server, no env leakage
for (existing_token in agent_channel_registry_list())
  agent_channel_registry_unbind(existing_token)
agent_channel_server_stop()
ok(
  ut_cmp_identical(length(agent_channel_registry_list()), 0L) &&
    ut_cmp_identical(agent_channel_server_info()$running, FALSE) &&
    ut_cmp_identical(Sys.getenv("OMICSVIEWER_AGENT_CHANNEL", unset = ""), ""),
  "suite leaves no bindings, no server, and no env leakage"
)
