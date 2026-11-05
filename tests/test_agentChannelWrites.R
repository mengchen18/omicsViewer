# Stage 3 server suite for the agent channel: the WRITE tier, the
# figure GET endpoint, and channel-tagged write logging
# (R/auxi_agentChannelServer.R) driven at three layers -
#
#   [W2] env knobs: OMICSVIEWER_AGENT_CHANNEL_WRITES truthiness; the
#        tier is decided at SERVER START and not re-read while running.
#   [W3] the PURE handler (agent_channel_handle) with tier
#        "read-write": initialize/tools/list/tools/call, read-tier
#        refusal of writes, rate cap counting writes, the bounded
#        write-receipt helpers, and run_control's settled-policy
#        receipt.
#   [W4] write logging (D5): channel_write_start / channel_write_end
#        with 16-hex session tag, duration, and a receipt summary that
#        carries counts + at most 5 key NAMES and never values.
#   [W5] attach hook: figures / figure_directory / figure_renderer
#        land in the registry entry; channel_bind logs the live tier.
#   [W6] loopback smoke, real server + real clients, each section
#        skipping gracefully when its client package is absent:
#        [W6-WS] write round-trips through the socket with STORE-state
#        assertions (store_snapshot) and create/get/update figure
#        round-trip (websocket section); [W6-HTTP] the figure GET
#        endpoint over plain HTTP - 200 PNG bytes, wrong token refused,
#        traversal refused, unknown ids 404 (curl section, hoisted OUT
#        of the websocket conditional - O2 - so figure-GET coverage
#        survives on machines with curl but no websocket; the literal
#        ../ traversal guarantee is asserted at the PARSER level - O3,
#        libcurl normalizes dot segments client-side).
#
# NEW FILE (not an extension of test_agentChannelServer.R): the
# Stage-2 suite stays a stable read-tier guard (its 68 assertions pin
# the Stage-2 wire contract); write-tier concerns get their own
# focused suite over the same harness pattern (the Stage-1 -> Stage-2
# duplication precedent, decisions I12-I17).
#
# HTTP client note: the client and the httpuv server share THIS
# process's single R event loop - a BLOCKING fetch (curl_fetch_memory
# / download.file) would deadlock the server, so the suite uses the
# curl package's ASYNC interface (curl_fetch_multi on a private pool
# + later::run_now pump), verified against a live httpuv server
# before this suite was written. Without the curl package the HTTP
# section skips gracefully.
#
# Run: Rscript tests/test_agentChannelWrites.R   (from the repo root)
#
# Env-var discipline: the channel is NEVER enabled outside a
# withr-style local scope (plain Sys.setenv with on.exit restore).

library(unittest, quietly = TRUE)

if (!requireNamespace("pkgload", quietly = TRUE) ||
    !requireNamespace("ellmer", quietly = TRUE) ||
    !requireNamespace("shinychat", quietly = TRUE)) {
  ok(TRUE, "agent channel write tests skipped because optional packages are unavailable")
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
  "agent_normalize_figure_spec", "agent_build_figure_data",
  "agent_build_figure_plot", "agent_render_figure",
  ".agent_figure_selection",
  ".agent_channel_wire_spec", ".agent_channel_token_tag",
  ".agent_channel_new_token", ".agent_channel_writes_enabled",
  ".agent_channel_write_tools", ".agent_channel_write_receipt",
  ".agent_channel_write_error", ".agent_channel_figure_parse",
  ".agent_channel_uri_bytes", ".agent_channel_figure_resolve"))
  assign(nm, get(nm, envir = asNamespace("omicsViewer")))

# ==========================================================================
# [W1] Registry harness (duplicate of the Stage-2 constructor; see
#      tests/test_agentChannelServer.R and STATUS decisions I12-I17)
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
figure_directory <- file.path(tempdir(), "agent-channel-writes")
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

# The figure renderer mirrors the module hook's Stage-3 closure
# (R/module_aiAssistant.R): re-render a stored normalized spec through
# the standard pipeline WITHOUT touching the figure registry/counter.
figure_renderer <- function(spec) {
  normalized <- agent_normalize_figure_spec(
    spec = spec,
    feature_data = isolate(feature_data()),
    sample_data = isolate(sample_data()),
    expression = isolate(expression_data()),
    selected_features = isolate(.agent_figure_selection(selected_features())),
    selected_samples = isolate(.agent_figure_selection(selected_samples())))
  figure_data <- agent_build_figure_data(
    spec = normalized,
    feature_data = isolate(feature_data()),
    sample_data = isolate(sample_data()),
    expression = isolate(expression_data()))
  agent_render_figure(
    agent_build_figure_plot(figure_data, normalized),
    directory = figure_directory,
    figure_id = "channel")
}

# ==========================================================================
# helpers
# ==========================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

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

stub_session <- function(token) {
  env <- new.env(parent = emptyenv())
  env$token <- token
  env$ended_callback <- NULL
  env$onEnded <- function(cb) env$ended_callback <- cb
  env
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

clean_slate <- function() {
  for (existing_token in agent_channel_registry_list())
    agent_channel_registry_unbind(existing_token)
  agent_channel_server_stop()
}

# ==========================================================================
# [W2] write-tier env knobs; tier is a server-start property
# ==========================================================================

ok(
  ut_cmp_identical(
    with_env(c(OMICSVIEWER_AGENT_CHANNEL_WRITES = "1"),
      .agent_channel_writes_enabled()), TRUE) &&
    ut_cmp_identical(
      with_env(c(OMICSVIEWER_AGENT_CHANNEL_WRITES = ""),
        .agent_channel_writes_enabled()), FALSE) &&
    ut_cmp_identical(
      with_env(c(OMICSVIEWER_AGENT_CHANNEL_WRITES = "off"),
        .agent_channel_writes_enabled()), FALSE),
  "writes gate: truthy values enable, empty/falsy words disable (off by default)"
)
ok(
  ut_cmp_identical(
    sort(.agent_channel_write_tools()),
    sort(c("set_omics_viewer_state", "set_scatter_view", "set_widgets",
           "run_control", "set_enrichment_parameters", "set_table_view",
           "create_figure", "update_figure"))),
  "the write-only tool set is exactly the eight write tools"
)

clean_slate()
ok(
  ut_cmp_identical(
    with_env(c(OMICSVIEWER_AGENT_CHANNEL = "1",
               OMICSVIEWER_AGENT_CHANNEL_PORT = as.character(find_free_port())),
      {
        started <- agent_channel_server_start()
        tier_no_writes <- agent_channel_server_info()$tier
        agent_channel_server_stop()
        c(started = started, tier_no_writes = tier_no_writes)
      }),
    c(started = TRUE, tier_no_writes = "read")),
  "server starts on the READ tier when OMICSVIEWER_AGENT_CHANNEL_WRITES is unset"
)
ok(
  ut_cmp_identical(
    with_env(c(OMICSVIEWER_AGENT_CHANNEL = "1",
               OMICSVIEWER_AGENT_CHANNEL_WRITES = "1",
               OMICSVIEWER_AGENT_CHANNEL_PORT = as.character(find_free_port())),
      {
        started <- agent_channel_server_start()
        tier_writes <- agent_channel_server_info()$tier
        # tier is a START property: clearing the env while running must
        # NOT downgrade the running server
        Sys.setenv(OMICSVIEWER_AGENT_CHANNEL_WRITES = "")
        tier_after_unset <- agent_channel_server_info()$tier
        agent_channel_server_stop()
        c(started = started, tier_writes = tier_writes,
          tier_after_unset = tier_after_unset)
      }),
    c(started = TRUE, tier_writes = "read-write",
      tier_after_unset = "read-write")),
  "server starts on the READ-WRITE tier when the writes env is truthy; the tier is decided at start and not re-read while running"
)
ok(
  ut_cmp_identical(Sys.getenv("OMICSVIEWER_AGENT_CHANNEL", unset = ""), "") &&
    ut_cmp_identical(
      Sys.getenv("OMICSVIEWER_AGENT_CHANNEL_WRITES", unset = ""), "") &&
    ut_cmp_identical(agent_channel_server_info()$running, FALSE),
  "with_env scopes restored both env vars and stopped the server"
)

# ==========================================================================
# [W3] pure handler on the read-write tier
# ==========================================================================

rw_state <- agent_channel_state_new(tools = registry$tools, logger = logger,
                                    session_token = domain$token,
                                    tier = "read-write")

## W3.1 initialize carries the live tier -----------------------------------

resp <- agent_channel_handle(rw_state, req_msg(1, "initialize"))
ok(
  ut_cmp_identical(err_code(resp), NULL) &&
    ut_cmp_identical(parse_response(resp)$result$tool_tier, "read-write"),
  "initialize reports tool_tier 'read-write' when the state serves writes"
)

## W3.2 tools/list exposes all 14 tools in tier order ----------------------

resp <- agent_channel_handle(rw_state, req_msg(2, "tools/list"))
listed <- parse_response(resp)$result
ok(
  ut_cmp_identical(
    vapply(listed$tools, function(spec) spec$name, character(1)),
    c("get_omics_viewer_state", "search_annotations", "summarize_annotation",
      "set_omics_viewer_state", "set_scatter_view", "find_controls",
      "describe_control", "set_widgets", "run_control",
      "set_enrichment_parameters", "set_table_view", "create_figure",
      "get_figure", "update_figure")),
  "tools/list exposes exactly the 14 read-write tier tools in tier order"
)

## W3.3 a write tool executes through the SAME dispatcher -------------------

snapshot_before <- store_snapshot(deps$store)
resp <- agent_channel_handle(rw_state, call_msg(3, "set_widgets",
  list(patch = '{"dataspace.expr_heatmap.heatmap_colors": "RdGy"}',
       `_intent` = "switch the palette")))
ok(
  ut_cmp_identical(err_code(resp), NULL),
  "tools/call set_widgets answers without error on the read-write tier"
)
snapshot_after <- store_snapshot(deps$store)
ok(
  ut_cmp_identical(
    snapshot_after$values[["dataspace.expr_heatmap.heatmap_colors"]], "RdGy") &&
    ut_cmp_identical(
      snapshot_before$values[["dataspace.expr_heatmap.heatmap_colors"]],
      "RdBu") &&
    ut_cmp_identical(
      snapshot_after$epochs[["dataspace.expr_heatmap.heatmap_colors"]],
      snapshot_before$epochs[["dataspace.expr_heatmap.heatmap_colors"]] + 1L),
  "set_widgets through the channel CHANGED the store: value flipped and epoch bumped (origin=agent pathway)"
)

## W3.4 read tier still refuses writes (single gate, tier name in error) ----

read_state <- agent_channel_state_new(tools = registry$tools,
                                      session_token = domain$token)
resp <- agent_channel_handle(read_state, call_msg(4, "set_widgets",
  list(patch = '{"dataspace.expr_heatmap.heatmap_colors": "BrBG"}',
       `_intent` = "try a write")))
ok(
  ut_cmp_identical(err_code(resp), -32602L) &&
    grepl("'read' tier", resp, fixed = TRUE),
  "the READ tier still refuses write tools (-32602, error names the tier)"
)
ok(
  ut_cmp_identical(
    store_snapshot(deps$store)$values[["dataspace.expr_heatmap.heatmap_colors"]],
    "RdGy"),
  "the refused write did not touch the store"
)

## W3.5 run_control: settled-policy rejection receipt ----------------------

resp <- agent_channel_handle(rw_state, call_msg(5, "run_control",
  list(action_id = "snapshot_save", `_intent` = "try a forbidden save")))
ok(
  ut_cmp_identical(err_code(resp), NULL) &&
    grepl("excluded_by_policy", resp, fixed = TRUE) &&
    grepl("permanently excluded", resp, fixed = TRUE),
  "run_control settled-policy rejection arrives as a self-describing receipt (ok result, excluded_by_policy)"
)
resp <- agent_channel_handle(rw_state, call_msg(6, "run_control",
  list(action_id = "action:dataspace.tab_feature.clear",
       `_intent` = "clear the table selection")))
snapshot_run <- store_snapshot(deps$store)
ok(
  ut_cmp_identical(err_code(resp), NULL) &&
    ut_cmp_identical(
      snapshot_run$actions[["dataspace.tab_feature.clear"]], 1L),
  "run_control ran the registered action and recorded its epoch through the channel"
)
resp <- agent_channel_handle(rw_state, call_msg(7, "run_control",
  list(action_id = "action:no.such_action", `_intent` = "typo")))
ok(
  ut_cmp_identical(err_code(resp), NULL) &&
    grepl("no.such_action", resp, fixed = TRUE) &&
    grepl("Unknown", resp, fixed = TRUE),
  "run_control with an unknown id returns the self-describing receipt (agent_action_apply catches into error text + suggestions)"
)
resp <- agent_channel_handle(rw_state, call_msg(12, "set_widgets",
  list(pat = "{}", `_intent` = "wrong arg name")))
ok(
  ut_cmp_identical(err_code(resp), -32602L) &&
    grepl("Unknown argument 'pat'", resp, fixed = TRUE),
  "a genuinely malformed WRITE call surfaces the dispatcher-layer error (-32602, exact-name validation I9)"
)

## W3.6 rate cap counts writes like reads (D6) ------------------------------

capped_rw <- agent_channel_state_new(tools = registry$tools,
                                     session_token = domain$token,
                                     tier = "read-write", max_calls = 2L)
ok(
  ut_cmp_identical(
    err_code(agent_channel_handle(capped_rw, call_msg(8, "find_controls"))),
    NULL) &&
    ut_cmp_identical(
      err_code(agent_channel_handle(capped_rw, call_msg(9, "set_widgets",
        list(patch = '{"dataspace.expr_heatmap.heatmap_colors": "PuOr"}',
             `_intent` = "second call")))),
      NULL) &&
    ut_cmp_identical(capped_rw$calls, 2L),
  "reads AND writes increment the same per-connection counter"
)
cap_write <- agent_channel_handle(capped_rw, call_msg(10, "set_widgets",
  list(patch = '{"dataspace.expr_heatmap.heatmap_colors": "BrBG"}',
       `_intent` = "third call")))
cap_read <- agent_channel_handle(capped_rw, call_msg(11, "find_controls"))
ok(
  ut_cmp_identical(err_code(cap_write), -32000L) &&
    grepl("Rate cap exceeded", cap_write, fixed = TRUE) &&
    ut_cmp_identical(err_code(cap_read), -32000L) &&
    grepl("Rate cap exceeded", cap_read, fixed = TRUE),
  "the cap hits writes exactly like reads (-32000, rate text)"
)
ok(
  ut_cmp_identical(
    store_snapshot(deps$store)$values[["dataspace.expr_heatmap.heatmap_colors"]],
    "PuOr"),
  "the cap-refused write never reached the store"
)

## W3.7 write-receipt helper unit probes -------------------------------------

synthetic <- ellmer::ContentToolResult(
  value = "x", extra = list(data = list(
    applied = paste0("key_", 1:8),
    unchanged = c("u1", "u2"),
    rejected = list(list(id = "bad1", reason = "nope"),
                    list(id = "bad2", reason = "nope")))))
receipt <- .agent_channel_write_receipt(synthetic)
ok(
  ut_cmp_identical(receipt$applied, 8L) &&
    ut_cmp_identical(receipt$unchanged, 2L) &&
    ut_cmp_identical(receipt$rejected, 2L) &&
    ut_cmp_identical(receipt$ran, 0L) &&
    ut_cmp_identical(receipt$keys, paste0("key_", 1:5)) &&
    ut_cmp_identical(receipt$applied_values, NULL),
  "receipt summary counts applied/unchanged/rejected and caps keys at 5 NAMES (values dropped)"
)
ok(
  ut_cmp_identical(.agent_channel_write_receipt(NULL),
                   list(applied = 0L, ran = 0L, unchanged = 0L,
                        rejected = 0L, keys = character())) &&
    ut_cmp_identical(
      .agent_channel_write_receipt(ellmer::ContentToolResult(value = "plain"))
        $applied, 0L),
  "receipt summary degrades to zero counts on absent/foreign payloads"
)
figure_receipt <- .agent_channel_write_receipt(ellmer::ContentToolResult(
  value = "x",
  extra = list(data = list(figure_id = "fig_9", row_count = 10L))))
ok(
  ut_cmp_identical(figure_receipt$figure_id, "fig_9"),
  "figure receipts carry the bounded figure_id (an id, never a value)"
)
ok(
  ut_cmp_identical(
    .agent_channel_write_error(paste(rep("x", 400L), collapse = "")),
    paste0(paste(rep("x", 200L), collapse = ""), "...")),
  "write-failure messages are bounded to 200 chars for the log event"
)

# ==========================================================================
# [W4] channel-tagged write logging (D5): bounded receipts, tags, no values
# ==========================================================================

log_dir <- file.path(tempdir(), "agent-channel-writes-logs")
unlink(log_dir, recursive = TRUE, force = TRUE)
log_domain <- shiny:::MockShinySession$new()
log_logger <- agent_logger_new(
  log_domain,
  list(enabled = TRUE, directory = log_dir, max_bytes = 10485760L))
agent_logger_set_enabled(log_logger, TRUE)
logged_state <- agent_channel_state_new(tools = registry$tools,
                                        logger = log_logger,
                                        session_token = log_domain$token,
                                        tier = "read-write")

invisible(agent_channel_handle(logged_state, req_msg(21, "initialize")))
log_mark <- length(readLines(log_logger$path))
invisible(agent_channel_handle(logged_state, call_msg(22, "set_widgets",
  list(patch = '{"dataspace.expr_heatmap.heatmap_colors": "RdYlBu"}',
       `_intent` = "a logged write"))))
write_lines <- readLines(log_logger$path)
start_lines <- grep("channel_write_start", write_lines, fixed = TRUE)
end_lines <- grep("channel_write_end", write_lines, fixed = TRUE)
ok(
  ut_cmp_identical(length(start_lines), 1L) &&
    ut_cmp_identical(length(end_lines), 1L),
  "one executed write emits exactly one channel_write_start and one channel_write_end"
)
end_event <- tryCatch(
  jsonlite::fromJSON(write_lines[end_lines], simplifyVector = FALSE)$details,
  error = function(e) NULL)
ok(
  !is.null(end_event) &&
    ut_cmp_identical(end_event$tool, "set_widgets") &&
    ut_cmp_identical(end_event$ok, TRUE) &&
    ut_cmp_identical(is.numeric(end_event$duration_ms), TRUE) &&
    ut_cmp_identical(nchar(end_event$session), 16L) &&
    ut_cmp_identical(end_event$session,
                     .agent_channel_token_tag(log_domain$token)),
  "channel_write_end carries tool, ok, duration, and the 16-hex session tag"
)
ok(
  !is.null(end_event) &&
    ut_cmp_identical(end_event$receipt$applied, 1L) &&
    ut_cmp_identical(end_event$receipt$keys,
                     "dataspace.expr_heatmap.heatmap_colors") &&
    ut_cmp_identical(end_event$receipt$applied_values, NULL) &&
    ut_cmp_identical(
      any(grepl("RdYlBu", write_lines, fixed = TRUE)), FALSE) &&
    ut_cmp_identical(
      any(grepl("applied_values", end_lines, fixed = TRUE)), FALSE),
  "the receipt summary is bounded: key NAME + counts only, the widget VALUE never enters the log"
)
ok(
  ut_cmp_identical(
    any(grepl(log_domain$token, write_lines, fixed = TRUE)), FALSE) &&
    ut_cmp_identical(
      any(grepl("channel_write", readLines(log_logger$path)[seq_len(log_mark)],
                fixed = TRUE)), FALSE),
  "write events carry the tag (never the session token) and only fire for write tools"
)
# a READ call on the same state emits no write events
invisible(agent_channel_handle(logged_state, call_msg(23, "find_controls")))
ok(
  ut_cmp_identical(
    length(grep("channel_write", readLines(log_logger$path), fixed = TRUE)),
    2L),
  "read-tool calls emit no write events (start/end pair only from the write)"
)
# a FAILED write logs the bounded error (malformed patch JSON -> the
# handler stops, the dispatcher answers ok = FALSE)
invisible(agent_channel_handle(logged_state, call_msg(24, "set_widgets",
  list(patch = "not json", `_intent` = "a failing write"))))
all_lines <- readLines(log_logger$path)
fail_lines <- grep("channel_write_end", all_lines, fixed = TRUE)
fail_event <- jsonlite::fromJSON(
  all_lines[[fail_lines[[length(fail_lines)]]]],
  simplifyVector = FALSE)$details
ok(
  ut_cmp_identical(fail_event$ok, FALSE) &&
    nzchar(fail_event$error) &&
    ut_cmp_identical(is.null(fail_event$receipt), TRUE),
  "a failed write logs ok=FALSE with a bounded error and no receipt"

)
# figure writes name the figure id
invisible(agent_channel_handle(logged_state, call_msg(25, "create_figure",
  list(template = "volcano", x = "ttest|RE_vs_ME|mean.diff",
       y = "ttest|RE_vs_ME|log.fdr", title = "Logged volcano",
       `_intent` = "a figure write"))))
fig_lines <- grep("channel_write_end", readLines(log_logger$path),
                  fixed = TRUE)
fig_event <- jsonlite::fromJSON(
  readLines(log_logger$path)[[fig_lines[[length(fig_lines)]]]],
  simplifyVector = FALSE)$details
ok(
  ut_cmp_identical(fig_event$receipt$figure_id, "fig_1"),
  "figure write receipts name the created figure id"
)

# ==========================================================================
# [W5] attach hook: figure fields land in the registry entry
# ==========================================================================

clean_slate()
ok(
  ut_cmp_identical(
    with_env(c(OMICSVIEWER_AGENT_CHANNEL = "1",
               OMICSVIEWER_AGENT_CHANNEL_WRITES = "1",
               OMICSVIEWER_AGENT_CHANNEL_PORT = as.character(find_free_port())),
      agent_channel_session_attach(
        session = stub_session("writes-attach-token"),
        tools = registry$tools, session_domain = domain,
        logger = log_logger,
        figures = figures, figure_directory = figure_directory,
        figure_renderer = figure_renderer)),
    TRUE),
  "attach binds with the Stage-3 figure fields when writes are enabled"
)
entry <- get("writes-attach-token", envir = get(
  ".agent_channel_registry", envir = asNamespace("omicsViewer")))
ok(
  ut_cmp_identical(identical(entry$figures, figures), TRUE) &&
    ut_cmp_identical(entry$figure_directory, figure_directory) &&
    ut_cmp_identical(is.function(entry$figure_renderer), TRUE) &&
    ut_cmp_identical(agent_channel_server_info()$tier, "read-write"),
  "the registry entry carries figures/figure_directory/figure_renderer and the server runs read-write"
)
log_lines_all <- readLines(log_logger$path)
bind_line <- grep('"channel_bind"', log_lines_all, fixed = TRUE, value = TRUE)
bind_event <- if (length(bind_line))
  jsonlite::fromJSON(bind_line[[length(bind_line)]],
                     simplifyVector = FALSE)$details else NULL
ok(
  !is.null(bind_event) && ut_cmp_identical(bind_event$tier, "read-write"),
  "channel_bind logs the live server tier (read-write), not a hard-coded one"
)
agent_channel_session_detach("writes-attach-token")
ok(
  ut_cmp_identical(agent_channel_server_info()$running, FALSE),
  "detach stops the server on last unbind (Stage-3 fields included)"
)

# ==========================================================================
# [W6] loopback smoke: socket writes + figure GET over plain HTTP
# ==========================================================================

if (!requireNamespace("websocket", quietly = TRUE)) {
  ok(TRUE, "socket smoke skipped: the websocket package is not installed")
} else {

with_env(c(OMICSVIEWER_AGENT_CHANNEL = "1",
           OMICSVIEWER_AGENT_CHANNEL_WRITES = "1",
           OMICSVIEWER_AGENT_CHANNEL_PORT =
             as.character(smoke_port <- find_free_port())), {
  clean_slate()
  agent_channel_registry_bind(
    "smoke-session-token",
    list(tools = registry$tools, session_domain = domain,
         logger = log_logger,
         figures = figures, figure_directory = figure_directory,
         figure_renderer = figure_renderer))
  ok(ut_cmp_identical(agent_channel_server_start(), TRUE),
     "server starts read-write on a random free port with a bound session")
  token <- agent_channel_server_info()$token
  base_url <- paste0("http://127.0.0.1:", smoke_port)

  wait_for <- function(probe, timeout = 20) {
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

  ## W6.1 tier + tools over the wire ----------------------------------------

  resp <- rpc(cli, req_msg(1, "initialize"))
  ok(
    ut_cmp_identical(err_code(resp), NULL) &&
      ut_cmp_identical(parse_response(resp)$result$tool_tier, "read-write"),
    "socket initialize reports the read-write tier"
  )
  resp <- rpc(cli, req_msg(2, "tools/list"))
  ok(
    ut_cmp_identical(err_code(resp), NULL) &&
      ut_cmp_identical(length(parse_response(resp)$result$tools), 14L),
    "socket tools/list returns all 14 read-write tier tools"
  )

  ## W6.2 set_widgets round trip changes the STORE --------------------------

  snapshot_ws_before <- store_snapshot(deps$store)
  resp <- rpc(cli, call_msg(3, "set_widgets",
    list(patch = '{"dataspace.expr_heatmap.heatmap_colors": "PiYG"}',
         `_intent` = "socket palette switch")))
  ok(ut_cmp_identical(err_code(resp), NULL),
     "socket set_widgets answers without error")
  snapshot_ws_after <- store_snapshot(deps$store)
  ok(
    ut_cmp_identical(
      snapshot_ws_after$values[["dataspace.expr_heatmap.heatmap_colors"]],
      "PiYG") &&
      ut_cmp_identical(
        snapshot_ws_after$epochs[["dataspace.expr_heatmap.heatmap_colors"]],
        snapshot_ws_before$epochs[["dataspace.expr_heatmap.heatmap_colors"]] +
          1L),
    "socket set_widgets APPLIED to the harness store (value + epoch changed)"
  )

  ## W6.3 run_control over the socket ---------------------------------------

  resp <- rpc(cli, call_msg(4, "run_control",
    list(action_id = "action:dataspace.tab_feature.clear",
         `_intent` = "clear from the socket")))
  ok(
    ut_cmp_identical(err_code(resp), NULL) &&
      grepl("ran", resp, fixed = TRUE),
    "socket run_control executed the registered action"
  )
  resp <- rpc(cli, call_msg(5, "run_control",
    list(action_id = "snapshot_save", `_intent` = "forbidden from socket")))
  ok(
    ut_cmp_identical(err_code(resp), NULL) &&
      grepl("excluded_by_policy", resp, fixed = TRUE) &&
      grepl("permanently excluded", resp, fixed = TRUE),
    "socket run_control settled-policy rejection arrives as the dispatcher receipt (snapshot lifecycle stays user-only)"
  )

  ## W6.4 figure round trip: create -> get -> update ------------------------

  resp <- rpc(cli, call_msg(6, "create_figure",
    list(template = "volcano", x = "ttest|RE_vs_ME|mean.diff",
         y = "ttest|RE_vs_ME|log.fdr", title = "Socket volcano",
         `_intent` = "volcano over the socket")))
  created_id <- tryCatch(
    jsonlite::fromJSON(
      parse_response(resp)$result$content[[1]]$text,
      simplifyVector = FALSE)$figure_id,
    error = function(e) NULL)
  ok(
    ut_cmp_identical(err_code(resp), NULL) &&
      ut_cmp_identical(grepl("^fig_[0-9]+$", created_id %||% ""), TRUE),
    paste0("socket create_figure rendered a figure (", created_id, ")")
  )
  resp <- rpc(cli, call_msg(7, "get_figure",
    list(figure_id = created_id, `_intent` = "read it back")))
  ok(
    ut_cmp_identical(err_code(resp), NULL) &&
      grepl("Socket volcano", resp, fixed = TRUE),
    "socket get_figure returns the stored spec"
  )
  resp <- rpc(cli, call_msg(8, "update_figure",
    list(figure_id = created_id,
         changes = list(labels = list(title = "Socket revision")),
         `_intent` = "retitle it")))
  revised_id <- tryCatch(
    jsonlite::fromJSON(
      parse_response(resp)$result$content[[1]]$text,
      simplifyVector = FALSE)$figure_id,
    error = function(e) NULL)
  ok(
    ut_cmp_identical(err_code(resp), NULL) &&
      ut_cmp_identical(grepl("^fig_[0-9]+$", revised_id %||% ""), TRUE) &&
      ut_cmp_identical(revised_id != created_id, TRUE) &&
      grepl("Socket revision", resp, fixed = TRUE),
    paste0("socket update_figure produced the patched revision (",
           revised_id, ")")
  )
  ok(
    ut_cmp_identical(length(isolate(figures())), 3L) &&
      ut_cmp_identical(isolate(figure_counter()), 3L) &&
      ut_cmp_identical(
        isolate(figures())[[revised_id]]$parent_id, created_id),
    "the figure registry tracked all revisions (incl. the W4 logged figure)"
  )

  ## W6.5 figure GET over plain HTTP ----------------------------------------
  # hoisted OUT of this websocket conditional (O2): the plain-HTTP
  # surface has no websocket dependency, so its coverage now runs in
  # the [W6-HTTP] section below whenever curl is available - even on
  # machines without the websocket package (which would otherwise
  # mask the whole subsection behind the skip above).

  ## W6.6 write logging through the socket path ------------------------------

  invisible(rpc(cli, call_msg(9, "set_widgets",
    list(patch = '{"dataspace.expr_heatmap.heatmap_colors": "BrBG"}',
         `_intent` = "a socket write to log"))))
  socket_log <- readLines(log_logger$path)
  ok(
    ut_cmp_identical(
      any(grepl("channel_write_start", socket_log, fixed = TRUE)) &&
        any(grepl("channel_write_end", socket_log, fixed = TRUE)), TRUE),
    "socket write calls emit channel_write_start/end events (D5)"
  )

  # cleanup
  try(cli$ws$close(), silent = TRUE)
  agent_channel_registry_unbind("smoke-session-token")
  agent_channel_server_stop()
  ok(ut_cmp_identical(agent_channel_server_info()$running, FALSE),
     "socket cleanup stopped the server")
})
}

# ==========================================================================
# [W6-HTTP] figure GET over plain HTTP (O2: hoisted out of the websocket
#           conditional above - runs whenever curl is available, even
#           without the websocket package; the plain-HTTP surface has
#           no WS dependency)
# ==========================================================================

if (!requireNamespace("curl", quietly = TRUE)) {
  ok(TRUE, "figure HTTP tests skipped: the curl package is not installed")
} else {

with_env(c(OMICSVIEWER_AGENT_CHANNEL = "1",
           OMICSVIEWER_AGENT_CHANNEL_WRITES = "1",
           OMICSVIEWER_AGENT_CHANNEL_PORT =
             as.character(smoke_port <- find_free_port())), {
  clean_slate()

  # Seed the figure registry through the PURE handler (the W4
  # pattern): no websocket needed - the dispatcher mutates the same
  # shared figures/figure_counter reactives the server resolves from.
  http_state <- agent_channel_state_new(
    tools = registry$tools, logger = logger,
    session_token = domain$token, tier = "read-write")
  resp <- agent_channel_handle(http_state, call_msg(31, "create_figure",
    list(template = "volcano", x = "ttest|RE_vs_ME|mean.diff",
         y = "ttest|RE_vs_ME|log.fdr", title = "HTTP volcano",
         `_intent` = "seed the HTTP endpoint")))
  http_created_id <- tryCatch(
    jsonlite::fromJSON(
      parse_response(resp)$result$content[[1]]$text,
      simplifyVector = FALSE)$figure_id,
    error = function(e) NULL)
  ok(
    ut_cmp_identical(err_code(resp), NULL) &&
      ut_cmp_identical(grepl("^fig_[0-9]+$", http_created_id %||% ""), TRUE),
    paste0("HTTP section: pure-handler create_figure seeded the figure registry (",
           http_created_id, ") without any websocket")
  )
  resp <- agent_channel_handle(http_state, call_msg(32, "update_figure",
    list(figure_id = http_created_id,
         changes = list(labels = list(title = "HTTP revision")),
         `_intent` = "seed a revision")))
  http_revised_id <- tryCatch(
    jsonlite::fromJSON(
      parse_response(resp)$result$content[[1]]$text,
      simplifyVector = FALSE)$figure_id,
    error = function(e) NULL)
  ok(
    ut_cmp_identical(err_code(resp), NULL) &&
      ut_cmp_identical(grepl("^fig_[0-9]+$", http_revised_id %||% ""), TRUE) &&
      ut_cmp_identical(http_revised_id != http_created_id, TRUE),
    "HTTP section: pure-handler update_figure seeded a distinct revision"
  )

  agent_channel_registry_bind(
    "http-figure-token",
    list(tools = registry$tools, session_domain = domain,
         logger = log_logger,
         figures = figures, figure_directory = figure_directory,
         figure_renderer = figure_renderer))
  ok(ut_cmp_identical(agent_channel_server_start(), TRUE),
     "HTTP section: server starts read-write on a random free port (websocket package not required)")
  token <- agent_channel_server_info()$token
  base_url <- paste0("http://127.0.0.1:", smoke_port)

  # async HTTP GET: client and server share this process's event loop,
  # so the fetch rides curl's multi pool while later services httpuv.
  http_get <- function(url) {
    resp <- NULL
    failed <- NULL
    pool <- curl::new_pool()
    curl::curl_fetch_multi(
      url,
      done = function(res) resp <<- res,
      fail = function(e) failed <<- e,
      pool = pool)
    deadline <- Sys.time() + 20
    while (is.null(resp) && is.null(failed) && Sys.time() < deadline) {
      later::run_now(0.05)
      curl::multi_run(timeout = 0, pool = pool)
    }
    if (!is.null(failed))
      structure(list(status = NA_integer_, content = raw(0)),
                class = "curl_fail", message = conditionMessage(failed))
    resp
  }

  get_png <- function(path) http_get(paste0(base_url, path))
  is_png <- function(bytes)
    length(bytes) >= 4L && identical(bytes[[1]], as.raw(0x89)) &&
    identical(bytes[[2]], as.raw(0x50))  # \x89P (NG)

  r <- get_png(paste0("/", token, "/figure/", http_created_id, "/preview.png"))
  ok(
    !is.null(r) && ut_cmp_identical(r$status, 200L) &&
      ut_cmp_identical(is_png(r$content), TRUE) &&
      ut_cmp_identical(length(r$content) > 500L, TRUE),
    "figure GET serves the preview PNG bytes (200, image/png, PNG magic)"
  )
  ok(
    ut_cmp_identical(
      grepl("image/png", rawToChar(r$headers), fixed = TRUE), TRUE),
    "figure GET answers with the image/png content type"
  )
  r <- get_png(paste0("/", token, "/figure/", http_created_id, "/full.png"))
  ok(
    ut_cmp_identical(r$status, 200L) &&
      ut_cmp_identical(is_png(r$content), TRUE) &&
      ut_cmp_identical(length(r$content) > 1000L, TRUE),
    "figure GET serves the full-resolution PNG (bigger than the preview)"
  )
  r_preview <- get_png(paste0("/", token, "/figure/", http_created_id,
                              "/preview.png"))
  ok(
    ut_cmp_identical(
      identical(r_preview$content, get_png(
        paste0("/", token, "/figure/", http_created_id, "/full.png"))$content),
      FALSE),
    "preview and full artifacts differ (renderer honors the size segment)"
  )
  r <- get_png(paste0("/", token, "/figure/", http_revised_id, "/preview.png"))
  ok(
    ut_cmp_identical(r$status, 200L) &&
      ut_cmp_identical(is_png(r$content), TRUE),
    "the revision figure is servable too"
  )

  # wrong token: the SAME silent refusal as WS - generic 404, no echo
  r <- get_png(paste0("/wrongtoken/figure/", http_created_id, "/preview.png"))
  ok(
    ut_cmp_identical(r$status, 404L) &&
      ut_cmp_identical(
        grepl(token, rawToChar(r$content), fixed = TRUE), FALSE),
    "wrong-token figure GET is refused with the generic 404 and no token echo"
  )
  # traversal attempts: encoded over the wire; literal at the PARSER
  r <- get_png(paste0("/", token, "/figure/..%2F..%2Fetc%2Fpasswd/preview.png"))
  ok(
    ut_cmp_identical(r$status, 404L) &&
      ut_cmp_identical(length(r$content) < 100L, TRUE),
    "percent-encoded traversal (..%2F..%2F) is refused (404, short generic body)"
  )
  r <- get_png(paste0("/", token, "/figure/fig_1%2Fsecrets/preview.png"))
  ok(
    ut_cmp_identical(r$status, 404L),
    "encoded slash inside the figure id is refused (404)"
  )
  # O3: a LITERAL-../ wire fetch was vacuous - libcurl normalizes dot
  # segments client-side, so the server never saw the traversal (the
  # 404 came back for an unrelated normalized path). The guarantee is
  # asserted where it actually lives: the pure parser, on literal
  # traversal shapes (a 5-segment ../ id fails the id regex; deeper
  # shapes fail the segment count), with a positive control.
  ok(
    ut_cmp_identical(
      .agent_channel_figure_parse(
        paste0("/", token, "/figure/", http_created_id, "/preview.png"),
        token),
      list(figure_id = http_created_id, size = "preview")),
    "parser control: the valid figure GET path parses (so the refusals below mean something)"
  )
  ok(
    is.null(.agent_channel_figure_parse(
      paste0("/", token, "/figure/../../../etc/passwd/preview.png"), token)) &&
      is.null(.agent_channel_figure_parse(
        paste0("/", token, "/figure/../preview.png"), token)) &&
      is.null(.agent_channel_figure_parse(
        paste0("/", token, "/figure/a/../../etc/passwd/preview.png"), token)),
    "literal traversal (../) is refused by the figure-path parser itself (O3: dot-segment ids fail the id regex, deeper shapes fail the segment count)"
  )
  # unknown figure and bad size
  r <- get_png(paste0("/", token, "/figure/no_such_figure/preview.png"))
  ok(
    ut_cmp_identical(r$status, 404L),
    "unknown figure id answers 404"
  )
  r <- get_png(paste0("/", token, "/figure/", http_created_id, "/medium.png"))
  ok(ut_cmp_identical(r$status, 404L), "non-png size answers 404")
  r <- get_png(paste0("/", token, "/figure/", http_created_id, "/full.PNG"))
  ok(ut_cmp_identical(r$status, 404L), "case-variant size answers 404")
  r <- get_png(paste0("/anything/else"))
  ok(
    ut_cmp_identical(r$status, 404L) &&
      ut_cmp_identical(rawToChar(r$content), "Not Found"),
    "non-figure plain HTTP keeps the generic 404 (no path echo)"
  )
  # the token never appears in any HTTP response body
  ok(
    ut_cmp_identical(
      any(grepl(token, vapply(
        list(get_png(paste0("/", token, "/figure/nope/preview.png")),
             get_png(paste0("/wrongtoken/figure/", http_created_id,
                            "/preview.png"))),
        function(x) rawToChar(x$content), character(1)),
        fixed = TRUE)), FALSE),
    "no HTTP response body ever contains the pairing token"
  )

  # cleanup
  agent_channel_registry_unbind("http-figure-token")
  agent_channel_server_stop()
  ok(ut_cmp_identical(agent_channel_server_info()$running, FALSE),
     "HTTP section: cleanup stopped the server")
})
}

# final hygiene: nothing bound, no server, no env leakage
clean_slate()
ok(
  ut_cmp_identical(length(agent_channel_registry_list()), 0L) &&
    ut_cmp_identical(agent_channel_server_info()$running, FALSE) &&
    ut_cmp_identical(Sys.getenv("OMICSVIEWER_AGENT_CHANNEL", unset = ""), "") &&
    ut_cmp_identical(
      Sys.getenv("OMICSVIEWER_AGENT_CHANNEL_WRITES", unset = ""), ""),
  "suite leaves no bindings, no server, and no env leakage"
)
