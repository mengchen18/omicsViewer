# Stage 1 live suite for the agent channel: a scripted "second agent"
# constructs the REAL tool registry (agent_tool_registry(...) over the
# demo dataset) and drives every tool through agent_channel_dispatch,
# then proves the effects identical to the existing test-hook pathway
# (module_agentTestHooks.R) on identically-seeded widget stores.
#
# Run: Rscript tests/test_agentChannelLive.R   (from the repo root)
#
# Test-only by design (AGENT_CHANNEL_PLAN.md Stage 1 / decision D11):
# no production file is loaded differently or edited here, and the full
# 337-case browser sweep is NOT replayed - that gate runs after Stages
# 2-3. This suite provides representative end-to-end coverage of the
# plug-and-play claim: registry -> dispatcher -> real effects.
#
# Compromises forced by full-fidelity registry construction are
# documented as decisions I11-I17 in AGENT_CHANNEL_STATUS.md.

library(unittest, quietly = TRUE)

if (!requireNamespace("pkgload", quietly = TRUE) ||
    !requireNamespace("ellmer", quietly = TRUE) ||
    !requireNamespace("shinychat", quietly = TRUE)) {
  ok(TRUE, "live agent channel tests skipped because optional packages are unavailable")
  quit(save = "no", status = 0)
}

# The package is NOT reinstalled for this suite; load the working tree
# (the agentChannel/snapshotRoundTrip convention, decision I8).
pkg_root <- if (file.exists("DESCRIPTION")) "." else
  if (file.exists(file.path("..", "DESCRIPTION"))) ".." else "."
pkgload::load_all(pkg_root, quiet = TRUE, attach_testthat = FALSE)
library(ellmer)

for (nm in c(
  "widget_store_new", "widget_store_child", "widget_binding", "widget_action",
  "store_register", "store_register_action", "store_seed", "store_apply",
  "store_read", "store_snapshot", "store_action_view",
  "agent_widget_apply", "agent_tool_registry", "agent_channel_dispatch",
  "agent_channel_toolspecs", "agent_normalize_state_update",
  "agent_normalize_scatter_view", "agent_normalize_enrichment_update",
  "agent_normalize_table_view_update", "agent_test_hooks_module",
  "agent_logger_new", "trisetter", "build_app_state", "agent_compact_state",
  "agent_annotation_catalog", "agent_figure_grammar"))
  assign(nm, get(nm, envir = asNamespace("omicsViewer")))

# ==========================================================================
# [1] Real dependencies from the demo dataset (2702 features x 60 samples)
# ==========================================================================

demo_path <- system.file("extdata", "demo.RDS", package = "omicsViewer")
if (!nzchar(demo_path))
  demo_path <- file.path(pkg_root, "inst", "extdata", "demo.RDS")
demo <- readRDS(demo_path)
fd <- Biobase::fData(demo)
pd <- Biobase::pData(demo)
mat <- Biobase::exprs(demo)

ok(
  ut_cmp_identical(nrow(fd), 2702L) && ut_cmp_identical(nrow(pd), 60L),
  "harness loads the real demo dataset (2702 features x 60 samples)"
)

# demo.RDS ships without gene-set annotations, but the enrichment tier
# requires attr(fData, "GS") (agent_normalize_enrichment_update); attach
# a small 3-column gene-set table (documented harness prep, decision
# I13). All other demo attributes (quickViews, columns) are untouched.
set.seed(1L)
fids <- rownames(fd)
gs_ids <- paste0("gs_set", c("A", "B", "C", "D"))
gs <- data.frame(
  featureId = factor(sample(fids, 300L, replace = TRUE), levels = fids),
  gsId = factor(sample(gs_ids, 300L, replace = TRUE), levels = gs_ids),
  weight = 1L)
attr(fd, "GS") <- gs
ok(
  ut_cmp_identical(length(unique(attr(fd, "GS")$gsId)), 4L),
  "harness attaches a synthetic GS attribute so the enrichment tier is exercisable"
)

# The real trisetter output feeds the store's cascade choices providers
# (the app's providers are the same expressions over a reactive triset;
# decision I12). The demo quickViews attributes ride along to the apply
# callbacks exactly as the app passes attr(v1(), "quickViews").
ftri <- trisetter(meta = fd, expr = mat, combine = "feature")
ptri <- trisetter(meta = pd, expr = mat, combine = "pheno")
quick_views <- list(feature = attr(fd, "quickViews"),
                    sample = attr(pd, "quickViews"))
ok(
  ut_cmp_identical(nrow(ftri) > 100L, TRUE) &&
    ut_cmp_identical(nrow(ptri) > 10L, TRUE) &&
    ut_cmp_identical(nrow(quick_views$feature), 6L) &&
    ut_cmp_identical(nrow(quick_views$sample), 1L),
  "real trisetter cascade tables and demo quick views are available"
)

# Tab vocabularies mirror L0's agent_data_tabs / agent_analysis_tabs
# computed for the demo dataset (GS attached -> ORA/fGSEA present).
DS_TABS <- c("Feature", "Feature table", "Sample", "Sample table", "Cor",
             "Heatmap", "Dynamic heatmap", "Expression", "GSList")
ANALYSIS_TABS <- c("Feature", "ORA", "fGSEA", "Geneshot", "Sample")

# Dataset reactives, exactly the contract ai_assistant_module receives.
feature_data <- reactive(fd)
sample_data <- reactive(pd)
expression_data <- reactive(mat)
selected_features <- reactiveVal(character())
selected_samples <- reactiveVal(character())

# ==========================================================================
# [2] Widget-store factory: a minimal-but-valid mirror of the app's own
#     registrations and seeding (blast-radius rule: new CALLER of the
#     store, no store internals touched). Bindings mirror module_meta_
#     scatter.R, L1_module_data_space.R, L1_module_result_space.R,
#     module_ora.R, module_dataTable.R and heatmapshinyApp.R.
# ==========================================================================

make_store <- function() {
  store <- widget_store_new()

  # dataspace navbar (L1_module_data_space.R: active_tab)
  sds <- widget_store_child(store, "dataspace")
  store_register(
    sds,
    widget_binding("active_tab", "navbar", label = "Data-space tab",
                   help = "Active tab of the data-space navbar",
                   values = DS_TABS))
  store_seed(sds, list(active_tab = "Feature"))

  # scatter axis cascades (module_meta_scatter.R): same kinds, same
  # depends_on topology, same choices-provider expressions; the providers
  # close over the STATIC demo triset (decision I12).
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
    # seed the demo's own default axes (the PCA quick view), split into
    # cascade parts exactly as the app seeds settled triples
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

  # resultspace navbar (L1_module_result_space.R: analyst_tab)
  srs <- widget_store_child(store, "resultspace")
  store_register(
    srs,
    widget_binding("analyst_tab", "navbar", label = "Analysis tab",
                   help = "Active tab of the analysis navbar",
                   choices_provider = function(v) ANALYSIS_TABS))
  store_seed(srs, list(analyst_tab = "Feature"))

  # feature table widgets + clear action (module_dataTable.R)
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

  # one heatmap widget (heatmapshinyApp.R: heatmap_colors) for set_widgets
  shm <- widget_store_child(store, "dataspace.expr_heatmap")
  store_register(
    shm,
    widget_binding("heatmap_colors", "select", label = "Heatmap color panel",
      help = "Diverging palette for the heatmap cells",
      values = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy",
                 "RdYlBu", "RdYlGn")))
  store_seed(shm, list(heatmap_colors = "RdBu"))

  # ORA collapse cascade + selected pathway row (module_ora.R)
  ko1 <- "resultspace.ora.xax_analysis"
  ko2 <- "resultspace.ora.xax_subset"
  sora <- widget_store_child(store, "resultspace.ora")
  store_register(
    sora,
    widget_binding("xax_analysis", "select", label = "Collapse category",
      help = "Category whose values features are collapsed on",
      choices_provider = function(v) unique(ftri[, 1])),
    widget_binding("xax_subset", "select", label = "Collapse subcategory",
      help = "Subcategory within the collapse category",
      depends_on = "xax_analysis",
      choices_provider = function(v) unique(ftri[ftri[, 1] %in% v[[ko1]], 2])),
    widget_binding("xax_variable", "select_cascaded", label = "Collapse variable",
      help = "Variable the features are collapsed on",
      depends_on = c("xax_analysis", "xax_subset"),
      choices_provider = function(v)
        ftri[ftri[, 1] %in% v[[ko1]] & ftri[, 2] %in% v[[ko2]], 3]),
    widget_binding("selected_row", "select", label = "Selected pathway",
      help = "Gene set whose results row is selected",
      choices_provider = function(v) gs_ids))
  store
}

ok(
  ut_cmp_identical(
    identical(store_snapshot(make_store()), store_snapshot(make_store())),
    TRUE),
  "the store factory produces identically-seeded stores (equivalence baseline)"
)

# ==========================================================================
# [3] L0-mirroring apply callbacks + compact-state builder
#
#     Thin validate + store_apply(origin = "agent") wrappers mirroring
#     L0_module_app.R's apply_agent_* (decision I14): the same
#     normalizers, the same canonical keys, the same receipts. Selections
#     go through a minimal recording selection bus (the app reads the
#     data-space module's ports; the ids/origin/mirror interface is
#     mirrored verbatim).
# ==========================================================================

make_deps <- function() {
  store <- make_store()
  sel_log <- new.env(parent = emptyenv())
  sel_log$feature <- NULL
  sel_log$sample <- NULL
  ports <- list(
    feature = list(apply = function(ids, origin, mirror)
      sel_log$feature <- ids),
    sample = list(apply = function(ids, origin, mirror)
      sel_log$sample <- ids))

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
    if (!is.null(validated$features))
      ports$feature$apply(ids = validated$features, origin = "system",
                          mirror = TRUE)
    if (!is.null(validated$samples))
      ports$sample$apply(ids = validated$samples, origin = "system",
                         mirror = TRUE)
    list(
      data_space_tab = validated$data_space_tab,
      analysis_space_tab = validated$analysis_space_tab,
      feature_count = length(validated$features %||% character()),
      sample_count = length(validated$samples %||% character()),
      example_features = utils::head(validated$features %||% character(), 20L),
      example_samples = utils::head(validated$samples %||% character(), 20L))
  }

  apply_scatter_view <- function(space, quick_view_id = NULL, x_axis = NULL,
                                 y_axis = NULL) {
    view <- agent_normalize_scatter_view(
      space = space, quick_view_id = quick_view_id, x_axis = x_axis,
      y_axis = y_axis, quick_views = quick_views,
      feature_columns = colnames(isolate(feature_data())),
      sample_columns = colnames(isolate(sample_data())))
    axis_data <- if (view$space == "feature") isolate(feature_data()) else
      isolate(sample_data())
    if (!any(vapply(c(view$x_axis, view$y_axis), function(nm)
      is.numeric(axis_data[[nm]]), logical(1))))
      stop("At least one scatter axis must contain numeric values.")
    prefix <- if (view$space == "feature") "dataspace.feature_space" else
      "dataspace.sample_space"
    split_axis <- function(axis) strsplit(axis, "|", fixed = TRUE)[[1]]
    patch <- stats::setNames(
      c(as.list(split_axis(view$x_axis)), as.list(split_axis(view$y_axis))),
      c(paste(prefix, c("x_analysis", "x_subset", "x_variable"), sep = "."),
        paste(prefix, c("y_analysis", "y_subset", "y_variable"), sep = ".")))
    shiny::isolate(store_apply(store, patch, origin = "agent", strict = TRUE))
    shiny::isolate(store_apply(
      store,
      list("dataspace.active_tab" = if (view$space == "feature")
        "Feature" else "Sample"),
      origin = "agent", strict = TRUE))
    view
  }

  apply_enrichment <- function(update) {
    validated <- agent_normalize_enrichment_update(update,
                                                   isolate(feature_data()))
    receipt <- shiny::isolate(store_apply(store, validated$patch,
                                          origin = "agent", strict = FALSE))
    list(method = validated$method, panel_tab = validated$tab,
         applied = receipt$applied, applied_values = receipt$diff,
         unchanged = receipt$skipped,
         rejected = receipt$rejected %||% list(), note = "panel note")
  }

  apply_table_view <- function(update) {
    validated <- agent_normalize_table_view_update(update)
    receipt <- shiny::isolate(store_apply(store, validated$patch,
                                          origin = "agent", strict = FALSE))
    list(table = validated$table, panel_tab = validated$tab,
         applied = receipt$applied, applied_values = receipt$diff,
         unchanged = receipt$skipped,
         rejected = receipt$rejected %||% list(), note = "panel note")
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
       apply_table_view = apply_table_view, state = state,
       sel_log = sel_log)
}

# ==========================================================================
# [4] The REAL registry: agent_tool_registry(...) called directly with
#     every formal sourced (module_aiAssistant.R builds it with exactly
#     these arguments). session_domain is a MockShinySession - the same
#     class shiny::testServer drives - so withReactiveDomain behaves as
#     in a session (decision I15).
# ==========================================================================

domain <- shiny:::MockShinySession$new()
logger <- agent_logger_new(domain)
figures <- reactiveVal(list())
figure_counter <- reactiveVal(0L)
figure_directory <- file.path(tempdir(), "agent-channel-live")
dir.create(figure_directory, showWarnings = FALSE, recursive = TRUE)

deps_live <- make_deps()
registry <- agent_tool_registry(
  state = deps_live$state,
  feature_data = feature_data,
  sample_data = sample_data,
  expression_data = expression_data,
  selected_features = selected_features,
  selected_samples = selected_samples,
  apply_state = deps_live$apply_state,
  apply_scatter_view = deps_live$apply_scatter_view,
  apply_enrichment = deps_live$apply_enrichment,
  apply_table_view = deps_live$apply_table_view,
  store = deps_live$store,
  figures = figures,
  figure_counter = figure_counter,
  figure_directory = figure_directory,
  session_domain = domain,
  logger = logger)

ok(
  ut_cmp_identical(length(registry$tools), 14L),
  "the real registry exposes all 14 tools (full tier: store + enrichment + table view)"
)

dispatch <- function(name, args)
  agent_channel_dispatch(registry$tools, name, args)

# ==========================================================================
# [5] toolspecs smoke over the live registry
# ==========================================================================

specs <- agent_channel_toolspecs(registry$tools)
ok(
  ut_cmp_identical(length(specs), 14L),
  "toolspecs project every registered tool"
)
ok(
  ut_cmp_identical(
    sort(names(specs)),
    sort(c("get_omics_viewer_state", "search_annotations",
           "summarize_annotation", "set_omics_viewer_state",
           "set_scatter_view", "find_controls", "describe_control",
           "set_widgets", "run_control", "set_enrichment_parameters",
           "set_table_view", "create_figure", "get_figure",
           "update_figure"))),
  "toolspecs are keyed by the wire names the dispatcher accepts"
)
ok(
  all(vapply(specs, function(spec)
    length(spec$description) == 1L && nzchar(spec$description), logical(1))),
  "every live spec carries a non-empty description"
)
ok(
  all(vapply(specs, function(spec)
    identical(spec$inputSchema$type, "object") &&
      "_intent" %in% spec$inputSchema$required, logical(1))),
  "every live inputSchema is an object requiring _intent"
)
ok(
  all(vapply(specs, function(spec)
    identical(spec$inputSchema$additionalProperties, FALSE), logical(1))),
  "every live inputSchema advertises additionalProperties: false (decision I5)"
)

# ==========================================================================
# [6] Dispatch all 14 tools through the channel with valid arguments.
#     Content assertions read value@extra$data - the structured original
#     for in-process consumers (decision I17); value@value is the
#     codec's JSON serialization and must never cross a wire (F3).
# ==========================================================================

## 6.1 get_omics_viewer_state (overview) ----------------------------------

r <- dispatch("get_omics_viewer_state", list(`_intent` = "read state"))
ok(ut_cmp_identical(r$ok, TRUE), "get_omics_viewer_state dispatches ok")
overview <- r$value@extra$data
ok(
  ut_cmp_identical(overview$dataset$dimensions[["features"]], 2702L) &&
    ut_cmp_identical(overview$dataset$dimensions[["samples"]], 60L),
  "state overview reports the real dataset dimensions"
)
ok(
  ut_cmp_identical(overview$scatter_view$feature$x$name,
                   "PCA|All|PC1(10.5%)") &&
    ut_cmp_identical(overview$scatter_view$feature$axis_mode, "quick"),
  "state overview anchors the seeded scatter view from the store"
)
ok(
  ut_cmp_identical(overview$capabilities$capability_count > 10L, TRUE),
  "state overview carries capability counts over the live store"
)

r <- dispatch("get_omics_viewer_state",
              list(sections = "quick_views", `_intent` = "list quick views"))
ok(
  ut_cmp_identical(r$ok, TRUE) &&
    ut_cmp_identical(
      vapply(r$value@extra$data$quick_views$feature,
             function(v) is.null(v$x) || is.character(v$x), logical(1)) |>
        all(), TRUE),
  "get_omics_viewer_state(sections = 'quick_views') returns full records"
)

## 6.2 search_annotations --------------------------------------------------

r <- dispatch("search_annotations",
              list(space = "feature", query = "", `_intent` = "catalog"))
ok(
  ut_cmp_identical(r$ok, TRUE) &&
    ut_cmp_identical(r$value@extra$data$feature$category_count > 3L, TRUE),
  "search_annotations with an empty query returns the annotation catalog"
)
r <- dispatch("search_annotations",
              list(space = "feature", query = "Gene.name",
                   max_results = 5L, `_intent` = "find the name column"))
ok(
  ut_cmp_identical(r$ok, TRUE) &&
    ut_cmp_identical(r$value@extra$data$matching_column_count >= 1L, TRUE),
  "search_annotations finds matching annotation columns on the demo data"
)

## 6.3 summarize_annotation ------------------------------------------------

r <- dispatch("summarize_annotation",
              list(space = "sample", column = "General|All|Origin",
                   `_intent` = "summarize origin groups"))
ok(
  ut_cmp_identical(r$ok, TRUE) &&
    ut_cmp_identical(r$value@extra$data$non_missing, 60L),
  "summarize_annotation summarizes a real sample column (60 observed)"
)

## 6.4 set_omics_viewer_state (write) --------------------------------------

r <- dispatch("set_omics_viewer_state",
              list(data_space_tab = "Sample table",
                   analysis_space_tab = "ORA",
                   features = utils::head(rownames(fd), 3L),
                   samples = utils::head(rownames(pd), 2L),
                   `_intent` = "switch tabs and selections"))
ok(ut_cmp_identical(r$ok, TRUE), "set_omics_viewer_state dispatches ok")
ok(
  ut_cmp_identical(
    store_read(deps_live$store, "dataspace.active_tab")[[1]], "Sample table") &&
    ut_cmp_identical(
      store_read(deps_live$store, "resultspace.analyst_tab")[[1]], "ORA"),
  "set_omics_viewer_state changed both tab keys through the store"
)
ok(
  ut_cmp_identical(length(deps_live$sel_log$feature), 3L) &&
    ut_cmp_identical(length(deps_live$sel_log$sample), 2L),
  "set_omics_viewer_state drove the selection bus with the validated ids"
)

## 6.5 set_scatter_view (write) --------------------------------------------

r <- dispatch("set_scatter_view",
              list(space = "feature",
                   x_axis = "ttest|RE_vs_ME|mean.diff",
                   y_axis = "ttest|RE_vs_ME|log.fdr",
                   `_intent` = "custom volcano axes"))
ok(ut_cmp_identical(r$ok, TRUE), "set_scatter_view dispatches ok")
ok(
  ut_cmp_identical(
    store_read(deps_live$store, "dataspace.feature_space.x_variable")[[1]],
    "mean.diff") &&
    ut_cmp_identical(
      store_read(deps_live$store, "dataspace.feature_space.y_variable")[[1]],
      "log.fdr") &&
    ut_cmp_identical(
      store_read(deps_live$store, "dataspace.active_tab")[[1]], "Feature"),
  "set_scatter_view wrote the validated axis cascade into the store"
)
ok(
  ut_cmp_identical(
    store_snapshot(deps_live$store)$epochs[["dataspace.feature_space.x_variable"]],
    2L),
  "the axis write bumped the key's epoch (seed epoch 1 -> 2)"
)
ok(
  ut_cmp_identical(r$value@extra$data$mode, "custom"),
  "set_scatter_view reports the custom axis mode"
)

## 6.6 find_controls -------------------------------------------------------

r <- dispatch("find_controls",
              list(query = "scatter color", `_intent` = "find color controls"))
ok(
  ut_cmp_identical(r$ok, TRUE) &&
    ut_cmp_identical(r$value@extra$data$match_count >= 1L, TRUE),
  "find_controls(query) matches live capability records"
)
r <- dispatch("find_controls",
              list(prefix = "dataspace.expr_heatmap", `_intent` = "list heatmap"))
ok(
  ut_cmp_identical(r$ok, TRUE) &&
    ut_cmp_identical(r$value@extra$data$widget_count >= 1L, TRUE),
  "find_controls(prefix) lists raw widgets under a canonical prefix"
)

## 6.7 describe_control ----------------------------------------------------

r <- dispatch("describe_control",
              list(id = "dataspace.feature_space.x_variable",
                   `_intent` = "describe the x axis control"))
ok(
  ut_cmp_identical(r$ok, TRUE) &&
    ut_cmp_identical(r$value@extra$data$kind, "select_cascaded") &&
    ut_cmp_identical(isTRUE(r$value@extra$data$current_value %in%
                              c("PC1(10.5%)", "mean.diff")), TRUE),
  "describe_control returns the live widget record with its current value"
)

## 6.8 set_widgets (write) -------------------------------------------------

r <- dispatch("set_widgets",
              list(patch = '{"dataspace.expr_heatmap.heatmap_colors": "RdGy"}',
                   `_intent` = "switch the palette"))
ok(ut_cmp_identical(r$ok, TRUE), "set_widgets dispatches ok")
ok(
  ut_cmp_identical(
    store_read(deps_live$store, "dataspace.expr_heatmap.heatmap_colors")[[1]],
    "RdGy"),
  "set_widgets applied the canonical widget patch to the store"
)

## 6.9 run_control ---------------------------------------------------------

r <- dispatch("run_control",
              list(action_id = "action:dataspace.tab_feature.clear",
                   `_intent` = "clear the table selection"))
ok(
  ut_cmp_identical(r$ok, TRUE) &&
    ut_cmp_identical(r$value@extra$data$ran, "dataspace.tab_feature.clear") &&
    ut_cmp_identical(r$value@extra$data$epoch, 1L),
  "run_control ran the registered action and recorded its epoch"
)
r <- dispatch("run_control",
              list(action_id = "snapshot_save", `_intent` = "try a save"))
ok(
  ut_cmp_identical(r$ok, TRUE) &&
    ut_cmp_identical(r$value@extra$data$excluded_by_policy, TRUE) &&
    grepl("permanently excluded", r$value@extra$data$error, fixed = TRUE),
  "run_control returns the settled-policy rejection receipt for snapshot lifecycle ids"
)

## 6.10 set_enrichment_parameters (write) ----------------------------------

r <- dispatch("set_enrichment_parameters",
              list(method = "ora", collapse = "ttest|RE_vs_ME|log.fdr",
                   selected_pathway = "gs_setA",
                   `_intent` = "configure ORA"))
ok(ut_cmp_identical(r$ok, TRUE), "set_enrichment_parameters dispatches ok")
ok(
  ut_cmp_identical(
    store_read(deps_live$store, "resultspace.ora.xax_analysis")[[1]], "ttest") &&
    ut_cmp_identical(
      store_read(deps_live$store, "resultspace.ora.xax_variable")[[1]],
      "log.fdr") &&
    ut_cmp_identical(
      store_read(deps_live$store, "resultspace.ora.selected_row")[[1]],
      "gs_setA") &&
    ut_cmp_identical(
      store_read(deps_live$store, "resultspace.analyst_tab")[[1]], "ORA"),
  "set_enrichment_parameters wrote the collapse cascade, pathway row and panel tab"
)

## 6.11 set_table_view (write) ---------------------------------------------

r <- dispatch("set_table_view",
              list(table = "feature_table",
                   columns = list("General|All|Gene.name",
                                  "PCA|All|PC1(10.5%)",
                                  "ttest|RE_vs_ME|mean.diff"),
                   multi_selection = TRUE,
                   column_filters = "{\"General|All|Gene.name\": \"MT\"}",
                   page = 2L,
                   `_intent` = "configure the feature table"))
ok(ut_cmp_identical(r$ok, TRUE), "set_table_view dispatches ok")
vals <- store_read(deps_live$store,
                   c("dataspace.tab_feature.columns",
                     "dataspace.tab_feature.multi_selection",
                     "dataspace.tab_feature.column_filters",
                     "dataspace.tab_feature.page",
                     "dataspace.active_tab"))
ok(
  ut_cmp_identical(
    vals[["dataspace.tab_feature.columns"]],
    c("General|All|Gene.name", "PCA|All|PC1(10.5%)",
      "ttest|RE_vs_ME|mean.diff")) &&
    ut_cmp_identical(vals[["dataspace.tab_feature.multi_selection"]], TRUE) &&
    ut_cmp_identical(
      vals[["dataspace.tab_feature.column_filters"]],
      c("General|All|Gene.name" = "MT")) &&
    ut_cmp_identical(vals[["dataspace.tab_feature.page"]], 2L) &&
    ut_cmp_identical(vals[["dataspace.active_tab"]], "Feature table"),
  "set_table_view applied columns, selection mode, filters, page and tab switch"
)

## 6.12 create_figure / get_figure / update_figure -------------------------

r <- dispatch("create_figure",
              list(template = "volcano",
                   x = "ttest|RE_vs_ME|mean.diff",
                   y = "ttest|RE_vs_ME|log.fdr",
                   label_top_n = 5L,
                   title = "Channel volcano",
                   `_intent` = "volcano of the demo comparison"))
ok(ut_cmp_identical(r$ok, TRUE), "create_figure dispatches ok")
fig <- r$value@extra$data
ok(
  ut_cmp_identical(fig$figure_id, "fig_1") &&
    ut_cmp_identical(fig$row_count, 2702L) &&
    ut_cmp_identical(fig$full_png_bytes > 1000L, TRUE),
  "create_figure rendered a real volcano over all 2702 features (PNG encoded)"
)

r <- dispatch("get_figure", list(figure_id = "fig_1", `_intent` = "read back"))
ok(
  ut_cmp_identical(r$ok, TRUE) &&
    ut_cmp_identical(r$value@extra$data$spec$labels$title, "Channel volcano"),
  "get_figure returns the stored compact spec"
)

r <- dispatch("update_figure",
              list(figure_id = "fig_1",
                   changes = list(labels = list(title = "Channel revision")),
                   `_intent` = "retitle the figure"))
ok(ut_cmp_identical(r$ok, TRUE), "update_figure dispatches ok")
rev <- r$value@extra$data
ok(
  ut_cmp_identical(rev$figure_id, "fig_2") &&
    ut_cmp_identical(rev$parent_figure_id, "fig_1") &&
    ut_cmp_identical(rev$spec$labels$title, "Channel revision"),
  "update_figure produced a new revision with the patched title"
)
ok(
  ut_cmp_identical(length(shiny::isolate(figures())), 2L) &&
    ut_cmp_identical(shiny::isolate(figure_counter()), 2L),
  "the figure registry and counter tracked both revisions"
)

## 6.13 channel boundary on the live surface -------------------------------

r <- dispatch("set_scatter_view",
              list(space = "feature", spa = "x", `_intent` = "typo probe"))
ok(
  ut_cmp_identical(r$ok, FALSE) &&
    grepl("Unknown argument 'spa'", r$error, fixed = TRUE),
  "the dispatcher rejects non-exact argument names on the live registry (I9)"
)

# ==========================================================================
# [7] Equivalence: the existing test-hook pathway vs the channel
#
#     Baseline: two identically-seeded dep sets (verified before the
#     ops). Hook pathway = the exact calls module_agentTestHooks.R's
#     run_hook makes (apply_scatter_view / agent_widget_apply /
#     apply_enrichment / apply_table_view / apply_state). Channel
#     pathway = the same logical op dispatched through the ToolDef.
#     Identical store snapshots (values AND epochs) prove identical
#     effects. A final round drives the REAL agent_test_hooks_module
#     through shiny::testServer (JSON payload -> agent_args_sanitize ->
#     apply) to cover the hook's own plumbing end to end.
# ==========================================================================

domain2 <- shiny:::MockShinySession$new()
figures2 <- reactiveVal(list())
figure_counter2 <- reactiveVal(0L)
figure_directory2 <- file.path(tempdir(), "agent-channel-equiv")
dir.create(figure_directory2, showWarnings = FALSE, recursive = TRUE)

deps_hook <- make_deps()
deps_chan <- make_deps()
registry2 <- agent_tool_registry(
  state = deps_chan$state,
  feature_data = feature_data,
  sample_data = sample_data,
  expression_data = expression_data,
  selected_features = selected_features,
  selected_samples = selected_samples,
  apply_state = deps_chan$apply_state,
  apply_scatter_view = deps_chan$apply_scatter_view,
  apply_enrichment = deps_chan$apply_enrichment,
  apply_table_view = deps_chan$apply_table_view,
  store = deps_chan$store,
  figures = figures2,
  figure_counter = figure_counter2,
  figure_directory = figure_directory2,
  session_domain = domain2,
  logger = agent_logger_new(domain2))
dispatch2 <- function(name, args)
  agent_channel_dispatch(registry2$tools, name, args)

equiv <- function(label) {
  snap_hook <- store_snapshot(deps_hook$store)
  snap_chan <- store_snapshot(deps_chan$store)
  ok(
    ut_cmp_identical(snap_hook$values, snap_chan$values) &&
      ut_cmp_identical(snap_hook$epochs, snap_chan$epochs) &&
      ut_cmp_identical(snap_hook$actions, snap_chan$actions),
    paste(label, "- hook pathway and channel dispatcher leave the stores identical")
  )
}

ok(
  ut_cmp_identical(
    identical(store_snapshot(deps_hook$store), store_snapshot(deps_chan$store)),
    TRUE),
  "equivalence stores start identically seeded"
)

## 7.1 set_scatter_view -----------------------------------------------------

deps_hook$apply_scatter_view(
  space = "feature", x_axis = "ttest|RE_vs_ME|mean.diff",
  y_axis = "ttest|RE_vs_ME|log.fdr")
r <- dispatch2("set_scatter_view",
               list(space = "feature",
                    x_axis = "ttest|RE_vs_ME|mean.diff",
                    y_axis = "ttest|RE_vs_ME|log.fdr",
                    `_intent` = "same op via the channel"))
ok(ut_cmp_identical(r$ok, TRUE), "equivalence set_scatter_view: channel ok")
equiv("set_scatter_view (custom axes)")

## 7.2 set_widgets ----------------------------------------------------------

patch_json <- '{"dataspace.expr_heatmap.heatmap_colors": "RdGy"}'
shiny::isolate(agent_widget_apply(deps_hook$store, patch_json))
r <- dispatch2("set_widgets", list(patch = patch_json,
                                   `_intent` = "same op via the channel"))
ok(ut_cmp_identical(r$ok, TRUE), "equivalence set_widgets: channel ok")
equiv("set_widgets")

## 7.3 set_table_view -------------------------------------------------------

tv <- list(table = "feature_table",
           columns = list("General|All|Gene.name", "PCA|All|PC1(10.5%)",
                          "ttest|RE_vs_ME|mean.diff"),
           multi_selection = TRUE,
           column_filters = "{\"General|All|Gene.name\": \"MT\"}",
           page = 2L)
deps_hook$apply_table_view(tv)
r <- dispatch2("set_table_view", c(tv, list(`_intent` = "same op via the channel")))
ok(ut_cmp_identical(r$ok, TRUE), "equivalence set_table_view: channel ok")
equiv("set_table_view")

## 7.4 set_enrichment_parameters -------------------------------------------

en <- list(method = "ora", collapse = "ttest|RE_vs_ME|log.fdr",
           selected_pathway = "gs_setA")
deps_hook$apply_enrichment(en)
r <- dispatch2("set_enrichment_parameters",
               c(en, list(`_intent` = "same op via the channel")))
ok(ut_cmp_identical(r$ok, TRUE), "equivalence set_enrichment_parameters: channel ok")
equiv("set_enrichment_parameters")

## 7.5 set_omics_viewer_state ----------------------------------------------

st <- list(data_space_tab = "Sample table", analysis_space_tab = "ORA",
           features = utils::head(rownames(fd), 3L),
           samples = utils::head(rownames(pd), 2L))
deps_hook$apply_state(st)
r <- dispatch2("set_omics_viewer_state",
               c(st, list(`_intent` = "same op via the channel")))
ok(ut_cmp_identical(r$ok, TRUE), "equivalence set_omics_viewer_state: channel ok")
equiv("set_omics_viewer_state")
ok(
  ut_cmp_identical(deps_hook$sel_log$feature, deps_chan$sel_log$feature) &&
    ut_cmp_identical(deps_hook$sel_log$sample, deps_chan$sel_log$sample),
  "both pathways drove the selection bus with identical ids"
)

## 7.6 full hook pipeline: the REAL module via testServer -------------------

qk_payload <- '{"space": "feature", "quick_view_id": "volcano_RE_vs_LE"}'
r <- dispatch2("set_scatter_view",
               list(space = "feature", quick_view_id = "volcano_RE_vs_LE",
                    `_intent` = "same op via the channel"))
ok(
  ut_cmp_identical(r$ok, TRUE) &&
    ut_cmp_identical(r$value@extra$data$mode, "quick") &&
    ut_cmp_identical(r$value@extra$data$quick_view_id, "volcano_RE_vs_LE"),
  "equivalence quick-view scatter: channel ok (quick mode resolved)"
)
hook_output <- NULL
shiny::testServer(
  agent_test_hooks_module,
  args = list(
    apply_state = deps_hook$apply_state,
    apply_scatter_view = deps_hook$apply_scatter_view,
    state = deps_hook$state,
    store = deps_hook$store,
    apply_enrichment = deps_hook$apply_enrichment,
    apply_table_view = deps_hook$apply_table_view,
    assistant = NULL),
  expr = {
    session$setInputs(op = "scatter", payload = qk_payload)
    session$setInputs(run = 1L)
    hook_output <<- output$result
  })
ok(
  !grepl("hook_error", hook_output, fixed = TRUE) &&
    grepl("volcano_RE_vs_LE", hook_output, fixed = TRUE),
  "the real agent_test_hooks_module ran the same op without error"
)
equiv("quick-view scatter via the real hook module (testServer)")
