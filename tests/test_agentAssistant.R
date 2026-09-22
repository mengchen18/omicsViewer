library(omicsViewer)
library(unittest, quietly = TRUE)

agent_environment_config <- omicsViewer:::agent_environment_config
agent_request_limit <- omicsViewer:::agent_request_limit
agent_validate_provider_config <- omicsViewer:::agent_validate_provider_config
agent_annotation_catalog <- omicsViewer:::agent_annotation_catalog
agent_search_annotations <- omicsViewer:::agent_search_annotations
agent_summarize_annotation <- omicsViewer:::agent_summarize_annotation
agent_normalize_state_update <- omicsViewer:::agent_normalize_state_update
agent_normalize_scatter_view <- omicsViewer:::agent_normalize_scatter_view
.agent_suggest <- omicsViewer:::`.agent_suggest`
agent_normalize_figure_spec <- omicsViewer:::agent_normalize_figure_spec
agent_build_figure_data <- omicsViewer:::agent_build_figure_data
agent_build_figure_plot <- omicsViewer:::agent_build_figure_plot

fd <- data.frame(
  `ttest|A_vs_B|log.fdr` = c(0.1, 0.4, 0.8),
  category = c("kinase", "phosphatase", "kinase"),
  check.names = FALSE,
  row.names = c("Gene1", "Gene2", "Gene3")
)
pd <- data.frame(
  group = c("WT", "WT", "KO", "KO"),
  batch = c(1, 2, 1, 2),
  row.names = c("S1", "S2", "S3", "S4")
)

config <- agent_validate_provider_config(
  provider = "openai", model = " test-model ", api_key = " secret ",
  base_url = " https://example.org/v1 "
)
ok(ut_cmp_identical(config$provider, "openai"), "provider config is normalized")
ok(ut_cmp_identical(config$model, "test-model"), "provider model is trimmed")
ok(ut_cmp_identical(config$base_url, "https://example.org/v1"), "provider base URL is trimmed")
ok(ut_cmp_identical(config$configured, TRUE), "provider config reports a key")
ok(
  ut_cmp_error(
    agent_validate_provider_config(provider = "openai", base_url = "ftp://example.org"),
    "API base URL must use HTTPS, or local HTTP for this machine only."
  ),
  "remote non-HTTPS provider base URLs are rejected"
)
ok(
  ut_cmp_identical(
    agent_validate_provider_config("openai", base_url = "http://127.0.0.1:8080/v1")$base_url,
    "http://127.0.0.1:8080/v1"
  ),
  "local HTTP provider endpoints are allowed"
)

old_request_limit <- Sys.getenv("OMICSVIEWER_LLM_MAX_REQUESTS")
Sys.setenv(OMICSVIEWER_LLM_MAX_REQUESTS = "1000")
on.exit(Sys.setenv(
  OMICSVIEWER_LLM_MAX_REQUESTS = if (identical(old_request_limit, "")) "" else old_request_limit
), add = TRUE)
ok(ut_cmp_identical(omicsViewer:::agent_request_limit(), 200L), "request limit is capped")
Sys.setenv(OMICSVIEWER_LLM_MAX_REQUESTS = "not-a-number")
ok(ut_cmp_identical(omicsViewer:::agent_request_limit(), 40L), "invalid request limit falls back to 40")
Sys.setenv(
  OMICSVIEWER_LLM_MAX_REQUESTS = if (identical(old_request_limit, "")) "" else old_request_limit
)

catalog <- agent_annotation_catalog(fd, pd)
ok(ut_cmp_identical(catalog$feature$rows, 3L), "feature catalog row count")
ok(
  ut_cmp_identical(catalog$feature$columns[["ttest|A_vs_B|log.fdr"]]$type, "numeric"),
  "feature catalog column type"
)
ok(
  ut_cmp_identical(catalog$sample$columns[["group"]]$unique_count, 2L),
  "sample catalog cardinality"
)

search <- agent_search_annotations(
  space = "feature", query = "gene1", feature_data = fd, sample_data = pd
)
ok(ut_cmp_identical(search$matching_ids, "Gene1"), "annotation search is case-insensitive")
ok(
  ut_cmp_identical(search$value_matches, list()),
  "annotation search does not invent value matches"
)
value_search <- agent_search_annotations(
  space = "sample", query = "ko", feature_data = fd, sample_data = pd
)
ok(
  ut_cmp_identical(value_search$value_matches[[1]]$matching_rows, 2L),
  "annotation value search counts matches"
)
ok(
  ut_cmp_identical(value_search$matching_id_count, 0L),
  "annotation value search does not confuse values with IDs"
)

summary <- agent_summarize_annotation(
  space = "feature", column = "ttest|A_vs_B|log.fdr", feature_data = fd, sample_data = pd
)
ok(ut_cmp_identical(summary$type, "numeric"), "numeric annotation summary type")
ok(ut_cmp_identical(summary$non_missing, 3L), "numeric annotation summary count")
ok(ut_cmp_identical(summary$summary[[1]], 0.1), "numeric annotation summary minimum")
ok(
  ut_cmp_error(
    agent_summarize_annotation("feature", "missing", fd, pd),
    "Unknown feature annotation column: missing"
  ),
  "unknown annotation columns are rejected"
)

update <- agent_normalize_state_update(
  update = list(
    data_space_tab = "Sample",
    features = c("Gene3", "Gene1", "Gene1"),
    samples = character()
  ),
  data_tabs = c("Feature", "Sample"),
  analysis_tabs = c("Feature", "ORA"),
  feature_ids = rownames(fd),
  sample_ids = rownames(pd)
)
ok(ut_cmp_identical(update$data_space_tab, "Sample"), "assistant tab update is validated")
ok(ut_cmp_identical(update$features, c("Gene3", "Gene1")), "assistant feature IDs are deduplicated")
ok(ut_cmp_identical(update$samples, character()), "assistant can explicitly clear samples")
unknown_id_error <- tryCatch(
  agent_normalize_state_update(
    list(features = "Missing"), c("Feature"), c("ORA"), rownames(fd), rownames(pd)
  ),
  error = function(e) conditionMessage(e)
)
ok(
  ut_cmp_identical(grepl("Unknown feature ID", unknown_id_error, fixed = TRUE), TRUE),
  "unknown assistant feature IDs are rejected"
)

# provider sentinel artifacts (glm flash serializes omitted optionals as
# literal "null"/"{}" strings): they must be treated as absent, not rejected
sentinel_update <- agent_normalize_state_update(
  update = list(
    data_space_tab = "Heatmap",
    analysis_space_tab = "null",
    features = "null",
    samples = "[]"
  ),
  data_tabs = c("Feature", "Sample", "Heatmap"),
  analysis_tabs = c("Feature", "ORA"),
  feature_ids = rownames(fd),
  sample_ids = rownames(pd)
)
ok(
  ut_cmp_identical(sentinel_update, list(data_space_tab = "Heatmap")),
  "literal 'null' tab/selection sentinels are treated as omitted"
)
ok(
  ut_cmp_error(
    agent_normalize_scatter_view(
      space = "feature", x_axis = "{}", y_axis = "ttest|A_vs_B|log.fdr",
      feature_columns = colnames(fd)
    ),
    "requires either quick_view_id or both"
  ),
  "scatter-view string sentinels are treated as omitted axes"
)
ok(
  ut_cmp_error(
    agent_normalize_state_update(
      list(data_space_tab = "{}"),
      c("Feature", "Sample"), c("Feature"), rownames(fd), rownames(pd)
    ),
    "contains no changes"
  ),
  "sentinel-only updates still report no changes"
)

views <- data.frame(
  id = "volcano", label = "Volcano", x = "ttest|A_vs_B|mean.diff",
  y = "ttest|A_vs_B|log.fdr", description = "volcano", source = "auto",
  stringsAsFactors = FALSE
)
view <- agent_normalize_scatter_view(
  space = "feature", quick_view_id = "volcano", quick_views = list(feature = views)
)
ok(ut_cmp_identical(view$x_axis, "ttest|A_vs_B|mean.diff"), "quick-view X axis is resolved")
ok(ut_cmp_identical(view$mode, "quick"), "quick-view mode is preserved")
ok(
  ut_cmp_error(
    agent_normalize_scatter_view(
      space = "feature", x_axis = "missing", y_axis = "ttest|A_vs_B|log.fdr",
      feature_columns = colnames(fd)
    ),
    "Unknown feature X-axis annotation: missing"
  ),
  "unknown custom scatter axes are rejected"
)

# WP2: did-you-mean suggestions make validation errors self-correcting
ok(
  ut_cmp_identical(
    .agent_suggest("Featue", c("Feature", "Feature table", "Sample")),
    "Feature"
  ),
  "suggest returns the case/prefix closest tab"
)
ok(
  ut_cmp_identical(
    grepl("Closest matches: Feature", 
          tryCatch(agent_normalize_state_update(
            list(data_space_tab = "Featur"), 
            c("Feature", "Sample"), c("Feature"), c("Gene1"), c("S1")),
            error = function(e) conditionMessage(e)), fixed = TRUE),
    TRUE),
  "invalid tab error carries closest matches"
)
ok(
  ut_cmp_identical(
    grepl("Closest matches: Gene1",
          tryCatch(agent_normalize_state_update(
            list(features = c("Gene10")),
            c("Feature"), c("Feature"), rownames(fd), rownames(pd)),
            error = function(e) conditionMessage(e)), fixed = TRUE),
    TRUE),
  "invalid feature ID error carries closest match and search hint"
)
ok(
  ut_cmp_identical(
    grepl("search_annotations", 
          tryCatch(agent_normalize_state_update(
            list(features = c("Gene10")),
            c("Feature"), c("Feature"), rownames(fd), rownames(pd)),
            error = function(e) conditionMessage(e)), fixed = TRUE),
    TRUE),
  "ID errors point at search_annotations"
)
ok(
  ut_cmp_identical(
    grepl("Closest matches: ttest|A_vs_B|log.fdr",
          tryCatch(agent_normalize_scatter_view(
            space = "feature",
            x_axis = "ttest|A_vs_B|log.fdr", y_axis = "ttest|A_vs_B|log.fd",
            quick_views = list(feature = data.frame(
              id = "v1", label = "v1", x = "a|b|c", y = "d|e|f",
              description = "", source = "auto", stringsAsFactors = FALSE)),
            feature_columns = colnames(fd), sample_columns = colnames(pd)),
            error = function(e) conditionMessage(e)), fixed = TRUE),
    TRUE),
  "invalid axis error suggests the closest annotation column"
)
ok(
  ut_cmp_identical(
    grepl("naming convention",
          tryCatch(agent_normalize_scatter_view(
            space = "feature", x_axis = "log.fdr", y_axis = "ttest|A_vs_B|log.fdr",
            quick_views = NULL,
            feature_columns = colnames(fd), sample_columns = colnames(pd)),
            error = function(e) conditionMessage(e)), fixed = TRUE),
    TRUE),
  "malformed axis error still explains the naming convention"
)
zero_hit <- agent_search_annotations("feature", "logg.fdrr", fd, pd)
ok(
  ut_cmp_identical(
    "ttest|A_vs_B|log.fdr" %in% zero_hit$suggestions,
    TRUE),
  "zero-hit searches return column suggestions instead of a dead end"
)
ok(
  ut_cmp_identical(
    grepl("Closest matches: sample__group",
          tryCatch({
            spec <- agent_normalize_figure_spec(list(
              layers = list(list(geom = "point", x = "sample__groupz", y = "__expression__")),
              data_source = "expression", features = rownames(fd)[1], 
              samples = rownames(pd)[1]),
              fd, pd, matrix(1, dimnames = list(rownames(fd)[1], rownames(pd)[1])),
              rownames(fd)[1], rownames(pd)[1])
            fig_data <- agent_build_figure_data(spec, fd, pd, 
              matrix(1, dimnames = list(rownames(fd)[1], rownames(pd)[1])))
            plot_spec <- agent_normalize_figure_spec(list(
              layers = list(list(geom = "point", x = "sample__groupz", y = "__expression__")),
              data_source = "expression", features = rownames(fd)[1], 
              samples = rownames(pd)[1]),
              fd, pd, matrix(1, dimnames = list(rownames(fd)[1], rownames(pd)[1])),
              rownames(fd)[1], rownames(pd)[1])
            agent_build_figure_plot(fig_data, plot_spec)
          }, error = function(e) conditionMessage(e)), fixed = TRUE),
    TRUE),
  "figure mapping errors suggest closest data columns"
)

# providers that serialize omitted optional strings as the literal "null"
ok(
  ut_cmp_identical(
    agent_normalize_scatter_view(
      space = "feature", quick_view_id = "null",
      x_axis = "ttest|A_vs_B|log.fdr", y_axis = "ttest|A_vs_B|log.fdr",
      quick_views = NULL,
      feature_columns = colnames(fd), sample_columns = colnames(pd))$mode,
    "custom"),
  "literal \"null\" quick_view_id is treated as omitted"
)

# ---- WP1: sections / progressive disclosure on agent_compact_state -----
agent_compact_state <- omicsViewer:::agent_compact_state
agent_scatter_view_from_store <- omicsViewer:::agent_scatter_view_from_store
build_app_state <- omicsViewer:::build_app_state
widget_store_new <- omicsViewer:::widget_store_new
widget_store_child <- omicsViewer:::widget_store_child
widget_binding <- omicsViewer:::widget_binding
store_register <- omicsViewer:::store_register
store_apply <- omicsViewer:::store_apply

full_state <- build_app_state(
  dataset = NULL,
  dataset_id = "unit.RDS",
  data_status = list(eset_active_tab = "Feature"),
  result_status = list(analyst_active_tab = "Feature"),
  selected_features = paste0("G", 1:25),
  selected_samples = c("S1", "S2"),
  label = "unit"
)
qv <- list(feature = views, sample = NULL)
state_args <- list(
  annotations = catalog,
  quick_views = qv,
  available_tabs = list(
    data_space = c("Feature", "Sample"),
    analysis_space = c("Feature", "ORA")
  ),
  figure_grammar = list(geoms = "point")
)

overview <- do.call(agent_compact_state,
                    c(list(state = full_state), state_args))
ok(is.null(overview$annotations),
   "WP1 overview omits the annotation catalog")
ok(is.null(overview$panels), "WP1 overview omits panel state")
ok(is.null(overview$figure_grammar), "WP1 overview omits the figure grammar")
ok(is.null(overview$scatter_view),
   "WP1 overview has no scatter_view without a widget store")
ok(ut_cmp_identical(overview$selection$features$count, 25L),
   "WP1 overview selection keeps the full count")
ok(ut_cmp_identical(length(overview$selection$features$ids), 20L),
   "WP1 overview selection example IDs are capped at 20")
ok(ut_cmp_identical(overview$selection$features$truncated, TRUE),
   "WP1 overview selection reports truncation")
ok(
  ut_cmp_identical(
    overview$quick_views$feature[[1]], list(id = "volcano", label = "Volcano")
  ),
  "WP1 overview quick views carry id+label only"
)
ok(
  ut_cmp_identical(
    overview$available_sections,
    c("annotations", "quick_views", "panels", "figure_grammar")
  ),
  "WP1 overview advertises the section menu"
)
ok(ut_cmp_identical(overview$active_tabs$data_space, "Feature"),
   "WP1 overview keeps active tabs")

mixed <- do.call(agent_compact_state,
                 c(list(state = full_state, sections = "annotations"), state_args))
ok(ut_cmp_identical(mixed$annotations, catalog),
   "requested annotations section returns the full catalog")
ok(is.null(mixed$panels), "unrequested sections stay absent")
ok(
  ut_cmp_identical(mixed$quick_views$feature[[1]], list(id = "volcano", label = "Volcano")),
  "quick views stay compact unless the section is requested"
)

full <- do.call(
  agent_compact_state,
  c(list(state = full_state,
         sections = c("quick_views", "panels", "figure_grammar")),
    state_args)
)
ok(
  ut_cmp_identical(full$quick_views$feature[[1]]$x, "ttest|A_vs_B|mean.diff"),
  "requested quick_views section returns full records"
)
ok(ut_cmp_identical(full$panels$data_space$eset_active_tab, "Feature"),
   "requested panels section returns bounded panel state")
ok(ut_cmp_identical(full$figure_grammar, list(geoms = "point")),
   "requested figure_grammar section is passed through")
ok(is.null(full$annotations),
   "sections compose: unrequested annotations stay absent")
ok(
  ut_cmp_error(
    do.call(agent_compact_state,
            c(list(state = full_state, sections = "annotaions"), state_args)),
    "Unknown state section"
  ),
  "unknown sections are rejected"
)
ok(
  ut_cmp_identical(
    grepl(
      "Closest matches: annotations",
      tryCatch(
        do.call(agent_compact_state,
                c(list(state = full_state, sections = "annotaion"), state_args)),
        error = function(e) conditionMessage(e)),
      fixed = TRUE),
    TRUE),
  "unknown-section errors carry closest-match suggestions"
)
ok(
  ut_cmp_identical(
    do.call(agent_compact_state,
            c(list(state = full_state, sections = "null"), state_args))$selection,
    overview$selection
  ),
  "literal 'null' sections sentinel is treated as omitted"
)
ok(
  ut_cmp_identical(
    do.call(agent_compact_state,
            c(list(state = full_state, sections = list("panels")), state_args))$panels,
    full$panels
  ),
  "list-shaped sections (fromJSON) are accepted"
)

# widget-store scatter anchors on the overview ground floor
st <- widget_store_new()
fs <- widget_store_child(st, "dataspace.feature_space")
store_register(
  fs,
  widget_binding("x_analysis", "select", values = c("ttest", "PCA")),
  widget_binding("x_subset", "select", values = c("A_vs_B", "All")),
  widget_binding("x_variable", "select", values = c("log.fdr", "PC1")),
  widget_binding("y_analysis", "select", values = c("ttest", "PCA")),
  widget_binding("y_subset", "select", values = c("A_vs_B", "All")),
  widget_binding("y_variable", "select", values = c("log.fc", "PC2")),
  widget_binding("axis_mode", "enum", values = c("quick", "custom"))
)
store_apply(
  fs,
  list(x_analysis = "ttest", x_subset = "A_vs_B", x_variable = "log.fdr",
       y_analysis = "ttest", y_subset = "A_vs_B", y_variable = "log.fc",
       axis_mode = "quick"),
  origin = "restore"
)
sv <- agent_scatter_view_from_store(st)
ok(
  ut_cmp_identical(
    sv$feature$x,
    list(analysis = "ttest", subset = "A_vs_B", variable = "log.fdr",
         name = "ttest|A_vs_B|log.fdr")
  ),
  "scatter_view reads axis triples from the widget store"
)
ok(ut_cmp_identical(sv$feature$axis_mode, "quick"),
   "scatter_view carries the axis mode")
ok(is.null(sv$sample), "unset sample scatter is NULL")
ok(
  ut_cmp_identical(
    do.call(agent_compact_state,
            c(list(state = full_state, store = st), state_args))$scatter_view,
    sv
  ),
  "overview embeds the store-backed scatter_view"
)
st2 <- widget_store_new()
fs2 <- widget_store_child(st2, "dataspace.feature_space")
store_register(fs2, widget_binding("axis_mode", "enum", values = c("quick", "custom")))
store_apply(fs2, list(axis_mode = "quick"), origin = "restore")
ok(is.null(agent_scatter_view_from_store(st2)),
   "axis mode without axes yields no scatter_view block")
