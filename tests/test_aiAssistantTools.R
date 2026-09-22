library(omicsViewer)
library(unittest, quietly = TRUE)

if (!requireNamespace("ellmer", quietly = TRUE) ||
    !requireNamespace("shinychat", quietly = TRUE)) {
  ok(TRUE, "AI assistant tool tests skipped because optional packages are unavailable")
  quit(save = "no", status = 0)
}

fd <- data.frame(
  score = c(1, 2, 3),
  category = c("kinase", "phosphatase", "kinase"),
  row.names = c("Gene1", "Gene2", "Gene3"),
  check.names = FALSE
)
pd <- data.frame(
  group = c("WT", "WT", "KO", "KO"),
  row.names = c("S1", "S2", "S3", "S4")
)

compact_state <- list(
  dataset = list(
    id = "demo.RDS",
    class = "ExpressionSet",
    dimensions = c(features = 3L, samples = 4L)
  ),
  active_tabs = list(data_space = "Feature", analysis_space = "Feature"),
  selection = list(
    features = list(count = 1L, ids = "Gene1", truncated = FALSE),
    samples = list(count = 0L, ids = character(), truncated = FALSE)
  ),
  available_tabs = list(
    data_space = c("Feature", "Sample"),
    analysis_space = c("Feature", "ORA")
  ),
  annotations = list(feature = list(rows = 3L), sample = list(rows = 4L)),
  quick_views = list(feature = list(), sample = list()),
  panels = list(data_space = list(), result_space = list()),
  state_policy = "widget-only"
)

captured <- list()
apply_state <- function(update) {
  captured$state <<- update
  list(
    data_space_tab = "Sample",
    analysis_space_tab = NULL,
    feature_count = 2L,
    sample_count = 0L,
    example_features = c("Gene1", "Gene3"),
    example_samples = character()
  )
}
apply_scatter_view <- function(space, quick_view_id = NULL,
                              x_axis = NULL, y_axis = NULL) {
  captured$scatter <<- list(
    space = space,
    quick_view_id = quick_view_id,
    x_axis = x_axis,
    y_axis = y_axis
  )
  list(space = space, mode = "custom", quick_view_id = NULL,
       x_axis = x_axis, y_axis = y_axis)
}

# canonical widget store for the S3 generic tier
widget_store_new <- omicsViewer:::widget_store_new
widget_binding <- omicsViewer:::widget_binding
store_register <- omicsViewer:::store_register
test_store <- widget_store_new()
store_register(test_store,
  widget_binding("dataspace.expr_heatmap.heatmap_colors", "select",
    label = "Heatmap color panel",
    help = "Diverging palette for the heatmap",
    values = c("BrBG", "PiYG", "RdBu", "RdGy", "RdYlBu")),
  widget_binding("dataspace.expr_heatmap.margin_bottom", "integer",
    label = "Bottom margin", help = "Bottom plot margin in lines",
    min = 1L, max = 20L))

# WP1: the state bridge is a builder function called as state(sections);
# mirror the real app wiring by building through agent_compact_state so the
# tool wrapper and the progressive-disclosure helper are tested together.
state_builder <- function(sections = NULL) {
  omicsViewer:::agent_compact_state(
    state = omicsViewer:::build_app_state(
      dataset = NULL,
      dataset_id = "demo.RDS",
      data_status = list(eset_active_tab = "Feature"),
      result_status = list(analyst_active_tab = "Feature"),
      selected_features = "Gene1",
      selected_samples = character(),
      label = "unit test"
    ),
    annotations = compact_state$annotations,
    quick_views = NULL,
    available_tabs = compact_state$available_tabs,
    sections = sections
  )
}

shiny::testServer(
  omicsViewer:::ai_assistant_module,
  args = list(
    state = state_builder,
    state_available = shiny::reactive(TRUE),
    feature_data = shiny::reactive(fd),
    sample_data = shiny::reactive(pd),
    expression_data = shiny::reactive(matrix(
      1:12, nrow = 3, dimnames = list(c("Gene1", "Gene2", "Gene3"), rownames(pd))
    )),
    selected_features = shiny::reactive(c("Gene1", "Gene2", "Gene3")),
    selected_samples = shiny::reactive(character()),
    apply_state = apply_state,
    apply_scatter_view = apply_scatter_view,
    store = test_store
  ),
  expr = {
    tools <- chat_object$client$get_tools()
    expected_tools <- c(
      "get_omics_viewer_state", "search_annotations", "summarize_annotation",
      "set_omics_viewer_state", "set_scatter_view", "create_figure", "update_figure",
      "list_widgets", "get_widget", "set_widgets"
    )
    ok(
      ut_cmp_identical(sort(names(tools)), sort(expected_tools)),
      "assistant exposes only the intended tool allowlist"
    )

    state_result <- tools$get_omics_viewer_state(`_intent` = "unit test")
    ok(
      ut_cmp_identical(state_result@value$dataset$id, "demo.RDS"),
      "state tool returns the compact current state"
    )
    ok(
      ut_cmp_identical(state_result@value$available_sections,
                       c("annotations", "quick_views", "panels", "figure_grammar")),
      "state tool overview advertises the section menu"
    )
    ok(
      is.null(state_result@value$annotations) && is.null(state_result@value$panels),
      "state tool omits full-detail sections by default"
    )
    state_section <- tools$get_omics_viewer_state(
      sections = c("annotations"), `_intent` = "unit test"
    )
    ok(
      ut_cmp_identical(
        state_section@value$annotations,
        list(feature = list(rows = 3L), sample = list(rows = 4L))
      ),
      "state tool returns requested sections in full"
    )
    search_result <- tools$search_annotations(
      space = "feature", query = "gene1", max_results = 2L,
      `_intent` = "unit test"
    )
    ok(
      ut_cmp_identical(search_result@value$matching_ids, "Gene1"),
      "search tool returns bounded matching IDs"
    )
    summary_result <- tools$summarize_annotation(
      space = "sample", column = "group", max_values = 2L,
      `_intent` = "unit test"
    )
    ok(
      ut_cmp_identical(summary_result@value$value_counts$WT, 2L),
      "summary tool returns categorical counts"
    )
    update_result <- tools$set_omics_viewer_state(
      data_space_tab = "Sample",
      features = c("Gene1", "Gene3"),
      samples = character(),
      `_intent` = "unit test"
    )
    ok(
      ut_cmp_identical(update_result@value$feature_count, 2L),
      "state-update tool returns a bounded receipt"
    )
    scatter_result <- tools$set_scatter_view(
      space = "feature",
      x_axis = "score",
      y_axis = "category",
      `_intent` = "unit test"
    )
    ok(
      ut_cmp_identical(scatter_result@value$x_axis, "score"),
      "scatter tool returns the validated axis update"
    )

    figure_spec <- list(
      data_source = "expression",
      layers = list(list(
        geom = "boxplot",
        x = "sample__group",
        y = "__expression__",
        fill = "sample__group",
        params = list(alpha = 0.65)
      )),
      labels = list(title = "Expression by group"),
      theme = "minimal",
      palette = "colorblind"
    )
    figure_result <- tools$create_figure(
      spec = figure_spec,
      `_intent` = "unit test"
    )
    ok(
      ut_cmp_identical(figure_result@value$data_source, "expression"),
      "figure tool renders expression data"
    )
    ok(
      ut_cmp_identical(figure_result@value$row_count, 12L),
      "figure tool reports bounded plotting rows"
    )
    ok(
      ut_cmp_identical(
        startsWith(figure_result@extra$display$html$children[[1]]$attribs$src, "data:image/png;base64,"),
        TRUE
      ),
      "figure tool returns an embedded chat preview"
    )

    updated_spec <- figure_spec
    updated_spec$labels$title <- "Expression by group, revised"
    updated_spec$theme <- "classic"
    updated_figure <- tools$update_figure(
      figure_id = figure_result@value$figure_id,
      spec = updated_spec,
      `_intent` = "unit test"
    )
    ok(
      ut_cmp_identical(
        updated_figure@value$parent_figure_id,
        figure_result@value$figure_id
      ),
      "figure update records its parent figure"
    )

    # ---- WP3: spec round-trip ---------------------------------------
    ok(
      !is.null(figure_result@value$spec) &&
        ut_cmp_identical(figure_result@value$spec$data_source, "expression"),
      "figure result carries the full normalized spec"
    )
    ok(
      ut_cmp_identical(figure_result@value$spec$layers[[1]]$y, "__expression__") &&
        ut_cmp_identical("mappings" %in% names(figure_result@value$spec$layers[[1]]), FALSE) &&
        ut_cmp_identical(figure_result@value$spec$layers[[1]]$params$alpha, 0.65) &&
        ut_cmp_identical(figure_result@value$spec$samples, rownames(pd)),
      "result spec is schema-shaped (flat aesthetics, defaulted params, resolved samples)"
    )
    roundtrip_spec <- figure_result@value$spec
    roundtrip_spec$labels$title <- "Expression by group, round-tripped"
    roundtrip_figure <- tools$update_figure(
      figure_id = updated_figure@value$figure_id,
      spec = roundtrip_spec,
      `_intent` = "unit test"
    )
    expected_roundtrip <- figure_result@value$spec
    expected_roundtrip$labels$title <- "Expression by group, round-tripped"
    ok(
      ut_cmp_identical(roundtrip_figure@value$spec, expected_roundtrip),
      "echoed spec re-submitted through update normalizes identically"
    )

    # ---- S3 generic widget tier ---------------------------------------
    widgets_result <- tools$list_widgets(section = "dataspace", `_intent` = "unit test")
    ok(
      ut_cmp_identical(widgets_result@value$widget_count, 2L) &&
        ut_cmp_identical(widgets_result@value$widgets$dataspace.expr_heatmap.heatmap_colors$kind,
                         "select"),
      "list_widgets returns section-filtered registry records"
    )
    widget_result <- tools$get_widget(
      id = "dataspace.expr_heatmap.heatmap_colors", `_intent` = "unit test")
    ok(
      ut_cmp_identical(widget_result@value$id, "dataspace.expr_heatmap.heatmap_colors") &&
        "RdGy" %in% widget_result@value$allowed_values,
      "get_widget describes one widget with allowed values"
    )
    set_widgets_result <- tools$set_widgets(
      patch = "{\"dataspace.expr_heatmap.heatmap_colors\": \"RdGy\", \"dataspace.expr_heatmap.margin_bottom\": 9}",
      `_intent` = "unit test"
    )
    ok(
      setequal(set_widgets_result@value$applied,
               c("dataspace.expr_heatmap.heatmap_colors",
                 "dataspace.expr_heatmap.margin_bottom")) &&
        ut_cmp_identical(
          set_widgets_result@value$applied_values$dataspace.expr_heatmap.margin_bottom,
          9L),
      "set_widgets applies a JSON patch with typed coercion"
    )
    rejected_result <- tools$set_widgets(
      patch = "{\"dataspace.expr_heatmap.heatmap_colors\": \"Spectral\"}",
      `_intent` = "unit test"
    )
    ok(
      length(rejected_result@value$applied) == 0L &&
        length(rejected_result@value$rejected) == 1L &&
        grepl("RdGy|RdYlBu", rejected_result@value$rejected[[1]]$reason),
      "set_widgets rejects invalid values per key with suggestions"
    )

    # ---- provider sentinel sweep ---------------------------------------
    # glm flash serializes omitted optionals as literal "null"/"{}"/"[]"
    # strings. Guards the whole boundary class: no tool whose normalization
    # happens inside the tool body may error or degrade on sentinels.
    swept_listing <- tools$list_widgets(
      section = "null", `_intent` = "sentinel sweep")
    ok(
      ut_cmp_identical(swept_listing@value$widget_count, 2L) &&
        is.null(swept_listing@value$section),
      "list_widgets treats a sentinel section as omitted"
    )
    swept_id_error <- tryCatch(
      tools$get_widget(id = "{}", `_intent` = "sentinel sweep"),
      error = function(e) conditionMessage(e)
    )
    ok(
      grepl("Unknown or not user-editable widget id", swept_id_error),
      "get_widget reports a sentinel id as unknown (required-arg honesty)"
    )
    swept_figures <- tools$create_figure(
      spec = list(
        data_source = "expression",
        layers = list(list(
          geom = "boxplot", x = "sample__group", y = "__expression__",
          fill = "sample__group"
        )),
        theme = "null", palette = "{}", facet_by = "null",
        labels = list(title = "null")
      ),
      `_intent` = "sentinel sweep"
    )
    ok(
      ut_cmp_identical(swept_figures@value$theme, "minimal") &&
        ut_cmp_identical(swept_figures@value$palette, "default") &&
        is.null(swept_figures@value$labels$title),
      "figure tool applies defaults for sentinel optionals"
    )
  }
)

ok(
  ut_cmp_identical(captured$state$features, c("Gene1", "Gene3")),
  "state tool forwards only its typed proposal"
)
ok(
  ut_cmp_identical(captured$scatter$space, "feature"),
  "scatter tool forwards only its typed proposal"
)
