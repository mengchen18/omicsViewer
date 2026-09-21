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

shiny::testServer(
  omicsViewer:::ai_assistant_module,
  args = list(
    state = shiny::reactive(compact_state),
    state_available = shiny::reactive(TRUE),
    feature_data = shiny::reactive(fd),
    sample_data = shiny::reactive(pd),
    expression_data = shiny::reactive(matrix(
      1:12, nrow = 3, dimnames = list(c("Gene1", "Gene2", "Gene3"), rownames(pd))
    )),
    selected_features = shiny::reactive(c("Gene1", "Gene2", "Gene3")),
    selected_samples = shiny::reactive(character()),
    apply_state = apply_state,
    apply_scatter_view = apply_scatter_view
  ),
  expr = {
    tools <- chat_object$client$get_tools()
    expected_tools <- c(
      "get_omics_viewer_state", "search_annotations", "summarize_annotation",
      "set_omics_viewer_state", "set_scatter_view", "create_figure", "update_figure"
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
