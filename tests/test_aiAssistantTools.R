library(omicsViewer)
library(unittest, quietly = TRUE)

if (!requireNamespace("ellmer", quietly = TRUE) ||
    !requireNamespace("shinychat", quietly = TRUE)) {
  ok(TRUE, "AI assistant tool tests skipped because optional packages are unavailable")
  quit(save = "no", status = 0)
}

fd <- data.frame(
  score = c(1, 2, 3),
  logFdr = c(5, 4, 3),
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

# canonical widget store for the S3 generic tier + the WP8 semantic tools
widget_store_new <- omicsViewer:::widget_store_new
widget_binding <- omicsViewer:::widget_binding
store_register <- omicsViewer:::store_register
store_read <- omicsViewer:::store_read
store_apply <- omicsViewer:::store_apply
test_store <- widget_store_new()
store_register(test_store,
  widget_binding("dataspace.expr_heatmap.heatmap_colors", "select",
    label = "Heatmap color panel",
    help = "Diverging palette for the heatmap",
    values = c("BrBG", "PiYG", "RdBu", "RdGy", "RdYlBu")),
  widget_binding("dataspace.expr_heatmap.margin_bottom", "integer",
    label = "Bottom margin", help = "Bottom plot margin in lines",
    min = 1L, max = 20L),
  widget_binding("dataspace.active_tab", "navbar", label = "Data-space tab",
    help = "Active tab of the data-space navbar",
    values = c("Feature", "Feature table", "Sample", "Sample table", "Expression")),
  widget_binding("resultspace.analyst_tab", "navbar", label = "Analysis tab",
    help = "Active tab of the analysis navbar",
    values = c("Feature", "ORA", "fGSEA")),
  widget_binding("dataspace.tab_pheno.columns", "multi_select",
    label = "Shown columns", help = "Columns displayed in the sample table",
    min = 1L,
    choices_provider = function(v) c("group", "batch")),
  widget_binding("dataspace.tab_pheno.multi_selection", "boolean",
    label = "Multiple selection", help = "Allow multi-row selection"),
  widget_binding("dataspace.tab_pheno.page", "integer", label = "Table page",
    help = "Visible page", min = 1L),
  widget_binding("dataspace.tab_pheno.column_filters", "mapping",
    label = "Column filters", help = "Per-column search patterns",
    choices_provider = function(v) c("group", "batch")),
  widget_binding("resultspace.ora.xax_analysis", "select",
    label = "Collapse category", help = "Collapse analysis",
    choices_provider = function(v) c("ttest")),
  widget_binding("resultspace.ora.xax_subset", "select",
    label = "Collapse subcategory", help = "Collapse subset",
    depends_on = "resultspace.ora.xax_analysis",
    choices_provider = function(v) c("A_vs_B")),
  widget_binding("resultspace.ora.xax_variable", "select_cascaded",
    label = "Collapse variable", help = "Collapse variable",
    depends_on = c("resultspace.ora.xax_analysis", "resultspace.ora.xax_subset"),
    choices_provider = function(v) c("log.fdr", "md")),
  widget_binding("resultspace.ora.selected_row", "select",
    label = "Selected pathway", help = "Gene set selected in the results table",
    choices_provider = function(v) c("gs1", "gs2")))

# WP8 apply callbacks mirroring the L0 wiring (validate + store_apply)
fd_gs <- fd
attr(fd_gs, "GS") <- data.frame(
  featureId = c("Gene1", "Gene2", "Gene3"),
  gsId = c("gs1", "gs1", "gs2")
)
fd_gs$`ttest|A_vs_B|log.fdr` <- c(1, 2, 3)
fd_gs$`ttest|A_vs_B|md` <- c(-1, 0, 1)
apply_enrichment <- function(update) {
  validated <- omicsViewer:::agent_normalize_enrichment_update(update, fd_gs)
  receipt <- store_apply(test_store, validated$patch, origin = "agent",
                         strict = FALSE)
  list(method = validated$method, panel_tab = validated$tab,
       applied = receipt$applied, applied_values = receipt$diff,
       unchanged = receipt$skipped, rejected = receipt$rejected %||% list(),
       note = "test")
}
apply_table_view <- function(update) {
  validated <- omicsViewer:::agent_normalize_table_view_update(update)
  receipt <- store_apply(test_store, validated$patch, origin = "agent",
                         strict = FALSE)
  list(table = validated$table, panel_tab = validated$tab,
       applied = receipt$applied, applied_values = receipt$diff,
       unchanged = receipt$skipped, rejected = receipt$rejected %||% list(),
       note = "test")
}

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
    apply_enrichment = apply_enrichment,
    apply_table_view = apply_table_view,
    store = test_store
  ),
  expr = {
    tools <- chat_object$client$get_tools()
    expected_tools <- c(
      "get_omics_viewer_state", "search_annotations", "summarize_annotation",
      "set_omics_viewer_state", "set_scatter_view", "create_figure", "update_figure",
      "list_widgets", "get_widget", "set_widgets",
      "set_enrichment_parameters", "set_table_view",
      "search_ui_capabilities", "get_ui_capability"
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

    # ---- WP8 semantic tools (enrichment + table view) -----------------
    enrichment_result <- tools$set_enrichment_parameters(
      method = "ora",
      collapse = "ttest|A_vs_B|log.fdr",
      selected_pathway = "gs1",
      `_intent` = "unit test"
    )
    ok(
      ut_cmp_identical(enrichment_result@value$panel_tab, "ORA") &&
        "resultspace.ora.xax_variable" %in% enrichment_result@value$applied &&
        "resultspace.analyst_tab" %in% enrichment_result@value$applied,
      "enrichment tool applies the ranking cascade and opens the panel"
    )
    ok(
      ut_cmp_identical(
        omicsViewer:::store_read(test_store)$resultspace.ora.xax_variable,
        "log.fdr"),
      "enrichment tool lands in the widget store"
    )
    enr_reject <- tools$set_enrichment_parameters(
      method = "ora", selected_pathway = "gs9", `_intent` = "unit test"
    )
    ok(
      "resultspace.ora.selected_row" %in%
        vapply(enr_reject@value$rejected, function(r) r$id, character(1)) &&
        grepl("gs1|gs2", enr_reject@value$rejected[[1]]$reason),
      "unknown pathway rows are rejected per key with suggestions"
    )
    enr_error <- tryCatch(
      tools$set_enrichment_parameters(
        method = "ora", collapse = "ttest|A_vs_B|logg.fdrr",
        `_intent` = "unit test"
      ),
      error = function(e) conditionMessage(e)
    )
    ok(
      grepl("Unknown feature annotation", enr_error) &&
        grepl("log.fdr", enr_error),
      "enrichment tool suggests closest columns on typos"
    )

    table_result <- tools$set_table_view(
      table = "sample_table",
      columns = c("group", "batch"),
      column_filters = '{"group": "KO"}',
      page = 2L,
      `_intent` = "unit test"
    )
    stored_after_table <- omicsViewer:::store_read(test_store)
    ok(
      ut_cmp_identical(table_result@value$panel_tab, "Sample table") &&
        ut_cmp_identical(stored_after_table$dataspace.tab_pheno.column_filters,
                         c(group = "KO")) &&
        ut_cmp_identical(stored_after_table$dataspace.tab_pheno.page, 2L) &&
        ut_cmp_identical(stored_after_table$dataspace.active_tab, "Sample table"),
      "table-view tool lands columns, filters, page, and the tab switch"
    )
    table_reject <- tools$set_table_view(
      table = "sample_table", columns = c("group", "nope"),
      `_intent` = "unit test"
    )
    ok(
      "dataspace.tab_pheno.columns" %in%
        vapply(table_reject@value$rejected, function(r) r$id, character(1)) &&
        "dataspace.active_tab" %in%
          c(table_reject@value$applied, table_reject@value$unchanged),
      "invalid table columns are rejected per key while the tab still opens"
    )
    table_error <- tryCatch(
      tools$set_table_view(table = "gene_table", `_intent` = "unit test"),
      error = function(e) conditionMessage(e)
    )
    ok(
      grepl("Unknown table", table_error),
      "unknown tables are rejected with suggestions"
    )

    # ---- WP9 discovery tools ------------------------------------------
    cap_search <- tools$search_ui_capabilities(
      query = "filter", `_intent` = "unit test"
    )
    ok(
      cap_search@value$match_count >= 2L &&
        any(vapply(cap_search@value$capabilities,
                   function(r) identical(r$id, "dataspace.tab_pheno.column_filters"),
                   logical(1))),
      "capability search matches widget capabilities by meaning"
    )
    cap_get <- tools$get_ui_capability(
      id = "set_table_view", `_intent` = "unit test"
    )
    ok(
      ut_cmp_identical(cap_get@value$id, "tool:set_table_view") &&
        ut_cmp_identical(cap_get@value$writable, TRUE),
      "capability get returns semantic tool records"
    )
    cap_error <- tryCatch(
      tools$get_ui_capability(id = "set_table_vie", `_intent` = "unit test"),
      error = function(e) conditionMessage(e)
    )
    ok(
      grepl("Closest matches", cap_error),
      "capability get suggests closest ids on typos"
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
      ut_cmp_identical(widgets_result@value$widget_count, 7L) &&
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
      ut_cmp_identical(swept_listing@value$widget_count, 12L) &&
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

    # ---- WP6: figure templates ---------------------------------------
    template_error <- tryCatch(
      tools$create_figure(
        template = "volcano", x = "score", y = "category",
        `_intent` = "unit test"
      ),
      error = function(e) conditionMessage(e)
    )
    ok(
      grepl("must be numeric", template_error),
      "template path validates through the server-side expansion"
    )
    conflict_result <- tools$create_figure(
      template = "volcano", x = "score", y = "score",
      spec = list(layers = list(list(geom = "point", x = "score", y = "score"))),
      `_intent` = "unit test"
    )
    ok(
      grepl(
        "Template arguments ignored",
        paste(conflict_result@value$warnings, collapse = " ")
      ) &&
        ut_cmp_identical(conflict_result@value$spec$layers[[1]]$x, "score"),
      "WP6b: spec takes precedence over template shorthand, with a warning"
    )
    missing_error <- tryCatch(
      tools$create_figure(`_intent` = "unit test"),
      error = function(e) conditionMessage(e)
    )
    ok(
      grepl("requires either a template", missing_error),
      "create_figure without template or spec is rejected"
    )
    unknown_column_error <- tryCatch(
      tools$create_figure(
        template = "volcano", x = "score", y = "log.fdr",
        `_intent` = "unit test"
      ),
      error = function(e) conditionMessage(e)
    )
    ok(
      grepl("Closest matches", unknown_column_error) &&
        grepl("logFdr", unknown_column_error),
      "template unknown columns suggest the closest match"
    )

    volcano_result <- tools$create_figure(
      template = "volcano",
      x = "score",
      y = "logFdr",
      label_top_n = 1L,
      title = "RE vs ME",
      `_intent` = "unit test"
    )
    ok(
      ut_cmp_identical(volcano_result@value$template, "volcano") &&
        ut_cmp_identical(volcano_result@value$data_source, "feature_annotation") &&
        ut_cmp_identical(volcano_result@value$layers[[2]]$geom, "vline") &&
        ut_cmp_identical(volcano_result@value$layers[[3]]$geom, "label") &&
        ut_cmp_identical(volcano_result@value$spec$features[[1]], "Gene1") &&
        ut_cmp_identical(volcano_result@value$labels$title, "RE vs ME"),
      "volcano template renders through the tool with metadata and labels"
    )

    boxplot_result <- tools$create_figure(
      template = "boxplot", x = "group",
      `_intent` = "unit test"
    )
    ok(
      ut_cmp_identical(boxplot_result@value$data_source, "expression") &&
        ut_cmp_identical(boxplot_result@value$template, "boxplot") &&
        ut_cmp_identical(boxplot_result@value$spec$layers[[1]]$x, "sample__group") &&
        ut_cmp_identical(boxplot_result@value$spec$layers[[1]]$fill, "sample__group") &&
        ut_cmp_identical(boxplot_result@value$spec$features, c("Gene1", "Gene2", "Gene3")),
      "boxplot expression mode uses the selected features and namespaced grouping"
    )

    swept_template <- tools$create_figure(
      template = "scatter", x = "score", y = "logFdr",
      color = "null", label_top_n = "null", title = "null", space = "{}",
      `_intent` = "sentinel sweep"
    )
    ok(
      ut_cmp_identical(swept_template@value$data_source, "feature_annotation") &&
        is.null(swept_template@value$spec$layers[[1]]$color) &&
        length(swept_template@value$spec$layers) == 1L &&
        grepl("score", swept_template@value$labels$title, fixed = TRUE),
      "template tool treats sentinel optionals exactly like omitted values"
    )

    template_revision <- volcano_result@value$spec
    template_revision$labels$title <- "Volcano, revised"
    template_revision$theme <- "classic"
    revised_template_figure <- tools$update_figure(
      figure_id = volcano_result@value$figure_id,
      spec = template_revision,
      `_intent` = "unit test"
    )
    ok(
      ut_cmp_identical(revised_template_figure@value$labels$title, "Volcano, revised") &&
        ut_cmp_identical(revised_template_figure@value$theme, "classic") &&
        is.null(revised_template_figure@value$template),
      "template figures revise through the echoed spec (WP3 round-trip)"
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
