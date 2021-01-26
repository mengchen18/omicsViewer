#' Application level 1 ui - data space
#' @param id id
#' @importFrom shinythemes shinytheme
#' @importFrom shinyWidgets tooltipOptions
#' 
L1_data_space_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    # tabsetPanel(
    navbarPage(
      "Eset",
      theme = shinytheme("spacelab"), 
      tabPanel('Features', meta_scatter_ui(ns('feature_space'))),
      tabPanel("Samples", meta_scatter_ui(ns("sample_space"))),
      tabPanel(
        "Heatmap",
        dropdownButton(
          inputId = "mydropdown",
          label = "Controls",
          circle = FALSE, status = "default", icon = icon("gear"),
          width = 700,
          tooltip = tooltipOptions(title = "Click to update heatmap and check legend!"),
          margin = "10px",
          tabsetPanel(
            tabPanel("Parameters", iheatmapInput(id = ns("heatmapViewer"))),
            tabPanel("Legend", iheatmapLegend(id = ns("heatmapViewer")))
          )
        ),
        iheatmapOutput(id = ns("heatmapViewer"))
      ),
      tabPanel("Sample-tab", dataTable_ui(ns("tab_pheno"))),
      tabPanel("Feature-tab", dataTable_ui(ns("tab_feature"))),
      tabPanel("Exprs", dataTable_ui(ns("tab_expr")))
    )
  )
}

#' Application level 1 module - data space
#' @param input input
#' @param output output
#' @param session session
#' @param expr reactive value; expression matrix
#' @param pdata reactive value; phenotype data
#' @param fdata reactive value; feature data
#' @param x1_s pre-selected x axis for sample space, treselector 1
#' @param x2_s pre-selected x axis for sample space, treselector 2
#' @param x3_s pre-selected x axis for sample space, treselector 3
#' @param y1_s pre-selected y axis for sample space, treselector 1
#' @param y2_s pre-selected y axis for sample space, treselector 2
#' @param y3_s pre-selected y axis for sample space, treselector 3
#' @param x1_f pre-selected x axis for feature space, treselector 1
#' @param x2_f pre-selected x axis for feature space, treselector 2
#' @param x3_f pre-selected x axis for feature space, treselector 3
#' @param y1_f pre-selected y axis for feature space, treselector 1
#' @param y2_f pre-selected y axis for feature space, treselector 2
#' @param y3_f pre-selected y axis for feature space, treselector 3
#' 
L1_data_space_module <- function(input, output, session, expr, pdata, fdata,
                                 x1_s = NULL, x2_s = NULL, x3_s = NULL,
                                 y1_s = NULL, y2_s = NULL, y3_s = NULL,
                                 x1_f = NULL, x2_f = NULL, x3_f = NULL,
                                 y1_f = NULL, y2_f = NULL, y3_f = NULL
                                 ) {
  
  s_heatmap <- callModule( iheatmapModule, 'heatmapViewer', mat = expr, pd = pdata, fd = fdata )
  
  # sample space
  s_sample_fig <- callModule(
    meta_scatter_module, id = "sample_space", reactive_meta = pdata, reactive_expr = expr, combine = "pheno", source = "scatter_meta_sample",
    reactive_x1 = reactive(x1_s), reactive_x2 = reactive(x2_s), reactive_x3 = reactive(x3_s),
    reactive_y1 = reactive(y1_s), reactive_y2 = reactive(y2_s), reactive_y3 = reactive(y3_s)
  )

  # feature space
  s_feature_fig <- callModule(
    meta_scatter_module, id = "feature_space", reactive_meta = fdata, reactive_expr = expr, combine = "feature", source = "scatter_meta_feature",
    reactive_x1 = reactive(x1_f), reactive_x2 = reactive(x2_f), reactive_x3 = reactive(x3_f),
    reactive_y1 = reactive(y1_f), reactive_y2 = reactive(y2_f), reactive_y3 = reactive(y3_f)
  )

  observe({
    print(x3_f)
    print(y3_f)
    })

  ## tables
  tab_pd <- callModule(
    dataTable_module, id = "tab_pheno",  reactive_data = pdata,
    reactiveSelectorMeta = s_sample_fig, reactiveSelectorHeatmap = s_heatmap, subset = "col"
    )
  tab_fd <- callModule(
    dataTable_module, id = "tab_feature",  reactive_data = fdata,
    reactiveSelectorMeta = s_feature_fig, reactiveSelectorHeatmap = s_heatmap, subset = "row"
  )
  tab_expr <- callModule(
    dataTable_module, id = "tab_expr",  reactive_data = reactive(
      cbind(data.frame(feature = rownames(expr()), expr()))
    ), reactiveSelectorMeta = s_feature_fig, reactiveSelectorHeatmap = s_heatmap, subset = "row", selector = FALSE)
  
  ### return selected feature and samples
  selectedFeatures <- reactiveVal()
  selectedSamples <- reactiveVal()

  observeEvent(s_heatmap(), {

    if (!is.null(s_heatmap()$brushed$row)) {
      selectedFeatures(s_heatmap()$brushed$row)
    } else if (!is.null(s_heatmap()$clicked))
      selectedFeatures(s_heatmap()$clicked["row"])

    if (!is.null(s_heatmap()$brushed$col)) {
      selectedSamples(s_heatmap()$brushed$col)
    } else if (!is.null(s_heatmap()$clicked))
      selectedSamples(s_heatmap()$clicked["col"])
  }
  )
  # 
  observeEvent(s_feature_fig(), {
    if (!is.null(s_feature_fig()$selected) && length(s_feature_fig()$selected) > 0) {
      selectedFeatures( s_feature_fig()$selected )
    } else if ( !is.null(s_feature_fig()$clicked) )
      selectedFeatures( s_feature_fig()$clicked )
  })

  observeEvent(s_sample_fig(), {
    if (!is.null(s_sample_fig()$selected) && length(s_sample_fig()$selected) > 0) {
      selectedSamples( s_sample_fig()$selected )
    } else if ( !is.null(s_sample_fig()$clicked) )
      selectedSamples( s_sample_fig()$clicked )
  })

  observe( selectedSamples(tab_pd()) )
  observe( selectedFeatures(tab_fd()) )
  observe( selectedFeatures(tab_expr()) )

  reactive({
    list(
      feature = selectedFeatures(),
      sample = selectedSamples()
    )
  })
}

