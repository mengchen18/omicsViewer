#' @description Application level 1 ui - data space
#' @param id id
#' @param activeTab one of "Feature", "Feature table", "Sample", "Sample table", "Heatmap"
#' @importFrom shinythemes shinytheme
#' @importFrom shinyWidgets tooltipOptions dropdown
#' 
L1_data_space_ui <- function(id, activeTab = "Feature") {
  ns <- NS(id)
  navbarPage(
    "Eset", id = ns("eset"),
    selected = activeTab,
    theme = shinytheme("spacelab"), 
    tabPanel('Feature', meta_scatter_ui(ns('feature_space'))),
    tabPanel("Feature table", dataTable_ui(ns("tab_feature"))),
    tabPanel("Sample", meta_scatter_ui(ns("sample_space"))),
    tabPanel("Sample table", dataTable_ui(ns("tab_pheno"))),
    tabPanel(
      "Heatmap",
      fluidRow(
        column(
          6,             
          dropdown(
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
          )
        ),
        column(
          6, align = "right",
          iheatmapClear(id = ns("heatmapViewer"))
        ),
        column(
          12, iheatmapOutput(id = ns("heatmapViewer"))
        )
      )
    ),
    tabPanel("Expression", dataTable_ui(ns("tab_expr")))
  )
}

#' @description Application level 1 module - data space
#' @param input input
#' @param output output
#' @param session session
#' @param expr reactive value; expression matrix
#' @param pdata reactive value; phenotype data
#' @param fdata reactive value; feature data
#' @param rowDendrogram row dendrogram list
#' @param reactive_x_s pre-selected x axis for sample space
#' @param reactive_y_s pre-selected y axis for sample space
#' @param reactive_x_f pre-selected x axis for feature space
#' @param reactive_y_f pre-selected y axis for feature space
#' @param status intial status

L1_data_space_module <- function(
  input, output, session, expr, pdata, fdata, 
  reactive_x_s = reactive(NULL), reactive_y_s = reactive(NULL),
  reactive_x_f = reactive(NULL), reactive_y_f = reactive(NULL),
  rowDendrogram = reactive(NULL), status = reactive(NULL)
) {
  
  ns <- session$ns

  # heatmap
  s_heatmap <- callModule( 
    iheatmapModule, 'heatmapViewer', mat = expr, pd = pdata, fd = fdata, 
    rowDendrogram = rowDendrogram, status = reactive(status()$eset_heatmap)
    )
  
  ### feature space
  r_feature_fig <- reactiveVal()
  s_feature_fig <- callModule(
    meta_scatter_module, id = "feature_space", reactive_meta = fdata, reactive_expr = expr, combine = "feature", source = "scatter_meta_feature",
    reactive_x = reactive_x_f, reactive_y = reactive_y_f, reactive_status = reactive(status()$eset_fdata_fig)
  )
  observe({ r_feature_fig ( s_feature_fig() ) })
  observeEvent(status(), { r_feature_fig ( list(
    clicked = status()$eset_fdata_fig$selection_clicked,
    selected = status()$eset_fdata_fig$selection_selected
   )) })  

  # sample space
  r_sample_fig <- reactiveVal()
  s_sample_fig <- callModule(
    meta_scatter_module, id = "sample_space", reactive_meta = pdata, reactive_expr = expr, combine = "pheno", source = "scatter_meta_sample",
    reactive_x = reactive_x_s, reactive_y = reactive_y_s, reactive_status = reactive(status()$eset_pdata_fig)
  )  
  observe({ r_sample_fig ( s_sample_fig() ) })
  observeEvent(status(), { r_sample_fig ( list(
    clicked = status()$eset_pdata_fig$selection_clicked,
    selected = status()$eset_pdata_fig$selection_selected
   )) })  
  
  ## tables
  tab_pd <- callModule(
    dataTable_module, id = "tab_pheno",  reactive_data = pdata, tab_status = reactive(status()$eset_pdata_tab),
    reactiveSelectorMeta = s_sample_fig, reactiveSelectorHeatmap = s_heatmap, subset = "col"
  )
  tab_fd <- callModule(
    dataTable_module, id = "tab_feature",  reactive_data = fdata, tab_status = reactive(status()$eset_fdata_tab),
    reactiveSelectorMeta = s_feature_fig, reactiveSelectorHeatmap = s_heatmap, subset = "row"
  )
  tab_expr <- callModule(
    dataTable_module, id = "tab_expr",  reactive_data = reactive(
      cbind(data.frame(feature = rownames(expr()), expr()))
    ), tab_status = reactive(status()$eset_exprs_tab),
    reactiveSelectorMeta = s_feature_fig, 
    reactiveSelectorHeatmap = s_heatmap, subset = "row", selector = FALSE)
  
  ### return selected feature and samples
  selectedFeatures <- reactiveVal()
  selectedSamples <- reactiveVal()
  
  observeEvent(s_heatmap(), {
    
    if (!is.null(s_heatmap()$brushed$row)) {
      selectedFeatures(s_heatmap()$brushed$row)
    } else if (!is.null(s_heatmap()$clicked))
      selectedFeatures(s_heatmap()$clicked["row"]) else
        selectedFeatures(character(0))
    
    if (!is.null(s_heatmap()$brushed$col)) {
      selectedSamples(s_heatmap()$brushed$col)
    } else if (!is.null(s_heatmap()$clicked))
      selectedSamples(s_heatmap()$clicked["col"]) else
        selectedSamples(character(0))
  })
  
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

  #### status for snapshot #####
  observe({
    if (!is.null(tb <- status()$eset_active_tab))
      updateNavbarPage(session = session, inputId = "eset", selected = tb)
    })
  ####

  reactive({
    l <- list(
      feature = selectedFeatures(),
      sample = selectedSamples()
    )
    attr(l, "status") <- list(
      eset_active_tab = input$eset,
      eset_pdata_tab = attr(tab_pd(), "status"), # -> pdata_tab
      eset_fdata_tab = attr(tab_fd(), "status"), # -> fdata_tab
      eset_exprs_tab = attr(tab_expr(), "status"), # -> exprs_tab
      eset_fdata_fig = attr(s_feature_fig(), "status"), # -> fdata_fig
      eset_pdata_fig = attr(s_sample_fig(), "status"), # -> pdata_fig
      eset_heatmap = attr(s_heatmap(), "status")
      )
    l
  })
}

