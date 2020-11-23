L1_data_space_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    # tabsetPanel(
    navbarPage(
      "Data analysis",
      theme = shinythemes::shinytheme("spacelab"), 
      tabPanel(
        "heatmap",
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
      tabPanel("Sample space", meta_scatter_ui(ns("sample_space"))),
      tabPanel('feature_space', meta_scatter_ui(ns('feature_space'))),
      tabPanel("phenotype_tab", dataTable_ui(ns("tab_pheno"))),
      tabPanel("feature_tab", dataTable_ui(ns("tab_feature"))),
      tabPanel("exprs", dataTable_ui(ns("tab_expr"), selector = FALSE))
    )
  )
}


L1_data_space_module <- function(input, output, session, expr, pdata, fdata) {
  
  # heatmap
  s_heatmap <- callModule( iheatmapModule, 'heatmapViewer', mat = expr, pd = pdata, fd = fdata )
  
  # sample space 
  s_sample_fig <- callModule(
    meta_scatter_module, id = "sample_space", reactive_meta = pdata, reactive_expr = expr, combine = "pheno", source = "scatter_meta_sample"
  )
  
  # feature space
  s_feature_fig <- callModule(
    meta_scatter_module, id = "feature_space", reactive_meta = fdata, reactive_expr = expr, combine = "feature", source = "scatter_meta_feature"
  )
  
  ## tables
  tab_pd <- callModule(dataTable_module, id = "tab_pheno",  reactive_data = pdata)
  tab_fd <- callModule(dataTable_module, id = "tab_feature",  reactive_data = fdata)
  tab_expr <- callModule(dataTable_module, id = "tab_expr",  reactive_data = reactive(
    cbind(data.frame(feature = rownames(expr()), expr()))
  ), selector = FALSE)
  
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
  
  observeEvent(s_feature_fig(), {
    if (!is.null(s_feature_fig()$selected) && length(s_feature_fig()$selected) > 0) {
      selectedFeatures( s_feature_fig()$selected )
    } else if ( !is.null(s_feature_fig()$clicked) ) 
      selectedFeatures( s_feature_fig()$clicked )
  })
  
  observeEvent(s_sample_fig(), {
    print( s_sample_fig() )
    if (!is.null(s_sample_fig()$selected) && length(s_sample_fig()$selected) > 0) {
      selectedSamples( s_sample_fig()$selected )
    } else if ( !is.null(s_sample_fig()$clicked) ) 
      selectedSamples( s_sample_fig()$clicked )
  })
  
  observe( selectedSamples(tab_pd()) )
  observe( selectedFeatures(tab_fd()) )
  observe( selectedFeatures(tab_expr()) )
  
  reactive({
    # print( selectedFeatures() )
    # print( selectedSamples() )
    list(
      feature = selectedFeatures(),
      sample = selectedSamples()
    )
  })

  
  # expression matrix
  
  
}

