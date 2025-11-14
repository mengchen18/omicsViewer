#' @description Utility enrichment analysis shiny ui
#' @param id id
#' @importFrom DT dataTableOutput
dose_response_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      7, plotOutput(ns("plot"), height = DOSE_RESPONSE_PLOT_HEIGHT)
    ),
    column(
      5, dataTableDownload_ui(ns("param")), style = paste0("margin-top: ", DOSE_RESPONSE_TABLE_MARGIN_TOP, ";")
    ),
    column(
      12, dataTableDownload_ui(ns("feature"))
    )
  )
}

#' @description Utility dose response shiny module
#' @param id module id
#' @param reactive_expr reactive expression matrix
#' @param reactive_phenoData reactive phenotype data
#' @param reactive_featureData reactive feature data
#' @param reactive_i reactive index of rows to be selected
#' @param reactive_attr_drc reactive dose response curve attributes
#' @importFrom fastmatch fmatch

dose_response_module <- function(
    id,
    reactive_expr,
    reactive_phenoData,
    reactive_featureData,
    reactive_i,
    reactive_attr_drc
) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns
  
  dr1 <- reactive({
    req(length(reactive_i()) == 1)
    req(reactive_featureData())
    req(reactive_attr_drc())
    v <- reactive_attr_drc()
    list(
      col_dose = paste("General", "All", v["dose_col"], sep = "|"),
      col_curveid = paste("General", "All", v["curveid_col"], sep = "|")
      )
  })
  
  tabs <- reactive({
    req(dr1())
    plotDCMat(
      expr = reactive_expr(), 
      pd = reactive_phenoData(), 
      fd = reactive_featureData(), 
      featid = reactive_i(),
      dose.var = dr1()$col_dose, 
      curve.var = dr1()$col_curveid,
      only.par = TRUE
    )
  })
  
  output$plot <- renderPlot({
    req(dr1())
    plotDCMat(
      expr = reactive_expr(), 
      pd = reactive_phenoData(), 
      fd = reactive_featureData(), 
      featid = reactive_i(),
      dose.var = dr1()$col_dose, 
      curve.var = dr1()$col_curveid
    )
  })
  
  dataTableDownload_module(
    "feature", reactive_table = reactive(tabs()$featInfo), prefix = "Response_featureInfo_"
  )
  dataTableDownload_module(
    "param", reactive_table = reactive(tabs()$par), prefix = "Response_par_"
  )

  }) # end moduleServer
}
