#' @description Utility enrichment analysis shiny ui
#' @param id id
#' @importFrom DT dataTableOutput
dose_response_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      7, plotOutput(ns("plot"), height = "550px")
    ),
    column(
      5, dataTableDownload_ui(ns("param")), style = "margin-top: 30px;"
    ),
    column(
      12, dataTableDownload_ui(ns("feature"))
    )
  )
}

#' @description Utility enrichment analysis shiny module
#' @param input input
#' @param output output
#' @param session session
#' @param reactive_i reactive index of rows to be selected (for ORA)
#' @param reactive_featureData reactive feature data
#' @importFrom fastmatch fmatch

dose_response_module <- function(
    input, output, session, 
    reactive_expr, 
    reactive_phenoData, 
    reactive_featureData,
    reactive_i,
    reactive_attr_drc
) {
  
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
  
  callModule(
    dataTableDownload_module, id = "feature", reactive_table = reactive(tabs()$featInfo), prefix = "Response_featureInfo_"
  )
  callModule(
    dataTableDownload_module, id = "param", reactive_table = reactive(tabs()$par), prefix = "Response_par_"
  )
}
