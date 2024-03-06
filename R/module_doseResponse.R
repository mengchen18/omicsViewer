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
#' @examples 
#' #' # source("Git/R/auxi_fgsea.R")
#' # source("Git/R/auxi_vectORA.R")
#' # source("Git/R/module_barplotGsea.R")
# dat <- readRDS("inst/extdata/demo.RDS")
# obj <- tallGS(dat)
# fd <- Biobase::fData(obj)
# fdgs <- attr(fd, "GS")
# selected_ids <- rownames(fd)[fd$`PCA|All|PC1(10.1%)` > 0.02]
# ui <- fluidPage(
#   enrichment_analysis_ui("ea")
# )
# server <- function(input, output, session) {
#   callModule(
#     enrichment_analysis_module, id = "ea",
#     reactive_featureData = reactive(fd), reactive_i = reactive(selected_ids)
#   )
# }
# shinyApp(ui, server)


dose_response_module <- function(
    input, output, session, 
    reactive_expr, 
    reactive_phenoData, 
    reactive_featureData,
    reactive_i
) {
  
  ns <- session$ns
  
  dr1 <- reactive({
    req(length(reactive_i()) == 1)
    req(reactive_featureData())
    attr(reactive_featureData(), "ResponseCurve")
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
  
  output$feature <- DT::renderDT(
    form
  )
  
  callModule(
    dataTableDownload_module, id = "feature", reactive_table = reactive(tabs()$featInfo), prefix = "Response_featureInfo_"
  )
  callModule(
    dataTableDownload_module, id = "param", reactive_table = reactive(tabs()$par), prefix = "Response_par_"
  )
}
