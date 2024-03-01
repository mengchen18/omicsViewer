#' @description Utility enrichment analysis shiny ui
#' @param id id
#' @importFrom DT dataTableOutput
dose_response_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # table
    # uiOutput(ns("error")),
    # DT::dataTableOutput(ns("stab")),
    # column(12, style = "margin-top: 0px;", triselector_ui(ns("tris_sample_general"), right_margin = "5"))
    # triselector_ui(ns("tris_ora"), right_margin = "5"),
    # dataTableDownload_ui(ns("stab")),
    # dataTableDownload_ui(ns("overlapTab"))
    # plotly barplot
    plotOutput(ns("plot"))
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
  input, output, session, reactive_featureData, reactive_i
) {
  
  ns <- session$ns
  
  reactive_dr <- reactive({
    attr(reactive_featureData(), "DoseResponse")
  })

  dr1 <- reactive({
    kk <<- reactive_dr()
    reactive_dr()[[reactive_i()]]
    })

  output$plot <- renderPlot({
    plotDC(dr1())
    })
}
