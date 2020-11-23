
# expr <- exprs(dat)
# pdata <- pData(dat)
# fdata <- fData(dat)
a <- installed.packages()
library(tidyr, lib.loc = "/usr/local/lib/R/site-library")
library(survminer)
library(survival)
library(fgsea)
library(Biobase)
library(shiny)
library(shinyBS)
library(shinydashboard)
library(networkD3)



dat <- readRDS("Dat/exampleEset.RDS")
.f <- list.files("Git/R", full.names = TRUE)
.f <- .f[-grep("/shinyApp.R", .f)]
for (.i in .f) source(.i)

ui <- fluidPage(
  uiOutput("summary"),
  br(),
  fluidRow(
    column(6, L1_data_space_ui('dataspace')),
    column(6, L1_result_space_ui("resultspace"))
  )
)

server <- function(input, output, session) {
  
  reactive_eset <- reactive( dat )
  
  expr <- reactive({
    req(reactive_eset())
    Biobase::exprs(reactive_eset())
  })
  
  output$summary <- renderUI({
    txt <- sprintf(
      '<h1 style="display:inline;">ExpressionSetViewer</h1> <h3 style="display:inline;">---- ExpressionSet with %s features and %s samples:</h3>', 
      nrow(expr()), ncol(expr())
      )
    HTML(txt)
  })
  
  pdata <-reactive({
    req(reactive_eset())
    Biobase::pData(reactive_eset())
  })
  
  fdata <-reactive({
    req(reactive_eset())
    Biobase::fData(reactive_eset())
  })
  
  v1 <- callModule(L1_data_space_module, id = "dataspace", expr = expr, pdata = pdata, fdata = fdata)
  
  callModule(L1_result_space_module, id = "resultspace", 
             reactive_expr = expr, 
             reactive_phenoData = pdata, 
             reactive_featureData = fdata, 
             reactive_i = reactive(v1()$feature), 
             reactive_highlight = reactive(v1()$sample)
  )
}
shinyApp(ui, server)
