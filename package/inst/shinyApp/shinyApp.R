library(fgsea)
library(tidyr, lib.loc = "/usr/local/lib/R/site-library")
library(survminer)
library(survival)
library(Biobase)
library(shiny)
library(shinyBS)
library(shinydashboard)
library(networkD3)
library(shinyWidgets)


# dat <- readRDS("Dat/exampleEset.RDS")
.f <- list.files("Git/R", full.names = TRUE)
.f <- .f[-grep("/shinyApp.R", .f)]
for (.i in .f) source(.i)

ui <- fluidPage(
  app_ui("app")
)

server <- function(input, output, session) {
  callModule(app_module, id = "app", dir = reactive("/media/share_baybioms/Projects/008_Bioinformatics/B032_ExpressionSetViewer/Dat/"))
}
shinyApp(ui, server)


