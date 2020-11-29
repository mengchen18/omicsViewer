library(ExpressionSetViewer)
source("landingPage.R")

server <- function(input, output, session) {
  
  ns <- session$ns
  
  v <- callModule(landingPage_module, id = "test")
  showLanding <- reactiveVal(TRUE)
  observe({
    showLanding ( length(v()) != 1 )
  })
  output$uis <- renderUI({
    req(showLanding())
    landingPage_ui("test") 
  })
  
  observe(print(v()))
  
  callModule(ExpressionSetViewer:::app_module, id = "app", dir = v)
  output$aout <- renderUI({
    req(!showLanding())
    ExpressionSetViewer:::app_ui("app")
  })
}


