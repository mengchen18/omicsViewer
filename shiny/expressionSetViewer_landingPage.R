landingPage_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$style("#landingpageid{margin: auto;}"),
    h1("Welcome to BayBioMS ExpressionSetViewer for proteomics!"),
    absolutePanel(id = "landingpageid", class = "myClass", fixed = TRUE
                  ,draggable = FALSE, left = 0, right = 0, top = "35%"
                  ,width = 500, height = "auto", 
                  wellPanel(
                    uiOutput(ns("error.ui")),
                    fluidRow(
                      column(10, textInput(ns("id"), "Your passcode", placeholder = "Enter your passcode here!")),
                      column(2, actionButton(ns("go"), "Go!"), style = "padding-top:25px")
                    )
                  )
    )
  )
}

landingPage_module <- function(input, output, session) {
  
  ns <- session$ns
  
  hashTab <- data.frame(
    rbind(
      c("demo", "/path/to/file"),
      c("6655457", "/some/other/path")
    ),
    stringsAsFactors = FALSE
  )
  
  m0 <- eventReactive(input$go, {
    id <- trimws(input$id)
    i <- which(hashTab[, 1] %in% id)
    if (length(i) != 1)
      return(NULL)
    i
  })
  
  output$error <- renderUI({
    HTML("There was an error with your passcode, please try again! <br><hr>")
  })
  
  output$error.ui <- renderUI({
    req(is.null(m0()))
    uiOutput(ns("error"))
  })
  
  reactive({
    req(m0())
    hashTab[m0(), 2]
  })
}

################ example use ######################
library(ExpressionSetViewer)

.f <- list.files("Git/R", full.names = TRUE)
.f <- .f[-grep("/shinyApp.R", .f)]
for (.i in .f) source(.i)

ui <- fluidPage(
  uiOutput("uis"),
  uiOutput("aout")
)

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
  
  callModule(app_module, id = "app", dir = v)
  output$aout <- renderUI({
    req(!showLanding())
    app_ui("app")
  })
}
options(error = traceback)
shinyApp(ui, server)
