#' Utility triselector ui
#' @param id id
#' 
triselector_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(1, offset = 0, style='padding-left:15px; padding-right:5px; padding-top:0px; padding-bottom:0px', 
             uiOutput(ns("groupLabel"))
      ),
      column(3, offset = 0, style='padding-left:15px; padding-right:5px; padding-top:0px; padding-bottom:0px',
             div(style="display: inline-block;vertical-align:top;", h5("Analysis:")),
             div(style="display: inline-block;vertical-align:top; width:65%;", 
                 uiOutput(ns("analysis.output")))),
      
      column(4, offset = 0, style='padding-left:15px; padding-right:5px; padding-top:0px; padding-bottom:0px',
             div(style="display: inline-block;vertical-align:top;", h5("Subset:")),
             div(style="display: inline-block;vertical-align:top; width:78%;", 
                 uiOutput(ns("subset.output")))),
      
      column(4, offset = 0, style='padding-left:15px; padding-right:5px; padding-top:0px; padding-bottom:0px',
             div(style="display: inline-block;vertical-align:top;", h5("Variable:")),
             div(style="display: inline-block;vertical-align:top; width:74%;", 
                 uiOutput(ns("variable.output"))))
    )
  )
}

#' Utility triselector ui
#' @param input input
#' @param output output
#' @param session session
#' @param reactive_x an nx3 matrix
#' @param reactive_selector1 default value for selector 1
#' @param reactive_selector2 default value for selector 2
#' @param reactive_selector3 default value for selector 3
#' @param label of the triselector
#' @param suspendWhenHidden seuspend when hidden
#' @examples 
#' # ### examples
#' # 
#' # dat <- readRDS("Dat/exampleEset.RDS")
#' # library(shiny)
#' # phenoData <- pData(dat)
#' # triset <- stringr::str_split_fixed(colnames(phenoData), '\\|', n= 3)
#' # 
#' # ui <- fluidPage(
#' #   triselector_ui("tres")
#' # )
#' # server <- function(input, output, session) {
#   v <- callModule(triselector_module, id = "tres", reactive_x = reactive(triset), label = #' "x-axis")
#' #   observe(print(v()))
#' # }
#' # 
#' # shinyApp(ui, server)
#' 
#' # ###
#' # 
#' # 
#' # dat <- readRDS("Dat/exampleEset.RDS")
#' # library(shiny)
#' # phenoData <- pData(dat)
#' # triset <- stringr::str_split_fixed(colnames(phenoData), '\\|', n= 3)
#' # 
#' # ui <- fluidPage(
#' #   triselector_ui("tres"),
#' #   triselector_ui("tres2")
#' # )
#' # server <- function(input, output, session) {
#' #   v1 <- callModule(triselector_module, id = "tres", x = reactive(triset))
#' #   v2 <- callModule(triselector_module, id = "tres2", x = reactive(triset), 
#' #                    reactive_selector1 = reactive(v1()$analysis), 
#' #                    reactive_selector2 = reactive(v1()$subset))
#' #   observe({
#' #     print("/////////////////////////")
#' #     print(v1())
#' #     })
#' # }
#' # 
#' # shinyApp(ui, server)
#' 
triselector_module <- function(input, output, session, 
                               reactive_x, 
                               reactive_selector1 = reactive(NULL), 
                               reactive_selector2 = reactive(NULL), 
                               reactive_selector3 = reactive(NULL),
                               label = "Group Label:",
                               suspendWhenHidden = TRUE) {
  
  ns <- session$ns
                
  output$groupLabel <- renderUI({
    h5(HTML(sprintf("<b>%s</b>", label)))
  })
  
  output$analysis.output <- renderUI({
    cc <- unique(reactive_x()[, 1])
    preselected <- NULL
    if (!is.null(reactive_selector1()))
      if (reactive_selector1() %in% cc)
        preselected <- reactive_selector1()
    selectInput(inputId = ns("analysis"), label = NULL, choices = cc, selectize = TRUE, selected = preselected)
  })
  output$subset.output <- renderUI(
    selectInput(inputId = ns("subset"), label = NULL, choices = NULL, selectize = TRUE)
  )
  output$variable.output <- renderUI(
    selectInput(inputId = ns("variable"), label = NULL, choices = NULL, selectize = TRUE)
  )
  outputOptions(output, "analysis.output", suspendWhenHidden = suspendWhenHidden )
  outputOptions(output, "subset.output", suspendWhenHidden = suspendWhenHidden )
  outputOptions(output, "variable.output", suspendWhenHidden = suspendWhenHidden )
  
  observe({
    preselected <- NULL
    cc <- unique(reactive_x()[reactive_x()[, 1] == input$analysis, 2])
    if (!is.null(reactive_selector2()))
      if (reactive_selector2() %in% cc)
        preselected <- reactive_selector2()
    updateSelectInput(session, inputId = "subset", choices = cc, selected = preselected)
  })
  
  observe({
    # req(input$analysis)
    # req(input$subset)
    cc <- reactive_x()[, 3][reactive_x()[, 1] == input$analysis & reactive_x()[, 2] == input$subset]
    if (length(cc) > 1)
      cc <- c("Select a variable!", cc)
    
    preselected <- NULL
    if (!is.null(reactive_selector3())) {
      preselected <- try(match.arg(reactive_selector3(), cc), silent = TRUE)
      if (inherits(preselected, "try-error"))
        preselected <- NULL
    }
    updateSelectInput(session, inputId = "variable", choices = cc, selected = preselected)
  })
  
  eventReactive(eventExpr = input$variable,
    list(analysis = input$analysis, subset = input$subset, variable = input$variable)
  )
}



