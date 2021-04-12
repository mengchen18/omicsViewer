#' The three-step selector - the ui function
#' @param id id
#' @export

triselector_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(2, offset = 0, align = "right", 
             style='padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px', 
             uiOutput(ns("groupLabel"))
      ),
      # column(3, offset = 0, style='padding-left:15px; padding-right:5px; padding-top:0px; padding-bottom:0px',
      #        div(style="display: inline-block;vertical-align:top;", h5("Analysis:")),
      #        div(style="display: inline-block;vertical-align:top; width:65%;", 
      #            uiOutput(ns("analysis.output")))),
      # 
      # column(4, offset = 0, style='padding-left:15px; padding-right:5px; padding-top:0px; padding-bottom:0px',
      #        div(style="display: inline-block;vertical-align:top;", h5("Subset:")),
      #        div(style="display: inline-block;vertical-align:top; width:78%;", 
      #            uiOutput(ns("subset.output")))),
      # 
      # column(4, offset = 0, style='padding-left:15px; padding-right:5px; padding-top:0px; padding-bottom:0px',
      #        div(style="display: inline-block;vertical-align:top;", h5("Variable:")),
      #        div(style="display: inline-block;vertical-align:top; width:74%;", 
      #            uiOutput(ns("variable.output"))))
      column(3, offset = 0, style='padding:2px;', uiOutput(ns("analysis.output"))),
      column(4, offset = 0, style='padding:2px;', uiOutput(ns("subset.output"))),
      column(3, offset = 0, style='padding:2px;', uiOutput(ns("variable.output")))
    )
  )
}

#' The three-step selector - the module function
#' @description The selector is used to select columns of phenotype and feature data. It maybe useful when
#'   develop modules. 
#' @param input input
#' @param output output
#' @param session session
#' @param reactive_x an nx3 matrix
#' @param reactive_selector1 default value for selector 1
#' @param reactive_selector2 default value for selector 2
#' @param reactive_selector3 default value for selector 3
#' @param label of the triselector
#' @param suspendWhenHidden seuspend when hidden
#' @export
#' @examples 
#' # dat <- readRDS("inst/extdata/demo.RDS")
#' # library(shiny)
#' # library(Biobase)
#' # fData <- fData(dat)
#' # triset <- stringr::str_split_fixed(colnames(fData), '\\|', n= 3)
#' # source("R/module_triselector.R")
#' # 
#' # ui <- fluidPage(
#' #   triselector_ui("tres"),
#' #   triselector_ui("tres2")
#' # )
#' # server <- function(input, output, session) {
#' #   v1 <- callModule(triselector_module, id = "tres", reactive_x = reactive(triset), 
#' #                    reactive_selector1 = reactive("ttest"),
#' #                    reactive_selector2 = reactive("RE_vs_ME"),
#' #                    reactive_selector3 = reactive("mean.diff")
#' #                    )
#' #   v2 <- callModule(triselector_module, id = "tres2", reactive_x = reactive(triset),
#' #                    reactive_selector1 = reactive("ttest"),
#' #                    reactive_selector2 = reactive("RE_vs_ME"),
#' #                    reactive_selector3 = reactive("log.fdr"))
#' #   observe({
#' #     print("/////////////////////////")
#' #     print(v1())
#' #     })
#' # }
#' # 
#' # shinyApp(ui, server)

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
  
  # init empty selectize input
  output$analysis.output <- renderUI({    
    cc <- unique(reactive_x()[, 1])
    # cc <- c("Select analysis", cc) ## placeholder
    selectInput(inputId = ns("analysis"), label = NULL, choices = cc, selectize = TRUE, selected = reactive_selector1())    
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
  
  # updat selectize input when reactive_x is given
  observe({
    req(input$analysis)
    cc <- unique(reactive_x()[reactive_x()[, 1] == input$analysis, 2])
    # cc <- c("Select subset", cc) # placeholder
    updateSelectInput(session, inputId = "subset", choices = cc, selected = reactive_selector2())
  })
    
  # observeEvent(list(input$subset, input$analysis), {    
  observe({    
    req(input$analysis)
    req(input$subset)
    
    cc <- reactive_x()[, 3][reactive_x()[, 1] == input$analysis & reactive_x()[, 2] == input$subset]
    cc <- c("Select a variable!", cc)
    preselected <- try(match.arg(reactive_selector3(), cc), silent = TRUE)
      if (inherits(preselected, "try-error"))
        preselected <- NULL
    
    updateSelectInput(session, inputId = "variable", choices = cc, selected = preselected)
  })
  
  # eventReactive(eventExpr = input$variable,
  reactive({
    req(input$variable)
    list(analysis = input$analysis, subset = input$subset, variable = input$variable)
  })
}



