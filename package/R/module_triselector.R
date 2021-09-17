#' The three-step selector - the ui function
#' @param id id
#' @export
#' @examples
#' if (interactive()) {
#'   library(shiny)
#'   library(Biobase)
#'   
#'   file <- system.file("extdata/demo.RDS", package = "ExpressionSetViewer")
#'   dat <- readRDS(file)
#'   fData <- fData(dat)
#'   triset <- stringr::str_split_fixed(colnames(fData), '\\|', n= 3)
#'   
#'   ui <- fluidPage(
#'     triselector_ui("tres"),
#'     triselector_ui("tres2")
#'   )
#'   server <- function(input, output, session) {
#'     v1 <- callModule(triselector_module, id = "tres", reactive_x = reactive(triset),
#'                      reactive_selector1 = reactive("ttest"),
#'                      reactive_selector2 = reactive("RE_vs_ME"),
#'                      reactive_selector3 = reactive("mean.diff")
#'     )
#'     v2 <- callModule(triselector_module, id = "tres2", reactive_x = reactive(triset),
#'                      reactive_selector1 = reactive("ttest"),
#'                      reactive_selector2 = reactive("RE_vs_ME"),
#'                      reactive_selector3 = reactive("log.fdr"))
#'     observe({
#'       print("/////////////////////////")
#'       print(v1())
#'     })
#'   }
#'   
#'   shinyApp(ui, server)
#' }
#' @return a tagList of UI components

triselector_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(2, offset = 0, align = "right", 
             style='padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px', 
             uiOutput(ns("groupLabel"))
      ),
      column(3, offset = 0, style='padding:2px;', selectInput(inputId = ns("analysis"), label = NULL, choices = NULL, selectize = TRUE)),
      column(4, offset = 0, style='padding:2px;', selectInput(inputId = ns("subset"), label = NULL, choices = NULL, selectize = TRUE)),
      column(3, offset = 0, style='padding:2px;', selectInput(inputId = ns("variable"), label = NULL, choices = NULL, selectize = TRUE))
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
#' @export
#' @examples 
#' if (interactive()) {
#'   library(shiny)
#'   library(Biobase)
#'   
#'   file <- system.file("extdata/demo.RDS", package = "ExpressionSetViewer")
#'   dat <- readRDS(file)
#'   fData <- fData(dat)
#'   triset <- stringr::str_split_fixed(colnames(fData), '\\|', n= 3)
#'   
#'   ui <- fluidPage(
#'     triselector_ui("tres"),
#'     triselector_ui("tres2")
#'   )
#'   server <- function(input, output, session) {
#'     v1 <- callModule(triselector_module, id = "tres", reactive_x = reactive(triset),
#'                      reactive_selector1 = reactive("ttest"),
#'                      reactive_selector2 = reactive("RE_vs_ME"),
#'                      reactive_selector3 = reactive("mean.diff")
#'     )
#'     v2 <- callModule(triselector_module, id = "tres2", reactive_x = reactive(triset),
#'                      reactive_selector1 = reactive("ttest"),
#'                      reactive_selector2 = reactive("RE_vs_ME"),
#'                      reactive_selector3 = reactive("log.fdr"))
#'     observe({
#'       print("/////////////////////////")
#'       print(v1())
#'     })
#'   }
#'   
#'   shinyApp(ui, server)
#' }
#' @return an reactive object containing the selected values

triselector_module <- function(input, output, session, 
                               reactive_x, 
                               reactive_selector1 = reactive(NULL), 
                               reactive_selector2 = reactive(NULL), 
                               reactive_selector3 = reactive(NULL),
                               label = "Group Label:") {
  
  ns <- session$ns
                
  output$groupLabel <- renderUI({
    h5(HTML(sprintf("<b>%s</b>", label)))
  })
  
  inte <- reactive({
    aa <- grepl("tris_feature_general", ns("x"))
    length(aa) > 0 && aa
  })
  
  observeEvent(list(reactive_selector1()), {
    if (length(names(input)) == 0)
      return(NULL)
    cc <- unique(reactive_x()[, 1])
    if (!is.null(reactive_selector1()))
      ss <- reactive_selector1() else
        ss <- cc[1]
    updateSelectInput(session, inputId = "analysis", choices = cc, selected = ss)
  })
  
  observeEvent(list(names(input), reactive_x()), {
    if (length(names(input)) == 0)
      return(NULL)
    cc <- unique(reactive_x()[, 1])
    if (input$analysis %in% cc)
      ss <- input$analysis else if (!is.null(reactive_selector1()))
        ss <- reactive_selector1() else
          ss <- cc[1]
    updateSelectInput(session, inputId = "analysis", choices = cc, selected = ss)
  })
  
  # bug fix
  observeEvent(input$analysis, {
    req (input$analysis == "")
    cc <- unique(reactive_x()[, 1])
    if (is.null(reactive_selector1()))
      ss <- cc[1] else
        ss <- reactive_selector1()    
      updateSelectInput(session, inputId = "analysis", choices = cc, selected = ss)
    })
  
  # updat selectize input when reactive_x is given
  observe({
    input$analysis
    req(input$analysis)
    cc <- unique(reactive_x()[reactive_x()[, 1] == input$analysis, 2])
    updateSelectInput(session, inputId = "subset", choices = cc, selected = reactive_selector2())
  })
  
  observe({
    input$analysis
    input$subset
    req(input$analysis)
    req(input$subset)
    
    cc <- reactive_x()[, 3][reactive_x()[, 1] == input$analysis & reactive_x()[, 2] == input$subset]
    cc <- c("Select a variable!", cc)
    preselected <- try(match.arg(reactive_selector3(), cc), silent = TRUE)
      if (inherits(preselected, "try-error"))
        preselected <- NULL
    updateSelectInput(session, inputId = "variable", choices = cc, selected = preselected)
  })
  
  reactive({
    req(input$variable)
    list(analysis = input$analysis, subset = input$subset, variable = input$variable)
  })
}



