#' Utility - extended figure control shiny ui
#' @param id id
#' @param circle circle icon for dropdown manu
#' @importFrom shinyWidgets dropdown
#' 
attr4selector_ui <- function(id, circle = TRUE) {
  ns <- NS(id)
  tagList(
    dropdown(
      margin = "25px",
      circle = circle, status = "default", icon = icon("gear"), width = "950px",
      tooltip = tooltipOptions(title = "Click to modify figure!"),
      br(),
      triselector_ui(ns("selectColorUI")),
      triselector_ui(ns("selectShapeUI")),
      triselector_ui(ns("selectSizeUI")),
      triselector_ui(ns("selectTooltipUI")),
      triselector_ui(ns("selectSearchCol")),
      fluidRow(
        column(1),
        column(9,  offset = 0, style='padding-left:15px; padding-right:5px; padding-top:0px; padding-bottom:0px',
               div(style="display: inline-block;vertical-align:top;", h5("Find and highlight:")),
               div(style="display: inline-block;vertical-align:top; width:70%;", 
                   selectInput(ns("searchon"), label = NULL, choices = NULL, multiple = TRUE, width = "800px"))
                   ),
        column(2, offset = 0, style='padding-right:21px; padding-top:0px; padding-bottom:0px',
               actionButton(ns("goSearch"), label = "Find", width = "120px"), align="right")
      )
    )
  )
}

#' Utility - extended figure control shiny module
#' @param input input
#' @param output output
#' @param session session
#' @param reactive_meta reactive meta info, usually phenotype data or feature data of ExpressoinSet
#' @param reactive_expr expression matrix
#' @param reactive_triset reactive value, a matrix of nx3 use for triselector
#' @examples 
#' #' # library(shiny)
#' # library(shinyjs)
#' # library(Biobase)
#' # library(shinyWidgets)
#' # source("Git/R/module_triselector.R")
#' # dat <- readRDS("Dat/exampleEset.RDS")
#' # pd <- pData(dat)
#' # fd <- fData(dat)
#' # expr <- exprs(dat)
#' # ts <- trisetter(meta = pd, expr = expr, combine = "pheno")
#' # # ts <- stringr::str_split_fixed(colnames(pd), pattern = "\\|", n = 3)
#' 
#' # ui <- fluidPage(
#' #   attr4selector_ui("a4test")
#' # )
#' 
#' # server <- function(input, output, session) {
#   k <- callModule(attr4selector_module, id = "a4test", reactive_meta=reactive(pd), #' reactive_expr=reactive(expr), reactive_triset = reactive(ts))
#' #   # observe(
#' #   #   print(k())
#' #   # )
#' # }#' 
#' # shinyApp(ui, server)
#' 
attr4selector_module <- function(
  input, output, session, reactive_meta=reactive(NULL), reactive_expr=reactive(NULL), reactive_triset = reactive(NULL)
) {
  
  selectColor <- callModule(triselector_module, id = "selectColorUI", reactive_x = reactive_triset, label = "Color", suspendWhenHidden = FALSE)
  selectShape <- callModule(triselector_module, id = "selectShapeUI", reactive_x = reactive_triset, label = "Shape", suspendWhenHidden = FALSE)
  selectSize <- callModule(triselector_module, id = "selectSizeUI", reactive_x = reactive_triset, label = "Size", suspendWhenHidden = FALSE)
  selectTooltip <- callModule(triselector_module, id = "selectTooltipUI", reactive_x = reactive_triset, label = "Tooltips", suspendWhenHidden = FALSE)
  searchOnCol <- callModule(triselector_module, id = "selectSearchCol", reactive_x = reactive_triset, label = "Search", suspendWhenHidden = FALSE)
  
  vv <- reactive( varSelector(searchOnCol(), expr = reactive_expr(), meta = reactive_meta()) )
  observe(
    updateSelectInput(session, "searchon", choices = vv())
  )
  
  params <- reactiveValues(highlight = NULL, highlightName = NULL, color = NULL, shape = NULL, size = NULL, tooltips = NULL)
  observeEvent(input$goSearch, {
    params$highlight <- which(vv() %in% input$searchon)
    params$highlightName <- searchOnCol()$variable
  })
  observe(
    params$color <- varSelector(selectColor(), reactive_expr(), reactive_meta(), alternative = selectShape()$variable)
  )
  observe(
    params$shape <- varSelector(selectShape(), reactive_expr(), reactive_meta(), alternative = selectColor()$variable)
  )
  observe(
    params$size <- varSelector(selectSize(), reactive_expr(), reactive_meta())
  )
  observe(
    params$tooltips <- varSelector(selectTooltip(), reactive_expr(), reactive_meta())
  )
  
  params
}

