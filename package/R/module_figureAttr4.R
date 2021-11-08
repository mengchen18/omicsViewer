#' @description Utility - extended figure control shiny ui
#' @param id id
#' @param circle circle icon for dropdown manu
#' @importFrom shinyWidgets dropdown textInputIcon updateTextInputIcon updateSwitchInput
#' 
attr4selector_ui <- function(id, circle = TRUE) {
  ns <- NS(id)
  dropdown(
    margin = "25px",
    circle = circle, status = "default", icon = icon("gear"), width = "700px",
    tooltip = tooltipOptions(title = "Click to modify figure!"),
    br(),
    triselector_ui(ns("selectColorUI")),
    triselector_ui(ns("selectShapeUI")),
    triselector_ui(ns("selectSizeUI")),
    triselector_ui(ns("selectTooltipUI")),
    triselector_ui(ns("selectSearchCol")),
    conditionalPanel(
      "1 == 2",
      checkboxInput(ns("showSearchBox"), label = "show", value = FALSE)
      ),
    conditionalPanel(
      "input.showSearchBox == true",
      ns = ns,    
      div(
        style='padding-left:100px; padding-right:0px; padding-top:0px; padding-bottom:0px',
        selectInput(ns("searchon"), label = NULL, choices = NULL, multiple = TRUE, width = "100%")
        )
      ),
    fluidRow(      
      column(2),
      column(
        4, offset = 0, style='padding:2px;', 
        textInputIcon(inputId = ns("xcut"), label = "Select points by x/y cutoffs", value = "log10(2)", placeholder = "e.g. -1 or -log10(2)", icon = list("x-cut"))),
      column(
        4, offset = 0, style='padding-left:5px; padding-right:2px; padding-top:27px; padding-bottom:2px;',         
        textInputIcon(inputId = ns("ycut"), label = NULL, value = "-log10(0.05)", placeholder = "e.g 2 or -log10(0.05)", icon = list("y-cut"))),
      column(
        2, offset = 0, style='padding-left:5px; padding-right:2px; padding-top:4px; padding-bottom:2px;', 
        selectInput(inputId = ns("scorner"), label = "Area", choices = "None", selectize = TRUE))
    )
  )
}

#' @description Utility - extended figure control shiny module
#' @param input input
#' @param output output
#' @param session session
#' @param reactive_meta reactive meta info, usually phenotype data or feature data of ExpressoinSet
#' @param reactive_expr expression matrix
#' @param reactive_triset reactive value, a matrix of nx3 use for triselector
#' @param pre_volcano logical; whether select areas using volcano cutoff
#' @param reactive_status the status of scatter plot, e.g. color variable, shape variable, etc. 
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
#' # ui <- fluidPage(
#' #   attr4selector_ui("a4test")
#' # ) 
#' # server <- function(input, output, session) {
#   k <- callModule(attr4selector_module, id = "a4test", reactive_meta=reactive(pd), #' reactive_expr=reactive(expr), reactive_triset = reactive(ts))
#' #   # observe(
#' #   #   print(k())
#' #   # )
#' # }#' 
#' # shinyApp(ui, server)
#' 
attr4selector_module <- function(
  input, output, session, reactive_meta=reactive(NULL), reactive_expr=reactive(NULL), 
  reactive_triset = reactive(NULL), pre_volcano = reactive(FALSE),
  reactive_status = reactive(NULL)
) {
  ns <- session$ns
  params <- reactiveValues(highlight = NULL, highlightName = NULL, color = NULL, shape = NULL, size = NULL, tooltips = NULL, cutoff = NULL)

  selectColor_s1 <- reactiveVal()
  selectColor_s2 <- reactiveVal()
  selectColor_s3 <- reactiveVal()
  
  selectShape_s1 <- reactiveVal()
  selectShape_s2 <- reactiveVal()
  selectShape_s3 <- reactiveVal()
  
  selectSize_s1 <- reactiveVal()
  selectSize_s2 <- reactiveVal()
  selectSize_s3 <- reactiveVal()
  
  selectTooltip_s1 <- reactiveVal()
  selectTooltip_s2 <- reactiveVal()
  selectTooltip_s3 <- reactiveVal()
  
  searchOnCol_s1 <- reactiveVal()
  searchOnCol_s2 <- reactiveVal()
  searchOnCol_s3 <- reactiveVal()

  selectColor <- callModule(triselector_module, id = "selectColorUI", reactive_x = reactive_triset, label = "Color",
    reactive_selector1 = selectColor_s1, reactive_selector2 = selectColor_s2, reactive_selector3 = selectColor_s3)#, suspendWhenHidden = FALSE)
  selectShape <- callModule(triselector_module, id = "selectShapeUI", reactive_x = reactive_triset, label = "Shape",
    reactive_selector1 = selectShape_s1, reactive_selector2 = selectShape_s2, reactive_selector3 = selectShape_s3)#, suspendWhenHidden = FALSE)
  selectSize <- callModule(triselector_module, id = "selectSizeUI", reactive_x = reactive_triset, label = "Size",
    reactive_selector1 = selectSize_s1, reactive_selector2 = selectSize_s2, reactive_selector3 = selectSize_s3)#, suspendWhenHidden = FALSE)
  selectTooltip <- callModule(triselector_module, id = "selectTooltipUI", reactive_x = reactive_triset, label = "Tooltips",
    reactive_selector1 = selectTooltip_s1, reactive_selector2 = selectTooltip_s2, reactive_selector3 = selectTooltip_s3)#, suspendWhenHidden = FALSE)
  searchOnCol <- callModule(triselector_module, id = "selectSearchCol", reactive_x = reactive_triset, label = "Search",
    reactive_selector1 = searchOnCol_s1, reactive_selector2 = searchOnCol_s2, reactive_selector3 = searchOnCol_s3)#, suspendWhenHidden = FALSE)

  vv <- reactive( varSelector(searchOnCol(), expr = reactive_expr(), meta = reactive_meta()) )

  pre_search <- reactiveVal()
  observe({    
    updateSelectInput(session, "searchon", choices = vv(), selected = pre_search())
  })
  observe(
    updateCheckboxInput(session, "showSearchBox", value = !is.null(vv()))
  )
  
  # val_xcut <- reactiveVal(NULL)
  # observe({    
  #   val_xcut( text2num(input$xcut) )
  # })
  # val_ycut <- reactiveVal(NULL)
  # observe({
  #   val_ycut( text2num(input$ycut) )
  # })

  val_xcut <- reactive({ text2num(input$xcut) })
  val_ycut <- reactive({ text2num(input$ycut) })

  observeEvent(list(val_xcut(), val_ycut()), {    
    if (is.numeric(val_xcut()) && is.null(val_ycut())) {
      ac <- c("None", "left", "right")
    } else if (is.null(val_xcut()) && is.numeric(val_ycut())) {
      ac <- c("None", "top", "bottom")
    } else if (is.numeric(val_xcut()) && is.numeric(val_ycut())) {
      ac <- c("None", "volcano", "left", "right", "top", "bottom", "topleft", "topright", "bottomleft", "bottomright")      
    } else 
      ac <- "None"
    ps <- ac[1]
    if (!is.null(input$scorner) && input$scorner %in% ac) 
      ps <- input$scorner
    updateSelectInput(session, inputId = "scorner", choices = ac, selected = ps)
  })  

  observeEvent(pre_volcano(), {    
    if (pre_volcano()) {      
      l <- list(x = val_xcut(), y = val_ycut(), corner = "volcano")
      attr(l, "seed") <- Sys.time()
      params$cutoff <- l
      updateSelectInput(session, inputId = "scorner", selected = "volcano")
    } else {      
      params$cutoff <- list(x = val_xcut(), y = val_ycut(), corner = "None")
      updateSelectInput(session, inputId = "scorner", selected = "None")
    }
  })
    
  searchValue <- reactiveVal()
  observe({
    foo <- function() searchValue(input$searchon)
    debounce(foo, 1000) 
    })  

  observe({    
    if (is.null(vv())) {
      updateSelectInput(session, "searchon", choices = NULL, selected = NULL)
      searchValue(NULL)
      pre_search(NULL)
    }
    })
  observe({    
    params$highlight <- which(vv() %in% searchValue())
    isolate( params$highlightName <- searchOnCol()$variable )
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

  acorner <- reactiveVal()    
  i_xcut <- reactiveVal()    
  i_ycut <- reactiveVal()    
  # observeEvent(input$actSelect, {
  # observeEvent(input$scorner, {    
  observe({
    req( !is.null(input$scorner) && nchar(input$scorner) != 0 )      
    acorner( input$scorner )
    i_xcut( input$xcut )
    i_ycut( input$ycut )
    l <- list(x = val_xcut(), y = val_ycut(), corner =  input$scorner)
    attr(l, "seed") <- Sys.time()
    params$cutoff <- l
  })

  observe(
    params$status <- list(
      selectColor = selectColor(),
      selectShape = selectShape(),
      selectSize = selectSize(),
      selectTooltip = selectTooltip(),
      searchOnCol = searchOnCol(),
      searchValue = searchValue(),
      xcut = i_xcut(),
      ycut = i_ycut(),
      acorner = acorner()
      # ,
      # nClickCorner = input$actSelect
    )
  )
  ############### restore status ##############
  observe({        
    if (is.null(s <- reactive_status()))
      return(NULL)

    selectColor_s1( s$selectColor[[1]] )
    selectColor_s2( s$selectColor[[2]] )
    selectColor_s3( s$selectColor[[3]] )
    })
  observe({        
    if (is.null(s <- reactive_status()))
      return(NULL)
    selectShape_s1( s$selectShape[[1]] )
    selectShape_s2( s$selectShape[[2]] )
    selectShape_s3( s$selectShape[[3]] )
  })

  observe({        
    if (is.null(s <- reactive_status()))
      return(NULL)
    selectSize_s1( s$selectSize[[1]] )
    selectSize_s2( s$selectSize[[2]] )
    selectSize_s3( s$selectSize[[3]] )
  })

  observe({        
    if (is.null(s <- reactive_status()))
      return(NULL)
    selectTooltip_s1( s$selectTooltip[[1]] )
    selectTooltip_s2( s$selectTooltip[[2]] )
    selectTooltip_s3( s$selectTooltip[[3]] )
  })

  observe({        
    if (is.null(s <- reactive_status()))
      return(NULL)
    searchOnCol_s1( s$searchOnCol[[1]] )
    searchOnCol_s2( s$searchOnCol[[2]] )
    searchOnCol_s3( s$searchOnCol[[3]] )
  })

  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return(NULL)
    updateTextInputIcon(session, "xcut", value = s$xcut)
    updateTextInputIcon(session, "ycut", value = s$ycut) 
    updateSelectInput(session, "scorner", selected = s$acorner)
  })
  
  observe({
    if (is.null(s <- reactive_status()))
      return(NULL)
    pre_search(s$searchValue)     
  })

  params
}

