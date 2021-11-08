#' @description Utility - scatter plot for meta shiny ui
#' @param id id
#' 
meta_scatter_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(1, attr4selector_ui(ns("a4selector"))), # style = "margin-top: 20px;", 
      column(11, 
             triselector_ui(ns("tris_main_scatter1")),
             triselector_ui(ns("tris_main_scatter2")))
    ),
    plotly_scatter_ui(ns("main_scatterOutput"), height = "666px"),
    actionButton(ns("clear"), "Clear selection and box")
  )
}

#' @description Utility - scatter plot for meta shiny module
#' @param input input
#' @param output output
#' @param session session
#' @param reactive_meta reactive meta data, phenotype data or feature data
#' @param reactive_expr reactive expression data
#' @param combine how to combine the expression and meta data, pheno or feature?
#' @param source source id for plotly object
#' @param reactive_x reactive value for pre-selected x-aixs
#' @param reactive_y reactive value for pre-selected y-aixs
#' @param reactive_status the status of scatter plot, e.g. x-, y-axis, color variable, shape variable, etc. 
#' #' # library(shiny)
#' # library(Biobase)
#' # dat <- readRDS("Dat/exampleEset.RDS")
#' # source("Git/R/module_triselector.R")
#' # source("Git/R/module_barplot.R")
#' # source("Git/R/module_scatter.R")
#' # source("Git/R/module_boxplot.R")
#' # 
#' # 
#' # ui <- fluidPage(
#' #   meta_scatter_ui("test_meta_scatter")
#' # )
#' # 
#' # server <- function(input, output, session) {
#' #   callModule(meta_scatter_module, id = "test_meta_scatter",
#' #              reactive_meta = reactive(pData(dat)),
#' #              combine = c("pheno", "feature")[1],
#' #              # reactive_meta=reactive(fData(dat)),
#' #              # combine = c("pheno", "feature")[2],
#' #              reactive_expr = reactive(exprs(dat))
#' #   )
#' # }
#' # 
#' # shinyApp(ui, server)
#' 
meta_scatter_module <- function(
  input, output, session, reactive_meta=reactive(NULL), reactive_expr=reactive(NULL), 
  combine = c("pheno", "feature"), source = "plotlyscattersource",
  reactive_x = reactive(NULL), reactive_y = reactive(NULL),
  reactive_status = reactive(NULL)
) {
  ns <- session$ns
  
  triset <- reactive( {
    ts <- trisetter(expr = reactive_expr(), meta = reactive_meta(), combine = combine[1])
    ts[ts[, 1] != "Surv", ]
  } )
  
  xax <- reactiveVal()
  observeEvent(reactive_x(), {    
    r <- list()
    if (!is.null(reactive_x())) {
      l <- strsplit(reactive_x(), "\\|")[[1]]
      r <- list(v1 = l[1], v2 = l[2], v3 = l[3])
    }
    xax(r)
  })
  
  yax <- reactiveVal()
  observe({
    r <- list()
    if (!is.null(reactive_y())) {
      l <- strsplit(reactive_y(), "\\|")[[1]]
      r <- list(v1 = l[1], v2 = l[2], v3 = l[3])
    } 
    yax(r)
  })  

  v1 <- callModule(triselector_module, id = "tris_main_scatter1", reactive_x = triset, label = "X-axis", 
                   reactive_selector1 = reactive(xax()$v1), 
                   reactive_selector2 = reactive(xax()$v2), 
                   reactive_selector3 = reactive(xax()$v3))
  v2 <- callModule(triselector_module, id = "tris_main_scatter2", reactive_x = triset, label = "Y-axis",
                   reactive_selector1 = reactive(yax()$v1), 
                   reactive_selector2 = reactive(yax()$v2), 
                   reactive_selector3 = reactive(yax()$v3))
  
  pre_vol <- reactiveVal(FALSE)
  # pre_vol <- reactive({
  observe({
    vv <- c("v1", "v2", "v3")
    if (all(vv %in% names(xax())) && all(vv %in% names(yax()))) {
      if (xax()$v1 == "ttest" && 
          yax()$v1 == "ttest" && 
          xax()$v3 == "mean.diff" && 
          yax()$v3 %in% c("log.fdr", "log.pvalue"))
        pre_vol(TRUE)
    }    
  })
  
  attr4select_status <- reactiveVal()
  attr4select <- callModule(
    attr4selector_module, id = "a4selector", reactive_meta = reactive_meta, reactive_expr = reactive_expr, 
    reactive_triset = triset, pre_volcano = pre_vol, reactive_status = attr4select_status
  )
  
  xycoord <- reactive({
    req(!v1()$variable %in% c("Select a variable!", ""))
    req(!v2()$variable %in% c("Select a variable!", ""))
    x <- varSelector(v1(), reactive_expr(), reactive_meta())
    y <- varSelector(v2(), reactive_expr(), reactive_meta())
    req(x)
    req(y)
    req(is.numeric(x) || is.numeric(y))
    req(length(x) == length(y))
    list( x = x, y = y )
  })

  rectval <- reactiveVal(NULL)
  observe({    
    req(xycoord())    
    rectval( line_rect(l = attr4select$cutoff, xycoord())$rect )
  })
  
  observeEvent(input$clear, {
    pre_vol(FALSE)
    rectval( NULL )
  })
  observeEvent(list(v1(), v2()), {
    if (is.null(attr4select$cutoff) || attr4select$cutoff$corner == "None") 
      rectval(NULL)
  })
  
  # ins <- reactiveVal(NA)
  scatter_vars <- reactive({
    req(l <- xycoord())
    l$source <- source
    l$xlab <- attr(l$x, "label")
    l$ylab <- attr(l$y, "label")
    l$color <- attr4select$color
    l$shape <- attr4select$shape
    l$size <- attr4select$size
    l$tooltips <- attr4select$tooltips
    l$highlight <- attr4select$highlight
    l$highlightName <- attr4select$highlightName
    l$rect <- rectval()
    l$inSelection <- NA
    l
  })
  
  # scatter plot:
  #  - single feature selected - numerical phenoData selected
  showRegLine <- reactiveVal(FALSE)
  htestV1 <- reactiveVal()
  htestV2 <- reactiveVal()
  v_scatter <- callModule(
    plotly_scatter_module, id = "main_scatterOutput", reactive_param_plotly_scatter = scatter_vars,
    reactive_regLine = showRegLine, htest_var1 = htestV1, htest_var2 = htestV2)
  observe({    
    if (!is.null(s <- reactive_status()))
      showRegLine(s$showRegLine) else
        showRegLine(v_scatter()$regline) 
    })
  
  selVal <- reactiveVal(
    list(
      clicked = character(0),
      selected = character(0)
    )
  )
  sbc <- reactiveVal(FALSE)
  
  observeEvent(list(input$clear, reactive_expr()), {
    selVal( list(
      clicked = character(0),
      selected = character(0)
    ) )
    sbc(FALSE)
  })
  
  clientSideSelection <- reactiveVal(character(0))
  observeEvent( v_scatter(), { 
    if (combine == "pheno") 
      l <- colnames(reactive_expr()) else
        l <- rownames(reactive_expr())
      u_c <- l[v_scatter()$clicked]
      u_s <- l[v_scatter()$selected]
      req( !identical( tmp <- c(u_c, u_s),  clientSideSelection() ) )
      clientSideSelection(tmp)
      selVal( list(
        clicked = u_c,
        selected = u_s
      ) )
      sbc(FALSE)
  })
  
  emptyValues <- function(...) {
    l <- list(...)
    j <- sapply(l, function(x) is.null(x) || length(x) == 0)
    all(j)
  }

  returnCornerSelection <- reactiveVal(TRUE)    
  observeEvent( rectval(), {          
    
    if (!returnCornerSelection())
      return(NULL)   

    rec <- rectval()
    if (is.null(rec)) {
      selVal(list(
        clicked = character(0),
        selected = character(0)
      )) 
      return(NULL)
    }
    req( cc <- xycoord() )
    if (combine == "pheno")
      l <- colnames(reactive_expr()) else
        l <- rownames(reactive_expr())
    
    i <- lapply(rec, function(r1) {
      which( cc$x > r1["x0"] & cc$x < r1["x1"] & cc$y > r1["y0"] & cc$y < r1["y1"] )
    })
    i <- sort(unique(unlist(i)))
    selVal( list(
      clicked = character(0),
      selected = l[ i ]
    ) )
    sbc(TRUE)
  } )

  ############## status save ############### 
  observe({
    sv <- selVal()
    attr(sv, "status") <- list(
      xax = v1(),
      yax = v2(), 
      showRegLine = showRegLine(),
      attr4 = attr4select$status,
      htestV1 = v_scatter()$htest_V1,
      htestV2 = v_scatter()$htest_V2,
      selection_clicked = sv$clicked,
      selection_selected = sv$selected,
      selectByCorner = sbc()
      )
    selVal(sv)
    })

  ############## status restore ###############
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return()
    xax(NULL)
    xax(list(v1 = s$xax[[1]], v2 = s$xax[[2]], v3 = s$xax[[3]]))
    })
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return()
    yax(NULL)
    yax(list(v1 = s$yax[[1]], v2 = s$yax[[2]], v3 = s$yax[[3]]))
    })
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return()
    attr4select_status(NULL)
    attr4select_status(s$attr4)    
    })

  # pure mimic - a bit strange
  # observeEvent(reactive_status(), { 
  #   req(xycoord())
  #   req(s <- reactive_status())
  #   l <- list(x = text2num(s$attr4$xcut), y =  text2num(s$attr4$ycut), corner = s$attr4$acorner)
  #   try(r <- line_rect(l = l, xycoord())$rect, silent = TRUE)
  #   if (inherits(r, "try-error")) return(NULL)
  #   j <- sapply(r, function(x) length(x) == 4 && is.numeric(x))
  #   if (any(!j)) return(NULL)
  #   printWithName(r$rect, "r$rect")
  #   rectval( r )
  # })

  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return()
    htestV1( s$htestV1 )
    htestV2( s$htestV2 )
    })  
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return()
    returnCornerSelection( s$selectByCorner )
    selVal( list(
      clicked = s$selection_clicked,
      selected = s$selection_selected
      ))
    })
  #############################################

  selVal
}
