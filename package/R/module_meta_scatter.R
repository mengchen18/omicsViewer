#' Utility - scatter plot for meta shiny ui
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

#' Utility - scatter plot for meta shiny module
#' @param input input
#' @param output output
#' @param session session
#' @param reactive_meta reactive meta data, phenotype data or feature data
#' @param reactive_expr reactive expression data
#' @param combine how to combine the expression and meta data, pheno or feature?
#' @param source source id for plotly object
#' @param reactive_x reactive value for pre-selected x-aixs
#' @param reactive_y reactive value for pre-selected y-aixs
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
  input, output, session, reactive_meta=reactive(NULL), reactive_expr=reactive(NULL), combine = c("pheno", "feature"), source = "plotlyscattersource",
  reactive_x = reactive(NULL), reactive_y = reactive(NULL)
) {
  ns <- session$ns
  
  triset <- reactive( {
    ts <- trisetter(expr = reactive_expr(), meta = reactive_meta(), combine = combine[1])
    ts[ts[, 1] != "Surv", ]
  } )
  
  xax <- reactive({    
    r <- list()
    if (!is.null(reactive_x())) {
      l <- strsplit(reactive_x(), "\\|")[[1]]
      r <- list(v1 = l[1], v2 = l[2], v3 = l[3])
    }
    r
  })
  
  yax <- reactive({
    r <- list()
    if (!is.null(reactive_y())) {
      l <- strsplit(reactive_y(), "\\|")[[1]]
      r <- list(v1 = l[1], v2 = l[2], v3 = l[3])
    } 
    r
  })
  
  v1 <- callModule(triselector_module, id = "tris_main_scatter1", reactive_x = triset, label = "X-axis", 
                   reactive_selector1 = reactive(xax()$v1), 
                   reactive_selector2 = reactive(xax()$v2), 
                   reactive_selector3 = reactive(xax()$v3))
  v2 <- callModule(triselector_module, id = "tris_main_scatter2", reactive_x = triset, label = "Y-axis",
                   reactive_selector1 = reactive(yax()$v1), 
                   reactive_selector2 = reactive(yax()$v2), 
                   reactive_selector3 = reactive(yax()$v3))
  
  pre_vol <- reactive({
    v <- FALSE
    vv <- c("v1", "v2", "v3")
    if (all(vv %in% names(xax())) && all(vv %in% names(yax()))) {
      if (xax()$v1 == "ttest" && 
          yax()$v1 == "ttest" && 
          xax()$v3 == "mean.diff" && 
          yax()$v3 %in% c("log.fdr", "log.pvalue"))
        v <- TRUE
    }
    v
  })
  
  attr4select <- callModule(
    attr4selector_module, id = "a4selector", reactive_meta = reactive_meta, reactive_expr = reactive_expr, reactive_triset = triset, pre_volcano = pre_vol
  )
  
  xycoord <- reactive({
    req(!v1()$variable %in% c("Select a variable!", ""))
    req(!v2()$variable %in% c("Select a variable!", ""))
    x <- varSelector(v1(), reactive_expr(), reactive_meta())
    y <- varSelector(v2(), reactive_expr(), reactive_meta())
    req(x)
    req(y)
    req(length(x) == length(y))
    list( x = x, y = y )
  })
  
  rectval <- reactiveVal(NULL)
  observe({
    req(xycoord())
    rectval( line_rect(l = attr4select$cutoff, xycoord())$rect )
  })
  observeEvent(input$clear, {
    rectval( NULL )
  })
  observeEvent(list(v1(), v2()), {
    if (is.null(attr4select$cutoff)) {
      rectval(NULL)
      return(NULL)
    }
    if (attr4select$cutoff$corner == "volcano")
      if (v1()[[1]] != "ttest" || v2()[[1]] != "ttest")
        rectval(NULL)
  })
  
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
    l
  })
  
  # scatter plot:
  #  - single feature selected - numerical phenoData selected
  showRegLine <- reactiveVal(FALSE)
  v_scatter <- callModule(
    plotly_scatter_module, id = "main_scatterOutput",
    reactive_param_plotly_scatter = scatter_vars,
    reactive_regLine = reactive( showRegLine())
  )
  observe( showRegLine(v_scatter()$regline) )
  
  selVal <- reactiveVal(
    list(
      clicked = character(0),
      selected = character(0)
    )
  )
  
  observeEvent(list(input$clear, reactive_expr()), {
    selVal( list(
      clicked = character(0),
      selected = character(0)
    ) )
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
  })
  
  observeEvent( rectval(), {
    req( rec <- rectval() )
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
  })
  
  selVal
}
