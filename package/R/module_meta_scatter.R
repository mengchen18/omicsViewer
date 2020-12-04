#' Utility - scatter plot for meta shiny ui
#' @param id id
#' 
meta_scatter_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(1, style = "margin-top: 20px;", attr4selector_ui(ns("a4selector"))),
      column(11, 
             triselector_ui(ns("tris_main_scatter1")),
             triselector_ui(ns("tris_main_scatter2")))
    ),
    plotly_scatter_ui(ns("main_scatterOutput"), height = "666px")
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
#' @param reactive_x1 reactive value for pre-selected x-aixs, triselect 1
#' @param reactive_x2 reactive value for pre-selected x-aixs, triselect 2
#' @param reactive_x3 reactive value for pre-selected x-aixs, triselect 3
#' @param reactive_y1 reactive value for pre-selected y-aixs, triselect 1
#' @param reactive_y2 reactive value for pre-selected y-aixs, triselect 2
#' @param reactive_y3 reactive value for pre-selected y-aixs, triselect 3
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
  reactive_x1 = reactive(NULL), reactive_x2 = reactive(NULL), reactive_x3 = reactive(NULL),
  reactive_y1 = reactive(NULL), reactive_y2 = reactive(NULL), reactive_y3 = reactive(NULL)
) {
  ns <- session$ns
  
  triset <- reactive( {
    ts <- trisetter(expr = reactive_expr(), meta = reactive_meta(), combine = combine[1])
    ts[ts[, 1] != "Surv", ]
  } )
  
  # 
  v1 <- callModule(triselector_module, id = "tris_main_scatter1", reactive_x = triset, label = "X-axis", 
                   reactive_selector1 = reactive_x1, 
                   reactive_selector2 = reactive_x2, 
                   reactive_selector3 = reactive_x3)
  v2 <- callModule(triselector_module, id = "tris_main_scatter2", reactive_x = triset, label = "Y-axis",
                   reactive_selector1 = reactive_y1, 
                   reactive_selector2 = reactive_y2, 
                   reactive_selector3 = reactive_y3)
  # v1 <- callModule(triselector_module, id = "tris_main_scatter1", reactive_x = triset, label = "X-axis")
  # v2 <- callModule(triselector_module, id = "tris_main_scatter2", reactive_x = triset, label = "Y-axis")
  attr4select <- callModule(
    attr4selector_module, id = "a4selector", reactive_meta = reactive_meta, reactive_expr = reactive_expr, reactive_triset = triset
  )
  
  scatter_vars <- reactive({
    req(!v1()$variable %in% c("Select a variable!", ""))
    req(!v2()$variable %in% c("Select a variable!", ""))
    
    l <- list(source = source)
    # x-axis
    l$x <- varSelector(v1(), reactive_expr(), reactive_meta())
    # y-axis
    l$y <- varSelector(v2(), reactive_expr(), reactive_meta())
    # xlab
    l$xlab <- attr(l$x, "label")
    # ylab
    l$ylab <- attr(l$y, "label")
    
    l$color <- attr4select$color
    l$shape <- attr4select$shape
    l$size <- attr4select$size
    l$tooltips <- attr4select$tooltips
    l$highlight <- attr4select$highlight
    l$highlightName <- attr4select$highlightName
    l
  })
  
  # scatter plot:
  #  - single feature selected - numerical phenoData selected
  showRegLine <- reactiveVal(FALSE)
  v_scatter <- callModule(plotly_scatter_module, id = "main_scatterOutput",
                          reactive_param_plotly_scatter = scatter_vars,
                          reactive_regLine = reactive( showRegLine()))
  observe(showRegLine(v_scatter()$regline))
  
  reactive({
    if (combine == "pheno") 
      l <- colnames(reactive_expr()) else
        l <- rownames(reactive_expr())
      
      # print(v_scatter())
      i1 <- v_scatter()$selected
      i2 <- v_scatter()$clicked
      list(
        clicked = l[i2],
        selected = l[i1]
      )
  })
}

