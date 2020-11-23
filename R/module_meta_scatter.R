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

meta_scatter_module <- function(
  input, output, session, reactive_meta=reactive(NULL), reactive_expr=reactive(NULL), combine = c("pheno", "feature"), source = "plotlyscattersource"
) {
  ns <- session$ns
  
  triset <- reactive( {
    ts <- trisetter(expr = reactive_expr(), meta = reactive_meta(), combine = combine[1])
    ts[ts[, 1] != "Surv", ]
  } )
  v1 <- callModule(triselector_module, id = "tris_main_scatter1", reactive_x = triset, label = "X-axis")
  v2 <- callModule(triselector_module, id = "tris_main_scatter2", reactive_x = triset, label = "Y-axis")
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
################# examples #####################

# # 
# library(shiny)
# library(Biobase)
# dat <- readRDS("Dat/exampleEset.RDS")
# source("shiny/module_triselector.R")
# source("shiny/module_barplot.R")
# source("shiny/module_scatter.R")
# source("shiny/module_boxplot.R")
# 
# 
# ui <- fluidPage(
#   meta_scatter_ui("test_meta_scatter")
# )
# 
# server <- function(input, output, session) {
#   callModule(meta_scatter_module, id = "test_meta_scatter",
#              reactive_meta = reactive(pData(dat)),
#              combine = c("pheno", "feature")[1],
#              # reactive_meta=reactive(fData(dat)),
#              # combine = c("pheno", "feature")[2],
#              reactive_expr = reactive(exprs(dat))
#   )
# }
# 
# shinyApp(ui, server)
