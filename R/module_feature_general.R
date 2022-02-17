#' @description utility - feature general ui
#' @param id id
#' @importFrom DT dataTableOutput

feature_general_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(1, style = "margin-top: 0px;", attr4selector_ui(ns("a4_gf"), circle = FALSE)),
      column(11, triselector_ui(ns("tris_feature_general")))
    ),
    uiOutput(ns("feature_general_plot")),
    dataTableDownload_ui(ns("mtab"))
  )
}

#' @description utility - feature general ui
#' @param input input
#' @param output output
#' @param session session
#' @param reactive_expr reactive expression matrix
#' @param reactive_i reactive row index to be highlighted
#' @param reactive_highlight reactive col index to be highlighted
#' @param reactive_phenoData reactive phenotype data
#' @param reactive_featureData reactive feature data
#' @param reactive_status saved status to restore
#' @importFrom DT renderDataTable
#' @importFrom reshape2 melt
#' @examples
#' #' # library(shiny)
#' # library(shinyBS)
#' # library(Biobase)
#' # dat <- readRDS("Dat/exampleEset.RDS")
#' # source("Git/R/module_triselector.R")
#' # source("Git/R/module_scatter.R")
#' # source("Git/R/module_contTableStats.R")
#' # source("Git/R/module_survival.R")
#' # source("Git/R/module_boxplot.R")
#' # source("Git/R/module_figureAttr4.R")
#' # source("Git/R/auxi_figureAttr4.R")
#' # 
#' # 
#' # ui <- fluidPage(
#' #   feature_general_ui("tres")
#' # )
#' # 
#' # server <- function(input, output, session) {
#' #   v <- callModule(feature_general_module, id = "tres",
#' #                   reactive_expr = reactive(exprs(dat)),
#' #                   # reactive_i = reactive(c(5, 6, 7)),
#' #                   reactive_i = reactive(c( 7)),
#' #                   reactive_highlight = reactive(c(3, 5, 10)),
#' #                   reactive_phenoData = reactive(pData(dat)),
#' #                   reactive_featureData = reactive(fData(dat)))
#' # }
#' # 
#' # shinyApp(ui, server)


feature_general_module <- function(input, output, session, 
                                   reactive_expr, reactive_i = reactive(NULL), 
                                   reactive_highlight = reactive(NULL),
                                   reactive_phenoData, 
                                   reactive_featureData,
                                   reactive_status = reactive(NULL)) {
  ns <- session$ns
  
  # selector
  triset <- reactive({
    ts <- trisetter(expr = reactive_expr(), meta = reactive_phenoData(), combine = "none")
    ts[ts[, 1] != "Surv", ]
  })
  
  xax <- reactiveVal()
  v1 <- callModule(
    triselector_module, id = "tris_feature_general", reactive_x = triset, label = 'Value', 
    reactive_selector1 = reactive(xax()$v1), 
    reactive_selector2 = reactive(xax()$v2), 
    reactive_selector3 = reactive(xax()$v3))

  attr4select_status <- reactiveVal()
  attr4select <- callModule(
    attr4selector_module, id = "a4_gf", reactive_meta = reactive_phenoData, reactive_expr = reactive_expr, 
    reactive_triset = triset, reactive_status = attr4select_status
  )
  
  # what to do
  pheno <- reactive({
    req(v1())
    cs <- do.call(paste, list(v1(), collapse = "|"))    
    if (!cs %in% colnames(reactive_phenoData()))
      return(NULL)
    reactive_phenoData()[, cs]
  })
  pheno_cat <- reactive({ is.factor(pheno()) || is.character(pheno()) })
  pheno_num <- reactive({ is.numeric(pheno()) })
  single_i <- reactive({ length(reactive_i()) == 1 })
  
  showBoxplot <- reactive(length(pheno()) == 0 || (pheno_num() && !single_i()) || 
    length(reactive_i()) == 0 || length(reactive_i()) >= 10)
  showBeeswarm <- reactive(pheno_cat() && length(reactive_i()) > 0 && length(reactive_i()) < 10)
  showScatter <- reactive( single_i () && pheno_num() )
  
  output$feature_general_plot <- renderUI({
    if (showBoxplot())
      return( plotly_boxplot_ui(ns("feature_general_boxplotly")) )
    if (showScatter())
      return( plotly_scatter_ui(ns("feature_general_scatter")) )
    if (showBeeswarm())
      return( plotly_scatter_ui(ns("feature_general_beeswarm")) )
  })
  
  # boxplot:
  #  - no phenoData selected 
  #  - multi feature selected - numerical phenoData selected - boxplot with external
  callModule(plotly_boxplot_module, id = "feature_general_boxplotly",
             reactive_param_plotly_boxplot = reactive({
               req(reactive_expr())
               ylab <- rownames(reactive_expr())[reactive_i()]
               if (length(ylab) > 1)
                 ylab <- "Abundance of selected features"
               ylab.extvar <- do.call(paste, list(v1(), collapse = "|"))
               list(
                 x = reactive_expr(), i = reactive_i(), 
                 highlight = match(reactive_highlight(), colnames(reactive_expr())), 
                 extvar = pheno(),
                 ylab = ylab, ylab.extvar = ylab.extvar)
             }), 
             reactive_checkpoint = showBoxplot
  )
  
  scatter_vars <- reactive({
    l <- list(source = "feature_general_module")
    l$color <- attr4select$color
    l$shape <- attr4select$shape
    l$size <- attr4select$size 
    l$tooltips <- attr4select$tooltips
    l$highlight <- attr4select$highlight
    l$highlightName <- attr4select$highlightName
    
    if (showScatter()) {
      l$x <- reactive_expr()[reactive_i(), ] 
      l$y <- pheno()
      l$xlab <- rownames(reactive_expr())[reactive_i()]
      l$ylab <- do.call(paste, list(v1(), collapse = "|"))
      if (is.null(l$tooltips))
        l$tooltips <- colnames(reactive_expr())
    }
    
    if (showBeeswarm()) {
      df <- melt(reactive_expr()[reactive_i(), , drop = FALSE])
      df$color <- rep(l$color, each = length(reactive_i()))
      df$pheno <- rep(pheno(), each = length(reactive_i()))
      xlab <- ""
      ylab <- rownames(reactive_expr())[reactive_i()]
      if (length(ylab) > 1)
        ylab <- "Abundance of multiple selected features"
      df <- na.omit(df)
      
      l$y <- df$value
      l$x <- df$pheno
      l$ylab <- ylab
      l$color <- df$color
      if (is.null(l$tooltips))
        l$tooltips <- sprintf("<b>Feature: </b>%s<br><b>Sample: </b>%s", df$Var1, df$Var2)
    }
    l
  })
  
  showRegLine <- reactiveVal(FALSE)
  htestV1 <- reactiveVal()
  htestV2 <- reactiveVal()
  v_scatter <- callModule(plotly_scatter_module, id = "feature_general_scatter",
                          reactive_param_plotly_scatter = scatter_vars,
                          reactive_checkpoint = showScatter,
                          reactive_regLine = reactive( showRegLine()))    
  observe({
    showRegLine(v_scatter()$regline) 
    })
  
  ## beeswarm:
  # - single feature selected - categorical phenoData selected
  # - multi feature selected - categorical phenoData selected
  v_beeswarm <- callModule(plotly_scatter_module, id = "feature_general_beeswarm",
                           reactive_param_plotly_scatter = scatter_vars, 
                           reactive_checkpoint = showBeeswarm,
                           htest_var1 = htestV1, htest_var2 = htestV2)
  
  metatab <- reactive({
    req(reactive_i())
    tab <- reactive_featureData()
    tab <- tab[, grep("^General\\|", colnames(tab)), drop = FALSE]
    tab <- tab[reactive_i(), , drop = FALSE]
    ic <- vapply(tab, is.numeric, logical(1)) & vapply(tab, is.integer, logical(1))
    tab[ic] <- lapply(tab[ic], signif, digits = 2)
    colnames(tab) <- sub('General\\|All\\|', "", colnames(tab))
    tab
  })
  
  callModule(
    dataTableDownload_module, id = "mtab", reactive_table = metatab, prefix = "FeatureTable_"
  )

  ## save and restore status
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return()
    xax(NULL)
    xax(list(v1 = s$xax[[1]], v2 = s$xax[[2]], v3 = s$xax[[3]]))
    })

  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return()
    attr4select_status(NULL)
    attr4select_status(s$attr4)    
    })

  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return()
    htestV1( s$htestV1 )
    htestV2( s$htestV2 )
    })  

  observeEvent(reactive_status(), {    
    if (!is.null(s <- reactive_status()))
      showRegLine(s$showRegLine) 
    })

  ## return status ##
  rv <- reactiveValues()
  observe( rv$xax <- v1() )
  observe( rv$showRegLine <- showRegLine() )
  observe( rv$attr4 <- attr4select$status )
  observe({
    rv$htestV1 <- v_beeswarm()$htestV1
    rv$htestV2 <- v_beeswarm()$htestV2
    })

  reactive({
    reactiveValuesToList(rv)
    })

}
