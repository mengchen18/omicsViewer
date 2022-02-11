#' @description Utility sample general ui
#' @param id id
#' 
sample_general_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(1, style = "margin-top: 0px;", attr4selector_ui(ns("a4_gp"), circle = FALSE)),
      column(11, triselector_ui(ns("tris_sample_general")))
    ),
    uiOutput(ns("sample_general_plot")),
    # DT::dataTableOutput(ns('mtab'))
    dataTableDownload_ui(ns("msatab"))
  )
}

#' @description Utility sample general module
#' @param input input
#' @param output output
#' @param session session
#' @param reactive_phenoData reactive phenotype data
#' @param reactive_j index for which row in phenotype data should be highlighted/selected
#' @param reactive_status saved status to restore
#' @examples 
#' #' # library(shiny)
#' # #
#' # source("Git/R/module_triselector.R")
#' # source("Git/R/module_scatter.R")
#' # source("Git/R/module_contTableStats.R")
#' # source("Git/R/module_survival.R")
#' # source("Git/R/module_figureAttr4.R")
#' # source("Git/R/auxi_figureAttr4.R")
#' # 
#' # dat <- readRDS("Dat/exampleEset.RDS")
#' # pd <- pData(dat)
#' # 
#' # ui <- fluidPage(
#' #   sample_general_ui("sample_general")
#' # )
#' # 
#' # server <- function(input, output, session) {
#   callModule(sample_general_module, id = "sample_general", reactive_phenoData = reactive(pd), #' reactive_j = reactive(sample(rownames(pd), size = 20)) )
#' # }
#' # 
#' # shinyApp(ui, server)
#'
sample_general_module <- function(input, output, session, reactive_phenoData, 
  reactive_j = reactive(NULL), reactive_status = reactive(NULL)) {
  
  ns <- session$ns
  
  triset <- reactive({
    trisetter(meta = reactive_phenoData(), combine = "none")
  })

  xax <- reactiveVal()
  v1 <- callModule(
    triselector_module, id = "tris_sample_general", reactive_x = triset, label = "Value",
    reactive_selector1 = reactive(xax()$v1), 
    reactive_selector2 = reactive(xax()$v2), 
    reactive_selector3 = reactive(xax()$v3))
  
  attr4select_status <- reactiveVal()
  attr4select <- callModule(
    attr4selector_module, id = "a4_gp", reactive_meta = reactive_phenoData, 
    reactive_triset = triset, reactive_status = attr4select_status
  )
  
  pheno <- reactive({
    req(v1()$variable)
    req(!v1()$variable %in% c("", "Select a variable!"))
    req(reactive_phenoData())
    cs <- do.call(paste, list(v1(), collapse = "|"))
    if (!cs %in% colnames(reactive_phenoData()))
      return(NULL)
    val <- reactive_phenoData()[, cs]
    
    if (v1()$analysis == "Surv") {
      type <- "surv"
    } else if (is.numeric(val)) {
      type <- "beeswarm"  
    } else if (is.character(val) || is.factor(val)) {
      type <- "table"
    } else {
      warnings("Unknown type of val: sample_general_module, return NULL!")
      return(NULL)
    }
    list(value = val, type = type)
  })
  
  select <- reactive({
    req(reactive_j())
    req(reactive_phenoData())
    select <- rep("Unselected", nrow(reactive_phenoData()))
    if (!is.null(reactive_j()))
      select[rownames(reactive_phenoData()) %in% reactive_j()] <- "selected"
    select
  })
  
  output$sample_general_plot <- renderUI({
    req(pheno()$type)
    if (pheno()$type == "beeswarm")
      return( plotly_scatter_ui(ns("sample_general_beeswarm")) )
    if (pheno()$type == "table")
      return( factorIndependency_ui(ns("sample_general_contab")) )
    if (pheno()$type == "surv")
      return( survival_ui(ns("sample_general_surv")) )
  })
  
  ## beeswarm
  # showRegLine <- reactiveVal(FALSE)  
  htestV1 <- reactiveVal()
  htestV2 <- reactiveVal()
  vs_scatter <- callModule(
    plotly_scatter_module, id = "sample_general_beeswarm", 
    reactive_param_plotly_scatter = reactive({
      req(reactive_j())
      req(pheno()$value)
      tooltips <- attr4select$tooltips
      if (is.null(tooltips))
        tooltips <- rownames(reactive_phenoData())      
      l <- list(
        x = select(), 
        y = pheno()$value,
        xlab = "", 
        ylab = do.call(paste, list(v1(), collapse = "|")),
        tooltips = tooltips
      )
      l$color <- attr4select$color
      l$shape <- attr4select$shape
      l$size <- attr4select$size
      l$highlight <- attr4select$highlight
      l$highlightName <- attr4select$highlightName
      l
    }), 
    reactive_regLine = reactive(FALSE), # showRegLine,
    reactive_checkpoint = reactive(pheno()$type == "beeswarm"),
    htest_var1 = htestV1, htest_var2 = htestV2)
  
  # cont table stats
  callModule(factorIndependency_module, id = "sample_general_contab", 
             x = select, y = reactive(pheno()$value),
             reactive_checkpoint = reactive(pheno()$type == "table")
  )
  
  ## survival
  callModule(survival_module, id = 'sample_general_surv', 
             reactive_resp = reactive(pheno()$value), reactive_strata = select,
             reactive_checkpoint = reactive(pheno()$type == "surv")
  )
  
  ## table
  metatab <- reactive({
    req(reactive_j())
    tab <- reactive_phenoData()
    tab <- tab[, grep("^General\\|", colnames(tab)), drop = FALSE]
    tab <- tab[reactive_j(), , drop = FALSE]
    ic <- vapply(tab, is.numeric, logical(1)) & vapply(tab, is.integer, logical(1))
    tab[ic] <- lapply(tab[ic], signif, digits = 2)
    colnames(tab) <- sub('General\\|All\\|', "", colnames(tab))
    tab
  })
  
  callModule(
    dataTableDownload_module, id = "msatab", reactive_table = metatab, prefix = "SampleTable_"
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

  ## return status ##
  rv <- reactiveValues()
  observe( rv$xax <- v1() )
  observe( rv$attr4 <- attr4select$status )
  observe({
    rv$htestV1 <- vs_scatter()$htestV1
    rv$htestV2 <- vs_scatter()$htestV2
    })

  reactive(
    reactiveValuesToList(rv)
    )
}
