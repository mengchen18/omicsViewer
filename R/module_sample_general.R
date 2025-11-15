#' @description Utility sample general ui
#' @param id id
#' 
sample_general_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Module description for AI browsers and screen readers
    div(class = "sr-only", id = ns("module-help"),
      tags$h4("About Sample Group Analysis"),
      tags$p("This module provides comprehensive statistical analysis to compare and visualize sample groups or patient cohorts. It offers multiple analysis types including beeswarm plots for continuous variables, contingency table analysis for categorical variables, and survival analysis for time-to-event data. This is the sample-centric complement to the feature analysis module."),
      tags$h4("When to use sample analysis"),
      tags$p("Use this module when you want to analyze relationships between sample metadata variables, compare patient groups, test associations between clinical variables, or perform survival analysis stratified by sample characteristics. This is essential for understanding sample clustering, identifying clinical correlations, and validating patient groupings."),
      tags$h4("How to interpret results"),
      tags$p("The visualization type depends on your selected variables: beeswarm plots show distributions of continuous variables across groups with statistical tests, contingency tables show relationships between categorical variables with chi-square or Fisher's exact test, and survival curves show time-to-event outcomes stratified by groups with log-rank test. The results table provides detailed statistical metrics appropriate for each analysis type.")
    ),
    fluidRow(
      shinydashboard::box(
        batch_comparison_ui(ns("batch_comp")),
        height = "490px",
        width = 12
      ),
      column(12, style = "margin-top: 0px;", triselector_ui(ns("tris_sample_general"), right_margin = "5")),
      shinydashboard::box(
        dataTableDownload_ui(ns("msatab")) ,
        height = "500px",
        width = 5
      ),
      shinydashboard::box(
        uiOutput(ns("sample_general_plot")),
        height = "500px",
        width = 7
      )
    )
  )
}

#' @description Utility sample general module
#' @param id module id
#' @param reactive_phenoData reactive phenotype data
#' @param reactive_expr reactive expression data
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
#   sample_general_module("sample_general", reactive_phenoData = reactive(pd), #' reactive_j = reactive(sample(rownames(pd), size = 20)) )
#' # }
#' #
#' # shinyApp(ui, server)
#'
sample_general_module <- function(id, reactive_phenoData, reactive_expr,
  reactive_j = reactive(NULL), reactive_status = reactive(NULL)) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns
  
  triset <- reactive({
    trisetter(meta = reactive_phenoData(), expr = reactive_expr(), combine = "pheno")
  })

  xax <- reactiveVal()
  v1 <- triselector_module(
    "tris_sample_general", reactive_x = triset, label = "Link selection to",
    reactive_selector1 = reactive(xax()$v1),
    reactive_selector2 = reactive(xax()$v2),
    reactive_selector3 = reactive(xax()$v3))

  attr4select_status <- reactiveVal()
  reactive_input <- reactive({
    req(reactive_phenoData())
    req(reactive_expr())
    ee <- t(reactive_expr())
    colnames(ee) <- paste0("Feature|Auto|", colnames(ee))
    cbind(reactive_phenoData(), ee)
  })
  attr4select <- attr4selector_module(
    "a4_gp", reactive_meta = reactive_input,
    reactive_triset = triset, reactive_status = attr4select_status
  )
  
  pheno <- reactive({
    req(v1()$variable)
    req(!v1()$variable %in% c("", "--select--"))
    req(reactive_input())
    cs <- do.call(paste, list(v1(), collapse = "|"))
    if (!cs %in% colnames(reactive_input()))
      return(NULL)
    val <- reactive_input()[, cs]
    
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
      r <-  plotly_scatter_ui(ns("sample_general_beeswarm")) 
    if (pheno()$type == "table")
      r <- factorIndependency_ui(ns("sample_general_contab"))
    if (pheno()$type == "surv")
      r <- survival_ui(ns("sample_general_surv")) 
    tagList(
      column(11, r),
      column(1, attr4selector_ui(ns("a4_gp"), circle = FALSE, right = TRUE))      
      )
  })
  
  ## beeswarm
  # showRegLine <- reactiveVal(FALSE)
  htestV1 <- reactiveVal()
  htestV2 <- reactiveVal()
  vs_scatter <- plotly_scatter_module(
    "sample_general_beeswarm",
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
  factorIndependency_module("sample_general_contab",
             x = select, y = reactive(pheno()$value),
             reactive_checkpoint = reactive(pheno()$type == "table")
  )

  ## survival
  survival_module('sample_general_surv',
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

  dataTableDownload_module(
    "msatab", reactive_table = metatab, prefix = "SampleTable_"
  )

  ## batch comparison module
  # Convert reactive_j (row names) to indices for batch_comparison_module
  reactive_i_samples <- reactive({
    req(reactive_j())
    req(reactive_phenoData())
    which(rownames(reactive_phenoData()) %in% reactive_j())
  })

  # Call batch comparison module
  batch_comp_selected <- batch_comparison_module(
    "batch_comp",
    reactive_expr = reactive_expr,
    reactive_phenoData = reactive_phenoData,
    reactive_featureData = reactive(NULL),  # No feature data available in this context
    reactive_i_samples = reactive_i_samples
  )

  # Update triselector when a row is selected from batch comparison
  observeEvent(batch_comp_selected(), {
    selected <- batch_comp_selected()
    req(!is.null(selected))

    if (selected$source == "phenotype") {
      # Parse variable_name to extract category|subcategory|variable
      var_name <- selected$data$variable_name
      parts <- strsplit(var_name, "\\|")[[1]]

      if (length(parts) >= 3) {
        xax(NULL)  # Clear first to trigger reactivity
        xax(list(v1 = parts[1], v2 = parts[2], v3 = parts[3]))
      }
    } else if (selected$source == "features") {
      # For features: set to "Feature" -> "Auto" -> [feature_name]
      feature_name <- selected$data$feature_name
      xax(NULL)  # Clear first to trigger reactivity
      xax(list(v1 = "Feature", v2 = "Auto", v3 = feature_name))
    }
  })

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

  }) # end moduleServer
}
