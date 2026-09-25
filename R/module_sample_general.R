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
#' @param store Optional child view of the canonical widget store
#'   (\code{\link{widget_store_child}}) for this module, e.g.
#'   \code{resultspace.sample_general}. When given, the module's
#'   user-editable surface (the link-variable cascade and the embedded
#'   attribute-4 panel) registers on the store; when NULL the legacy
#'   status-restore path is kept.
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
  reactive_j = reactive(NULL), reactive_status = reactive(NULL),
  store = NULL) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns
  
  triset <- reactive({
    trisetter(meta = reactive_phenoData(), expr = reactive_expr(), combine = "pheno")
  })

  xax <- reactiveVal()

  # ------------------------------------------------------------------
  # Canonical widget-store bindings (control plane, plan section 6, S4):
  # the link-variable cascade (beeswarm / contingency / survival switch).
  # The attr4 panel registers under <prefix>.attr4.* inside
  # attr4selector_module. Choices providers mirror triset() WITHOUT req()
  # (a shiny.validation condition inside a provider aborts the whole
  # store_apply transaction) and include the Surv category, which decides
  # the survival view.
  # ------------------------------------------------------------------
  if (!is.null(store)) {
    .sg_ts <- function() {
      if (is.null(reactive_phenoData())) return(NULL)
      ts <- tryCatch(
        trisetter(meta = reactive_phenoData(), expr = reactive_expr(),
                  combine = "pheno"),
        shiny.silent.error = function(e) NULL, error = function(e) NULL)
      if (is.null(ts) || !is.matrix(ts) || nrow(ts) == 0L) NULL else ts
    }
    .sg_ts1 <- function() {
      ts <- .sg_ts()
      if (is.null(ts)) character(0) else unique(ts[, 1])
    }
    .sg_ts2 <- function(a) {
      ts <- .sg_ts()
      if (is.null(ts) || is.null(a) || !nzchar(a)) character(0)
      else unique(ts[ts[, 1] %in% a, 2])
    }
    .sg_ts3 <- function(a, b) {
      ts <- .sg_ts()
      if (is.null(ts)) character(0)
      else {
        i <- rep(TRUE, nrow(ts))
        if (!is.null(a) && nzchar(a)) i <- i & ts[, 1] %in% a
        if (!is.null(b) && nzchar(b)) i <- i & ts[, 2] %in% b
        unique(ts[i, 3])
      }
    }
    kx1 <- paste0(store$prefix, ".xax_analysis")
    kx2 <- paste0(store$prefix, ".xax_subset")
    store_register(
      store,
      widget_binding("xax_analysis", "select",
        label = "Link variable category",
        help = paste("Annotation category of the variable the selected",
                     "samples are analyzed against"),
        choices_provider = function(v) .sg_ts1()),
      widget_binding("xax_subset", "select", label = "Link variable subcategory",
        help = "Subcategory within the link-variable category",
        depends_on = "xax_analysis",
        choices_provider = function(v) .sg_ts2(v[[kx1]])),
      widget_binding("xax_variable", "select_cascaded", label = "Link variable",
        help = paste("Variable driving the analysis view: numeric gives a",
                     "beeswarm, categorical a contingency table, Surv a",
                     "Kaplan-Meier curve"),
        depends_on = c("xax_analysis", "xax_subset"),
        choices_provider = function(v) .sg_ts3(v[[kx1]], v[[kx2]]))
    )
    .sg_store_observers <- list()
    .sg_keep <- function(obs) {
      .sg_store_observers[[length(.sg_store_observers) + 1L]] <<- obs
      invisible(obs)
    }
    # programmatic triple writes (status restore, batch link) share one
    # helper; user picks sync back through the v1() observer below
    .sg_apply_triple <- function(p1, p2, p3, origin = "system") {
      patch <- stats::setNames(list(p1, p2, p3),
                               c("xax_analysis", "xax_subset", "xax_variable"))
      tryCatch(store_apply(store, patch, origin = origin, strict = FALSE),
               error = function(e) NULL)
    }
    v1 <- triselector_module(
      "tris_sample_general", reactive_x = triset, label = "Link selection to",
      reactive_selector1 = store_watch(store, "xax_analysis"),
      reactive_selector2 = store_watch(store, "xax_subset"),
      reactive_selector3 = store_watch(store, "xax_variable"),
      reactive_axis_request = store_epoch(store))
  } else {
    .sg_apply_triple <- function(p1, p2, p3, origin = "system") {
      xax(NULL)
      xax(list(v1 = p1, v2 = p2, v3 = p3))
    }
    v1 <- triselector_module(
      "tris_sample_general", reactive_x = triset, label = "Link selection to",
      reactive_selector1 = reactive(xax()$v1),
      reactive_selector2 = reactive(xax()$v2),
      reactive_selector3 = reactive(xax()$v3))
  }

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
    reactive_triset = triset, reactive_status = attr4select_status,
    store = store
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
             reactive_checkpoint = reactive(pheno()$type == "surv"),
             store = store
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
    reactive_i_samples = reactive_i_samples,
    store = store
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
        .sg_apply_triple(parts[1], parts[2], parts[3])
      }
    } else if (selected$source == "features") {
      # For features: set to "Feature" -> "Auto" -> [feature_name]
      feature_name <- selected$data$feature_name
      .sg_apply_triple("Feature", "Auto", feature_name)
    }
  })

  # Store glue: UI -> store sync for the cascade. The binding fires on
  # settled triples only (WP2 store_bind_triselector): mid-cascade echoes
  # never enter the store, and store pushes are acknowledged through the
  # same path once the triselector confirms. The cascade has no meaningful
  # default, so unlike the scalar widgets it is not seeded and fills on the
  # first pick.
  if (!is.null(store)) {
    store_bind_triselector(store,
      keys = c(analysis = "xax_analysis", subset = "xax_subset",
               variable = "xax_variable"),
      sel = v1, keep = .sg_keep)
  }

  ## save and restore status
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return()
    if (identical(length(s$xax), 3L)) {
      tr <- lapply(s$xax, function(x)
        if (is.null(x) || !nzchar(x) || identical(x, "--select--")) NULL else x)
      if (!any(vapply(tr, is.null, logical(1))))
        .sg_apply_triple(tr[[1]], tr[[2]], tr[[3]], origin = "restore")
    }
    })

  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return()
    attr4select_status(NULL)
    attr4select_status(s$attr4)    
    })

  ## return status ##
  rv <- reactiveValues()
  observe( rv$xax <- v1() )
  observe( rv$attr4 <- attr4select$status )
  reactive(
    reactiveValuesToList(rv)
    )

  }) # end moduleServer
}
