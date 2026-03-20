#' @description utility - feature general ui
#' @param id id
#' @importFrom DT dataTableOutput

feature_general_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Module description for AI browsers and screen readers
    div(class = "sr-only", id = ns("module-help"),
      tags$h4("About Feature Statistical Analysis"),
      tags$p("This module performs comprehensive statistical analysis to compare selected features (genes/proteins) across different sample groups. It provides multiple visualization options including beeswarm boxplots for distribution comparison and ROC/PR curves for classification performance evaluation."),
      tags$h4("When to use this analysis"),
      tags$p("Use this module when you want to test if a feature shows significant differences between groups (e.g., disease vs. control, treatment vs. untreated). The beeswarm plot is ideal for visualizing distributions and identifying outliers, while ROC/PR curves help assess how well a feature can distinguish between two classes."),
      tags$h4("How to interpret results"),
      tags$p("For boxplots: The box shows the interquartile range (25th to 75th percentile), the line inside is the median, and individual points show data distributions. P-values indicate statistical significance of group differences. For ROC curves: Area Under Curve (AUC) near 1.0 indicates excellent classification, 0.5 is random chance. The curve shows trade-off between sensitivity and specificity at different thresholds.")
    ),
    tags$h3("Variable Selection and Plot Type", class = "sr-only", `aria-label` = "Controls for selecting statistical test results and choosing between beeswarm boxplot or ROC curve visualization"),
    fluidRow(
      column(12, style = "margin-top: 0px;", triselector_ui(ns("tris_feature_general"), right_margin = "5")),
      column(11, uiOutput(ns("feature_general_plot"))),
      column(
        1,
        attr4selector_ui(ns("a4_gf"), circle = FALSE, right = TRUE),
        radioGroupButtons(
              inputId = ns("internal_radio"),
              label = " ",
              size = "xs",
              choices = c("Bees", "Curve"),
              direction = "vertical",
              status = "success"
            ) %>%
          tagAppendAttributes(`data-testid` = paste0(id, "-plot-type-selector"))
        )
    ),
    tags$h3("Statistical Results Table", class = "sr-only", `aria-label` = "Table showing statistical test results including p-values, effect sizes, and group comparisons for selected features"),
    dataTableDownload_ui(ns("mtab")),
    # Hidden text summary for AI browsers and screen readers
    div(class = "sr-only", `aria-live` = "polite", `aria-atomic` = "true",
        uiOutput(ns("plotSummary")))
  )
}

#' @description utility - feature general module
#' @param id module id
#' @param reactive_expr reactive expression matrix
#' @param reactive_i reactive row index to be highlighted
#' @param reactive_highlight reactive col index to be highlighted
#' @param reactive_phenoData reactive phenotype data
#' @param reactive_featureData reactive feature data
#' @param reactive_status saved status to restore
#' @importFrom DT renderDataTable
#' @importFrom reshape2 melt
#' @importFrom shinyWidgets radioGroupButtons updateRadioGroupButtons
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
#' #   v <- feature_general_module("tres",
#' #                   reactive_expr = reactive(exprs(dat)),
#' #                   # reactive_i = reactive(c(5, 6, 7)),
#' #                   reactive_i = reactive(c( 7)),
#' #                   reactive_highlight = reactive(c(3, 5, 10)),
#' #                   reactive_phenoData = reactive(pData(dat)),
#' #                   reactive_featureData = reactive(fData(dat)))
#' # }
#' #
#' # shinyApp(ui, server)


feature_general_module <- function(id,
                                   reactive_expr, reactive_i = reactive(NULL),
                                   reactive_highlight = reactive(NULL),
                                   reactive_phenoData,
                                   reactive_featureData,
                                   reactive_status = reactive(NULL)) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns
  
  # selector
  triset <- reactive({
    ts <- trisetter(expr = reactive_expr(), meta = reactive_phenoData(), combine = "pheno")
    ts[ts[, 1] != "Surv", ]
  })
  
  xax <- reactiveVal()
  v1 <- triselector_module(
    "tris_feature_general", reactive_x = triset, label = "Link to variable",
    reactive_selector1 = reactive(xax()$v1),
    reactive_selector2 = reactive(xax()$v2),
    reactive_selector3 = reactive(xax()$v3))

  attr4select_status <- reactiveVal()
  attr4select <- attr4selector_module(
    "a4_gf", reactive_meta = reactive_phenoData, reactive_expr = reactive_expr,
    reactive_triset = triset, reactive_status = attr4select_status
  )
  
  reactive_input <- reactive({
    req(reactive_expr())
    req(reactive_phenoData())
    e <- t(reactive_expr())
    colnames(e) <- paste0("Feature|Auto|", colnames(e))
    cbind(reactive_phenoData(), e)
    })

  # what to do
  pheno <- reactive({
    req(v1())
    cs <- do.call(paste, list(v1(), collapse = "|"))    
    if (!cs %in% colnames(reactive_input()))
      return(NULL)
    reactive_input()[, cs]
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
    if (showBeeswarm()) {
      if (input$internal_radio == "Bees")
        r <- plotly_scatter_ui(ns("feature_general_beeswarm")) else
          r <- plot_roc_pr_ui(ns("feature_general_roc_pr")) 
      r          
    }
  })

  
  rh <- reactiveVal()
  observeEvent(reactive_highlight(), {
    r <- reactive_highlight()
    if ( is.null(r) || is.logical(r) )
      return(NULL)
    
    if ( length(r) == 0 && !is.null(rh()) && length(rh()) > 0) {
      rh(integer(0))
    } else if (length(r) > 0) {
      rh( match(r, colnames(reactive_expr())) ) 
    }    
    })

  # boxplot:
  #  - no phenoData selected
  #  - multi feature selected - numerical phenoData selected - boxplot with external
  plotly_boxplot_module("feature_general_boxplotly",
             reactive_param_plotly_boxplot = reactive({
               req(reactive_expr())
               ylab <- rownames(reactive_expr())[reactive_i()]
               if (length(ylab) > 1)
                ylab <- "Abundance of selected features" else if (length(ylab) == 0)
                  ylab <- "Relative abundance" else if (is.na(ylab))
                    ylab <- "Relative abundance"
               ylab.extvar <- do.call(paste, list(v1(), collapse = "|"))
               list(
                 x = reactive_expr(), i = reactive_i(),
                 highlight = rh(),
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
    
    req(all(reactive_i() %in% rownames(reactive_expr())) || 
      all(reactive_i() <= nrow(reactive_expr())))
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
       ylab <- "Abundance of selected features" else if (length(ylab) == 0)
         ylab <- "Relative abundance" else if (is.na(ylab))
           ylab <- "Relative abundance" 
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
  v_scatter <- plotly_scatter_module("feature_general_scatter",
                          reactive_param_plotly_scatter = scatter_vars,
                          reactive_checkpoint = showScatter,
                          reactive_regLine = reactive( showRegLine()))
  observe({
    showRegLine(v_scatter()$regline)
    })

  ## beeswarm:
  # - single feature selected - categorical phenoData selected
  # - multi feature selected - categorical phenoData selected
  v_beeswarm <- plotly_scatter_module("feature_general_beeswarm",
                           reactive_param_plotly_scatter = scatter_vars,
                           reactive_checkpoint = showBeeswarm,
                           htest_var1 = htestV1, htest_var2 = htestV2)

  plot_roc_pr_module("feature_general_roc_pr",
    reactive_param = scatter_vars, reactive_checkpoint = showBeeswarm)
  
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

  dataTableDownload_module(
    "mtab", reactive_table = metatab, prefix = "FeatureTable_"
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

  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return()
    if (!is.null(s$plotType))
      updateRadioGroupButtons(session, "internal_radio", selected = s$plotType)
    })

  ## Hidden plot summary for AI browsers ##
  output$plotSummary <- renderUI({
    plot_type <- input$internal_radio

    if (is.null(plot_type)) {
      return(NULL)
    }

    summary_text <- if (plot_type == "Bees") {
      # Beeswarm/boxplot summary
      n_features <- length(reactive_i())
      n_samples <- ncol(reactive_expr())
      n_highlighted <- if (!is.null(reactive_highlight())) length(reactive_highlight()) else 0

      var_name <- if (!is.null(v1()) && !is.null(v1()$variable)) {
        v1()$variable
      } else {
        "selected variable"
      }

      parts <- c(
        sprintf("Boxplot showing expression distribution across %d samples.", n_samples)
      )

      if (n_features > 0) {
        parts <- c(parts, sprintf("%d feature(s) highlighted in the plot.", n_features))
      }

      if (n_highlighted > 0) {
        parts <- c(parts, sprintf("%d sample(s) highlighted.", n_highlighted))
      }

      if (!is.null(pheno()) && length(pheno()) > 0) {
        parts <- c(parts, sprintf("External variable %s overlaid on plot.", var_name))
      }

      paste(parts, collapse = " ")

    } else {
      # ROC/PR curve summary
      "ROC and precision-recall curves showing classification performance for binary outcome."
    }

    tags$p(summary_text)
  })

  ## return status ##
  rv <- reactiveValues()
  observe( rv$xax <- v1() )
  observe( rv$showRegLine <- showRegLine() )
  observe( rv$attr4 <- attr4select$status )
  observe( rv$plotType <- input$internal_radio )
  observe({
    rv$htestV1 <- v_beeswarm()$htestV1
    rv$htestV2 <- v_beeswarm()$htestV2
    })

  reactive({
    reactiveValuesToList(rv)
    })

  }) # end moduleServer
}
