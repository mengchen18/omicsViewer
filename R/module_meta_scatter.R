#' Meta Scatter Plot UI Function
#'
#' @description
#' Creates the user interface for the metadata scatter plot visualization module.
#' Provides interactive scatter plots with advanced selection tools including
#' corner selection for volcano plots.
#'
#' @param id Character. Namespace ID for the Shiny module. Must match the ID
#'   used in \code{\link{meta_scatter_module}}.
#'
#' @return
#' A \code{tagList} containing:
#' \itemize{
#'   \item Figure attribute selector (color, shape, size controls)
#'   \item Clear selection button
#'   \item X-axis and Y-axis variable selectors
#'   \item Interactive plotly scatter plot with lasso/box selection
#' }
#'
#' @family visualization modules
#' @seealso
#' \code{\link{meta_scatter_module}} for the corresponding server logic.
#' \code{\link{plotly_scatter_module}} for the scatter plot implementation.
#'
#' @keywords internal
#' @importFrom shinyWidgets actionBttn
#'
meta_scatter_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$h3("Plot Controls and Variable Selection", class = "sr-only", `aria-label` = "Controls for customizing scatter plot appearance including color, shape, size mapping and selecting X and Y axis variables"),
    fluidRow(
      column(1,
        attr4selector_ui(ns("a4selector")),
        actionBttn(ns("clear"), "Clear figure selection", style = "minimal", color = "primary", size = "xs") %>%
          tagAppendAttributes(`data-testid` = paste0(id, "-clear-selection-button"))
      ), # style = "margin-top: 20px;",
      column(11,
             triselector_ui(ns("tris_main_scatter1")),
             triselector_ui(ns("tris_main_scatter2")))
    ),
    tags$h3("Interactive Scatter Plot Visualization", class = "sr-only", `aria-label` = "Scatter plot with lasso and box selection tools, regression line option, and corner selection for volcano plots"),
    plotly_scatter_ui(ns("main_scatterOutput"), height = META_SCATTER_PLOT_HEIGHT),
    # Hidden text summary for AI browsers and screen readers
    div(class = "sr-only", `aria-live` = "polite", `aria-atomic` = "true",
        uiOutput(ns("plotSummary")))
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
  id, reactive_meta=reactive(NULL), reactive_expr=reactive(NULL),
  combine = c("pheno", "feature"), source = "plotlyscattersource",
  reactive_x = reactive(NULL), reactive_y = reactive(NULL),
  reactive_status = reactive(NULL)
) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns

  # Helper: Get feature/sample names based on combine mode
  get_names <- function() {
    if (combine == "pheno") {
      colnames(reactive_expr())
    } else {
      rownames(reactive_expr())
    }
  }

  triset <- reactive( {
    ts <- trisetter(expr = reactive_expr(), meta = reactive_meta(), combine = combine[1])
    ts[ts[, 1] != "Surv", ]
  } )

  # Axis config: Parse pipe-separated strings into list(v1, v2, v3)
  # Use reactiveVal pattern (not pure reactive) to prevent multiple invalidations
  xax <- reactiveVal()
  observe({
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

  v1 <- triselector_module("tris_main_scatter1", reactive_x = triset, label = "X-axis",
                   reactive_selector1 = reactive(xax()$v1),
                   reactive_selector2 = reactive(xax()$v2),
                   reactive_selector3 = reactive(xax()$v3))
  v2 <- triselector_module("tris_main_scatter2", reactive_x = triset, label = "Y-axis",
                   reactive_selector1 = reactive(yax()$v1),
                   reactive_selector2 = reactive(yax()$v2),
                   reactive_selector3 = reactive(yax()$v3))

  # Detect volcano plot: x=mean.diff, y=log.fdr/log.pvalue (both from ttest)
  pre_vol <- reactive({
    # Check if selections are valid
    if (is.null(v1()) || is.null(v2())) return(FALSE)
    if (is.null(v1()$analysis) || is.null(v2()$analysis)) return(FALSE)
    if (is.null(v1()$variable) || is.null(v2()$variable)) return(FALSE)

    # Check for volcano plot pattern
    v1()$analysis == "ttest" &&
      v2()$analysis == "ttest" &&
      v1()$variable == "mean.diff" &&
      v2()$variable %in% c("log.fdr", "log.pvalue")
  })

  attr4select_status <- reactiveVal()
  attr4select <- attr4selector_module(
    "a4selector", reactive_meta = reactive_meta, reactive_expr = reactive_expr,
    reactive_triset = triset, pre_volcano = pre_vol, reactive_status = attr4select_status
  )

  xycoord <- reactive({
    req(v1()$variable)
    req(v2()$variable)
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

  # Track clear button clicks
  clear_counter <- reactiveVal(0)
  observeEvent(input$clear, {
    clear_counter(clear_counter() + 1)
  })

  # Rectangle for corner selection (volcano plot)
  rectval <- reactive({
    # Recalculate when clear button clicked
    clear_counter()

    # Return NULL if no cutoff selected or "None" corner
    if (is.null(attr4select$cutoff) || attr4select$cutoff$corner == "None") {
      return(NULL)
    }
    
    # Force dependency on axis changes
    v1()
    v2()

    # Calculate rectangle based on coordinates and cutoff
    # Note: xycoord() already depends on v1() and v2(), so we don't need
    # to touch them explicitly - reactive graph handles transitive dependencies
    coords <- xycoord()
    if (is.null(coords)) return(NULL)

    line_rect(l = attr4select$cutoff, coords)$rect
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
    l$inSelection <- NA
    l
  })

  showRegLine <- reactiveVal(FALSE)
  htestV1 <- reactiveVal()
  htestV2 <- reactiveVal()
  v_scatter <- plotly_scatter_module(
    "main_scatterOutput", reactive_param_plotly_scatter = scatter_vars,
    reactive_regLine = showRegLine, htest_var1 = htestV1, htest_var2 = htestV2)
  observe({
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
  
  # Workaround: Track previous selection to prevent redundant updates
  # Plotly events can fire even when selection hasn't actually changed,
  # causing unnecessary reactive chain invalidations. We store the previous
  # selection and only update selVal when it truly changes.
  clientSideSelection <- reactiveVal(character(0))
  observeEvent( v_scatter(), {
    l <- get_names()
    u_c <- l[v_scatter()$clicked]
    u_s <- l[v_scatter()$selected]

    # Only update if selection actually changed
    req( !identical( tmp <- c(u_c, u_s),  clientSideSelection() ) )

    clientSideSelection(tmp)
    selVal( list(
      clicked = u_c,
      selected = u_s
    ) )
    sbc(FALSE)
  })

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
    l <- get_names()

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
  # Consolidate all status restoration into single observer
  observeEvent(reactive_status(), {
    s <- reactive_status()
    if (is.null(s)) return()

    # Restore axis selections
    xax(list(v1 = s$xax[[1]], v2 = s$xax[[2]], v3 = s$xax[[3]]))
    yax(list(v1 = s$yax[[1]], v2 = s$yax[[2]], v3 = s$yax[[3]]))

    # Restore attribute selector status
    attr4select_status(NULL)
    attr4select_status(s$attr4)

    # Restore regression line setting
    showRegLine(s$showRegLine)

    # Restore hypothesis test variables
    htestV1(s$htestV1)
    htestV2(s$htestV2)

    # Restore selection state
    returnCornerSelection(s$selectByCorner)
    selVal(list(
      clicked = s$selection_clicked,
      selected = s$selection_selected
    ))
  })
  #############################################

  # Generate hidden text summary for AI browsers and screen readers
  output$plotSummary <- renderUI({
    req(scatter_vars())
    vars <- scatter_vars()

    # Calculate basic statistics
    x_vals <- vars$x
    y_vals <- vars$y
    n_points <- length(x_vals)

    # Build summary text
    summary_parts <- c(
      sprintf("Scatter plot visualization with %d data points.", n_points),
      sprintf("X-axis: %s.", vars$xlab %||% "Variable"),
      sprintf("Y-axis: %s.", vars$ylab %||% "Variable")
    )

    # Add numeric range info for numeric axes
    if (is.numeric(x_vals)) {
      summary_parts <- c(summary_parts,
        sprintf("X-axis range: %.3g to %.3g.", min(x_vals, na.rm = TRUE), max(x_vals, na.rm = TRUE)))
    }
    if (is.numeric(y_vals)) {
      summary_parts <- c(summary_parts,
        sprintf("Y-axis range: %.3g to %.3g.", min(y_vals, na.rm = TRUE), max(y_vals, na.rm = TRUE)))
    }

    # Add correlation for numeric-numeric plots
    if (is.numeric(x_vals) && is.numeric(y_vals)) {
      cor_result <- tryCatch({
        cor.test(x_vals, y_vals, use = "complete.obs")
      }, error = function(e) NULL)

      if (!is.null(cor_result)) {
        summary_parts <- c(summary_parts,
          sprintf("Pearson correlation: r = %.3f, p-value = %.3g.",
                  cor_result$estimate, cor_result$p.value))
      }
    }

    # Add selection info
    sel <- selVal()
    if (length(sel$selected) > 0) {
      summary_parts <- c(summary_parts,
        sprintf("%d points currently selected.", length(sel$selected)))
    }

    tags$p(paste(summary_parts, collapse = " "))
  })

  selVal

  }) # end moduleServer
}
