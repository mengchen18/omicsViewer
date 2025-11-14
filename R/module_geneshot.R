#' Geneshot Search UI Function
#'
#' @description
#' Creates the user interface for the Geneshot gene search module. Allows users
#' to search for genes associated with arbitrary text queries using the Geneshot
#' API and visualize results by publication metrics.
#'
#' @param id Character. Namespace ID for the Shiny module. Must match the ID
#'   used in \code{\link{geneshot_module}}.
#'
#' @return
#' A \code{tagList} containing:
#' \itemize{
#'   \item Text input for search terms (semicolon-separated)
#'   \item Submit button to trigger search
#'   \item Gene ID mapper selector
#'   \item Interactive scatter plot of search results
#'   \item Results table with download functionality
#' }
#'
#' @family search modules
#' @seealso
#' \code{\link{geneshot_module}} for the corresponding server logic.
#' \code{\link{getAutoRIF}} for the underlying API function.
#'
#' @references
#' Lachmann et al. (2019) Geneshot: search engine for ranking genes from
#' arbitrary text queries. Nucleic Acids Research 47(W1):W571-W577.
#'
#' @keywords internal
#' @importFrom shinybusy show_modal_spinner remove_modal_spinner
#'
geneshot_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        width = 8, 
        textInput(ns("term"), label = NULL, value = "", width = "100%", 
          placeholder = "Search any term here, multiple separated by (;) ...")
        ),
      column(
        width = 4, align = "right",
        actionButton(ns("submit"), label = "Search related genes!", width = "100%")
        ),
      column(
        width = 11, triselector_ui(ns("geneNameCol")))
      ),    
    plotly::plotlyOutput(ns("plt")),
    tags$br(),
    dataTableDownload_ui(ns("autorif"))
  )
}

#' Geneshot Search Server Function
#'
#' @description
#' Server logic for the Geneshot gene search module. Queries the Geneshot API
#' to find genes associated with user-provided text terms, ranks results by
#' publication metrics, and highlights overlap with selected features.
#'
#' @param id Character. Namespace ID for the Shiny module. Must match the ID
#'   used in \code{\link{geneshot_ui}}.
#'
#' @param pdata Reactive expression. Returns a data.frame of sample metadata
#'   (phenotype data). Currently unused but required for module consistency.
#'
#' @param fdata Reactive expression. Returns a data.frame of feature metadata
#'   with features as rows. Used for ID mapping and overlap detection.
#'
#' @param expr Reactive expression. Returns a numeric matrix with features as
#'   rows and samples as columns. Currently unused but required for module
#'   consistency.
#'
#' @param feature_selected Reactive expression. Returns a character vector or
#'   integer vector of currently selected feature IDs/indices. Used to highlight
#'   overlap with search results.
#'
#' @param sample_selected Reactive expression. Returns selected sample IDs.
#'   Currently unused but required for module consistency.
#'
#' @param object Reactive expression. Returns the full ExpressionSet or
#'   SummarizedExperiment object. Currently unused but required for module
#'   consistency.
#'
#' @param reactive_status Reactive expression. Returns a list containing saved
#'   state for session restoration. Optional. Default: \code{reactive(NULL)}.
#'
#' @details
#' The module performs the following steps:
#' \enumerate{
#'   \item Accepts semicolon-separated search terms from user
#'   \item Queries Geneshot API via \code{\link{getAutoRIF}}
#'   \item Filters results to top 20 genes (GENESHOT_TOP_RESULTS)
#'   \item Highlights genes that overlap with currently selected features
#'   \item Displays results in interactive scatter plot and table
#' }
#'
#' Handles API errors gracefully with user notifications.
#'
#' @return
#' A reactive expression returning a list with module state for session saving.
#'
#' @family search modules
#' @seealso
#' \code{\link{geneshot_ui}} for the corresponding UI function.
#' \code{\link{getAutoRIF}} for the API call function.
#'
#' @keywords internal
#'
geneshot_module <- function(
  id, pdata, fdata, expr, feature_selected, sample_selected, object,
  reactive_status = reactive(NULL)
) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns

  triset <- reactive({
    fd <- fdata()
    cn <- colnames(fd)[!vapply(fd, is.numeric, logical(1)) & !grepl("^GS\\|", colnames(fd))]
    str_split_fixed(cn, "\\|", n = 3)
  })

  xax <- reactiveVal()
  v1 <- triselector_module(
    "geneNameCol", reactive_x = triset, label = "Map ID",
    reactive_selector1 = reactive(xax()$v1),
    reactive_selector2 = reactive(xax()$v2),
    reactive_selector3 = reactive(xax()$v3)
  )
  

  
  rif <- reactiveVal()
  observeEvent(input$submit, {
    show_modal_spinner(text = "Querying Geneshot database ...")
    res <- getAutoRIF(trimws(strsplit(input$term, ";")[[1]]), filter = TRUE)
    remove_modal_spinner()

    # Check for errors or empty results
    if (is.null(res)) {
      showNotification(
        "Geneshot query failed or returned no results. Please check your search terms and internet connection.",
        type = "error",
        duration = 10
      )
      return()
    }

    if (nrow(res) == 0) {
      showNotification(
        "No genes identified. Please try different search terms or check the spelling.",
        type = "warning",
        duration = 10
      )
      return()
    }

    # Process results
    res$selected <- ""
    res <- res[order(res$rank, decreasing = TRUE), ]
    dft <- res[seq_len(min(nrow(res), GENESHOT_TOP_RESULTS)), ]
    outliers <- list(
      x = dft$n,
      y = dft$perc,
      text = dft$gene,
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      ax = 10,
      ay = -20
    )

    rif(
      list(
      tab = res,
      outliers = outliers
    ))

    # Show success notification
    showNotification(
      sprintf("Found %d genes related to your search.", nrow(res)),
      type = "message",
      duration = 3
    )
  })

  outliersLabs <- reactive( rif()$outliers )
  
  gtab <- reactive({
    req(df <- rif()$tab)
    cid <- paste(unlist(v1()), collapse = "|")
    if ( cid %in% colnames(fdata()) )
      df$selected[ df$gene %in% fdata()[feature_selected(), cid] ] <- "+"    
    df
    })  

  rifRow <- dataTableDownload_module(
    "autorif", reactive_table = reactive({req(gtab()); gtab()}), prefix = "autoRIF", pageLength = DEFAULT_TABLE_PAGE_LENGTH)

  ##
  output$plt <- plotly::renderPlotly({
    req(df <- gtab())
    fig <- plotly_scatter(
      x = df$n, y = df$perc, xlab = "# of publication", ylab = "Publications with Search Term(s) / Total Publications",
      color = df$selected, size = 10, tooltips=df$gene, shape = "select"
    )
    fig <- plotly::layout(fig$fig, annotations = outliersLabs()) #rif()$outliers)
    plotly::config(
      fig,
      toImageButtonOptions = list(
        format = "svg",
        filename = "omicsViewerPlot",
        width = 700,
        height = 700
        )
      )
  })

  ### status save and restore

  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return()
    xax(NULL)
    xax(list(v1 = s$xax[[1]], v2 = s$xax[[2]], v3 = s$xax[[3]]))
    rif(s$rif)
    updateTextInput(session, inputId = "term", value = s$term)
    })

  rv <- reactiveValues()
  observe(rv$xax <- v1())
  observe(rv$rif <- rif())
  observe(rv$term <- input$term)

  reactive({
    reactiveValuesToList(rv)
    })

  }) # end moduleServer
}