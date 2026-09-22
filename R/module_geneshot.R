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
    # Module description for AI browsers and screen readers
    div(class = "sr-only", id = ns("module-help"),
      tags$h4("About Geneshot Literature Search"),
      tags$p("Geneshot is a search engine that finds genes and proteins associated with arbitrary text queries by mining the scientific literature. It ranks genes by how frequently they co-occur with your search terms in PubMed abstracts and other biomedical texts."),
      tags$h4("When to use Geneshot"),
      tags$p("Use Geneshot when you want to discover genes related to a disease, phenotype, biological process, or any concept of interest. It's particularly useful for literature-based hypothesis generation, finding candidate genes for follow-up studies, or understanding which genes are most studied in relation to your topic of interest."),
      tags$h4("How to interpret results"),
      tags$p("Genes are ranked by their association score with your search terms. The scatter plot shows publication count (x-axis) versus percentage of publications mentioning your search terms (y-axis). Genes in the upper right are both well-studied and highly associated with your query. The plus symbol (+) in the 'selected' column indicates genes that overlap with your currently selected features in the app.")
    ),
    fluidRow(
      column(
        width = 8,
        textInput(ns("term"), label = NULL, value = "", width = "100%",
          placeholder = "Search any term here, multiple separated by (;) ...") %>%
          tagAppendAttributes(`data-testid` = paste0(id, "-search-input"))
        ),
      column(
        width = 4, align = "right",
        actionButton(ns("submit"), label = "Search related genes!", width = "100%") %>%
          tagAppendAttributes(`data-testid` = paste0(id, "-search-button"))
        ),
      column(
        width = 11, triselector_ui(ns("geneNameCol")))
      ),
    plotly::plotlyOutput(ns("plt")),
    tags$br(),
    dataTableDownload_ui(ns("autorif")),
    # Hidden text summary for AI browsers
    div(class = "sr-only", `aria-live` = "polite", `aria-atomic` = "true",
        uiOutput(ns("searchSummary")))
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
#' @param store Optional child view of the canonical widget store
#'   (\code{\link{widget_store_child}}) for this module, e.g.
#'   \code{resultspace.geneshot}. When given, the search-term input and
#'   the ID-mapper cascade register on the store; when NULL the legacy
#'   status-restore path is kept. The Search button is deliberately not
#'   registered: it is a stateless command, not a widget value.
#'
geneshot_module <- function(
  id, pdata, fdata, expr, feature_selected, sample_selected, object,
  reactive_status = reactive(NULL), store = NULL
) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns

  triset <- reactive({
    fd <- fdata()
    cn <- colnames(fd)[!vapply(fd, is.numeric, logical(1)) & !grepl("^GS\\|", colnames(fd))]
    str_split_fixed(cn, "\\|", n = 3)
  })

  # ------------------------------------------------------------------
  # Canonical widget-store bindings (control plane, plan section 6, S4).
  # The search term and the ID-mapper cascade register on the store; the
  # cascade is driven by store_watch selectors (meta_scatter pattern).
  # ------------------------------------------------------------------
  .gs_store_observers <- list()
  .gs_keep <- function(obs) {
    .gs_store_observers[[length(.gs_store_observers) + 1L]] <<- obs
    invisible(obs)
  }
  xax <- reactiveVal()
  if (!is.null(store)) {
    .gs_ts <- function() {
      fd <- tryCatch(fdata(), shiny.silent.error = function(e) NULL,
                     error = function(e) NULL)
      if (is.null(fd)) return(NULL)
      cn <- colnames(fd)[!vapply(fd, is.numeric, logical(1)) &
                           !grepl("^GS\\|", colnames(fd))]
      if (length(cn) == 0L) return(NULL)
      str_split_fixed(cn, "\\|", n = 3)
    }
    .gs_ts1 <- function() {
      ts <- .gs_ts()
      if (is.null(ts)) character(0) else unique(ts[, 1])
    }
    .gs_ts2 <- function(a) {
      ts <- .gs_ts()
      if (is.null(ts) || is.null(a) || !nzchar(a)) character(0)
      else unique(ts[ts[, 1] %in% a, 2])
    }
    .gs_ts3 <- function(a, b) {
      ts <- .gs_ts()
      if (is.null(ts)) character(0)
      else {
        i <- rep(TRUE, nrow(ts))
        if (!is.null(a) && nzchar(a)) i <- i & ts[, 1] %in% a
        if (!is.null(b) && nzchar(b)) i <- i & ts[, 2] %in% b
        unique(ts[i, 3])
      }
    }
    kg1 <- paste0(store$prefix, ".xax_analysis")
    kg2 <- paste0(store$prefix, ".xax_subset")
    store_register(
      store,
      widget_binding("term", "string", label = "Search terms",
        help = paste("Semicolon-separated literature search terms for the",
                     "Geneshot AutoRIF query")),
      widget_binding("xax_analysis", "select", label = "ID mapper category",
        help = "Annotation category used to map Geneshot genes to features",
        choices_provider = function(v) .gs_ts1()),
      widget_binding("xax_subset", "select", label = "ID mapper subcategory",
        help = "Subcategory within the ID-mapper category",
        depends_on = "xax_analysis",
        choices_provider = function(v) .gs_ts2(v[[kg1]])),
      widget_binding("xax_variable", "select_cascaded", label = "ID mapper column",
        help = paste("Feature annotation column Geneshot genes are matched",
                     "against for the overlap highlight"),
        depends_on = c("xax_analysis", "xax_subset"),
        choices_provider = function(v) .gs_ts3(v[[kg1]], v[[kg2]]))
    )
    v1 <- triselector_module(
      "geneNameCol", reactive_x = triset, label = "Map ID",
      reactive_selector1 = store_watch(store, "xax_analysis"),
      reactive_selector2 = store_watch(store, "xax_subset"),
      reactive_selector3 = store_watch(store, "xax_variable"),
      reactive_axis_request = store_epoch(store))
    .gs_read_tris <- function(sel)
      tryCatch(sel(), shiny.silent.error = function(e) NULL,
               error = function(e) NULL)
    .gs_component_set <- function(sel)
      !is.null(sel) &&
        nzchar(sel$analysis %||% "") && !identical(sel$analysis, "--select--") &&
        nzchar(sel$subset %||% "") && !identical(sel$subset, "--select--") &&
        nzchar(sel$variable %||% "") && !identical(sel$variable, "--select--")
    .gs_keep(observe({
      xv <- .gs_read_tris(v1)
      if (.gs_component_set(xv)) {
        store_sync_from_ui(store, "xax_analysis", xv$analysis)
        store_sync_from_ui(store, "xax_subset", xv$subset)
        store_sync_from_ui(store, "xax_variable", xv$variable)
      }
    }))
    .gs_keep(observeEvent(input$term, {
      if (!is.null(input$term) && nzchar(input$term))
        store_sync_from_ui(store, "term", input$term)
    }, ignoreInit = TRUE))
    # seed the term once it exists and is non-empty; the cascade has no
    # meaningful default and fills on the first pick
    .gs_seeded <- FALSE
    .gs_keep(observe({
      if (.gs_seeded) return(NULL)
      if (is.null(input$term)) return(NULL)
      .gs_seeded <<- TRUE
      held <- store_read(store, "term")
      if (is.null(held[[paste0(store$prefix, ".term")]]) &&
          nzchar(input$term))
        tryCatch(store_apply(store, list(term = input$term),
                             origin = "system", strict = FALSE),
                 error = function(e) NULL)
    }))
    # store -> UI push for external writes only (pending entries mark them)
    .gs_epoch <- store_epoch(store)
    .gs_root_store <- if (is.null(store$parent)) store else store$parent
    .gs_keep(observe({
      .gs_epoch()
      tm <- store_read(store, "term")[[1]]
      if (!is.null(tm) &&
          !is.null(.gs_root_store$pending[[paste0(store$prefix, ".term")]]))
        updateTextInput(session, "term", value = tm)
    }))
  } else {
    v1 <- triselector_module(
      "geneNameCol", reactive_x = triset, label = "Map ID",
      reactive_selector1 = reactive(xax()$v1),
      reactive_selector2 = reactive(xax()$v2),
      reactive_selector3 = reactive(xax()$v3)
    )
  }

  
  rif <- reactiveVal()
  observeEvent(input$submit, {
    # Validate input terms
    if (is.null(input$term) || nchar(trimws(input$term)) == 0) {
      showNotification(
        "Please enter at least one search term",
        type = "warning",
        duration = 5
      )
      return()
    }

    terms <- trimws(strsplit(input$term, ";")[[1]])
    terms <- terms[nchar(terms) > 0]  # Remove empty terms

    if (length(terms) == 0) {
      showNotification(
        "No valid search terms provided",
        type = "warning",
        duration = 5
      )
      return()
    }

    show_modal_spinner(text = "Querying Geneshot database ...")
    res <- getAutoRIF(terms, filter = TRUE)
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
    if (!is.null(store)) {
      # single transactional path (per-key resilient, meta_scatter style)
      patch <- list()
      if (!is.null(s$term) && nzchar(s$term)) patch$term <- s$term
      if (identical(length(s$xax), 3L)) {
        tr <- lapply(s$xax, function(x)
          if (is.null(x) || !nzchar(x) || identical(x, "--select--")) NULL else x)
        if (!is.null(tr[[1]])) patch$xax_analysis <- tr[[1]]
        if (!is.null(tr[[2]])) patch$xax_subset <- tr[[2]]
        if (!is.null(tr[[3]])) patch$xax_variable <- tr[[3]]
      }
      if (length(patch))
        tryCatch(store_apply(store, patch, origin = "restore", strict = FALSE),
                 error = function(e) NULL)
    } else {
      xax(NULL)
      xax(list(v1 = s$xax[[1]], v2 = s$xax[[2]], v3 = s$xax[[3]]))
      updateTextInput(session, inputId = "term", value = s$term)
    }
    })

  rv <- reactiveValues()
  observe(rv$xax <- v1())
  observe(rv$term <- input$term)

  reactive({
    reactiveValuesToList(rv)
    })

  }) # end moduleServer
}