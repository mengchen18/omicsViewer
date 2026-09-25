#' STRING Network Analysis UI Function
#'
#' @description
#' Creates the user interface for the STRING protein-protein interaction network
#' analysis module. Displays network visualization, enrichment results, and
#' taxonomy configuration.
#'
#' @param id Character. Namespace ID for the Shiny module. Must match the ID
#'   used in \code{\link{string_module}}.
#'
#' @return
#' A \code{tagList} containing:
#' \itemize{
#'   \item Taxonomy code input
#'   \item Feature count display with max limit warning
#'   \item Action button to trigger analysis
#'   \item Enrichment results table with download
#'   \item Interactive network visualization with label toggle
#' }
#'
#' @family network modules
#' @seealso
#' \code{\link{string_module}} for the corresponding server logic.
#' \code{\link{stringNetwork}} for network API calls.
#' \code{\link{stringGSA}} for enrichment analysis.
#'
#' @keywords internal
#' @importFrom shinybusy show_modal_spinner remove_modal_spinner
#'
string_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Module description for AI browsers and screen readers
    div(class = "sr-only", id = ns("module-help"),
      tags$h4("About STRING Protein-Protein Interaction Network"),
      tags$p("STRING (Search Tool for the Retrieval of Interacting Genes/Proteins) is a database of known and predicted protein-protein interactions. This module queries the STRING database to build an interaction network for your selected proteins, showing how they might work together in biological processes."),
      tags$h4("When to use STRING"),
      tags$p("Use this analysis when you want to understand how your proteins of interest interact with each other, identify protein complexes, or discover hub proteins that connect many others. This is especially useful for interpreting proteomics data or understanding the molecular machinery involved in a biological process."),
      tags$h4("How to interpret results"),
      tags$p("The network visualization shows proteins as nodes and interactions as edges connecting them. The interaction table lists all protein pairs with confidence scores (0-1, higher is better) and evidence types (experimental, database, text mining, etc.). Connected clusters of proteins often represent functional modules or complexes. You can also see pathway enrichment results to understand what biological processes the network is involved in.")
    ),
    tags$h3("Query Configuration", class = "sr-only", `aria-label` = "Input for organism taxonomy code and run button to query STRING database. Maximum 300 proteins allowed"),
    fluidRow(
      column(
        4, offset = 0, style='padding-left:15px; padding-right:2px; padding-top:0px; padding-bottom:0px',
        textInputIcon(inputId = ns("tax"), label = NULL, value = "9606", icon = list("Taxomony Code")) %>%
          tagAppendAttributes(`data-testid` = paste0(id, "-taxonomy-input"))),
      column(
        6, offset = 0, style='padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px',
        verbatimTextOutput(ns("error.msg")) %>%
          tagAppendAttributes(`aria-live` = "assertive")),
      column(
        2, offset = 0, style='padding-left:15px; padding-right:2px; padding-top:0px; padding-bottom:0px',
        actionButton(ns("run"), "Run!") %>%
          tagAppendAttributes(`data-testid` = paste0(id, "-run-analysis-button")))
    ),
    uiOutput(ns("noresRet")),
    tags$h3("Interaction Results Table", class = "sr-only", `aria-label` = "Table of protein-protein interactions from STRING database with confidence scores and evidence types"),
    dataTableDownload_ui(ns("strtab")),
    tags$h3("Network Visualization", class = "sr-only", `aria-label` = "Interactive force-directed network graph showing protein interactions with node labels toggle"),
    checkboxInput(ns("showLabel"), label = "Show labels", value = FALSE) %>%
      tagAppendAttributes(`data-testid` = paste0(id, "-show-labels-toggle")),
    forceNetworkOutput(ns("network")),
    # Hidden text summary for AI browsers
    div(class = "sr-only", `aria-live` = "polite", `aria-atomic` = "true",
        uiOutput(ns("networkSummary")))
  )
}

#' STRING Network Analysis Server Function
#'
#' @description
#' Server logic for the STRING protein-protein interaction network analysis
#' module. Queries the STRING database API to retrieve protein interactions
#' and perform functional enrichment analysis on selected features.
#'
#' @param id Character. Namespace ID for the Shiny module. Must match the ID
#'   used in \code{\link{string_ui}}.
#'
#' @param reactive_ids Reactive expression. Returns a character vector of
#'   protein/gene identifiers to query. Maximum of 300 features allowed
#'   (enforced by STRING API limits).
#'
#' @param reactive_status Reactive expression. Returns a list containing
#'   saved state for session restoration (tax ID, label settings). Optional.
#'   Default: 
#'
#' @param store Optional child view of the canonical widget store
#'   (\code{\link{widget_store_child}}) for this module, e.g.
#'   \code{resultspace.stringdb}. When given, the taxonomy input, the
#'   labels checkbox, and the enrichment-table row selection (which
#'   highlights network nodes) register on the store; when NULL the legacy
#'   status-restore path is kept. The Run button is deliberately not
#'   registered: it is a stateless command (an API query trigger), not a
#'   widget value.\code{reactive(NULL)}.
#'
#' @param active Reactive expression. Returns a logical indicating whether
#'   the module should auto-run on initialization. Used for session restoration.
#'   Default: \code{reactive(FALSE)}.
#'
#' @details
#' ## API Interaction
#' The module makes two separate API calls to STRING database:
#' \enumerate{
#'   \item Network query: Retrieves protein-protein interactions
#'   \item Enrichment query: Performs gene set enrichment analysis
#' }
#'
#' Both queries are rate-limited and may fail with network errors. The module
#' handles errors gracefully with user notifications.
#'
#' ## Feature Limits
#' - Maximum 300 genes (STRING_MAX_GENES constant)
#' - Maximum 999 network edges displayed (STRING_MAX_NETWORK_EDGES constant)
#'
#' @return
#' A reactive expression returning a list with:
#' \itemize{
#'   \item tax: Taxonomy ID used
#'   \item showLabel: Logical, whether labels are shown in network
#' }
#'
#' @family network modules
#' @seealso
#' \code{\link{string_ui}} for the corresponding UI function.
#' \code{\link{stringNetwork}} for network retrieval.
#' \code{\link{stringGSA}} for enrichment analysis.
#' \code{\link{stringD3Net}} for network visualization.
#'
#' @keywords internal
#' @importFrom networkD3 renderForceNetwork forceNetworkOutput
#'
string_module <- function(
  id, reactive_ids, reactive_status = reactive(NULL), active = reactive(FALSE),
  store = NULL
) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns
  
  # ------------------------------------------------------------------
  # Canonical widget-store bindings (control plane, plan section 6, S4).
  # Taxonomy code, network-label toggle, and the enrichment row selection
  # (drives the network highlight) register on the store. The Run button
  # stays unregistered (stateless command). Observer retention mandatory.
  # ------------------------------------------------------------------
  .str_store_observers <- list()
  .str_keep <- function(obs) {
    .str_store_observers[[length(.str_store_observers) + 1L]] <<- obs
    invisible(obs)
  }
  if (!is.null(store)) {
    store_register(
      store,
      widget_binding("taxonomy", "string", label = "Taxonomy code",
        help = paste("NCBI taxonomy identifier of the organism queried in",
                     "the STRING database (9606 = human)")),
      widget_binding("show_labels", "boolean", label = "Show network labels",
        help = "Draw gene labels on the interaction network")
    )
    .str_root_store <- if (is.null(store$parent)) store else store$parent
    .str_keep(observeEvent(input$tax, {
      if (!is.null(input$tax) && nzchar(input$tax))
        store_sync_from_ui(store, "taxonomy", input$tax)
    }, ignoreInit = TRUE))
    .str_keep(observeEvent(input$showLabel, {
      if (!is.null(input$showLabel))
        store_sync_from_ui(store, "show_labels", input$showLabel)
    }, ignoreInit = TRUE))
    # seed once the inputs exist (restore-first-wins)
    .str_seeded <- FALSE
    .str_keep(observe({
      if (.str_seeded) return(NULL)
      if (is.null(input$tax) || is.null(input$showLabel)) return(NULL)
      .str_seeded <<- TRUE
      store_seed(store, list(taxonomy = input$tax,
                             show_labels = input$showLabel))
    }))
    # store -> UI push for external writes only (pending entries mark them)
    .str_epoch <- store_epoch(store)
    .str_keep(observe({
      .str_epoch()
      tx <- store_read(store, "taxonomy")[[1]]
      if (!is.null(tx) &&
          !is.null(.str_root_store$pending[[paste0(store$prefix, ".taxonomy")]]))
        updateTextInputIcon(session, "tax", value = tx)
      lb <- store_read(store, "show_labels")[[1]]
      if (!is.null(lb) &&
          !is.null(.str_root_store$pending[[paste0(store$prefix, ".show_labels")]]))
        updateCheckboxInput(session, "showLabel", value = lb)
    }))
  }
  
  overflow <- reactive({
    length(reactive_ids()) > STRING_MAX_GENES
  })

  output$error.msg <- renderText({
    sprintf("%s features selected [MAX %d FEATURES ALLOWED!]",
            length(reactive_ids()), STRING_MAX_GENES)
  })
  
  nk <- reactiveVal()
  observeEvent( input$run, {
    # Validate inputs before API call
    ids <- reactive_ids()
    validation <- validate_character_vector(
      ids,
      name = "selected features",
      min_length = 1,
      max_length = STRING_MAX_GENES,
      allow_empty = FALSE
    )

    if (!validation$valid) {
      showNotification(
        validation$message,
        type = "error",
        duration = 10
      )
      return()
    }

    # Validate taxonomy ID
    if (is.null(input$tax) || nchar(input$tax) == 0) {
      showNotification(
        "Taxonomy ID is required",
        type = "error",
        duration = 5
      )
      return()
    }

    show_modal_spinner(text = "Querying STRING network ...")
    r <- stringNetwork(genes = ids, taxid = input$tax)
    remove_modal_spinner()

    # Check for API errors
    if (inherits(r, "character")) {
      # Error message returned
      showNotification(
        paste("STRING Network Error:", r),
        type = "error",
        duration = 10
      )
      nk(r)
      return()
    }

    if (is.data.frame(r)) {
      if (nrow(r) > STRING_MAX_NETWORK_EDGES) {
        r <- r[order(r$score, decreasing = TRUE), ]
        r <- r[seq_len(STRING_MAX_NETWORK_EDGES), ]
        showNotification(
          sprintf("Showing top %d network edges (from %d total)",
                  STRING_MAX_NETWORK_EDGES, nrow(r)),
          type = "warning",
          duration = 5
        )
      }
    }
    nk(r)
  })
  
  gs <- reactiveVal()
  gs <- eventReactive( input$run, {
    req(!overflow())
    show_modal_spinner(text = "Querying STRING enrichment database ...")
    tab <- stringGSA(genes = reactive_ids(), taxid = input$tax)
    remove_modal_spinner()

    # Check for API errors
    if (inherits(tab, "character")) {
      # Error message returned
      showNotification(
        paste("STRING Enrichment Error:", tab),
        type = "error",
        duration = 10
      )
      return(tab)
    }

    if (!is.data.frame(tab)) {
      showNotification(
        "STRING Enrichment returned unexpected data format.",
        type = "error",
        duration = 10
      )
      return("Error: Unexpected data format")
    }

    colnames(tab) <- c(
      "category", "term", "gene number", "background number",
      "TaxonId", "inputGenes", "preferredNames", "p value", "fdr", "description")
    tab
  })

  nores <- reactive( {
    !is.data.frame(nk()) || !is.data.frame(gs()) 
  })
  
  output$nores.msg <- renderText({
    nores()
    req(nores())
    c(gs(), nk())[c(!is.data.frame(gs()), !is.data.frame(nk()))]
  })
  
  output$noresRet <- renderUI({    
    verbatimTextOutput(ns("nores.msg"))
  })
  
  highlightP <- reactiveVal(1)
  
  output$network <- renderForceNetwork({
    req(!nores())
    req(nrow(nk()) > 0)
    stringD3Net(ntwk = nk(), gsa = gs(), i = highlightP(), label = input$showLabel)
  })
  
  strtab_df <- eventReactive(gs(), {
    req(!nores())
    req(!overflow())
    req(nrow(gs()) > 0)
    gs()[, c("category", "term", "gene number", "background number", "p value", "fdr", "description")]
  })
  
  tt <- dataTableDownload_module(
    "strtab", reactive_table = strtab_df, prefix = "FeatureTable_",
    reactive_row_ids = reactive({
      df <- tryCatch(strtab_df(), shiny.silent.error = function(e) NULL,
                     error = function(e) NULL)
      if (is.data.frame(df) && "term" %in% colnames(df))
        as.character(df$term) else NULL
    }),
    store = store, store_key = "selected_row",
    store_label = "Selected enrichment term",
    store_help = paste("STRING enrichment term whose row is selected;",
                       "highlights its genes in the network")
  )

  observe({
    req(tt())
    highlightP(tt())
  })

  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return()
    if (!is.null(store)) {
      # single transactional path (per-key resilient, meta_scatter style)
      patch <- list()
      if (!is.null(s$tax) && nzchar(s$tax)) patch$taxonomy <- s$tax
      if (!is.null(s$showLabel)) patch$show_labels <- s$showLabel
      if (length(patch))
        tryCatch(store_apply(store, patch, origin = "restore", strict = FALSE),
                 error = function(e) NULL)
    } else {
      updateTextInputIcon(session, "tax", value = s$tax)
      updateCheckboxInput(session, "showLabel", value = s$showLabel)
    }
    if (active())
      shinyjs::click("run")
    })

  reactive({
    list(tax = input$tax, showLabel = input$showLabel)
  })

  }) # end moduleServer
}


