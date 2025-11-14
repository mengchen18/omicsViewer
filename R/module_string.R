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
    fluidRow(
      column(
        4, offset = 0, style='padding-left:15px; padding-right:2px; padding-top:0px; padding-bottom:0px',
        textInputIcon(inputId = ns("tax"), label = NULL, value = "9606", icon = list("Taxomony Code"))),
      column(
        6, offset = 0, style='padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px',
        verbatimTextOutput(ns("error.msg"))),
      column(
        2, offset = 0, style='padding-left:15px; padding-right:2px; padding-top:0px; padding-bottom:0px',
        actionButton(ns("run"), "Run!"))
    ),
    uiOutput(ns("noresRet")),
    dataTableDownload_ui(ns("strtab")),
    checkboxInput(ns("showLabel"), label = "Show labels", value = FALSE),
    forceNetworkOutput(ns("network"))    
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
#'   Default: \code{reactive(NULL)}.
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
  id, reactive_ids, reactive_status = reactive(NULL), active = reactive(FALSE)
) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns
  
  overflow <- reactive({
    length(reactive_ids()) > STRING_MAX_GENES
  })

  output$error.msg <- renderText({
    sprintf("%s features selected [MAX %d FEATURES ALLOWED!]",
            length(reactive_ids()), STRING_MAX_GENES)
  })
  
  nk <- reactiveVal()
  # nk <- eventReactive( input$run, {
  observeEvent( input$run, {
    show_modal_spinner(text = "Querying STRING network ...")
    r <- stringNetwork(genes = reactive_ids(), taxid = input$tax)# reactive_taxid())
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
  # observeEvent(input$run, {
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
    # gs(tab)
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
  
  tt <- dataTableDownload_module(
    "strtab", reactive_table = eventReactive(gs(), {
      req(!nores())
      req(!overflow())
      req(nrow(gs()) > 0)
      gs()[, c("category", "term", "gene number", "background number", "p value", "fdr", "description")]
    }), prefix = "FeatureTable_"
  )

  observe({
    req(tt())
    highlightP(tt())
  })

  # return and restore session
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return()
    # nk(s$nk)
    # gs(s$gs)
    updateTextInputIcon(session, "tax", value = s$tax)
    updateCheckboxInput(session, "showLabel", value = s$showLabel)
    if (active())
      shinyjs::click("run")
    })

  reactive({
    list( tax = input$tax, showLabel = input$showLabel ) # nk = nk(), gs = gs(),
    })

  }) # end moduleServer
}


