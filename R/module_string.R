#' @description Utility string ui
#' @param id id
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

#' @description Utility string module
#' @param id module id
#' @param reactive_ids ids passed to stringDB
#' @param reactive_status reactive status for restoring
#' @param active reactive logical indicating if active
#' @importFrom networkD3 renderForceNetwork forceNetworkOutput
#' @examples
#' #' # # # ####################
#' # library(shiny)
#' # source("Git/R/auxi_queryStringdb.R")
#' # dat <- readRDS("Dat/exampleEset.RDS")
#' # fd <- Biobase::fData(dat)
# ids <- fd$`General|All|Protein ID`[which(fd$`t-test|RE_BR|pval` < 0.05 & fd$`t-test|RE_BR|md` > #' 0.5)]
#' #
#' # ui <- fluidPage(
#' #   string_ui("str")
#' # )
#' #
#' # server <- function(input, output, session) {
#' #   string_module("str", reactive_ids = reactive(ids))
#' # }
#' #
#' # shinyApp(ui, server)
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
    r <- stringNetwork(genes = reactive_ids(), taxid = input$tax)# reactive_taxid())
    if (is.data.frame(r)) {
      if (nrow(r) > STRING_MAX_NETWORK_EDGES) {
        r <- r[order(r$score, decreasing = TRUE), ]
        r <- r[seq_len(STRING_MAX_NETWORK_EDGES), ]
      }
    }
    nk(r)
  })
  
  gs <- reactiveVal()
  gs <- eventReactive( input$run, {
  # observeEvent(input$run, {
    req(!overflow())
    show_modal_spinner(text = "Querying database ...")
    tab <- stringGSA(genes = reactive_ids(), taxid = input$tax) 
    remove_modal_spinner()
    if (inherits(tab, "character")) {
      # gs(tab)
      return(tab)
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


