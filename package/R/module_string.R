#' Utility string ui
#' @param id id
#' @importFrom shinybusy show_modal_spinner remove_modal_spinner
#' 
string_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        4, offset = 0, style='padding-left:15px; padding-right:2px; padding-top:0px; padding-bottom:0px',
        textInputAddon(inputId = ns("tax"), label = NULL, value = "9606", addon = "Taxomony Code")),
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

#' Utility string module
#' @param input input
#' @param output output
#' @param session session
#' @param reactive_ids ids passed to stringDB
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
#' #   callModule(string_module, id = "str", reactive_ids = reactive(ids))
#' # }
#' # 
#' # shinyApp(ui, server)
#' 
string_module <- function(
  input, output, session, reactive_ids
) {
  
  ns <- session$ns
  
  overflow <- reactive({
    length(reactive_ids()) > 300
  })
  
  output$error.msg <- renderText({
    sprintf("%s features selected [MAX 300 FEATURES ALLOWED!]", length(reactive_ids()))
  })
  
  nk <- eventReactive( input$run, {
    r <- stringNetwork(genes = reactive_ids(), taxid = input$tax)# reactive_taxid()) 
    if (is.data.frame(r)) {
      if (nrow(r) > 999) {
        r <- r[order(r$score, decreasing = TRUE), ]
        r <- r[1:999, ]
      }
    }
    r
  })
  
  gs <- eventReactive(input$run, {
    req(!overflow())
    show_modal_spinner(text = "Querying database ...")
    tab <- stringGSA(genes = reactive_ids(), taxid = input$tax) 
    remove_modal_spinner()
    if (inherits(tab, "character")) 
      return(tab)
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
  
  highlight <- reactiveVal(1)
  
  output$network <- renderForceNetwork({
    req(!nores())
    req(nrow(nk()) > 0)
    stringD3Net(ntwk = nk(), gsa = gs(), i = highlight(), label = input$showLabel)
  })
  
  tt <- callModule(
    dataTableDownload_module, id = "strtab", reactive_table = eventReactive(gs(), {
      req(!nores())
      req(!overflow())
      req(nrow(gs()) > 0)
      gs()[, c("category", "term", "gene number", "background number", "p value", "fdr", "description")]
    }), prefix = "FeatureTable_"
  )
  
  observe({
    req(tt())
    highlight(tt())
  })
}


