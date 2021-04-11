#' Utility string ui
#' @param id id
#' @importFrom shinybusy show_modal_spinner remove_modal_spinner
#' 
string_ui <- function(id) {
  ns <- NS(id)
  tagList(
      fluidRow(
        column(3, offset = 0, style='padding-left:15px; padding-right:5px; padding-top:0px; padding-bottom:0px',
               div(style="display: inline-block;vertical-align:top;", h5("Tax ID: ")),
               div(style="display: inline-block;vertical-align:top; width:50%", 
                   textInput(ns("tax"), label = NULL, value = "9606"))),
        column(6, offset = 0, style='padding-left:15px; padding-right:5px; padding-top:0px; padding-bottom:0px',
               div(style="display: inline-block;vertical-align:top;", 
                   uiOutput(ns("errorOrRun")))),
        column(3, offset = 0, style='padding-left:15px; padding-right:5px; padding-top:0px; padding-bottom:0px',
               align = "right", uiOutput(ns("showButton")) )
    ),
    uiOutput(ns("noresRet")),
    DT::dataTableOutput(ns("strtab")),    
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
  
  output$errorOrRun <- renderUI({
    tagList(
      verbatimTextOutput(ns("error.msg")),
      if (!overflow())
        actionButton(ns("run"), "Run!")
    )
  })
  
  overflow <- reactive({
    length(reactive_ids()) > 300
  })
  
  output$error.msg <- renderText({
    req(overflow())
    sprintf("%s features selected, allow max 300 input features!", length(reactive_ids()))
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
    show_modal_spinner(text = "Querying database ...")
    tab <- stringGSA(genes = reactive_ids(), taxid = input$tax) 
    colnames(tab) <- c("category", "term", "gene number", "background number", "TaxonId", "inputGenes", "preferredNames", "p value", "fdr", "description")
    remove_modal_spinner()
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
    req(!overflow())
    stringD3Net(ntwk = nk(), gsa = gs(), i = highlight(), label = input$showLabel)
  })
  
  output$strtab <- DT::renderDataTable({
    req(!nores())
    req(!overflow())
    req(nrow(gs()) > 0)
    tab <- gs()[, c("category", "term", "gene number", "background number", "p value", "fdr", "description")]
    DT::datatable(data = tab, options = list(scrollX = TRUE), rownames = FALSE, selection = "single")
  })
  
  output$showButton <- renderUI({
    # print( nrow(gs()) > 0 )
    # req(nrow(gs()) > 0)
    downloadButton(ns("downloadData"), ".tsv")
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("StringORA_", Sys.time(), ".tsv")
    },
    content = function(file) {
      tab <- gs()
      ic <- which(sapply(tab, is.list))
      if (length(ic) > 0) {
        for (ii in ic) {
          tab[, ii] <- sapply(tab[, ii], paste, collapse = ";")
        }
      }
      write.table(tab, file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    }
  )
  
  observe({
    req(input$strtab_rows_selected)
    highlight(input$strtab_rows_selected)
  })
}


