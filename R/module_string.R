
string_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(4, offset = 0, style='padding-left:15px; padding-right:5px; padding-top:0px; padding-bottom:0px',
             div(style="display: inline-block;vertical-align:top;", h5("Taxonomy code: ")),
             div(style="display: inline-block;vertical-align:top; width:50%", 
                 textInput(ns("tax"), label = NULL, value = "9606"))),
      column(3, offset = 0, style='padding-left:15px; padding-right:5px; padding-top:0px; padding-bottom:0px',
             div(style="display: inline-block;vertical-align:top;", 
                 uiOutput(ns("errorOrRun"))))
    ),
    uiOutput(ns("noresRet")),
    shinycssloaders::withSpinner(
      color="#0dc5c1",
      DT::dataTableOutput(ns("strtab"))
    ),
    wellPanel(
      checkboxInput(ns("showLabel"), label = "Show labels", value = FALSE),
      forceNetworkOutput(ns("network"))
    )
  )
}

string_module <- function(
  input, output, session, reactive_ids# , reactive_taxid = reactive(9606)
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
  
  nk <- eventReactive(
    input$run, {
      r <- stringNetwork(genes = reactive_ids(), taxid = input$tax)# reactive_taxid()) 
      if (is.data.frame(r)) {
        if (nrow(r) > 999) {
          r <- r[order(r$score, decreasing = TRUE), ]
          r <- r[1:999, ]
        }
      }
      r
    })
  gs <- eventReactive(input$run, stringGSA(genes = reactive_ids(), taxid = input$tax) )
  
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
    tab <- gs()[, c("category", "term", "description", "number_of_genes", "number_of_genes_in_background", "p_value", "fdr")]
    colnames(tab) <- c("category", "term", "desc", "genes", "background", "p value", "fdr")
    DT::datatable(data = tab, options = list(scrollX = TRUE), rownames = FALSE, selection = "single")
  })
  
  observe({
    req(input$strtab_rows_selected)
    highlight(input$strtab_rows_selected)
  })
}


# # # ####################
# library(shiny)
# library(shinycssloaders)
# source("Git/R/auxi_queryStringdb.R")
# dat <- readRDS("Dat/exampleEset.RDS")
# fd <- Biobase::fData(dat)
# ids <- fd$`General|All|Protein ID`[which(fd$`t-test|RE_BR|pval` < 0.05 & fd$`t-test|RE_BR|md` > 0.5)]
# 
# ui <- fluidPage(
#   string_ui("str")
# )
# 
# server <- function(input, output, session) {
#   callModule(string_module, id = "str", reactive_ids = reactive(ids))
# }
# 
# shinyApp(ui, server)

