
string_ui <- function(id) {
  ns <- NS(id)
  tagList(
    actionButton(ns("run"), "Run!"),
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
  input, output, session, reactive_ids, reactive_taxid = reactive(9606)
) {
  
  ns <- session$ns
  
  nk <- eventReactive(input$run, stringNetwork(genes = reactive_ids(), taxid = reactive_taxid()) )
  gs <- eventReactive(input$run, stringGSA(genes = reactive_ids(), taxid = reactive_taxid()) )
  
  highlight <- reactiveVal(1)
  
  output$network <- renderForceNetwork({
    req(nrow(nk()) > 0)
    stringD3Net(ntwk = nk(), gsa = gs(), i = highlight(), label = input$showLabel)
  })
  
  output$strtab <- DT::renderDataTable({
    print(nrow(gs()))
    req(nrow(gs()) > 0)
    tab <- gs()[, c("category", "term", "description", "number_of_genes", "number_of_genes_in_background", "p_value", "fdr")]
    colnames(tab) <- c("category", "term", "desc", "genes", "background", "p value", "fdr")
    print(dim(tab))
    DT::datatable(data = tab, options = list(scrollX = TRUE), rownames = FALSE, selection = "single")
  })
  
  observe({
    req(input$strtab_rows_selected)
    highlight(input$strtab_rows_selected)
  })
}

# 
# #################### 
library(shiny)
library(shinycssloaders)
source("Git/R/auxi_queryStringdb.R")
dat <- readRDS("Dat/exampleEset.RDS")
fd <- Biobase::fData(dat)
ids <- fd$`General|All|Protein ID`[which(fd$`t-test|RE_BR|pval` < 0.05 & fd$`t-test|RE_BR|md` > 0.5)]

ui <- fluidPage(
  string_ui("str")
)

server <- function(input, output, session) {
  callModule(string_module, id = "str", reactive_ids = reactive(ids))
}

shinyApp(ui, server)

