#' utility - dataTable shiny UI
#' @param id id
#' @param selector where to show the selector
#' 
dataTable_ui <- function(id, selector = TRUE) {
  ns <- NS(id)
  tagList(
    if (selector)
      triselector_ui(ns("select")),
    DT::dataTableOutput(ns("table"))
  )
}

#' utility - dataTable shiny module
#' @param input input
#' @param output output
#' @param session session
#' @param reactive_data the data to be shown, a tabular objet
#' @param selector whether a selector should be added to the output
#' @importFrom stringr str_split_fixed
#' @examples 
#' # library(shiny)
#' # source("Git/R/module_triselector.R")
#' # 
#' # dat <- readRDS("Dat/exampleEset.RDS")
#' # expr <- exprs(dat)
#' # pdata <- pData(dat)
#' # fdata <- fData(dat)
#' # 
#' # ui <- fluidPage(
#' #   dataTable_ui("dttest")
#' # )
#' # 
#' # server <- function(input, output, session) {
#' #   callModule(dataTable_module, id = "dttest",  reactive_data = reactive(pdata))
#' # }
#' # 
#' # shinyApp(ui, server)
#' # 
#' # 
#' # ####
#' # ui <- fluidPage(
#' #   dataTable_ui("dttest")
#' # )
#' # 
#' # server <- function(input, output, session) {
#' #   callModule(dataTable_module, id = "dttest",  reactive_data = reactive(fdata))
#' # }
#' # 
#' # shinyApp(ui, server)
#' # 
#' # ###
#' # ui <- fluidPage(
#' #   dataTable_ui("dttest", selector = FALSE)
#' # )
#' # 
#' # server <- function(input, output, session) {
#' #   callModule(dataTable_module, id = "dttest",  reactive_data = reactive(expr), selector = FALSE)
#' # }
#' # 
#' # shinyApp(ui, server)
#' 
dataTable_module <- function(input, output, session, reactive_data, selector = TRUE) {
  
  rdd <- reactive({
    if (is.matrix(reactive_data())) {
      x <- as.data.frame(reactive_data()) 
    } else if (is.data.frame(reactive_data()))
      x <- reactive_data() else
        stop('reactive_data shold be either a matrix or data.frame')
    x
  })
  
  cols <- reactive({
    cn <- colnames(rdd())
    if (selector) {
      i <- grepl("^General\\|", cn)
      req(any(i))
      opt <- str_split_fixed(cn[!i], pattern = "\\|", n = 3)
      cn <- cn[i]
    } else 
      opt <- NULL
    list(shown = cn, opt = opt)
  })
  
  addcols <- callModule(triselector_module, id = "select", reactive_x = reactive({
    req(nrow(cols()$opt) > 0)
    cols()$opt
  }), label = "Add")
  
  scn <- reactiveVal(NULL)
  observe(
    scn(cols()$shown)
  )
  observeEvent(addcols(), {
    req(!addcols()$variable %in% c("", "Select a variable!"))
    oc <- scn()
    nc <- unique(c(oc, paste(addcols(), collapse = "|")))
    scn(nc)
  })
  
  formatTab <- function(tab) {
    numc <- which(sapply(tab, function(x) is.numeric(x) && !is.integer(x)))
    dt <- DT::datatable( 
      tab,
      selection =  "multiple",
      rownames = FALSE,
      filter = "top",
      class="table-bordered compact",
      options = list(scrollX = TRUE, scrollY = "800px", paging = FALSE, dom = 't')
    )
    dt <- DT::formatRound(dt, columns = numc, digits = 3)
    DT::formatStyle(dt, columns = 1:ncol(tab), fontSize = '90%')
  }

  output$table <- DT::renderDataTable({
    tab <- rdd()[, scn()]
    formatTab(tab)
    # i <- sapply(tab, is.numeric)
    # if (any(i))
    #   tab[i] <- lapply(tab[i], signif, digits = 3)
    # DT::datatable(tab, options = list(
    #   scrollX = TRUE, pageLength = 20
    #   ), 
    #   rownames = FALSE, selection = "multiple"
    #   )
  })
  
  reactive({
    rownames(reactive_data())[input$table_rows_selected]
  })
}



