#' utility - dataTable shiny UI
#' @param id id
#' @importFrom shinyWidgets switchInput
#' 
dataTable_ui <- function(id) {
  ns <- NS(id)
  tagList(
    shinyWidgets::switchInput( inputId = ns("multisel"), label = "Multiple selection" , labelWidth = "120px"),
    uiOutput(ns("selector")),    
    DT::dataTableOutput(ns("table"))
  )
}

#' utility - dataTable shiny module
#' @param input input
#' @param output output
#' @param session session
#' @param reactive_data the data to be shown, a tabular objet
#' @param selector whether a selector should be added to the output
#' @param columns columns to show
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
dataTable_module <- function(input, output, session, reactive_data, selector = TRUE, columns = NULL) {
  
  ns <- session$ns

  rdd <- reactive({
    if (is.matrix(reactive_data())) {
      x <- as.data.frame(reactive_data()) 
    } else if (is.data.frame(reactive_data()))
      x <- reactive_data() else
        stop('reactive_data shold be either a matrix or data.frame')
    x
  })

  cols <- reactive({    
    
    cn <- intersect(columns, colnames(rdd()))    
    if (length(cn) == 0)
      cn <- grep("^General\\|", colnames(rdd()), ignore.case = TRUE, value = TRUE)
    if (length(cn) == 0)
      cn <- colnames(rdd())

    opt <- NULL
    if (selector) {
      optx <- setdiff(colnames(rdd()), cn)
      if (length(optx) > 0)
        opt <- str_split_fixed(optx, pattern = "\\|", n = 3)
    }      
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

  output$selector <- renderUI({
    req(nrow(cols()$opt) > 0)
    triselector_ui(ns("select"))
    })
  
  formatTab <- function(tab, sel) {    
    dt <- DT::datatable( 
      tab,
      selection =  c("single", "multiple")[as.integer(sel)+1],
      rownames = FALSE,
      filter = "top",
      class="table-bordered compact",
      options = list(scrollX = TRUE, pageLength = 25, dom = 'tip')
    )
    DT::formatStyle(dt, columns = 1:ncol(tab), fontSize = '90%')
  }

  output$table <- DT::renderDataTable({
    tab <- rdd()[, scn()]    
    i <- which(sapply(tab, function(x) is.numeric(x) && !is.integer(x)))
    if (any(i))
      tab[i] <- lapply(tab[i], round, digits = 4)
    formatTab(tab, sel = input$multisel)
  })
  
  reactive({
    rownames(reactive_data())[input$table_rows_selected]
  })
}



