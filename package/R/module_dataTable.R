#' utility - dataTable shiny UI
#' @param id id
#' @importFrom shinyWidgets switchInput
#' 
dataTable_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(3, actionButton(ns("clear"), "Show all")),
      column(6, align = "center", 
        shinyWidgets::switchInput( inputId = ns("multisel"), label = "Multiple selection" , labelWidth = "120px")
        ),
      column(3, dataTableDownload_ui(ns("downloadTable"), showTable = FALSE), align="right")
    ),
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
#' @param reactiveSelectorMeta object returned by pdata or fdata
#' @param reactiveSelectorHeatmap object returned by heatmap
#' @param subset subset row or column for heatmap
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
dataTable_module <- function(
  input, output, session, reactive_data, selector = TRUE, columns = NULL, 
  reactiveSelectorMeta = reactive(NULL), reactiveSelectorHeatmap = reactive(NULL), subset = c("row", "col", "none")[1]) {
  
  ns <- session$ns
  
  # 
  selectedRowOrCol <- reactiveVal(TRUE)
  
  notNullAndPosLength <- function(x) !is.null(x) && length(x) > 0
  
  if (subset %in% c('row', 'col')) {
    observeEvent(reactiveSelectorHeatmap(), {
      # req(reactiveSelectorHeatmap())
      if (subset == "row") {
        if (notNullAndPosLength(reactiveSelectorHeatmap()$brushed$row)) {
          selectedRowOrCol(reactiveSelectorHeatmap()$brushed$row)
        } else if (notNullAndPosLength(reactiveSelectorHeatmap()$clicked)) {
          selectedRowOrCol(reactiveSelectorHeatmap()$clicked["row"])
        } else
          selectedRowOrCol(TRUE)
      } else if (subset == "col") {
        if (notNullAndPosLength(reactiveSelectorHeatmap()$brushed$col)) {
          selectedRowOrCol(reactiveSelectorHeatmap()$brushed$col)
        } else if (notNullAndPosLength(reactiveSelectorHeatmap()$clicked)) {
          selectedRowOrCol(reactiveSelectorHeatmap()$clicked["col"])
        } else
          selectedRowOrCol(TRUE)
      } else 
        stop("Unknown subset, should be either 'row' or 'col'!")
    })
    
    observeEvent(reactiveSelectorMeta(), {
      # req(reactiveSelectorMeta())
      if (notNullAndPosLength(reactiveSelectorMeta()$selected)) {
        selectedRowOrCol( reactiveSelectorMeta()$selected )
      } else if ( notNullAndPosLength(reactiveSelectorMeta()$clicked) ) {
        selectedRowOrCol( reactiveSelectorMeta()$clicked )
      } else
        selectedRowOrCol(TRUE)
    })
  }

  observeEvent( reactive_data(), selectedRowOrCol(TRUE) )
  observeEvent( input$clear, selectedRowOrCol(TRUE) )

  rdd <- reactive({
    if (is.matrix(reactive_data())) {
      x <- as.data.frame(reactive_data()) 
    } else if (is.data.frame(reactive_data()))
      x <- reactive_data() else
        stop('reactive_data shold be either a matrix or data.frame')
    x <- x[selectedRowOrCol(), ]
  })
  
  callModule(
    dataTableDownload_module, id = "downloadTable", reactive_table = rdd, prefix = "viewerTable_"
  )

  cols <- eventReactive(reactive_data(), {    
    
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
      class="table-bordered compact nowrap",
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
  
  # selVal <- reactiveVal(character(0))
  # observeEvent(input$clear, {
  #   selVal(character(0))
  #   })
  # observeEvent(list(input$table_rows_selected, input$clear), {
  #   if (notNullAndPosLength(input$table_rows_selected))
  #     selVal( rownames(rdd())[input$table_rows_selected] ) else
  #      selVal( character(0) )
  #   })
  # reactive( selVal() )
  eventReactive( input$table_rows_selected, {
    req (notNullAndPosLength(input$table_rows_selected))
    rownames(rdd())[input$table_rows_selected] 
    } )
}



