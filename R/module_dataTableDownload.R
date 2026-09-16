#' @description utility - dataTable for download shiny UI
#' @param id id
#' @param showTable logical, if the table should be shown
#' 
dataTableDownload_ui <- function(id, showTable = TRUE) {
  ns <- NS(id)
  if (showTable) {
    r <- tagList(
      uiOutput(ns("showButton")),
      DT::dataTableOutput(ns("table"))      
      )
    
  } else
    r <- uiOutput(ns("showButton"))
  r
}

#' @description utility - dataTable for download shiny module
#' @description A subset of columns can be shown, specified by reactive_cols.
#'   The entire table will be downloaded.
#' @param id module id
#' @param reactive_table table to show
#' @param reactive_cols columns to be shown
#' @param prefix file name prefix
#' @param pageLength how many row per page
#' @param sortBy sort by column (name)
#' @param decreasing logical; sort by decreasing or not
#' @param tab_status table initial status
#' @importFrom utils write.table
#' @examples
#' # source("R/module_triselector.R")
#' # library(shiny)
#' # library(stringr)
#' #
#' # dat <- readRDS("../Dat/exampleEset.RDS")
#' # pdata <- pData(dat)
#' #
#' # ui <- fluidPage(
#' #   dataTableDownload_ui("dtd")
#' # )
#' # server <- function(input, output, session) {
#' #   dataTableDownload_module("dtd", reactive_table = reactive(pdata), reactive_cols = reactive(1:6), prefix = "testdownload")
#' # }
#' # shinyApp(ui, server)
#'
dataTableDownload_module <- function(id, reactive_table, tab_status = reactive(NULL),
  reactive_cols=reactive(NULL), prefix = "", pageLength = 10, sortBy = NULL,
  decreasing = TRUE, reactive_row_ids = reactive(NULL)) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns
  notNullAndPositiveLength <- function(x) !is.null(x) && length(x) > 0

  rtab <- reactive({
    req(tt <- reactive_table())
    if (is.matrix(tt))
      tt <- as.data.frame(tt, stringsAsFactors = FALSE)
    tt
    })

  # Some callers still pass a plain list while newer callers pass a reactive.
  tabStatus <- reactive({
    if (is.function(tab_status)) tab_status() else tab_status
  })

  rowIds <- reactive({
    if (is.null(reactive_row_ids()))
      return(NULL)
    as.character(reactive_row_ids())
  })

  selectedRows <- reactive({
    st <- tabStatus()
    if (is.null(st))
      return(NULL)
    if (is.null(st$selected_rows) || length(st$selected_rows) == 0)
      return(st$rows_selected)
    ids <- rowIds()
    if (is.null(ids))
      return(st$rows_selected)
    # tabsort()$index maps displayed rows back to rows in reactive_table().
    match(as.character(st$selected_rows), ids[tabsort()$index])
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0(prefix, Sys.time(), ".tsv")
    },
    content = function(file) {
      tab <- rtab()
      ic <- which(vapply(tab, is.list, logical(1)))
      if (length(ic) > 0) {
        for (ii in ic) {
          tab[, ii] <- vapply(tab[, ii], paste, collapse = ";", FUN.VALUE = character(1))
        }
      }
      write.table(tab, file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    }
  )
  
  output$showButton <- renderUI({
    req(rtab())
    downloadLink(ns("downloadData"), "Save table")
  })
  
  formatTab <- function(tab, sel = 0, pageLength = pageLength) {
    dt <- DT::datatable(
      tab,
      selection = list(
        mode = c("single", "multiple")[as.integer(sel) + 1],
        selected = selectedRows(),
        target = "row"
      ),
      rownames = FALSE,
      filter = "top",
      class="table-bordered compact nowrap",
      # Server snapshots are authoritative; DataTable browser-local state is
      # deliberately not enabled because it is machine-specific.
      options = list(scrollX = TRUE, dom = 'tip',
        searchCols = getSearchCols(tabStatus()), order = getOrderCols(tabStatus()),
        displayStart = tabStatus()$start,
        pageLength = restore_table_page_length(tabStatus()$length, pageLength = pageLength)
        )
    )
    DT::formatStyle(dt, columns = seq_len(ncol(tab)), fontSize = '90%')
  }

  tabsort <- reactive({
    req(tab <- rtab())
    index <- seq_len( nrow(tab) )
    if (!is.null(sortBy)) {
      if (sortBy %in% colnames(tab)) {
        o <- order(tab[, sortBy], decreasing = decreasing)
        tab <- tab[o, ]
        index <- index[o]
      }
    }
    if (!is.null( reactive_cols() ))
      tab <- tab[, reactive_cols()]
    ic <- which(vapply(tab, function(x) is.numeric(x) & !is.integer(x), logical(1)))
    if ( length(ic) > 0 )
      tab[, ic] <- lapply(tab[, ic, drop = FALSE], signif, digits = 3)

    list(tab = tab, index = index)
  })    

  output$table <- DT::renderDataTable(    
    formatTab(tabsort()$tab, pageLength = pageLength)
  )
  
  reactive({
    ii <- tabsort()$index[input$table_rows_selected]
    sta <- data_table_widget_state(input$table_state)
    if (is.null(sta))
      sta <- list()
    sta$rows_selected <- input$table_rows_selected
    ids <- rowIds()
    if (!is.null(ids) && notNullAndPositiveLength(ii))
      sta$selected_rows <- unique(ids[ii])
    attr(ii, "status") <- sta
    ii
    })

  }) # end moduleServer
}

