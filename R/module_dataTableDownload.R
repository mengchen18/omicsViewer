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
#' @param reactive_row_ids stable row identifiers (same length/order as
#'   \code{reactive_table()}); used for snapshot row-selection persistence
#'   and as the store choices when \code{store_key} is given
#' @param store Optional child view of the canonical widget store
#'   (\code{\link{widget_store_child}}). No-op unless \code{store_key} is
#'   also given.
#' @param store_key Leaf id under \code{store} for the row-selection
#'   binding. Register ONLY where row selection drives a downstream view
#'   (ORA/fGSEA results, STRING enrichment, batch links); see the settled
#'   S4 decision in AGENT_ACCURACY_PLAN.md section 6.2. Requires
#'   \code{reactive_row_ids} yielding stable semantic ids.
#' @param store_label label for the registered binding
#' @param store_help help text for the registered binding
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
  decreasing = TRUE, reactive_row_ids = reactive(NULL),
  store = NULL, store_key = NULL, store_label = "Selected row",
  store_help = NULL) {

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

  # ------------------------------------------------------------------
  # Canonical widget-store bindings (control plane, plan section 6, S4).
  # Row selection registers ONLY at call sites where the click drives a
  # downstream view with no other agent interface (see the settled
  # decision in AGENT_ACCURACY_PLAN.md 6.2). The desired row id lives in a
  # reactiveVal feeding formatTab's selection$selected: pushes re-render
  # the table with the row preselected and DT's report of the selection
  # acknowledges through the sync observer. Tables without store_key keep
  # the legacy tab_status path untouched.
  # ------------------------------------------------------------------
  storeTab <- NULL
  .dtd_root_store <- NULL
  .dtd_full_key <- NULL
  if (!is.null(store) && !is.null(store_key)) {
    storeTab <- store
    .dtd_root_store <- if (is.null(store$parent)) store else store$parent
    .dtd_full_key <- paste0(store$prefix, ".", store_key)
    .dtd_ids <- function() {
      ids <- tryCatch(rowIds(),
                      shiny.silent.error = function(e) NULL,
                      error = function(e) NULL)
      if (is.null(ids)) character(0) else as.character(ids)
    }
    store_register(
      storeTab,
      widget_binding(
        store_key, "select", label = store_label,
        help = store_help %||% paste("Row currently selected in the table;",
                                     "drives the connected view"),
        choices_provider = function(v) .dtd_ids())
    )
  }
  # desired row id (store pushes and user clicks converge here)
  .dtd_sel_id <- reactiveVal(NULL)

  selectedRows <- reactive({
    if (!is.null(storeTab)) {
      # store-backed tables: the desired row id is canonical (survives
      # re-renders; cleared when the id is gone from the current table)
      sid <- .dtd_sel_id()
      if (length(sid) == 1L && nzchar(sid)) {
        ids <- rowIds()
        ts <- tabsort()
        if (!is.null(ids) && nrow(ts$tab) > 0) {
          p <- match(sid, ids)
          p <- p[!is.na(p)]
          w <- which(ts$index %in% p)
          return(if (length(w)) w else integer(0))
        }
      }
      return(NULL)
    }
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

  # ------------------------------------------------------------------
  # Store glue for the row-selection binding (acknowledgement-aware:
  # pushes re-render the table and ack through the same observer that
  # mirrors user clicks). Observer retention is mandatory.
  # ------------------------------------------------------------------
  .dtd_keep_list <- list()
  .dtd_keep <- function(obs) {
    .dtd_keep_list[[length(.dtd_keep_list) + 1L]] <<- obs
    invisible(obs)
  }
  .dtd_current_id <- reactive({
    rows <- input$table_rows_selected
    if (is.null(rows) || length(rows) == 0)
      return(NULL)
    ids <- rowIds()
    if (is.null(ids))
      return(NULL)
    ii <- tabsort()$index[rows]
    ii <- ii[!is.na(ii)]
    id <- unique(ids[ii])
    if (length(id) == 1L) id else NULL
  })
  if (!is.null(storeTab)) {
    # snapshot of the input at push time: while a push is in flight, the
    # input still reports the PRE-render selection; that stale report is
    # not a user edit (skipped). Any actual input CHANGE is processed -
    # the user always wins over an in-flight push.
    .dtd_push_input <- reactiveVal(NULL)
    .dtd_keep(observe({
      id <- .dtd_current_id()
      pend <- .dtd_root_store$pending[[.dtd_full_key]]
      stale_report <- !is.null(pend) && !identical(pend$value, id) &&
        identical(input$table_rows_selected, .dtd_push_input())
      if (stale_report)
        return(NULL)
      if (!identical(id, .dtd_sel_id()))
        .dtd_sel_id(id)
      store_sync_from_ui(storeTab, store_key, id)
    }))
    # one-shot seeding once the table has rendered a selection
    # (restore-first-wins: held keys are never overwritten)
    .dtd_seeded <- FALSE
    .dtd_keep(observe({
      if (.dtd_seeded) return(NULL)
      if (is.null(input$table_rows_selected)) return(NULL)
      .dtd_seeded <<- TRUE
      store_seed(storeTab,
                 stats::setNames(list(isolate(.dtd_current_id())), store_key))
    }))
    # store -> UI push: re-render with the pushed row preselected; DT
    # reports the selection and the sync observer above acknowledges
    .dtd_epoch <- store_epoch(storeTab)
    .dtd_keep(observe({
      .dtd_epoch()
      v <- store_read(storeTab, store_key)[[1]]
      if (!is.null(v) &&
          !is.null(.dtd_root_store$pending[[.dtd_full_key]])) {
        .dtd_push_input(input$table_rows_selected)
        .dtd_sel_id(v)
      }
    }))
  }
  
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

