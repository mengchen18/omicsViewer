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

  # One-shot restore slot (todo 2.8): the saved DT state (filters, order,
  # page, pre-selected rows) applies on the FIRST render after a status
  # change only; later re-renders (new ranking, selection changes, page
  # flips) must not re-apply it - a restore used to keep pre-selecting the
  # snapshot's rows at their new positions on every re-render.
  .dtd_restore_slot <- reactiveVal(NULL)
  observeEvent(tabStatus(), {
    st <- tabStatus()
    if (is.list(st))
      .dtd_restore_slot(st)
  })
  .dtd_consume_restore <- function() {
    st <- .dtd_restore_slot()  # reactive read: a restore re-renders the table
    if (!is.null(st)) {
      # R-M1: defer while the table has no data yet - consuming here would
      # validate st$selected_rows against an empty/pre-restore table and
      # silently drop the saved row selection; the slot stays armed and is
      # consumed by the first render that actually has data
      ts <- tryCatch(tabsort(), error = function(e) NULL)
      if (is.null(ts) || is.null(ts$tab) || nrow(ts$tab) == 0)
        return(NULL)
      isolate(.dtd_restore_slot(NULL))
    }
    st
  }

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

  selectedRows <- function(st = NULL) {
    if (!is.null(storeTab)) {
      # store-backed tables: the desired row id is canonical (survives
      # re-renders; cleared when the id is gone from the current table).
      # isolate (R-M2): reading the id inside renderDataTable used to make
      # every selection change re-render the WHOLE table (page reset,
      # filters/order lost); live selection updates go through the DT
      # proxy observer below instead
      sid <- isolate(.dtd_sel_id())
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
    if (is.null(st))
      return(NULL)
    if (is.null(st$selected_rows) || length(st$selected_rows) == 0)
      return(st$rows_selected)
    ids <- rowIds()
    if (is.null(ids))
      return(st$rows_selected)
    # tabsort()$index maps displayed rows back to rows in reactive_table().
    match(as.character(st$selected_rows), ids[tabsort()$index])
  }
  
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
  
  formatTab <- function(tab, sel = 0, pageLength = pageLength, st = NULL) {
    dt <- DT::datatable(
      tab,
      selection = list(
        mode = c("single", "multiple")[as.integer(sel) + 1],
        selected = selectedRows(st),
        target = "row"
      ),
      rownames = FALSE,
      filter = "top",
      class="table-bordered compact nowrap",
      # Server snapshots are authoritative; DataTable browser-local state is
      # deliberately not enabled because it is machine-specific.
      options = list(scrollX = TRUE, dom = 'tip',
        searchCols = getSearchCols(st), order = getOrderCols(st),
        displayStart = st$start,
        pageLength = restore_table_page_length(st$length, pageLength = pageLength)
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

  # R-M2 fallback epoch: a push that lands BEFORE the table has ever
  # rendered cannot go through the DT proxy (no table object yet, e.g.
  # the first paint after a tab switch); bumping the epoch re-renders
  # with the selection in the widget payload instead. Gated on
  # table_state (NULL only before the first draw), not on the selection
  # input (which is legitimately NULL after a deselect).
  .dtd_render_epoch <- reactiveVal(0L)
  output$table <- DT::renderDataTable({
    .dtd_render_epoch()
    formatTab(tabsort()$tab, pageLength = pageLength, st = .dtd_consume_restore())
  })

  # ------------------------------------------------------------------
  # Store glue for the row-selection binding (acknowledgement-aware:
  # pushes select the row through the DT proxy without re-rendering;
  # user clicks and pushes acknowledge through the same sync observer,
  # which is guarded against stale positional reports across table-data
  # changes (R-H5). Observer retention is mandatory.
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
  # table-data epoch (R-H5): increments on every tabsort recomputation
  # (a reactive, so every reader sees the SAME value within a flush - no
  # observer-ordering race). A positional report derived from an input
  # that was set under an OLDER table must not be re-mapped against the
  # new table: the derived id would name a DIFFERENT row (the ORA results
  # table used to silently switch the selected pathway this way)
  .dtd_epoch_n <- 0L
  .dtd_tab_epoch <- reactive({
    tabsort()
    .dtd_epoch_n <<- .dtd_epoch_n + 1L
    .dtd_epoch_n
  })
  .dtd_input_epoch <- NULL  # epoch under which the last selection report was accepted
  if (!is.null(storeTab)) {
    # snapshot of the input at push time: while a push is in flight, the
    # input still reports the PRE-render selection; that stale report is
    # not a user edit (skipped). Any actual input CHANGE is processed -
    # the user always wins over an in-flight push.
    .dtd_push_input <- reactiveVal(NULL)
    .dtd_keep(observe({
      id <- .dtd_current_id()
      epoch <- .dtd_tab_epoch()
      if (!identical(.dtd_input_epoch, epoch)) {
        # R-H5: the selection input predates the current table data -
        # skip this (stale positional) report and consume the transition;
        # the re-render/proxy re-establishes the id-based selection and
        # the browser re-reports, which lands with matched epochs
        .dtd_input_epoch <<- epoch
        return(NULL)
      }
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
    # store -> UI push: set the desired id; the proxy observer below
    # selects the row WITHOUT re-rendering the table, and DT's report of
    # the selection acknowledges through the sync observer above
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
    # live selection updates via the DT proxy (R-M2): pushes and user
    # clicks converge on .dtd_sel_id; applying them through the proxy
    # preserves page/sort/filter where a re-render reset them. A table
    # that has never rendered falls back to the re-render path
    .dtd_keep(observe({
      .dtd_sel_id()
      if (is.null(input$table_state)) {
        .dtd_render_epoch(isolate(.dtd_render_epoch()) + 1L)
        return(NULL)
      }
      rows <- isolate(selectedRows())
      # dataTableProxy applies session$ns() itself - pass the module-local
      # id (ns("table") double-prefixes and the message lands on a
      # non-existent table)
      DT::selectRows(
        DT::dataTableProxy("table", session = session,
                           deferUntilFlush = FALSE),
        if (notNullAndPositiveLength(rows)) rows else NULL)
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

