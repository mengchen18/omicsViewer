#' @description utility - dataTable shiny UI
#' @param id id
#' @importFrom shinyWidgets switchInput
#' 
dataTable_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      # column(2, dropdown(
      #   margin = "25px", status = "default", icon = icon("cog"), width = "700px",
      #   tooltip = tooltipOptions(title = "Add more columns to the table!"),
      #    )),
      column(3,
        actionButton(ns("clear"), "Show all") %>%
          tagAppendAttributes(`data-testid` = paste0(id, "-clear-filter-button"))
      ),
      column(6, align = "center",
        shinyWidgets::switchInput( inputId = ns("multisel"), label = "Multiple_selection" , labelWidth = "125px") %>%
          tagAppendAttributes(`data-testid` = paste0(id, "-multiselect-toggle"))
      ),
      column(3, dataTableDownload_ui(ns("downloadTable"), showTable = FALSE), align="right")
    ),
    uiOutput(ns("selector")),
    DT::dataTableOutput(ns("table"))
  )
}

#' @description utility - dataTable shiny module
#' @param id module id
#' @param reactive_data the data to be shown, a tabular objet
#' @param selector whether a selector should be added to the output
#' @param columns columns to show
#' @param tab_status table initial status, reactive object
#' @param tab_rows rows to be shown
#' @param store Optional child view of the canonical widget store
#'   (\code{\link{widget_store_child}}) for this table instance. When given,
#'   the multi-row-selection switch and the set of shown columns are
#'   registered as agent-controllable bindings (control plane, plan
#'   section 6, S4). DataTable-internal browser state (search, ordering,
#'   pagination) stays on the \code{tab_status} snapshot path.
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
#' #   dataTable_module("dttest",  reactive_data = reactive(pdata))
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
#' #   dataTable_module("dttest",  reactive_data = reactive(fdata))
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
#' #   dataTable_module("dttest",  reactive_data = reactive(expr), selector = FALSE)
#' # }
#' #
#' # shinyApp(ui, server)
#'
dataTable_module <- function(
  id, reactive_data, selector = TRUE, columns = NULL,
  tab_status = reactive(NULL), tab_rows = reactive(NULL),
  store = NULL
  ) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns
  # 
  selectedRowOrCol <- reactiveVal(TRUE)
  
  notNullAndPosLength <- function(x) !is.null(x) && length(x) > 0

  observeEvent( reactive_data(), selectedRowOrCol(TRUE) )
  observeEvent( input$clear, selectedRowOrCol(TRUE) )
  observe( selectedRowOrCol( tab_rows() ) )

  rdd <- reactive({
    req(reactive_data())
    if (is.matrix(reactive_data())) {
      x <- as.data.frame(reactive_data()) 
    } else if (is.data.frame(reactive_data()))
      x <- reactive_data() else
        stop('reactive_data shold be either a matrix or data.frame')
    x[selectedRowOrCol(), , drop = FALSE]
  })
  
  dataTableDownload_module(
    "downloadTable", reactive_table = rdd, prefix = "viewerTable_"
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

  addcols <- triselector_module("select", reactive_x = reactive({
    req(cols()$opt)
    req(nrow(cols()$opt) > 0)
    cols()$opt
  }), label = "Add column")
  
  scn <- reactiveVal(NULL)
  observe(
    scn(cols()$shown)
  )
  observeEvent(addcols(), {
    req(!addcols()$variable %in% c("", "--select--"))
    oc <- scn()
    nc <- unique(c(oc, paste(addcols(), collapse = "|")))
    nc <- intersect(nc, colnames(rdd()))
    req(nc)
    scn(nc)
  })

  # ------------------------------------------------------------------
  # Canonical widget-store bindings (control plane, plan section 6, S4).
  # The table's user-editable surface: the multi-row-selection switch and
  # the set of shown columns (built through the Add-column selector;
  # snapshots restore it wholesale, so the store models it as a
  # multi_select). DataTable-internal browser state (search, ordering,
  # pagination) stays on the tab_status path.
  # ------------------------------------------------------------------
  if (!is.null(store)) {
    # Defensive choices source (NO req(): a shiny.validation condition
    # inside a choices_provider would abort the whole store_apply
    # transaction; see the heatmap S4 notes)
    .dt_col_choices <- function() {
      d <- reactive_data()
      if (is.matrix(d) || is.data.frame(d)) colnames(d) else character(0)
    }
    store_register(
      store,
      widget_binding("multi_selection", "boolean",
        label = "Multiple selection",
        help = "Allow selecting more than one row in the table"),
      widget_binding("columns", "multi_select", label = "Shown columns",
        help = paste("Table columns currently displayed; the Add-column",
                     "selector appends entries; at least one column must",
                     "remain selected"),
        min = 1L,
        choices_provider = function(v) .dt_col_choices())
    )
    .dt_store_root <- if (is.null(store$parent)) store else store$parent
    .dt_keys <- c(multi_selection = "multisel", columns = "scn")
    # Create the epoch reactive ONCE and keep strong references to every
    # store-glue observer (observer-GC rule, see heatmap/meta_scatter).
    .dt_epoch <- store_epoch(store)
    .dt_observers <- list()
    .dt_keep <- function(obs) {
      .dt_observers[[length(.dt_observers) + 1L]] <<- obs
      invisible(obs)
    }

    # UI -> store: switch edits and column-set changes (triselector adds,
    # status restores, store pushes) all mirror through the same sync;
    # it is acknowledgement-aware, so pushed values ack their pending
    # entries instead of counting as user overrides.
    .dt_keep(observeEvent(input$multisel, {
      store_sync_from_ui(store, "multi_selection", input$multisel)
    }, ignoreInit = TRUE))
    .dt_keep(observe({
      if (is.null(scn())) return(NULL)
      store_sync_from_ui(store, "columns", scn())
    }))

    # Seed unset keys with the widget defaults once the inputs exist, so
    # discovery tools report real current values from the start; restores
    # and agent applies that land first win (seeding skips held keys).
    .dt_seeded <- FALSE
    .dt_keep(observe({
      if (.dt_seeded) return(NULL)
      if (is.null(input$multisel) || is.null(scn())) return(NULL)
      .dt_seeded <<- TRUE
      held <- store_read(store, names(.dt_keys))
      patch <- list()
      if (is.null(held[[paste0(store$prefix, ".multi_selection")]]))
        patch$multi_selection <- input$multisel
      if (is.null(held[[paste0(store$prefix, ".columns")]]))
        patch$columns <- scn()
      if (length(patch))
        tryCatch(store_apply(store, patch, origin = "system", strict = FALSE),
                 error = function(e) NULL)
    }))

    # Store -> UI push for external writes only (pending entries mark
    # them). Columns push defensively: values that no longer exist in the
    # data are dropped, and an empty intersection is skipped rather than
    # blanking the table.
    .dt_keep(observe({
      .dt_epoch()
      vals <- store_read(store, names(.dt_keys))
      pending <- .dt_store_root$pending
      invisible(lapply(names(.dt_keys), function(key) {
        full <- paste0(store$prefix, ".", key)
        value <- vals[[full]]
        if (is.null(value) || is.null(pending[[full]]))
          return(NULL)
        if (identical(key, "multi_selection")) {
          updateSwitchInput(session, "multisel", value = value)
        } else {
          okcols <- intersect(value, .dt_col_choices())
          if (length(okcols))
            scn(okcols)
        }
      }))
    }))
  }

  output$selector <- renderUI({
    req(cols()$opt)
    req(nrow(cols()$opt) > 0)
    triselector_ui(ns("select"))
    })
  
  tabStatus <- reactive({
    st <- tab_status()
    if (!is.list(st)) return(NULL)
    st
  })

  selectedRows <- reactive({
    st <- tabStatus()
    if (is.null(st))
      return(NULL)
    if (is.null(st$selected_rows) || length(st$selected_rows) == 0)
      return(st$rows_selected)
    rn <- rownames(rdd())
    if (is.null(rn))
      return(st$rows_selected)
    match(as.character(st$selected_rows), rn)
  })

  formatTab <- function(tab, sel) {    
    ci <- unname(which(vapply(tab, inherits, c('factor', "character"), FUN.VALUE = logical(1))))
    if (length(ci) > 0)
    tab[ci] <- lapply(tab[ci], function(x) {
      x[is.na(x)] <- ""
      x
    })
    dt <- DT::datatable(
      tab,
      selection = list(mode = c("single", "multiple")[as.integer(sel)+1], selected = selectedRows(), target = "row"),
      rownames = FALSE,
      filter = "top",
      class="table-bordered compact nowrap",
      caption = htmltools::tags$caption(
        style = "caption-side: top; text-align: left; font-weight: normal; padding: 5px 0;",
        sprintf("Data table with %d rows and %d columns. Use column filters to search and sort data.",
                nrow(tab), ncol(tab))
      ),
      options = list(
        scrollX = TRUE, dom = 'tip',
        columnDefs = list(list(
          targets = ci-1,
          render = DT::JS(
            "function(data, type, row, meta) {",
            "return type === 'display' && data.length > 50 ?",
            "'<span title=\"' + data + '\">' + data.substr(0, 50) + '...</span>' : data;",
            "}")
        )),
        # Server snapshots are authoritative; do not restore browser-local state.
        searchCols = getSearchCols(tabStatus()), order = getOrderCols(tabStatus()),
        displayStart = tabStatus()$start,
        pageLength = restore_table_page_length(tabStatus()$length, pageLength = DEFAULT_TABLE_PAGE_LENGTH_LARGE)
        )
    )
    DT::formatStyle(dt, columns = seq_len(ncol(tab)), fontSize = '90%')
  }

  observeEvent(tabStatus(), {
    if (is.null(tabStatus()))
      return(NULL)
    if (!is.null(store)) {
      # single transactional path (per-key resilient): the column set is
      # intersected with the live colnames by the multi_select validator
      patch <- list()
      if (!is.null(tabStatus()$showColumns))
        patch$columns <- tabStatus()$showColumns
      if (!is.null(tabStatus()$multiSelection))
        patch$multi_selection <- tabStatus()$multiSelection
      if (length(patch))
        tryCatch(store_apply(store, patch, origin = "restore", strict = FALSE),
                 error = function(e) NULL)
    } else {
      if (!is.null(i <- tabStatus()$showColumns))
        scn(i)
      updateSwitchInput(session, "multisel", value = tabStatus()$multiSelection)
    }
  })
  
  output$table <- DT::renderDataTable({
    req(scn())
    tab <- rdd()[, scn(), drop = FALSE]
    i <- which(vapply(tab, function(x) is.numeric(x) && !is.integer(x), logical(1)))
    if (any(i))
      tab[i] <- lapply(tab[i], round, digits = 4)
    formatTab(tab, sel = input$multisel)
  })

  # outputOptions(output, "table", suspendWhenHidden = FALSE)
  tabproxy <- dataTableProxy(ns("table"))
  
  eventReactive( list(input$table_rows_selected, input$table_state), {
    r <- character(0)
    if (!is.null(tab_rows()))
      r <- tab_rows()
    if (notNullAndPosLength(input$table_rows_selected))
      r <- rownames(rdd())[input$table_rows_selected]
    sta <- data_table_widget_state(input$table_state)
    if (is.null(sta))
      sta <- list()
    sta$showColumns <- scn()
    sta$multiSelection <- input$multisel
    sta$rows_selected <- input$table_rows_selected
    rn <- rownames(rdd())
    if (!is.null(rn) && notNullAndPosLength(input$table_rows_selected))
      sta$selected_rows <- unique(rn[input$table_rows_selected])
    attr(r, "status") <- sta
    r
    })

  }) # end moduleServer
}



