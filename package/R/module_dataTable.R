#' @description utility - dataTable shiny UI
#' @param id id
#' @importFrom shinyWidgets switchInput
#' 
dataTable_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      # column(2, dropdown( 
      #   margin = "25px", status = "default", icon = icon("gear"), width = "700px",
      #   tooltip = tooltipOptions(title = "Add more columns to the table!"),
      #    )),
      column(3, actionButton(ns("clear"), "Show all")),
      column(6, align = "center", 
             shinyWidgets::switchInput( inputId = ns("multisel"), label = "Multiple_selection" , labelWidth = "125px")
      ),
      column(3, dataTableDownload_ui(ns("downloadTable"), showTable = FALSE), align="right")
    ),
    uiOutput(ns("selector")),
    DT::dataTableOutput(ns("table"))
  )
}

#' @description utility - dataTable shiny module
#' @param input input
#' @param output output
#' @param session session
#' @param reactive_data the data to be shown, a tabular objet
#' @param selector whether a selector should be added to the output
#' @param columns columns to show
#' @param tab_status table initial status, reactive object
#' @param tab_rows rows to be shown
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
  tab_status = reactive(NULL), tab_rows = reactive(NULL) 
  ) {
  
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
    x[selectedRowOrCol(), ]
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
    req(cols()$opt)
    req(nrow(cols()$opt) > 0)
    cols()$opt
  }), label = "Add column")
  
  scn <- reactiveVal(NULL)
  observe(
    scn(cols()$shown)
  )
  observeEvent(addcols(), {
    req(!addcols()$variable %in% c("", "Select a variable!"))
    oc <- scn()
    nc <- unique(c(oc, paste(addcols(), collapse = "|")))
    nc <- intersect(nc, colnames(rdd()))
    req(nc)
    scn(nc)
  })

  output$selector <- renderUI({
    req(cols()$opt)
    req(nrow(cols()$opt) > 0)
    triselector_ui(ns("select"))
    })
  
  formatTab <- function(tab, sel) {    
    ci <- unname(which(sapply(tab, inherits, c('factor', "character"))))
    if (length(ci) > 0)
    tab[ci] <- lapply(tab[ci], function(x) {
      x[is.na(x)] <- ""
      x
    })
    dt <- DT::datatable( 
      tab,
      selection =  list(mode = c("single", "multiple")[as.integer(sel)+1], selected = tab_status()$rows_selected, target = "row"),
      rownames = FALSE,
      filter = "top",
      class="table-bordered compact nowrap",
      options = list(
        scrollX = TRUE, pageLength = 25, dom = 'tip',
        columnDefs = list(list(
          targets = ci-1,
          render = DT::JS(
            "function(data, type, row, meta) {",
            "return type === 'display' && data.length > 50 ?",
            "'<span title=\"' + data + '\">' + data.substr(0, 50) + '...</span>' : data;",
            "}")
        )),
        stateSave = TRUE,  stateDuration = -1,
        searchCols = getSearchCols(tab_status()), order = getOrderCols(tab_status()),
        displayStart = tab_status()$start
        )
    )
    DT::formatStyle(dt, columns = 1:ncol(tab), fontSize = '90%')
  }

  observeEvent(tab_status(), {
    if (!is.null( i <- tab_status()$showColumns ))
      scn (i)
    updateSwitchInput(session, "multisel", value = tab_status()$multiSelection)    
    }) 
  
  output$table <- DT::renderDataTable({
    req(scn())
    tab <- rdd()[, scn()]
    i <- which(sapply(tab, function(x) is.numeric(x) && !is.integer(x)))
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
    sta <- input$table_state
    sta$showColumns <- scn()
    sta$multiSelection <- input$multisel    
    sta$rows_selected <- input$table_rows_selected    
    attr(r, "status") <- sta
    r
    })
}



