#' @importFrom shinybusy show_modal_spinner remove_modal_spinner
geneshot_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # select variable
    triselector_ui(ns("geneNameCol")),
    fluidRow(
      column(
        width = 10, 
        textInput(ns("term"), label = NULL, value = "", width = "100%", placeholder = "Search any term here, multiple separated by (;) ...")
      ),
      column(
        width = 2, align = "right", 
        actionButton(ns("submit"), label = "Submit!")
      )
    ),
    plotly::plotlyOutput(ns("plt")),
    tags$br(),
    dataTableDownload_ui(ns("autorif"))
  )
}

geneshot_module <- function(
  input, output, session, pdata, fdata, expr, feature_selected, sample_selected, object,
  reactive_status = reactive(NULL)
) {
  ns <- session$ns
  # select stats from feature data
  triset <- reactive({
    fd <- fdata()
    cn <- colnames(fd)[!sapply(fd, is.numeric) & !grepl("^GS\\|", colnames(fd))]
    str_split_fixed(cn, "\\|", n = 3)
  })

  xax <- reactiveVal()
  v1 <- callModule(
    triselector_module, id = "geneNameCol", reactive_x = triset, label = "Gene name",
    reactive_selector1 = reactive(xax()$v1), 
    reactive_selector2 = reactive(xax()$v2), 
    reactive_selector3 = reactive(xax()$v3)
  )
  

  
  rif <- reactiveVal()
  observeEvent(input$submit, {
    show_modal_spinner(text = "Querying database ...")
    res <- getAutoRIF(trimws(strsplit(input$term, ";")[[1]]), filter = TRUE)    
    if (!is.null(res) && nrow(res) > 0) {
      # req(nrow(res) > 0)
      res$selected <- ""
      res <- res[order(res$rank, decreasing = TRUE), ]
      dft <- res[1:min(nrow(res), 20), ]
      outliers <- list(
        x = dft$n,
        y = dft$perc,
        text = dft$gene,
        xref = "x",
        yref = "y",
        showarrow = TRUE,
        ax = 10,
        ay = -20
      )
    } 
    remove_modal_spinner()
    req(!is.null(res) && nrow(res) > 0)    
    rif(
      list(
      tab = res,
      outliers = outliers
    ))    
  })

  outliersLabs <- reactive( rif()$outliers )
  
  gtab <- reactive({
    req(df <- rif()$tab)
    cid <- paste(unlist(v1()), collapse = "|")
    if ( cid %in% colnames(fdata()) )
      df$selected[ df$gene %in% fdata()[feature_selected(), cid] ] <- "+"    
    df
    })  

  rifRow <- callModule(
    dataTableDownload_module,    
    id = "autorif", reactive_table = reactive({req(gtab()); gtab()}), prefix = "autoRIF", pageLength = 10)
  
  ##
  output$plt <- plotly::renderPlotly({
    req(df <- gtab())
    fig <- plotly_scatter(
      x = df$n, y = df$perc, xlab = "# of publication", ylab = "Publications with Search Term(s) / Total Publications", 
      color = df$selected, size = 10, tooltips=df$gene, shape = "select"
    )
    fig <- plotly::layout(fig$fig, annotations = outliersLabs()) #rif()$outliers)
    plotly::config(
      fig,
      toImageButtonOptions = list(
        format = "svg",
        filename = "ExpressionSetViewerPlot",
        width = 700,
        height = 700
        )
      )
  })

  ### status save and restore

  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return()
    xax(NULL)
    xax(list(v1 = s$xax[[1]], v2 = s$xax[[2]], v3 = s$xax[[3]]))
    # gtab(s$table)
    # outliersLabs(s$outliers)
    rif( s$rif )
    updateTextInput(session, inputId = "term", value = s$term)
    })

  ##
  rv <- reactiveValues()
  observe( rv$xax <- v1() )
  observe( rv$rif <- rif() )
  # observe( rv$outliers <- rif()$outliers )
  observe( rv$term <- input$term )

  reactive({
    reactiveValuesToList(rv)
    })

}