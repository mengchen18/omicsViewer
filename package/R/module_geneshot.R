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
      ),
      plotly::plotlyOutput(ns("plt")),
      tags$br(),
      dataTableDownload_ui(ns("autorif"))
    )
  )
}

geneshot_module <- function(
  input, output, session, pdata, fdata, expr, feature_selected, sample_selected, object
) {
  ns <- session$ns
  # select stats from feature data
  triset <- reactive({
    fd <- fdata()
    cn <- colnames(fd)[!sapply(fd, is.numeric) & !grepl("^GS\\|", colnames(fd))]
    str_split_fixed(cn, "\\|", n = 3)
  })
  v1 <- callModule(
    triselector_module, id = "geneNameCol", reactive_x = triset, label = "Gene name"
  )
  
  rif <- eventReactive(input$submit, {
    res <- getAutoRIF(trimws(strsplit(input$term, ";")[[1]]), filter = TRUE)    
    req(nrow(res) > 0)
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
    list(
      tab = res,
      outliers = outliers
    )
  })
  
  gtab <- reactive({
    df <- rif()$tab
    cid <- paste(unlist(v1()), collapse = "|")
    if ( cid %in% colnames(fdata()) )
      df$selected[ df$gene %in% fdata()[feature_selected(), cid] ] <- "+"
    df
    })

  rifRow <- callModule(
    dataTableDownload_module,    
    id = "autorif", reactive_table = gtab, prefix = "autoRIF", pageLength = 10)
  
  ##
  output$plt <- plotly::renderPlotly({
    df <- gtab()
    fig <- plotly_scatter(
      x = df$n, y = df$perc, xlab = "# of publication", ylab = "Publications with Search Term(s) / Total Publications", 
      color = df$selected, size = 10, tooltips=df$gene, shape = "select"
    )
    plotly::layout(fig$fig, annotations = rif()$outliers)
  })
}