app_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    uiOutput(ns("summary")),
    br(),
    absolutePanel(
      top = 8, right = 150, style = "z-index: 9999;",
      selectInput(inputId = ns("selectFile"), label = NULL, choices = NULL, selectize = TRUE)
    ),
    absolutePanel(
      top = 5, right = 20, style = "z-index: 9999;",
      downloadButton(outputId = ns("download"), label = "Download", class = NULL)
    ),
    fluidRow(
      column(6, L1_data_space_ui(ns('dataspace'))),
      column(6, L1_result_space_ui(ns("resultspace")))
    )
  )
}

app_module <- function(input, output, session, dir) {
  
  ns <- session$ns
  
  observe({
    req(dir())
    ll <- list.files(dir(), pattern = ".RDS$", ignore.case = TRUE)
    updateSelectInput(session = session, inputId = "selectFile", choices = ll, selected = "")
  })
  
  reactive_eset <- reactive({
    # req(dir())
    # print(input$selectFile)
    req(input$selectFile)
    flink <- file.path(dir(), input$selectFile)
    readRDS(flink)
  })
  
  expr <- reactive({
    req(reactive_eset())
    Biobase::exprs(reactive_eset())
  })
  
  pdata <-reactive({
    req(reactive_eset())
    Biobase::pData(reactive_eset())
  })
  
  fdata <-reactive({
    req(reactive_eset())
    Biobase::fData(reactive_eset())
  })
  
  output$summary <- renderUI({
    if (is.null(input$selectFile) || input$selectFile == "")
      return(HTML('<h1 style="display:inline;">ExpressionSetViewer</h1>'))
    txt <- sprintf(
      '<h1 style="display:inline;">ExpressionSetViewer</h1> <h3 style="display:inline;">  --   %s features and %s samples:</h3>', 
      nrow(expr()), ncol(expr())
    )
    HTML(txt)
  })
  
  v1 <- callModule(L1_data_space_module, id = "dataspace", expr = expr, pdata = pdata, fdata = fdata)
  
  callModule(L1_result_space_module, id = "resultspace",
             reactive_expr = expr,
             reactive_phenoData = pdata,
             reactive_featureData = fdata,
             reactive_i = reactive(v1()$feature),
             reactive_highlight = reactive(v1()$sample)
  )
}



