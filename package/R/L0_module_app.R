#' application level 0 UI
#' @param id id
#' 
app_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    uiOutput(ns("summary")),
    br(),
    absolutePanel(
      top = 8, right = 120, style = "z-index: 9999;",
      selectizeInput(inputId = ns("selectFile"), label = NULL, choices = NULL, options = list(placeholder = "Select a dataset here") )
    ),
    absolutePanel(
      top = 5, right = 20, style = "z-index: 9999;",
      downloadButton(outputId = ns("download"), label = "xlsx", class = NULL)
    ),
    fluidRow(
      column(6, L1_data_space_ui(ns('dataspace'))),
      column(6, L1_result_space_ui(ns("resultspace")))
    )
  )
}

#' Application level 0 module
#' @param input input
#' @param output output
#' @param session session
#' @param dir reactive; directory containing the .RDS file of ExpressionSet
#' @importFrom Biobase exprs pData fData
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics abline axis barplot image mtext par plot text
#' @importFrom stats 
#'  as.dendrogram
#'  as.dist
#'  chisq.test
#'  cor.test
#'  fisher.test
#'  hclust
#'  lm
#'  na.omit
#'  p.adjust
#'  predict
#'  quantile
#'  t.test
#'  uniroot
#'  wilcox.test
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#' 
#' 
app_module <- function(input, output, session, dir) {
  
  ns <- session$ns
  
  observe({
    req(dir())
    ll <- list.files(dir(), pattern = ".RDS$", ignore.case = TRUE)
    updateSelectizeInput(session = session, inputId = "selectFile", choices = ll, selected = "")
  })
  
  reactive_eset <- reactive({
    req(input$selectFile)
    flink <- file.path(dir(), input$selectFile)
    readRDS(flink)
  })
  
  expr <- reactive({
    req(reactive_eset())
    exprs(reactive_eset())
  })
  
  pdata <-reactive({
    req(reactive_eset())
    pData(reactive_eset())
  })
  
  fdata <-reactive({
    req(reactive_eset())
    fData(reactive_eset())
  })
  
  output$download <- downloadHandler(
    filename = function() {
      paste0("ExpressenSet", Sys.time(), ".xlsx")
    },
    content = function(file) {
      td <- function(tab) {
        ic <- which(sapply(tab, is.list))
        if (length(ic) > 0) {
          for (ii in ic) {
            tab[, ii] <- sapply(tab[, ii], paste, collapse = ";")
          }
        }
        tab
      }
      withProgress(message = 'Writing table', value = 0, {
        wb <- createWorkbook(creator = "BayBioMS")
        addWorksheet(wb, sheetName = "Phenotype info")
        addWorksheet(wb, sheetName = "Feature info")
        addWorksheet(wb, sheetName = "Expression")
        incProgress(1/4, detail = "expression matrix")
        id <- paste0("ID", 1:nrow(expr()))
        writeData(wb, sheet = "Expression", data.frame(ID = id, expr()))
        incProgress(1/4, detail = "feature table")
        writeData(wb, sheet = "Feature info", td(cbind(ID = id, fdata())))
        incProgress(1/4, detail = "phenotype table")
        writeData(wb, sheet = "Phenotype info", td(pdata()))
        incProgress(1/4, detail = "Saving table")
        saveWorkbook(wb, file = file, overwrite = TRUE)
      })
      
    }
  )
  
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



