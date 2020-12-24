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
#' @param additionalTabs additional tabs added to "Analyst" panel
#' @param esetLoader function to load the eset object, if an RDS file, should be "readRDS"
#' @param exprsGetter function to get the expression matrix from eset
#' @param pDataGetter function to get the phenotype data from eset
#' @param fDataGetter function to get the feature data from eset
#' @param defaultAxisGetter function to get the default axes to be visualized. It should be a function with two 
#'   arguments: x - the object loaded to the viewer; what - one of "sx", "sy", "fx" and "fy", representing the 
#'   sample space x-axis, sample space y-axis, feature space x-axis and feature space y-axis respectively. 
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
app_module <- function(
  input, output, session, dir, additionalTabs = NULL, 
  esetLoader = readRDS, exprsGetter = exprs, pDataGetter = pData, fDataGetter = fData, 
  defaultAxisGetter = function(x, what=c("sx", "sy", "fx", "fy")[1]) attr(x, what)
) {
  
  ns <- session$ns
  
  observe({
    req(dir())
    ll <- list.files(dir(), pattern = ".RDS$", ignore.case = TRUE)
    updateSelectizeInput(session = session, inputId = "selectFile", choices = ll, selected = "")
  })
  
  reactive_eset <- reactive({
    req(input$selectFile)
    flink <- file.path(dir(), input$selectFile)
    esetLoader(flink)
  })
  
  expr <- reactive({
    req(reactive_eset())
    exprsGetter(reactive_eset())
  })
  
  pdata <-reactive({
    req(reactive_eset())
    pDataGetter(reactive_eset())
  })
  
  fdata <-reactive({
    req(reactive_eset())
    fDataGetter(reactive_eset())
  })
  
  ########################
  
  ps <- reactive({
    req(eset <- reactive_eset())
    l <- list()
    if (!is.null(t0 <- defaultAxisGetter(eset, "sx"))){
      t0 <- strsplit(t0, "\\|")[[1]]
      l$x1_s <- t0[1]
      l$x2_s <- t0[2]
      l$x3_s <- t0[3]
    }
    if (!is.null(t0 <- defaultAxisGetter(eset, "sy"))) {
      t0 <- strsplit(t0, "\\|")[[1]]
      l$y1_s <- t0[1]
      l$y2_s <- t0[2]
      l$y3_s <- t0[3]
    }
    if (!is.null(t0 <- defaultAxisGetter(eset, "fx"))){ 
      t0 <- strsplit(t0, "\\|")[[1]]
      l$x1_f <- t0[1]
      l$x2_f <- t0[2]
      l$x3_f <- t0[3]
    }
    if (!is.null(t0 <- defaultAxisGetter(eset, "fy"))){ 
      t0 <- strsplit(t0, "\\|")[[1]]
      l$y1_f <- t0[1]
      l$y2_f <- t0[2]
      l$y3_f <- t0[3]
    }
    l
  })
  #####################
  
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
  
  v1 <- callModule(L1_data_space_module, id = "dataspace", expr = expr, pdata = pdata, fdata = fdata,
                   x1_s = ps()$x1_s, x2_s = ps()$x2_s, x3_s = ps()$x3_s,
                   y1_s = ps()$y1_s, y2_s = ps()$y2_s, y3_s = ps()$y3_s,
                   x1_f = ps()$x1_f, x2_f = ps()$x2_f, x3_f = ps()$x3_f,
                   y1_f = ps()$y1_f, y2_f = ps()$y2_f, y3_f = ps()$y3_f
                   )
  
  callModule(L1_result_space_module, id = "resultspace",
             reactive_expr = expr,
             reactive_phenoData = pdata,
             reactive_featureData = fdata,
             reactive_i = reactive(v1()$feature),
             reactive_highlight = reactive(v1()$sample),
             additionalTabs = additionalTabs,
             object = reactive_eset
  )
}



