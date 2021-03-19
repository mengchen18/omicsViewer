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
      selectizeInput(inputId = ns("selectFile"), label = NULL, choices = NULL, width = "550px", options = list(placeholder = "Select a dataset here") )
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
#' @param filePattern file pattern to be displayed.
#' @param additionalTabs additional tabs added to "Analyst" panel
#' @param esetLoader function to load the eset object, if an RDS file, should be "readRDS"
#' @param exprsGetter function to get the expression matrix from eset
#' @param pDataGetter function to get the phenotype data from eset
#' @param fDataGetter function to get the feature data from eset
#' @param defaultAxisGetter function to get the default axes to be visualized. It should be a function with two 
#'   arguments: x - the object loaded to the viewer; what - one of "sx", "sy", "fx" and "fy", representing the 
#'   sample space x-axis, sample space y-axis, feature space x-axis and feature space y-axis respectively. 
#' @param appName name of the application
#' @param appVersion version of the application
#' @importFrom Biobase exprs pData fData
#' @importFrom utils packageVersion
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
  input, output, session, dir, filePattern = ".RDS$", additionalTabs = NULL, 
  esetLoader = readRDS, exprsGetter = exprs, pDataGetter = pData, fDataGetter = fData, 
  defaultAxisGetter = function(x, what=c("sx", "sy", "fx", "fy", "dendrogram")[1]) attr(x, what),
  appName = "ExpressionSetViewer", appVersion = packageVersion("ExpressionSetViewer")
) {
  
  ns <- session$ns
  
  observe({
    req(dir())
    ll <- list.files(dir(), pattern = filePattern, ignore.case = TRUE)
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

  validEset <- function(expr, pd, fd) {
    i1 <- identical(rownames(expr), rownames(fd))
    i2 <- identical(colnames(expr), rownames(pd))
    if (!(i1 && i2))
      return(
        list(
          FALSE, "The rownames/colnames of exprs not matched to row names of feature data/phenotype data!"
          )
        )
    TRUE
  }

  observe({
    x <- validEset(expr = expr(), pd = pdata(), fd = fdata())
    if (!x[[1]]) {
      showModal(modalDialog(
        title = "Problem in data!",
        x[[2]]
      ))
    }
  })  
  ########################  
  d_s_x <- reactive( {
    req(eset <- reactive_eset())
    defaultAxisGetter(eset, "sx") 
    })
  d_s_y <- reactive( {
    req(eset <- reactive_eset())
    defaultAxisGetter(eset, "sy") 
    })
  d_f_x <- reactive( {
    req(eset <- reactive_eset())
    defaultAxisGetter(eset, "fx") 
    })
  d_f_y <- reactive( {
    req(eset <- reactive_eset())
    defaultAxisGetter(eset, "fy") 
    })
  reactive_rdg <- reactive({
    req(eset <- reactive_eset())
    defaultAxisGetter(eset, "dendrogram")
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
      return(HTML(sprintf('<h1 style="display:inline;">%s</h1>', appName)))
    txt <- sprintf(
      '<h1 style="display:inline;">%s</h1> <h3 style="display:inline;"><sup>%s</sup>  --   %s features and %s samples:</h3>', 
      appName, paste0("v", appVersion), nrow(expr()), ncol(expr())
    )
    HTML(txt)
  })
  
  v1 <- callModule(L1_data_space_module, id = "dataspace", expr = expr, pdata = pdata, fdata = fdata,
                    reactive_x_s = d_s_x, reactive_y_s = d_s_y, reactive_x_f = d_f_x, reactive_y_f = d_f_y,
                    rowDendrogram = reactive_rdg
                   )


  ri <- reactiveVal()
  observeEvent( v1(), ri(v1()$feature) )
  observeEvent( expr(), ri(NULL) )

  rh <- reactiveVal()
  observeEvent( v1(), rh(v1()$sample) )
  observeEvent( expr(), rh(NULL) )
  
  callModule(L1_result_space_module, id = "resultspace",
             reactive_expr = expr,
             reactive_phenoData = pdata,
             reactive_featureData = fdata,
             reactive_i = reactive(ri()), # reactive(v1()$feature),
             reactive_highlight = reactive(rh()), # reactive(v1()$sample),
             additionalTabs = additionalTabs,
             object = reactive_eset
  )
}



