#' application level 0 UI
#' @param id id
#' @param showDropList logical; whether to show the dropdown list to select RDS file, if the 
#'   ESVObj is given, this should be set to "FALSE"
#' @param activeTab one of "Feature", "Feature table", "Sample", "Sample table", "Heatmap"
#' @export
#' @examples
#' if (interactive()) {
#'   dir <- system.file("extdata", package = "ExpressionSetViewer")
#'   server <- function(input, output, session) {
#'     callModule(app_module, id = "app", dir = reactive(dir))
#'   }
#'   ui <- fluidPage(
#'     app_ui("app")
#'   )
#'   shinyApp(ui = ui, server = server)
#' }
#' @return a list of UI components

app_ui <- function(id, showDropList = TRUE, activeTab = "Feature") {
  ns <- NS(id)

  comp <- list(
    useShinyjs(),
    style = "background:white;",
    absolutePanel(
      top = 5, right = 20, style = "z-index: 9999;", width = 115, 
      downloadButton(outputId = ns("download"), label = "xlsx", class = NULL),
      actionButton(ns("snapshot"), label = NULL, icon = icon("camera-retro"))
    ),
    column(6, L1_data_space_ui(ns('dataspace'), activeTab = activeTab)),
    column(6, L1_result_space_ui(ns("resultspace"))))

  if (showDropList) {
    l2 <- list(
      uiOutput(ns("summary")),
      br(),
      absolutePanel(
        top = 8, right = 140, style = "z-index: 9999;",
        selectizeInput( inputId = ns("selectFile"), label = NULL, choices = NULL, 
          width = "500px", options = list(placeholder = "Select a dataset here") )

      ))
    comp <- c(l2, comp)
    }
  do.call(fluidRow, comp)
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
#' @param ESVObj the ESV object
#'   given, the drop down list should be disable in the "ui" component.
#' @importFrom Biobase exprs pData fData
#' @importFrom utils packageVersion
#' @importFrom DT renderDT DTOutput
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
#' @export
#' @examples
#' if (interactive()) {
#'   dir <- system.file("extdata", package = "ExpressionSetViewer")
#'   server <- function(input, output, session) {
#'     callModule(app_module, id = "app", dir = reactive(dir))
#'   }
#'   ui <- fluidPage(
#'     app_ui("app")
#'   )
#'   shinyApp(ui = ui, server = server)
#' }
#' @return do not return any values

app_module <- function(
  input, output, session, dir, filePattern = ".RDS$", additionalTabs = NULL, ESVObj = reactive(NULL),
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
    # try to get global object first
    if (!is.null(ESVObj())) 
      return( tallGS(ESVObj()) )

    # otherwise load from disk
    req(input$selectFile)
    flink <- file.path(dir(), input$selectFile)
    sss <- file.size(flink)
    if (sss > 1e7)
      show_modal_spinner(text = "Loading data ...")
    v <- esetLoader(flink)
    if (sss > 1e7)
      remove_modal_spinner()
    v
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
        id <- rownames(tab)
        if (is.null(id))
          id <- paste0("ID", 1:nrow(tab))
        data.frame(ID = id, tab)
      }
      withProgress(message = 'Writing table', value = 0, {
        wb <- createWorkbook(creator = "BayBioMS")
        addWorksheet(wb, sheetName = "Phenotype info")
        addWorksheet(wb, sheetName = "Feature info")
        addWorksheet(wb, sheetName = "Expression")
        addWorksheet(wb, sheetName = "Geneset annot")
        incProgress(1/5, detail = "expression matrix")
        writeData(wb, sheet = "Expression", td(expr()))
        incProgress(1/5, detail = "feature table")
        writeData(wb, sheet = "Feature info", td(fdata()))
        incProgress(1/5, detail = "phenotype table")
        writeData(wb, sheet = "Phenotype info", td(pdata()))
        incProgress(1/5, detail = "writing geneset annotation")
        writeData(wb, sheet = "Geneset annot", attr(fdata(), "GS"))
        incProgress(1/5, detail = "Saving table")
        saveWorkbook(wb, file = file, overwrite = TRUE)
      })
    }
  )
  
  output$summary <- renderUI({
    req(is.null(ESVObj()))
    if (is.null(input$selectFile) || input$selectFile == "" || !is.numeric(nrow(expr())))
      return(HTML(sprintf('<h1 style="display:inline;">%s</h1>', appName)))
    txt <- sprintf(
      '<h1 style="display:inline;">%s</h1> <h3 style="display:inline;"><sup>%s</sup>  --   %s features and %s samples:</h3>', 
      appName, paste0("v", appVersion), nrow(expr()), ncol(expr())
    )
    HTML(txt)
  })
  
  v1 <- callModule(
    L1_data_space_module, id = "dataspace", expr = expr, pdata = pdata, fdata = fdata,
    reactive_x_s = d_s_x, reactive_y_s = d_s_y, reactive_x_f = d_f_x, reactive_y_f = d_f_y,
    rowDendrogram = reactive_rdg, status = esv_status
  )
  
  sameValues <- function(a, b) {
    if (is.null(a) || is.null(b)) 
      return(FALSE)    
    all(sort(a) == sort(b)) 
    }
  ri <- reactiveVal()
  observeEvent( v1(), {     
    req(!is.na(f <- v1()$feature))    
    oldvalues <- ri()    
    if (sameValues(oldvalues, f)) return(NULL)    
    ri( f ) 
    })
  observeEvent( expr(), ri(NULL) )
  
  rh <- reactiveVal()
  observeEvent( v1(), {
    req(!is.na(s <- v1()$sample))    
    oldvalues <- rh()    
    if (sameValues(oldvalues, s)) return(NULL)    
    rh( s )
    })
  observeEvent( expr(), rh(NULL) )

  # observe({
  #   ri ( esv_status()$active_feature )
  #   rh ( esv_status()$active_sample )
  #   })
  
  
  v2 <- callModule(L1_result_space_module, id = "resultspace",
                   reactive_expr = expr,
                   reactive_phenoData = pdata,
                   reactive_featureData = fdata,
                   reactive_i = ri, # reactive(v1()$feature),
                   reactive_highlight = rh, # reactive(v1()$sample),
                   additionalTabs = additionalTabs,
                   object = reactive_eset,
                   status = esv_status
  )

  # =======================================================
  # =======================================================
  # ================= snapshot function ===================
  # =======================================================
  # =======================================================

  savedSS <- reactiveVal()
  observe({
    fl <- paste0("ESVSnapshot_", input$selectFile, "_")
    ff <- list.files(dir(), pattern = fl)
    if (length(ff) == 0)
      return(NULL)
    r <- sub(fl, "", ff)
    r <- sub(".ESS$", "", r)
    df <- data.frame("name" = r, link = ff, stringsAsFactors = FALSE, check.names = FALSE)    
    savedSS(df)
    })
  
  shinyInput <- function(FUN, len, id, ...) {
    inputs <- c()
    for (i in len) {
      inputs <- c(inputs, as.character(FUN(paste0(id, i), ...)))
    }
    inputs
  }

  output$tab_saveSS <- renderDT({
    req( nrow(dt <- savedSS()) > 0)    
    dt$delete <- shinyInput(
      actionButton, dt$name, 'deletess_', label = "Delete", onclick = sprintf('Shiny.setInputValue(\"%s\",  this.id)', ns("deletess_button")) 
      )
    DT::datatable(
      dt[, c(1, 3), drop = FALSE], rownames = FALSE, colnames = c(NULL, NULL, NULL), 
      selection = list(mode = "single", target = "cell", selectable = 1), escape = FALSE,
      options = list(
        dom = "t", autoWidth = FALSE, style="compact-hover", scrollY = "450px", 
        paging = FALSE, columns = list(list(width = "85%"), list(width = "15%"))
        )
      )
    })

  selectedSS <- reactiveVal()
  observe({    
    ss <- input$tab_saveSS_cells_selected
    if (length(ss) == 0 || ss[2] > 0)
      return(NULL)
    selectedSS( ss[1])
    })

  observeEvent(list(v1(), v2()), {
    selectedSS( NULL )
    })

  deleteSS <- reactiveVal()
  observeEvent(input$deletess_button, {
    selectedRow <- sub("deletess_", "", input$deletess_button)
    deleteSS(selectedRow)
  })

  observeEvent(input$snapshot, {
    showModal(
      modalDialog(
        title = NULL,
        fluidRow(
          column(9, textInput(ns("snapshot_name"), label = "Save new snapshot", placeholder = "snapshot name", width = "100%")),
          column(3, style = "padding-top:25px", actionButton(ns("snapshot_save"), label = "Save")),
          ),        
        hr(),
        strong("Load saved snapshots:"),
        DTOutput(ns("tab_saveSS")),
        footer = NULL,
        easyClose = TRUE
        )
      )
    })

  # observe({
  #   print(ri())
  #   print(rh())
  #   })

  observeEvent(input$snapshot_save, {
    req(input$selectFile)
    df <- savedSS()
    if (input$snapshot_name %in% df$name) {      
      showModal(modalDialog(
        title = "FAILED!",  
        "Snapshot with this name already exists, please give a different name."
        ))
      return(NULL)
    }
    obj <- c(attr(v1(), "status"), v2(), active_feature = list(ri()), active_sample = list(rh()))
    flink <- file.path(dir(), paste0("ESVSnapshot_", input$selectFile, "_", input$snapshot_name, ".ESS"))
    saveRDS(obj, flink)
    df <- rbind(df, data.frame(name = input$snapshot_name, link = basename(flink)), stringsAsFactors = FALSE)
    dt <- df[order(df$name), ]
    savedSS(dt)
    removeModal()
    })

  esv_status <- reactiveVal()
  observeEvent(selectedSS(), {
    req(nrow(df <- savedSS()) > 0)    
    if (length(i <- selectedSS()) == 0)
      return(NULL)
    removeModal()
    esv_status( readRDS(file.path(dir(), df[i, 2])) )    
    })  

  observeEvent(deleteSS(), {
    req(nrow( df <- savedSS() ) > 0 )
    req( i <- match( deleteSS(), df$name ))
    df <- savedSS()
    unlink(file.path(dir(), df[i, 2]))
    df <- df[-i, , drop = FALSE]
    savedSS(df)    
    })

  # observe(
  #   print(esv_status())
  #   )

  # 
  # eset - active tab
  # eset - feature - x, y, color, size, ...
  # eset - sample - x, y, color, size, ...
  # returnData <- reactiveVal()
  # observe( returnData( v1()$data ) )
  # # observe( returnData( v2()$data ) )
  # observe(print(returnData()))
}



