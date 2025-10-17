#' omicsViewer Application UI (Level 0)
#'
#' @description
#' Generates the user interface for the main omicsViewer application. This function creates
#' a responsive layout with data exploration panels, snapshot functionality, and data export
#' capabilities. Primarily intended for developers extending the application.
#'
#' @param id Character. Namespace ID for the Shiny module. Must match the ID used in
#'   \code{\link{app_module}}.
#' @param showDropList Logical. Whether to display the file selection dropdown menu.
#'   Set to FALSE when providing data directly via \code{ESVObj} parameter in
#'   \code{\link{app_module}}. Default: TRUE.
#' @param activeTab Character. Initial tab to display when data is loaded. Options:
#'   \itemize{
#'     \item "Feature" - Feature space scatter plot
#'     \item "Feature table" - Feature metadata table
#'     \item "Sample" - Sample space scatter plot
#'     \item "Sample table" - Sample metadata table
#'     \item "Cor" - Correlation heatmap
#'     \item "Heatmap" - Expression heatmap
#'     \item "Dynamic heatmap" - Interactive heatmap with selection
#'     \item "Expression" - Expression matrix table
#'     \item "GSList" - Gene set membership table
#'   }
#'   Default: "Feature".
#'
#' @return
#' A \code{fluidRow} containing the complete UI structure, including:
#' \itemize{
#'   \item File selection dropdown (if \code{showDropList = TRUE})
#'   \item Data summary display
#'   \item Export and snapshot buttons
#'   \item Two-column layout with data space (left) and analysis space (right)
#' }
#'
#' @export
#' @importFrom shinyjs useShinyjs hidden
#'
#' @seealso
#' \code{\link{app_module}} for the corresponding server logic.
#' \code{\link{omicsViewer}} for the high-level application launcher.
#'
#' @examples
#' if (interactive()) {
#'   dir <- system.file("extdata", package = "omicsViewer")
#'   ui <- fluidPage(
#'     app_ui("app", showDropList = TRUE, activeTab = "Feature")
#'   )
#'   server <- function(input, output, session) {
#'     app_module("app", .dir = reactive(dir))
#'   }
#'   shinyApp(ui = ui, server = server)
#' }
#' @keywords internal

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
    shinyjs::hidden(
      div(id = ns("contents"),
        column(6, L1_data_space_ui(ns('dataspace'), activeTab = activeTab)),
        column(6, L1_result_space_ui(ns("resultspace")))
        )
      )    
    )

  if (showDropList) {
    l2 <- list(
      shinycssloaders::withSpinner(
        uiOutput(ns("summary")), hide.ui = FALSE, type = 8, color = "green"
        ),      
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
#' @description Function should only be used for the developers
#' @param id module id
#' @param .dir reactive; directory containing the .RDS file of \code{ExpressionSet} or \code{SummarizedExperiment}
#' @param filePattern file pattern to be displayed.
#' @param additionalTabs additional tabs added to "Analyst" panel
#' @param esetLoader function to load the eset object, if an RDS file, should be "readRDS"
#' @param exprsGetter function to get the expression matrix from eset
#' @param imputeGetter function to get the imputed expression matrix from eset, only used when exporting imputed data to excel
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
#' @importFrom DT renderDT DTOutput dataTableProxy
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics abline axis barplot image mtext par plot text
#' @importFrom stats
#'  as.dendrogram
#'  as.dist
#'  as.hclust
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
#'   dir <- system.file("extdata", package = "omicsViewer")
#'   server <- function(input, output, session) {
#'     app_module("app", .dir = reactive(dir))
#'   }
#'   ui <- fluidPage(
#'     app_ui("app")
#'   )
#'   shinyApp(ui = ui, server = server)
#' }
#' @return do not return any values

app_module <- function(
  id, .dir, filePattern = ".(RDS|db|sqlite|sqlite3)$", additionalTabs = NULL, ESVObj = reactive(NULL),
  esetLoader = readESVObj, exprsGetter = getExprs, pDataGetter = getPData, fDataGetter = getFData,
  imputeGetter = getExprsImpute, defaultAxisGetter = getAx,
  appName = "omicsViewer", appVersion = packageVersion("omicsViewer")
) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns

  ll <- reactive({
    req(.dir())
    list.files(.dir(), pattern = filePattern, ignore.case = TRUE)
    })

  observe({
    req(ll())
    updateSelectizeInput(session = session, inputId = "selectFile", choices = ll(), selected = "")
  })
  
  reactive_eset <- reactive({
    # try to get global object first
    if (!is.null(ESVObj())) {
      updateSelectizeInput(session, "selectFile", choices = c("ESVObj.RDS", ll()), selected = "ESVObj.RDS")
      return( tallGS(ESVObj()) )
    }
    # otherwise load from disk
    req(input$selectFile)

    # Comprehensive input validation for file loading
    # 1. Validate file path - prevent directory traversal
    if (grepl("\\.\\.", input$selectFile) || grepl("/", input$selectFile) || grepl("\\\\", input$selectFile)) {
      showNotification("Invalid file name - path traversal not allowed", type = "error", duration = 10)
      return(NULL)
    }

    flink <- file.path(.dir(), input$selectFile)

    # 2. Check file existence
    if (!file.exists(flink)) {
      showNotification("File not found", type = "error", duration = 5)
      return(NULL)
    }

    # 3. Validate file extension
    allowed_ext <- c(".RDS", ".rds", ".db", ".sqlite", ".sqlite3")
    file_ext <- tolower(tools::file_ext(flink))
    if (!paste0(".", file_ext) %in% tolower(allowed_ext)) {
      showNotification(
        sprintf("Invalid file type. Allowed: %s", paste(allowed_ext, collapse = ", ")),
        type = "error",
        duration = 10
      )
      return(NULL)
    }

    # 4. Check file size and warn if too large
    sss <- file.size(flink)
    max_size <- 2e9  # 2GB limit
    if (sss > max_size) {
      showNotification(
        sprintf("File too large (%.1f GB). Maximum size is %.1f GB. Consider using database format.",
                sss/1e9, max_size/1e9),
        type = "error",
        duration = NULL
      )
      return(NULL)
    }

    # 5. Show progress for large files
    if (sss > 1e7)
      show_modal_spinner(text = "Loading data ...")

    # 6. Load with error handling
    v <- tryCatch(
      esetLoader(flink),
      error = function(e) {
        if (sss > 1e7) remove_modal_spinner()
        showNotification(
          sprintf("Error loading file: %s", e$message),
          type = "error",
          duration = NULL
        )
        return(NULL)
      },
      warning = function(w) {
        message("Warning during file loading: ", w$message)
      }
    )

    if (sss > 1e7)
      remove_modal_spinner()

    # 7. Validate loaded object
    if (is.null(v)) {
      showNotification("Failed to load data - file may be corrupted", type = "error", duration = 10)
      return(NULL)
    }

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
    i1 <- all(rownames(expr) == rownames(fd))
    i2 <- all(colnames(expr) == rownames(pd))
    if (!(i1 && i2))
      return(
        list(
          FALSE, "The rownames/colnames of exprs not matched to row names of feature data/phenotype data!"
        )
      )
    TRUE
  }
  
  vEset <- reactiveVal(FALSE)
  observe({    
    req(expr())
    req(pdata())
    req(fdata())
    x <- validEset(expr = expr(), pd = pdata(), fd = fdata())
    if (!x[[1]]) {
      showModal(modalDialog(
        title = "Problem in data!",
        x[[2]]
      ))
    } else {
      vEset( TRUE )
      shinyjs::show("contents")
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
  cormat <- reactive( {
    req(eset <- reactive_eset())
    attr(eset, "cormat") 
  })


  #####################
  
  output$download <- downloadHandler(
    filename = function() {
      paste0("ExpressenSet", Sys.time(), ".xlsx")
    },
    content = function(file) {
      td <- function(tab) {
        ic <- which(vapply(tab, is.list, logical(1)))
        if (length(ic) > 0) {
          for (ii in ic) {
            tab[, ii] <- vapply(tab[, ii], paste, collapse = ";", FUN.VALUE = character(1))
          }
        }
        id <- rownames(tab)
        if (is.null(id))
          id <- paste0("ID", seq_len(nrow(tab)))
        data.frame(ID = id, tab)
      }
      
      ig <- imputeGetter(reactive_eset())
      withProgress(message = 'Writing table', value = 0, {
        wb <- createWorkbook(creator = "BayBioMS")
        addWorksheet(wb, sheetName = "Phenotype info")
        addWorksheet(wb, sheetName = "Feature info")
        addWorksheet(wb, sheetName = "Expression")
        addWorksheet(wb, sheetName = "Geneset annot")
        incProgress(1/5, detail = "expression matrix")
        writeData(wb, sheet = "Expression", td(expr()))
        if (!is.null(ig)) {
          addWorksheet(wb, sheetName = "Expression_imputed")
          writeData(wb, sheet = "Expression_imputed", td(ig))
        }
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
    if (! vEset()) {
      txt <- sprintf(
      '<h1 style="display:inline;">%s</h1> <h3 style="display:inline;"><sup>%s</sup></h3>', 
      appName, paste0("v", appVersion))
    } else {
    txt <- sprintf(
      '<h1 style="display:inline;">%s</h1> <h3 style="display:inline;"><sup>%s</sup>  --   %s features and %s samples:</h3>', 
      appName, paste0("v", appVersion), nrow(expr()), ncol(expr()))
    }
    HTML(txt)
  })

  # v1 <- reactiveVal()
  v1 <- L1_data_space_module(
    "dataspace", expr = expr, pdata = pdata, fdata = fdata,
    reactive_x_s = d_s_x, reactive_y_s = d_s_y, reactive_x_f = d_f_x, reactive_y_f = d_f_y,
    status = esv_status, cormat = cormat
  )

  sameValues <- function(a, b) {
    if (is.null(a) || is.null(b))
      return(FALSE)
    all(sort(a) == sort(b))
    }
  ri <- reactiveVal()
  observeEvent( v1(), {
    ri( c(v1()$feature) )
    })
  observeEvent( expr(), ri(NULL) )

  rh <- reactiveVal()
  observeEvent( v1(), {
    rh( c( v1()$sample ) )
    })
  observeEvent( expr(), rh(NULL) )

  v2 <- L1_result_space_module("resultspace",
                   reactive_expr = expr,
                   reactive_phenoData = pdata,
                   reactive_featureData = fdata,
                   reactive_i = ri,
                   reactive_highlight = rh,
                   additionalTabs = additionalTabs,
                   object = reactive_eset,
                   status = esv_status)

  # =======================================================
  # =======================================================
  # ================= snapshot function ===================
  # =======================================================
  # =======================================================

  dir <- reactiveVal()
  observe({
    dd <- getwd()
    if (!is.null(.dir()))
      dd <- .dir()
    dir(dd)
    })

  savedSS <- reactiveVal()
  observe({
    # req(input$selectFile)
    req(.dir())
    if (is.null(input$selectFile) || nchar(input$selectFile) == 0)
      fs <- "ESVObj.RDS" else
        fs <- input$selectFile

    fl <- paste0("ESVSnapshot_", fs, "_")
    ff <- list.files(.dir(), pattern = fl)
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
    req( nrow(dt <- savedSS()) > 0 )    
    dt$delete <- shinyInput(
      actionButton, dt$name, 'deletess_', label = "Delete", onclick = sprintf('Shiny.setInputValue(\"%s\",  this.id)', ns("deletess_button")) 
      )
    DT::datatable(
      dt[, c(1, 3), drop = FALSE], rownames = FALSE, colnames = c(NULL, NULL, NULL), 
      selection = list(mode = "single", target = "cell", selectable = -cbind(seq_len(nrow(dt)), 1)), escape = FALSE,
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

  observeEvent(input$snapshot_save, {    
    # req(input$selectFile)
    if (is.null(input$selectFile) || nchar(input$selectFile) == 0)
      fs <- "ESVObj.RDS" else
        fs <- input$selectFile    

    df <- savedSS()
    if (!is.null(df)) {
      if (input$snapshot_name %in% df$name) {              
        showModal(modalDialog(
          title = "FAILED!",  
          "Snapshot with this name already exists, please give a different name."
          ))
        return(NULL)
      }
    }
    obj <- c(attr(v1(), "status"), v2(), active_feature = list(ri()), active_sample = list(rh()))
    flink <- file.path(.dir(), paste0("ESVSnapshot_", fs, "_", input$snapshot_name, ".ESS"))
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
    esv_status(NULL)  
    esv_status( readRDS(file.path(.dir(), df[i, 2])) )    
    })

  observeEvent(deleteSS(), {
    req(nrow( df <- savedSS() ) > 0 )
    req( i <- match( deleteSS(), df$name ))
    df <- savedSS()
    unlink(file.path(.dir(), df[i, 2]))
    df <- df[-i, , drop = FALSE]
    savedSS(df)
    })

  }) # end moduleServer
}



