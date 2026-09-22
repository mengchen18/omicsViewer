#' omicsViewer Application UI (Level 0)
#'
#' @description
#' Generates the user interface for the main omicsViewer application. This function creates
#' a responsive layout with data exploration panels, snapshot functionality, data export,
#' and an optional AI assistant. Primarily intended for developers extending the application.
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
#'   \item Optional floating AI assistant launcher and settings/chat drawer
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
    # Skip navigation links for keyboard users and AI browsers
    tags$a(
      href = paste0("#", ns("main-content")),
      class = "sr-only sr-only-focusable",
      `data-testid` = "skip-to-main",
      "Skip to main content"
    ),
    tags$a(
      href = paste0("#", ns("data-panel")),
      class = "sr-only sr-only-focusable",
      "Skip to data exploration"
    ),
    tags$a(
      href = paste0("#", ns("analysis-panel")),
      class = "sr-only sr-only-focusable",
      "Skip to analysis tools"
    ),
    # CSS for hiding sr-only elements (screen reader and AI browser only)
    tags$head(
      tags$style(HTML("
        .sr-only {
          position: absolute !important;
          width: 1px !important;
          height: 1px !important;
          padding: 0 !important;
          margin: -1px !important;
          overflow: hidden !important;
          clip: rect(0, 0, 0, 0) !important;
          white-space: nowrap !important;
          border: 0 !important;
        }
        .sr-only-focusable:focus {
          position: static !important;
          width: auto !important;
          height: auto !important;
          overflow: visible !important;
          clip: auto !important;
          white-space: normal !important;
        }
      "))
    ),
    # JSON-LD schema for AI browsers and machine readability
    tags$head(
      tags$script(
        type = "application/ld+json",
        HTML('{
          "@context": "https://schema.org",
          "@type": "WebApplication",
          "name": "omicsViewer",
          "description": "Interactive visualization and analysis platform for omics data (proteomics, transcriptomics, genomics). Supports multi-dimensional data exploration, statistical testing, pathway enrichment, and network analysis.",
          "applicationCategory": "Bioinformatics",
          "applicationSubCategory": "Omics Data Analysis",
          "operatingSystem": "Web browser",
          "offers": {
            "@type": "Offer",
            "price": "0",
            "priceCurrency": "USD"
          },
          "featureList": [
            "2D scatter plots with correlation analysis and regression lines for feature and sample metadata",
            "Interactive boxplots with statistical tests (t-test, ANOVA, Kruskal-Wallis) for group comparisons",
            "Gene set over-representation analysis (ORA) using hypergeometric test",
            "Fast gene set enrichment analysis (fGSEA) with leading edge identification",
            "Protein-protein interaction network visualization via STRING database integration",
            "Literature-based gene association discovery through Geneshot API",
            "Post-translational modification (PTM) motif enrichment analysis",
            "Dose-response curve fitting with EC50/IC50 estimation using 4-parameter logistic model",
            "Kaplan-Meier survival analysis with log-rank test",
            "ROC and precision-recall curve generation for binary classification",
            "Correlation heatmaps with hierarchical clustering",
            "Expression heatmaps with dendrogram and sample grouping",
            "Dynamic heatmap with interactive row selection and subsetting",
            "Contingency table analysis with chi-square and Fisher exact tests",
            "Searchable data tables for features, samples, expression values, and gene sets",
            "State snapshot management for reproducible analysis workflows",
            "Optional session-local AI assistant with bounded state inspection and validated view updates",
            "AI-generated declarative ggplot2 figures with in-chat previews and high-resolution downloads",
            "Data export to Excel format with all annotations"
          ],
          "softwareRequirements": "Modern web browser with JavaScript enabled",
          "permissions": "No special permissions required"
        }')
      )
    ),
    style = "background:white;",
    absolutePanel(
      top = 5, right = 20, style = "z-index: 9999;", width = 115,
      downloadButton(outputId = ns("download"), label = "xlsx", class = NULL) %>%
        tagAppendAttributes(`data-testid` = "app-download-dataset-button"),
      actionButton(ns("snapshot"), label = NULL, icon = icon("camera-retro")) %>%
        tagAppendAttributes(`data-testid` = "app-snapshot-button",
                           title = "Manage snapshots")
    ),
    ai_assistant_ui(ns("assistant")),
    # Headless-test-only bridge hooks; never rendered unless explicitly
    # enabled via OMICSVIEWER_TEST_HOOKS for the assistant UI-effect suite.
    if (agent_test_hooks_enabled())
      agent_test_hooks_ui(ns("agentTestHooks")),
    # Main content area with semantic HTML
    tags$main(
      role = "main",
      id = ns("main-content"),
      `aria-label` = "Main application content",
      shinyjs::hidden(
        div(id = ns("contents"),
          tags$section(
            `aria-label` = "Data exploration panel",
            id = ns("data-panel"),
            column(6, L1_data_space_ui(ns('dataspace'), activeTab = activeTab))
          ),
          tags$section(
            `aria-label` = "Analysis panel",
            id = ns("analysis-panel"),
            column(6, L1_result_space_ui(ns("resultspace")))
          )
        )
      )
    )
    )

  if (showDropList) {
    l2 <- list(
      shinycssloaders::withSpinner(
        uiOutput(ns("summary")), hide.ui = FALSE, type = 8, color = "green"
        ),
      # Aria-live region for loading status announcements
      div(
        class = "sr-only",
        `aria-live` = "polite",
        `aria-atomic` = "true",
        uiOutput(ns("loadingStatus"))
      ),
      br(),
      # Navigation element for dataset selection
      tags$nav(
        `aria-label` = "Dataset selection",
        absolutePanel(
          top = 8, right = 140, style = "z-index: 9999;",
          selectizeInput( inputId = ns("selectFile"), label = NULL, choices = NULL,
            width = "500px", options = list(placeholder = "Select a dataset here") ) %>%
            tagAppendAttributes(`data-testid` = "app-dataset-selector")
        )
      ))
    comp <- c(l2, comp)
    }
  do.call(fluidRow, comp)
}

#' omicsViewer Application Server Logic (Level 0)
#'
#' @description
#' Implements the main server-side logic for the omicsViewer Shiny application. Handles data
#' loading, validation, state management, snapshot functionality, and orchestrates communication
#' between sub-modules. Uses modern Shiny module pattern with \code{moduleServer}.
#' Primarily intended for developers extending the application.
#'
#' @param id Character. Namespace ID for the Shiny module. Must match the ID used in
#'   \code{\link{app_ui}}.
#' @param .dir Reactive expression. Returns the directory path containing data files
#'   (ExpressionSet or SummarizedExperiment .RDS files).
#' @param filePattern Character. Regular expression to filter displayed files.
#'   Default: \code{".(RDS|db|sqlite|sqlite3)$"} (case-insensitive).
#' @param additionalTabs List or NULL. Custom analysis modules to add to the application.
#'   Each element should contain: \code{tabName}, \code{moduleName}, \code{moduleUi}, and
#'   \code{moduleServer}. Default: NULL (no additional tabs).
#' @param ESVObj Reactive expression. Returns a pre-loaded ExpressionSet or SummarizedExperiment
#'   object, bypassing file loading. Default: \code{reactive(NULL)}.
#' @param esetLoader Function. Loads data objects from disk. Takes file path as input,
#'   returns ExpressionSet or SummarizedExperiment. Default: \code{readESVObj}.
#' @param exprsGetter Function. Extracts expression matrix from loaded object.
#'   Default: \code{getExprs}.
#' @param imputeGetter Function. Extracts imputed expression matrix (if available) for
#'   Excel export. Should return NULL if no imputed data. Default: \code{getExprsImpute}.
#' @param pDataGetter Function. Extracts sample/phenotype metadata. Default: \code{getPData}.
#' @param fDataGetter Function. Extracts feature metadata. Default: \code{getFData}.
#' @param defaultAxisGetter Function. Determines default plot axes. Takes object and
#'   \code{what} ("sx", "sy", "fx", "fy") as arguments. Default: \code{getAx}.
#' @param appName Character. Application name displayed in UI. Default: "omicsViewer".
#' @param appVersion Character or package_version. Version shown in UI.
#'   Default: current package version.
#'
#' @details
#' The module coordinates several key functionalities:
#' \itemize{
#'   \item \strong{Data Loading}: Validates file paths, checks file sizes, loads with error handling
#'   \item \strong{Data Validation}: Ensures rownames/colnames consistency across expression and metadata
#'   \item \strong{State Management}: Tracks selected features/samples across sub-modules
#'   \item \strong{Snapshots}: Save and restore analysis states to disk (.ESS files)
#'   \item \strong{Data Export}: Generate Excel files with expression data, metadata, and gene sets
#'   \item \strong{Module Coordination}: Manages data space (L1_data_space_module) and
#'         result space (L1_result_space_module) interactions
#'   \item \strong{AI Assistant}: Optionally connects a session-local ellmer chat to
#'         bounded application-state and annotation tools
#' }
#'
#' Security features include path traversal prevention, file type validation,
#' and size limits (2GB maximum).
#'
#' @return
#' NULL (invisibly). The module manages reactive state internally and communicates
#' with child modules. No explicit return value.
#'
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
#'
#' @seealso
#' \code{\link{app_ui}} for the corresponding UI function.
#' \code{\link{L1_data_space_module}}, \code{\link{L1_result_space_module}} for sub-modules.
#' \code{\link{omicsViewer}} for the high-level launcher.
#'
#' @examples
#' if (interactive()) {
#'   dir <- system.file("extdata", package = "omicsViewer")
#'   ui <- fluidPage(app_ui("app"))
#'   server <- function(input, output, session) {
#'     app_module("app", .dir = reactive(dir))
#'   }
#'   shinyApp(ui = ui, server = server)
#' }
#' @keywords internal

app_module <- function(
  id, .dir, filePattern = ".(RDS|db|sqlite|sqlite3)$", additionalTabs = NULL, ESVObj = reactive(NULL),
  esetLoader = readESVObj, exprsGetter = getExprs, pDataGetter = getPData, fDataGetter = getFData,
  imputeGetter = getExprsImpute, defaultAxisGetter = getAx,
  appName = "omicsViewer", appVersion = packageVersion("omicsViewer"),
  store = NULL
) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns

  # Canonical widget store: created here by default; callers may inject
  # their own (embedding contexts and headless tests observe it directly).
  app_store <- store %||% widget_store_new()

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

  # Loading status for screen readers and AI browsers
  output$loadingStatus <- renderUI({
    if (!is.null(reactive_eset()) && vEset()) {
      tags$span(sprintf(
        "Dataset loaded successfully: %s features and %s samples",
        nrow(expr()), ncol(expr())
      ))
    } else if (!is.null(input$selectFile) && nchar(input$selectFile) > 0) {
      tags$span("Loading dataset, please wait...")
    } else {
      tags$span("No dataset selected. Please select a dataset to begin.")
    }
  })

  v1 <- L1_data_space_module(
    "dataspace", expr = expr, pdata = pdata, fdata = fdata,
    reactive_x_s = d_s_x, reactive_y_s = d_s_y, reactive_x_f = d_f_x, reactive_y_f = d_f_y,
    status = reactive(esv_status()$panels$data_space), cormat = cormat,
    store = app_store
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
                   status = reactive(esv_status()$panels$result_space),
                   store = app_store)

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

  # =====================================================================
  # Semantic, versioned snapshot state
  # =====================================================================
  esv_status <- reactiveVal(NULL)

  current_dataset_id <- reactive({
    if (!is.null(ESVObj()))
      return("ESVObj.RDS")
    if (is.null(input$selectFile) || !nzchar(input$selectFile))
      return("ESVObj.RDS")
    input$selectFile
  })

  # Dataset changes invalidate all previous panel state. Child modules receive
  # NULL first and the restored object second, making restoration transactional
  # at the top-level state boundary.
  dataset_signature <- reactive({
    req(reactive_eset())
    paste(current_dataset_id(), dataset_fingerprint(reactive_eset(), id = current_dataset_id()))
  })
  observeEvent(dataset_signature(), {
    esv_status(NULL)
    ri(NULL)
    rh(NULL)
  })

  # =====================================================================
  # Optional session-local AI assistant
  #
  # The state bridge reuses the semantic snapshot representation, but stays
  # outside persisted .ESS files. Model tools receive this compact reactive
  # snapshot and can propose only narrowly validated UI transitions.
  # =====================================================================
  agent_data_tabs <- c(
    "Feature", "Feature table", "Sample", "Sample table", "Cor",
    "Heatmap", "Dynamic heatmap", "Expression", "GSList"
  )

  agent_analysis_tabs <- reactive({
    req(fdata())
    tabs <- "Feature"
    if (!is.null(attr(fdata(), "GS")))
      tabs <- c(tabs, "ORA", "fGSEA")
    if (any(grepl("^ResponseCurve\\|", colnames(fdata()))))
      tabs <- c(tabs, "Response")
    if (any(grepl("^StringDB\\|", colnames(fdata()))))
      tabs <- c(tabs, "StringDB")
    if (any(grepl("^SeqLogo\\|", colnames(fdata()))))
      tabs <- c(tabs, "SeqLogo")
    if (length(additionalTabs) > 0)
      tabs <- c(tabs, vapply(additionalTabs, function(x) x$tabName, character(1)))
    unique(c(tabs, "Geneshot", "Sample"))
  })

  agent_full_state <- reactive({
    req(vEset())
    req(reactive_eset())

    data_status <- tryCatch(
      attr(v1(), "status"),
      shiny.silent.error = function(e) list(),
      error = function(e) stop(e)
    )
    result_status <- tryCatch(
      v2(),
      shiny.silent.error = function(e) list(),
      error = function(e) stop(e)
    )
    build_app_state(
      dataset = reactive_eset(),
      dataset_id = current_dataset_id(),
      data_status = data_status,
      result_status = result_status,
      selected_features = ri(),
      selected_samples = rh(),
      label = "AI assistant current state"
    )
  })

  agent_state_available <- reactive({
    req(reactive_eset())
    req(vEset())
    TRUE
  })

  agent_state <- reactive({
    full_state <- agent_full_state()
    agent_compact_state(
      state = full_state,
      annotations = agent_annotation_catalog(fdata(), pdata()),
      quick_views = attr(v1(), "quickViews"),
      available_tabs = list(
        data_space = agent_data_tabs,
        analysis_space = agent_analysis_tabs()
      ),
      figure_grammar = agent_figure_grammar()
    )
  })

  apply_agent_state <- function(update) {
    full_state <- isolate(agent_full_state())
    if (is.null(full_state))
      stop("No dataset is currently available.")

    validated <- agent_normalize_state_update(
      update = update,
      data_tabs = agent_data_tabs,
      analysis_tabs = isolate(agent_analysis_tabs()),
      feature_ids = rownames(isolate(expr())),
      sample_ids = colnames(isolate(expr()))
    )

    if (!is.null(validated$data_space_tab)) {
      full_state$app$data_active_tab <- validated$data_space_tab
      full_state$panels$data_space$eset_active_tab <- validated$data_space_tab
    }
    if (!is.null(validated$analysis_space_tab)) {
      full_state$app$analysis_active_tab <- validated$analysis_space_tab
      full_state$panels$result_space$analyst_active_tab <- validated$analysis_space_tab
    }
    if (!is.null(validated$features))
      full_state$selection$features <- validated$features
    if (!is.null(validated$samples))
      full_state$selection$samples <- validated$samples

    # Use the same transactional boundary as snapshot restoration. Child
    # modules distinguish NULL from a state object and safely fill gaps.
    esv_status(NULL)
    esv_status(full_state)
    ri(full_state$selection$features)
    rh(full_state$selection$samples)

    list(
      data_space_tab = if (is.null(validated$data_space_tab)) NULL else full_state$app$data_active_tab,
      analysis_space_tab = if (is.null(validated$analysis_space_tab)) NULL else full_state$app$analysis_active_tab,
      feature_count = length(full_state$selection$features),
      sample_count = length(full_state$selection$samples),
      example_features = utils::head(full_state$selection$features, 20L),
      example_samples = utils::head(full_state$selection$samples, 20L)
    )
  }

  apply_agent_scatter_view <- function(space, quick_view_id = NULL,
                                       x_axis = NULL, y_axis = NULL) {
    full_state <- isolate(agent_full_state())
    if (is.null(full_state))
      stop("No dataset is currently available.")

    view <- agent_normalize_scatter_view(
      space = space,
      quick_view_id = quick_view_id,
      x_axis = x_axis,
      y_axis = y_axis,
      quick_views = isolate(attr(v1(), "quickViews")),
      feature_columns = colnames(isolate(fdata())),
      sample_columns = colnames(isolate(pdata()))
    )

    state_key <- if (view$space == "feature") "eset_fdata_fig" else "eset_pdata_fig"
    axis_data <- if (view$space == "feature") isolate(fdata()) else isolate(pdata())
    if (!any(vapply(
      c(view$x_axis, view$y_axis),
      function(nm) is.numeric(axis_data[[nm]]),
      logical(1)
    ))) {
      stop("At least one scatter axis must contain numeric values.")
    }

    split_axis <- function(axis) {
      parts <- strsplit(axis, "|", fixed = TRUE)[[1]]
      as.list(stats::setNames(parts, c("v1", "v2", "v3")))
    }
    # Deliberately do NOT write axisMode: a scatter-view change updates the
    # axes only, and the display mode (quick badges vs custom triselectors)
    # is left exactly as the user set it (plan section 6, diff-only writes).
    full_state$panels$data_space[[state_key]]$xax <- split_axis(view$x_axis)
    full_state$panels$data_space[[state_key]]$yax <- split_axis(view$y_axis)
    full_state$app$data_active_tab <- if (view$space == "feature") "Feature" else "Sample"
    full_state$panels$data_space$eset_active_tab <- full_state$app$data_active_tab

    esv_status(NULL)
    esv_status(full_state)
    ri(full_state$selection$features)
    rh(full_state$selection$samples)

    view
  }

  ai_assistant_module(
    "assistant",
    state = agent_state,
    state_available = agent_state_available,
    feature_data = fdata,
    sample_data = pdata,
    expression_data = expr,
    selected_features = ri,
    selected_samples = rh,
    apply_state = apply_agent_state,
    apply_scatter_view = apply_agent_scatter_view,
    store = app_store
  )

  # Test-only: drive the exact agent apply callbacks from headless tests.
  if (agent_test_hooks_enabled())
    agent_test_hooks_module(
      "agentTestHooks",
      apply_state = apply_agent_state,
      apply_scatter_view = apply_agent_scatter_view,
      state = agent_state,
      store = app_store
    )

  savedSS <- reactiveVal(
    data.frame(name = character(), link = character(), schema = integer(),
               created_at = character(), package_version = character(),
               stringsAsFactors = FALSE)
  )
  snapshot_refresh <- reactiveVal(0L)

  observe({
    req(.dir())
    snapshot_refresh()
    dsid <- sanitize_snapshot_name(current_dataset_id(), fallback = "ESVObj.RDS")
    prefix <- paste0("ESVSnapshot_", dsid, "_")
    ff <- list.files(.dir(), pattern = "\\.ESS$", ignore.case = TRUE)
    ff <- ff[startsWith(ff, prefix)]

    if (length(ff) == 0) {
      savedSS(data.frame(name = character(), link = character(), schema = integer(),
                         created_at = character(), package_version = character(),
                         stringsAsFactors = FALSE))
      return(NULL)
    }

    # Read compact metadata only; large panel payloads stay on disk.
    meta <- lapply(ff, function(f) tryCatch({
      x <- readRDS(file.path(.dir(), f))
      list(
        schema = if (is.null(x$schema_version)) NA_integer_ else as.integer(x$schema_version),
        created_at = x$created_at %.or_default% NA_character_,
        package_version = x$package_version %.or_default% NA_character_
      )
    }, error = function(e) list(schema = NA_integer_, created_at = NA_character_,
                                package_version = NA_character_)))

    savedSS(data.frame(
      name = sub("\\.ESS$", "", sub(prefix, "", ff)),
      link = ff,
      schema = vapply(meta, function(x) x$schema, integer(1)),
      created_at = vapply(meta, function(x) x$created_at, character(1)),
      package_version = vapply(meta, function(x) x$package_version, character(1)),
      stringsAsFactors = FALSE
    ))
  })

  shinyInput <- function(FUN, len, id, ...) {
    inputs <- c()
    for (i in len) {
      inputs <- c(inputs, as.character(FUN(paste0(id, i), ...)))
    }
    inputs
  }

  output$tab_saveSS <- renderDT({
    req(nrow(dt <- savedSS()) > 0)
    dt$delete <- shinyInput(
      actionButton, dt$name, "deletess_", label = "Delete",
      onclick = sprintf('Shiny.setInputValue("%s", this.id)', ns("deletess_button"))
    )
    dt$info <- ifelse(
      is.na(dt$schema),
      "legacy",
      paste0(
        "v", dt$schema,
        ifelse(is.na(dt$created_at), "", paste0(" | ", dt$created_at)),
        ifelse(is.na(dt$package_version), "", paste0(" | ", dt$package_version))
      )
    )
    DT::datatable(
      dt[, c("name", "info", "delete"), drop = FALSE],
      rownames = FALSE, colnames = c("Name", "Metadata", ""),
      selection = list(mode = "single", target = "cell",
                       selectable = -cbind(seq_len(nrow(dt)), 3)),
      escape = FALSE,
      options = list(
        dom = "t", autoWidth = FALSE, style = "compact-hover", scrollY = "450px",
        paging = FALSE,
        columns = list(list(width = "40%"), list(width = "42%"), list(width = "18%"))
      )
    )
  })

  selectedSS <- reactiveVal()
  observe({
    ss <- input$tab_saveSS_cells_selected
    if (length(ss) == 0 || ss[2] > 1)
      return(NULL)
    selectedSS(ss[1])
  })

  observeEvent(list(v1(), v2()), {
    selectedSS(NULL)
  })

  deleteSS <- reactiveVal()
  observeEvent(input$deletess_button, {
    selectedRow <- sub("deletess_", "", input$deletess_button, fixed = TRUE)
    deleteSS(selectedRow)
  })

  observeEvent(deleteSS(), {
    req(nrow(df <- savedSS()) > 0)
    req(i <- match(deleteSS(), df$name))
    showModal(modalDialog(
      title = "Delete snapshot",
      sprintf("Delete snapshot %s? This cannot be undone.", df$name[i]),
      footer = tagList(
        actionButton(ns("snapshot_delete_cancel"), "Cancel"),
        actionButton(ns("snapshot_delete_confirm"), "Delete", class = "btn-danger")
      ),
      easyClose = TRUE
    ))
  })

  observeEvent(input$snapshot_delete_cancel, {
    removeModal()
    deleteSS(NULL)
  })

  observeEvent(input$snapshot_delete_confirm, {
    req(nrow(df <- savedSS()) > 0)
    req(i <- match(deleteSS(), df$name))
    unlink(file.path(.dir(), df$link[i]))
    removeModal()
    deleteSS(NULL)
    snapshot_refresh(snapshot_refresh() + 1L)
  })

  observeEvent(input$snapshot, {
    showModal(
      modalDialog(
        title = NULL,
        fluidRow(
          column(9, textInput(ns("snapshot_name"), label = "Save new snapshot", placeholder = "snapshot name", width = "100%")),
          column(3, style = "padding-top:25px", actionButton(ns("snapshot_save"), label = "Save"))
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
    req(vEset())
    req(reactive_eset())
    name <- sanitize_snapshot_name(input$snapshot_name, fallback = paste0("snapshot-", format(Sys.time(), "%Y%m%d-%H%M%S")))

    df <- savedSS()
    if (name %in% df$name) {
      showNotification(sprintf("Snapshot name %s is already in use.", name), type = "error")
      return(NULL)
    }

    flink <- file.path(.dir(), snapshot_file_name(name, dataset_id = current_dataset_id(), fallback = name))
    if (file.exists(flink)) {
      showNotification("A snapshot file with this name already exists.", type = "error")
      return(NULL)
    }

    # Child state reactives can contain unmet req() conditions while optional
    # panels initialize (notably in testServer/headless sessions). Treat those
    # silent validation results as empty panel state rather than losing the
    # snapshot; genuine errors still abort below.
    data_status <- tryCatch(
      attr(v1(), "status"),
      shiny.silent.error = function(e) list(),
      error = function(e) stop(e)
    )
    result_status <- tryCatch(
      v2(),
      shiny.silent.error = function(e) list(),
      error = function(e) stop(e)
    )
    obj <- build_app_state(
      dataset = reactive_eset(),
      dataset_id = current_dataset_id(),
      data_status = data_status,
      result_status = result_status,
      selected_features = ri(),
      selected_samples = rh(),
      label = name
    )
    # Canonical widget-store state rides along (S4 start): keeps every
    # registered widget's desired value in one authoritative snapshot.
    obj$widget_store <- store_snapshot(app_store)
    write_app_state(obj, flink)
    snapshot_refresh(snapshot_refresh() + 1L)
    removeModal()
    showNotification(sprintf("Snapshot %s saved.", name), type = "message", duration = 3)
  })

  observeEvent(selectedSS(), {
    req(vEset())
    req(nrow(df <- savedSS()) > 0)
    if (length(i <- selectedSS()) == 0)
      return(NULL)

    removeModal()
    ss <- tryCatch(readRDS(file.path(.dir(), df$link[i])), error = function(e) {
      showNotification("Could not read the selected snapshot.", type = "error")
      NULL
    })
    req(ss)

    ss <- tryCatch(
      validate_app_state(ss, dataset = reactive_eset(), dataset_id = current_dataset_id()),
      error = function(e) {
        showNotification(paste("Invalid snapshot:", conditionMessage(e)), type = "error")
        NULL
      },
      warning = function(w) {
        showNotification(paste("Snapshot restored with warnings:", conditionMessage(w)), type = "warning", duration = 10)
        suppressWarnings(validate_app_state(ss, dataset = reactive_eset(), dataset_id = current_dataset_id()))
      }
    )
    req(ss)

    # Reset before restore. Child modules distinguish this boundary and can
    # safely restore defaults for fields absent from an older snapshot.
    esv_status(NULL)
    esv_status(ss)

    selection <- normalize_selection(ss$selection)
    ri(selection$features)
    rh(selection$samples)

    # Canonical widget-store restore (S4 start): applies per-key-resilient
    # through the same transactional protocol the agent uses. Panel-status
    # restoration above stays authoritative for not-yet-migrated modules;
    # diff-only writes make the overlap idempotent. Older snapshots without
    # widget_store skip this.
    if (!is.null(ss$widget_store) && is.list(ss$widget_store$values)) {
      tryCatch(
        store_restore(app_store, ss$widget_store),
        error = function(e)
          warning("Widget-store snapshot restore failed: ", conditionMessage(e))
      )
    }
  })


  }) # end moduleServer
}



