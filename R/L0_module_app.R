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
    # omicsViewer theme: namespaced CSS (see auxi_uiTheme.R) plus a small
    # SVG favicon for standalone sessions. Everything is scoped through the
    # .omicsviewer-app root class set on the wrapper below, so embedding the
    # viewer into a host app neither overrides nor leaks these rules.
    tags$head(
      tags$link(
        rel = "icon",
        href = paste0(
          "data:image/svg+xml,",
          "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 64 64'>",
          "<rect width='64' height='64' rx='14' fill='%230e7490'/>",
          "<text x='32' y='44' font-family='Arial' font-size='34' ",
          "font-weight='bold' fill='white' text-anchor='middle'>oV</text>",
          "</svg>"
        )
      ),
      omicsviewer_ui_theme()
    ),
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
      top = 12, right = 20, style = "z-index: 9999;", width = 115,
      downloadButton(outputId = ns("download"), label = tagList(icon("file-excel"), "xlsx"),
        class = "omicsviewer-tool-btn", icon = NULL,
        title = "Download the complete dataset (expression matrix, feature and sample annotations, gene sets) as an Excel workbook") %>%
        tagAppendAttributes(`data-testid` = "app-download-dataset-button"),
      actionButton(ns("snapshot"), label = NULL, icon = icon("camera-retro"),
        class = "omicsviewer-tool-btn") %>%
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
  # scoping root for the namespaced theme (auxi_uiTheme.R); wrapping the
  # row keeps every rule inside this subtree when embedded in host apps
  tags$div(class = "omicsviewer-app", do.call(fluidRow, comp))
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
      # M4: convert SummarizedExperiment first (tallGS reads fData/pData,
      # which only exist on ExpressionSet) - documented ESVObj = se used to
      # fail here
      return( tallGS(asEsetWithAttr(ESVObj())) )
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

    # 6. Load with error handling. Warnings are surfaced as console
    # messages but must NOT abort the load: a tryCatch warning *handler*
    # would swallow the warning and return NULL, turning any coercion or
    # reshape warning into "file may be corrupted" (todo 1.4).
    v <- tryCatch(
      withCallingHandlers(
        esetLoader(flink),
        warning = function(w) {
          message("Warning during file loading: ", conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        if (sss > 1e7) remove_modal_spinner()
        showNotification(
          sprintf("Error loading file: %s", e$message),
          type = "error",
          duration = NULL
        )
        return(NULL)
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
    # L8: identical() with non-NULL checks - all(rownames(x) == ...) recycles
    # (a NULL side passes as all-TRUE) and a length-1 rowname side never fails
    i1 <- !is.null(rownames(expr)) && !is.null(rownames(fd)) &&
      identical(rownames(expr), rownames(fd))
    i2 <- !is.null(colnames(expr)) && !is.null(rownames(pd)) &&
      identical(colnames(expr), rownames(pd))
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
      paste0("ExpressionSet_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
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
      '<div class="omicsviewer-titleblock"><h1 class="omicsviewer-app-title"><i class="fa fa-dna" aria-hidden="true"></i>%s</h1> <h3 class="omicsviewer-app-sub"><sup>%s</sup></h3></div>',
      appName, paste0("v", appVersion))
    } else {
    txt <- sprintf(
      '<div class="omicsviewer-titleblock"><h1 class="omicsviewer-app-title"><i class="fa fa-dna" aria-hidden="true"></i>%s</h1> <h3 class="omicsviewer-app-sub"><sup>%s</sup>  --   %s features and %s samples:</h3></div>',
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

  # =====================================================================
  # Semantic, versioned snapshot state
  #
  # Defined BEFORE the L1 modules on purpose (todo 2.5): the dataset-switch
  # reset observer below must run BEFORE the module seed observers in the
  # same flush, so clearing the dataset-scoped store keys lands first and
  # the seeds (which skip held keys) refill with the NEW dataset's
  # defaults. Registered after the modules, the reset ran last and the
  # already-seeded keys were cleared with their one-shot gates consumed.
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
  .dsig_last <- reactiveVal(NULL)
  observeEvent(dataset_signature(), {
    sig <- dataset_signature()
    prev <- isolate(.dsig_last())
    if (identical(sig, prev)) return(NULL)
    .dsig_last(sig)
    esv_status(NULL)
    ri(NULL)
    rh(NULL)
    # Dataset switch (todo 2.5): reset the dataset-scoped widget-store keys
    # so dataset A's axes/params cannot leak into dataset B. The FIRST load
    # must not reset - the widget-default seeds may already have fired with
    # their one-shot gates consumed, and there is nothing to clear anyway.
    # After the reset the natural repair paths refill what has a
    # dataset-derived value (scatter axes re-seed from the new dataset's
    # configured defaults, table columns re-derive from the data); user
    # choices become legitimately unset until made again.
    if (!is.null(prev))
      tryCatch(store_reset(app_store, c("dataspace", "resultspace")),
               error = function(e) NULL)
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
  .expr_reset_last <- reactiveVal(NULL)
  observeEvent( expr(), {
    # same-value recomputes of the expression reactive must not clear the
    # live selection (observeEvent does not dedupe); dims + boundary ids
    # are the dataset-change granularity this reset means. ONE observer
    # for both spaces (N1): rh used to have a bare reset without the
    # signature gate, and a shared marker across two observers would let
    # the first consume the transition and spare the second
    e <- expr()
    sig <- paste(nrow(e), ncol(e), head(rownames(e), 1), tail(rownames(e), 1),
                 head(colnames(e), 1), tail(colnames(e), 1), sep = "|")
    if (identical(sig, isolate(.expr_reset_last()))) return(NULL)
    .expr_reset_last(sig)
    ri(NULL)
    rh(NULL)
  } )

  rh <- reactiveVal()
  observeEvent( v1(), {
    rh( c( v1()$sample) )
    })

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

  # N1: a `dir` reactiveVal here was written on every flush and never read
  # (the `.dir()` closure is the live accessor) - removed (todo 4.6/N1)

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

  # WP1 progressive disclosure: the compact state is built on demand for
  # exactly the requested sections. Expensive parts (annotation catalog,
  # figure grammar) are only computed when their section is asked for;
  # quick-view id/label lists and the widget-store scatter view always ride
  # along on the overview (plan section 3, settled decision). The whole
  # builder runs inside isolate() so it is safe to call from ellmer tool
  # contexts and test hooks outside the reactive graph.
  agent_state <- function(sections = NULL) {
    shiny::isolate({
      full_state <- agent_full_state()
      if (is.null(full_state))
        return(NULL)
      wanted <- sections %||% character()
      agent_compact_state(
        state = full_state,
        annotations = if ("annotations" %in% wanted)
          agent_annotation_catalog(fdata(), pdata()),
        quick_views = attr(v1(), "quickViews"),
        available_tabs = list(
          data_space = agent_data_tabs,
          analysis_space = agent_analysis_tabs()
        ),
        figure_grammar = if ("figure_grammar" %in% wanted)
          agent_figure_grammar(),
        sections = sections,
        store = app_store
      )
    })
  }

  apply_agent_state <- function(update) {
    # todo 4.1: the one write plane - tabs go through store_apply (the
    # same validation/diffing/user-wins path every other agent write
    # uses) and selections through the selection bus; the former
    # full-snapshot esv_status replay is retired.
    if (is.null(isolate(agent_state_available())))
      stop("No dataset is currently available.")

    validated <- agent_normalize_state_update(
      update = update,
      data_tabs = agent_data_tabs,
      analysis_tabs = isolate(agent_analysis_tabs()),
      feature_ids = rownames(isolate(expr())),
      sample_ids = colnames(isolate(expr()))
    )

    tab_patch <- list()
    if (!is.null(validated$data_space_tab))
      tab_patch[["dataspace.active_tab"]] <- validated$data_space_tab
    if (!is.null(validated$analysis_space_tab))
      tab_patch[["resultspace.analyst_tab"]] <- validated$analysis_space_tab
    if (length(tab_patch))
      store_apply(app_store, tab_patch, origin = "agent", strict = TRUE)

    # Selections: the bus is the single canonical layer, so an agent
    # selection is durable (the old status replay was transient - module
    # echoes reverted it within a flush). Mirror semantics: the agent's
    # action is a semantic selection, NOT a table interaction - figure and
    # heatmap origins mirror their ids into the table ROW FILTER, but the
    # historical agent contract (and the Tier A expectations) keep the
    # tables unfiltered (mirror TRUE) so table paging/sorting views
    # survive an agent selection.
    ports <- tryCatch(isolate(v1())$selection,
                      shiny.silent.error = function(e) NULL,
                      error = function(e) NULL)
    if (is.null(ports))
      stop("Data-space panel is not ready; retry in a moment.")
    if (!is.null(validated$features))
      ports$feature$apply(ids = validated$features, origin = "system",
                          mirror = TRUE)
    if (!is.null(validated$samples))
      ports$sample$apply(ids = validated$samples, origin = "system",
                         mirror = TRUE)

    list(
      data_space_tab = validated$data_space_tab,
      analysis_space_tab = validated$analysis_space_tab,
      feature_count = length(validated$features %||% character()),
      sample_count = length(validated$samples %||% character()),
      example_features = utils::head(validated$features %||% character(), 20L),
      example_samples = utils::head(validated$samples %||% character(), 20L)
    )
  }

  apply_agent_scatter_view <- function(space, quick_view_id = NULL,
                                       x_axis = NULL, y_axis = NULL) {
    if (is.null(isolate(agent_state_available())))
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

    axis_data <- if (view$space == "feature") isolate(fdata()) else isolate(pdata())
    if (!any(vapply(
      c(view$x_axis, view$y_axis),
      function(nm) is.numeric(axis_data[[nm]]),
      logical(1)
    ))) {
      stop("At least one scatter axis must contain numeric values.")
    }

    # todo 4.1: axes are WRITTEN through the canonical store keys (the
    # same keys get_omics_viewer_state READS its scatter_view anchors
    # from), never through the legacy panel-status transport. The
    # cascaded triple validates jointly; the display mode (quick badges
    # vs custom triselectors) is deliberately left untouched.
    prefix <- if (view$space == "feature") "dataspace.feature_space" else
      "dataspace.sample_space"
    split_axis <- function(axis) {
      parts <- strsplit(axis, "|", fixed = TRUE)[[1]]
      parts
    }
    patch <- stats::setNames(
      c(as.list(split_axis(view$x_axis)), as.list(split_axis(view$y_axis))),
      c(paste(prefix, c("x_analysis", "x_subset", "x_variable"), sep = "."),
        paste(prefix, c("y_analysis", "y_subset", "y_variable"), sep = ".")))
    store_apply(app_store, patch, origin = "agent", strict = TRUE)
    store_apply(
      app_store,
      list("dataspace.active_tab" = if (view$space == "feature") "Feature" else "Sample"),
      origin = "agent", strict = TRUE)

    view
  }

  # WP8: semantic tier-1 tools as thin validate + store_apply wrappers over
  # the canonical widget store (plan section 6.3). The normalizers in
  # auxi_agentCapabilities.R build full canonical-id patches; per-key
  # resilience keeps one invalid value (e.g. a pathway row that only exists
  # after the ranking recomputes) from vetoing the rest, with the
  # suggestion-bearing rejection returned to the model for self-correction.
  apply_agent_enrichment <- function(update) {
    validated <- agent_normalize_enrichment_update(update, isolate(fdata()))
    receipt <- shiny::isolate(
      store_apply(app_store, validated$patch, origin = "agent", strict = FALSE)
    )
    list(
      method = validated$method,
      panel_tab = validated$tab,
      applied = receipt$applied,
      applied_values = receipt$diff,
      unchanged = receipt$skipped,
      rejected = receipt$rejected %||% list(),
      note = paste(
        "The panel opens and the enrichment recomputes from the currently",
        "selected features (ORA) or the full ranking (fGSEA); re-check the",
        "results through the app tables or get_omics_viewer_state."
      )
    )
  }

  apply_agent_table_view <- function(update) {
    validated <- agent_normalize_table_view_update(update)
    receipt <- shiny::isolate(
      store_apply(app_store, validated$patch, origin = "agent", strict = FALSE)
    )
    list(
      table = validated$table,
      panel_tab = validated$tab,
      applied = receipt$applied,
      applied_values = receipt$diff,
      unchanged = receipt$skipped,
      rejected = receipt$rejected %||% list(),
      note = paste(
        "The table's tab opens with the requested view; filters match",
        "substring (case-insensitive) per column and paging starts at 1."
      )
    )
  }

  assistant_api <- ai_assistant_module(
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
    apply_enrichment = apply_agent_enrichment,
    apply_table_view = apply_agent_table_view,
    store = app_store
  )

  # Test-only: drive the exact agent apply callbacks from headless tests.
  if (agent_test_hooks_enabled())
    agent_test_hooks_module(
      "agentTestHooks",
      apply_state = apply_agent_state,
      apply_scatter_view = apply_agent_scatter_view,
      apply_enrichment = apply_agent_enrichment,
      apply_table_view = apply_agent_table_view,
      state = agent_state,
      store = app_store,
      assistant = assistant_api
    )

  # todo 4.4(c): snapshot save/restore extracted to snapshot_module()
  # (R/L0_module_snapshot.R); the panel-status delivery, restore
  # controller and receipt flow are unchanged.
  snapshot_module(
    input = input, output = output, session = session,
    .dir = .dir,
    dataset_id = current_dataset_id,
    dataset_ready = vEset,
    eset = reactive_eset,
    data_space = v1,
    result_status = v2,
    selection_features = ri,
    selection_samples = rh,
    esv_status = esv_status,
    store = app_store,
    assistant = assistant_api
  )



  }) # end moduleServer
}



