#' @description Application level 1 ui - data space
#' @param id id
#' @param activeTab one of "Feature", "Feature table", "Sample", "Sample table", "Heatmap"
#' @importFrom shinythemes shinytheme
#' @importFrom shinyWidgets tooltipOptions dropdown
#' 
L1_data_space_ui <- function(id, activeTab = "Feature") {
  ns <- NS(id)
  navbarPage(
    "Data", id = ns("eset"),
    selected = activeTab,
    theme = shinytheme("spacelab"),
    tabPanel('Feature',
      tags$h2("Feature Scatter Plot", class = "sr-only", `aria-label` = "2D scatter plot visualization of feature metadata with correlation analysis and interactive selection"),
      div(class = "sr-only",
        tags$h4("About Feature Scatter Plot"),
        tags$p("Feature scatter plots visualize relationships between two feature-level variables (e.g., fold change vs. p-value for volcano plots, or PCA coordinates). Each point represents one gene/protein. This view helps identify interesting features, assess data quality, and select features for downstream analysis."),
        tags$h4("When to use this view"),
        tags$p("Use feature scatter plots to explore statistical test results, visualize dimensionality reduction (PCA, t-SNE), create volcano plots, or examine correlations between feature annotations. Click and drag to select features of interest for analysis in the right panel."),
        tags$h4("How to interpret"),
        tags$p("Points represent individual features. Colors and shapes can encode additional variables. Use the lasso or box selection tools to select features. Regression lines show trends. Correlation coefficient indicates strength of linear relationship.")
      ),
      meta_scatter_ui(ns('feature_space'))),
    tabPanel("Feature table",
      tags$h2("Feature Metadata Table", class = "sr-only", `aria-label` = "Searchable table of feature annotations including gene names, protein IDs, statistical test results, and functional classifications"),
      div(class = "sr-only",
        tags$h4("About Feature Metadata Table"),
        tags$p("The feature table displays all annotations for genes/proteins in your dataset. Each row represents one feature with columns for identifiers, statistical results, functional annotations, and gene set memberships. This is your primary reference for feature information."),
        tags$h4("When to use this view"),
        tags$p("Use the feature table to search for specific genes, sort by statistical significance, filter by fold change thresholds, or review detailed annotations. Export functionality allows you to save filtered results for reports or further analysis."),
        tags$h4("How to navigate"),
        tags$p("Use the search box to find specific genes. Click column headers to sort. Use filters to subset data. The table is fully interactive and updates based on your selections in other views.")
      ),
      dataTable_ui(ns("tab_feature"))),
    tabPanel("Sample",
      tags$h2("Sample Scatter Plot", class = "sr-only", `aria-label` = "2D scatter plot visualization of sample metadata with group comparisons and interactive selection"),
      div(class = "sr-only",
        tags$h4("About Sample Scatter Plot"),
        tags$p("Sample scatter plots visualize relationships between sample-level variables or show sample distributions in reduced dimensions (PCA, t-SNE, UMAP). Each point represents one sample/patient. This view reveals sample clustering, batch effects, and relationships between clinical variables."),
        tags$h4("When to use this view"),
        tags$p("Use sample scatter plots to assess data quality, identify outliers, visualize experimental groups, explore PCA results, or examine correlations between clinical variables. Select samples to highlight them across all visualizations."),
        tags$h4("How to interpret"),
        tags$p("Points represent individual samples. Clustering indicates similarity. Colors typically show experimental groups or phenotypes. Distance between points reflects dissimilarity in the plotted variables.")
      ),
      meta_scatter_ui(ns("sample_space"))),
    tabPanel("Sample table",
      tags$h2("Sample Metadata Table", class = "sr-only", `aria-label` = "Searchable table of sample annotations including experimental conditions, phenotypes, and clinical variables"),
      div(class = "sr-only",
        tags$h4("About Sample Metadata Table"),
        tags$p("The sample table displays all clinical and experimental annotations for samples in your dataset. Each row represents one sample with columns for identifiers, experimental conditions, phenotypes, and other metadata. This is your primary reference for sample information."),
        tags$h4("When to use this view"),
        tags$p("Use the sample table to review experimental design, verify sample annotations, identify samples with specific characteristics, or export metadata for external analysis. This is essential for validating groupings and understanding your cohort."),
        tags$h4("How to navigate"),
        tags$p("Search for specific samples or conditions. Sort by any column. Filter to find samples matching criteria. All tables are interactive and synchronize with other views.")
      ),
      dataTable_ui(ns("tab_pheno"))),
    tabPanel(
      "Cor",
      tags$h2("Correlation Heatmap", class = "sr-only", `aria-label` = "Interactive correlation matrix heatmap showing pairwise correlations between samples with hierarchical clustering dendrogram"),
      div(class = "sr-only",
        tags$h4("About Correlation Heatmap"),
        tags$p("The correlation heatmap displays pairwise correlations between all samples in your dataset. Hierarchical clustering organizes samples with similar expression profiles together, revealing sample relationships, experimental groupings, and potential batch effects. This is a key quality control and exploratory visualization."),
        tags$h4("When to use this view"),
        tags$p("Use the correlation heatmap to assess overall data quality, identify outlier samples, validate experimental groupings, detect batch effects, or explore sample relationships. Strong blocks in the heatmap indicate groups of similar samples."),
        tags$h4("How to interpret"),
        tags$p("Color intensity indicates correlation strength (darker = stronger). The dendrogram shows hierarchical relationships. Samples cluster based on expression similarity. Well-designed experiments show clear blocks corresponding to biological groups.")
      ),
      tags$h3("Heatmap Controls", class = "sr-only"),
      fluidRow(
        column(
          6,
          dropdown(
            inputId = "mydropdown2",
            label = "Controls",
            circle = FALSE, status = "default", icon = icon("cog"),
            width = 700,
            tooltip = tooltipOptions(title = "Click to update heatmap and check legend!"),
            margin = "10px",
            tabsetPanel(
              tabPanel("Parameters", iheatmapInput(id = ns("corheatmapViewer"), scaleOn = "none")),
              tabPanel("Legend", iheatmapLegend(id = ns("corheatmapViewer")))
            )
          )
        ),
        column(
          6, align = "right",
          iheatmapClear(id = ns("corheatmapViewer"))
        ),
        tags$h3("Heatmap Visualization", class = "sr-only"),
        column(
          12, iheatmapOutput(id = ns("corheatmapViewer"))
        )
      )
    ),
    tabPanel(
      "Heatmap",
      tags$h2("Expression Heatmap", class = "sr-only", `aria-label` = "Interactive expression heatmap showing feature abundance across samples with hierarchical clustering and customizable color scales"),
      div(class = "sr-only",
        tags$h4("About Expression Heatmap"),
        tags$p("The expression heatmap displays abundance values for features (rows) across samples (columns). Hierarchical clustering groups similar features and samples together, revealing expression patterns, co-regulated genes, and sample groupings. This is essential for exploring large-scale expression patterns."),
        tags$h4("When to use this view"),
        tags$p("Use the expression heatmap to visualize patterns in your selected features, identify co-expressed gene clusters, validate sample groupings, or create publication-quality figures. Customize colors, clustering, and annotations through the controls panel."),
        tags$h4("How to interpret"),
        tags$p("Each row is a feature, each column is a sample. Color indicates expression level (red = high, blue = low by default). Dendrograms show clustering relationships. Blocks of similar color indicate coordinated expression patterns across features or samples.")
      ),
      tags$h3("Heatmap Controls", class = "sr-only"),
      fluidRow(
        column(
          6,
          dropdown(
            inputId = "mydropdown",
            label = "Controls",
            circle = FALSE, status = "default", icon = icon("cog"),
            width = 700,
            tooltip = tooltipOptions(title = "Click to update heatmap and check legend!"),
            margin = "10px",
            tabsetPanel(
              tabPanel("Parameters", iheatmapInput(id = ns("heatmapViewer"))),
              tabPanel("Legend", iheatmapLegend(id = ns("heatmapViewer")))
            )
          )
        ),
        column(
          6, align = "right",
          iheatmapClear(id = ns("heatmapViewer"))
        ),
        tags$h3("Heatmap Visualization", class = "sr-only"),
        column(
          12, iheatmapOutput(id = ns("heatmapViewer"))
        )
      )
    ),
    tabPanel(
      "Dynamic heatmap",
      tags$h2("Dynamic Heatmap with Selection", class = "sr-only", `aria-label` = "Interactive expression heatmap with row selection capabilities for subsetting features based on expression patterns"),
      div(class = "sr-only",
        tags$h4("About Dynamic Heatmap"),
        tags$p("The dynamic heatmap allows interactive row selection to subset features based on expression patterns. Click individual rows or drag to select multiple features showing interesting patterns. Selected features are automatically updated in the analysis panel for downstream investigation."),
        tags$h4("When to use this view"),
        tags$p("Use the dynamic heatmap when you want to identify and select features with specific expression patterns, such as genes upregulated in certain samples, co-expressed gene clusters, or features showing gradual changes across conditions. This enables pattern-based feature selection for detailed analysis."),
        tags$h4("How to interact"),
        tags$p("Click rows to select individual features. Shift-click or drag to select multiple rows. Selected features are highlighted and become the active feature set for all downstream analyses. Use this to drill down into specific expression patterns.")
      ),
      tags$h3("Heatmap Controls", class = "sr-only"),
      fluidRow(
        column(
          6,
          dropdown(
            inputId = "mydropdown3",
            label = "Controls",
            circle = FALSE, status = "default", icon = icon("cog"),
            width = 700,
            tooltip = tooltipOptions(title = "Click to update heatmap and check legend!"),
            margin = "10px",
            tabsetPanel(
              tabPanel("Parameters", iheatmapInput(id = ns("dynheatmapViewer"))),
              tabPanel("Legend", iheatmapLegend(id = ns("dynheatmapViewer")))
            )
          )
        ),
        column(
          6, align = "right",
          iheatmapClear(id = ns("dynheatmapViewer"))
        ),
        tags$h3("Heatmap Visualization", class = "sr-only"),
        column(
          12, iheatmapOutput(id = ns("dynheatmapViewer"))
        )
      )
    ),
    tabPanel("Expression",
      tags$h2("Expression Matrix Table", class = "sr-only", `aria-label` = "Searchable table showing quantitative expression values for all features across all samples"),
      div(class = "sr-only",
        tags$h4("About Expression Matrix Table"),
        tags$p("The expression matrix table displays raw quantitative values for all features across all samples. Each row represents a feature (gene/protein) and each column represents a sample. This is the numerical data underlying all visualizations and analyses in the application."),
        tags$h4("When to use this view"),
        tags$p("Use the expression matrix when you need to examine specific numerical values, verify data for particular features or samples, export raw data for external analysis, or validate visualization results. This provides direct access to your quantitative measurements."),
        tags$h4("How to navigate"),
        tags$p("Search for specific features or samples. Sort columns to find highest/lowest values. Export functionality allows you to download the data matrix. Values shown are as imported, with any transformations or normalizations applied.")
      ),
      dataTable_ui(ns("tab_expr"))),
    tabPanel("GSList",
      tags$h2("Gene Set Membership Table", class = "sr-only", `aria-label` = "Table showing which features belong to which gene sets or functional categories for pathway enrichment analysis"),
      div(class = "sr-only",
        tags$h4("About Gene Set Membership Table"),
        tags$p("The gene set list shows which features belong to which functional categories, pathways, or custom gene sets. This table is used by enrichment analysis tools (ORA, fGSEA) to test for over-representation or enrichment. Each row links a feature to its gene set memberships."),
        tags$h4("When to use this view"),
        tags$p("Use the gene set list to understand which pathways or functional categories your features are annotated to, verify gene set assignments before enrichment analysis, or explore the breadth of functional annotations available in your dataset."),
        tags$h4("How to interpret"),
        tags$p("Each row shows a feature and its associated gene sets. Features may belong to multiple gene sets. The format is typically 'GeneSet|Database|Term'. This annotation enables pathway enrichment analysis to identify biological themes in your data.")
      ),
      gslist_ui(ns("gsList")))
  )
}

#' @description Application level 1 module - data space
#' @param id module id
#' @param expr reactive value; expression matrix
#' @param pdata reactive value; phenotype data
#' @param fdata reactive value; feature data
#' @param cormat reactive value; correlation matrix. if not given, calculated on the fly.
#' @param reactive_x_s pre-selected x axis for sample space
#' @param reactive_y_s pre-selected y axis for sample space
#' @param reactive_x_f pre-selected x axis for feature space
#' @param reactive_y_f pre-selected y axis for feature space
#' @param status intial status

L1_data_space_module <- function(
  id, expr, pdata, fdata,
  reactive_x_s = reactive(NULL), reactive_y_s = reactive(NULL),
  reactive_x_f = reactive(NULL), reactive_y_f = reactive(NULL),
  cormat = reactive(NULL), status = reactive(NULL)
) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns
  
  # ====== calculate correlation matrix and related component =====

  cmat <- reactive({
    req(expr())
    if (ncol(expr()) <= 3)
      return(NULL)
    if (is.null(cormat())) {      
      cc <- cor(expr(), use = "pairwise.complete.obs")
      diag(cc) <- NA
      hcl <- hclust(as.dist(1-cor(t(cc), use= "pair")), method = "ward.D")
      dend <- as.dendrogram(hcl)
      dl <- list( pearson_ward.D = list(ord = hcl$order, hcl = dend) )
      attr(cc, "rowDendrogram") <- dl
      attr(cc, "colDendrogram") <- dl
      return(cc)
    } else if (nrow(cormat()) == ncol(cormat()) && nrow(cormat()) == ncol(expr())) {
      return(cormat())
    } else
      stop("Incorrect dimension of cormat!")    
    })

  s_cor_heatmap <- iheatmapModule(
    'corheatmapViewer', mat = cmat, pd = pdata, fd = pdata,
    status = reactive(status()$eset_cor_heatmap), fill.NA = FALSE
  )
  
  # # heatmap
  s_heatmap <- iheatmapModule(
    'heatmapViewer', mat = expr, pd = pdata, fd = fdata,
    status = reactive(status()$eset_heatmap)
  )

  # ============ feature space - scatter plot ===========
  # r_feature_fig <- reactiveVal()
  s_feature_fig <- meta_scatter_module(
    "feature_space", reactive_meta = fdata, reactive_expr = expr,
    combine = "feature", source = "scatter_meta_feature", reactive_x = reactive_x_f,
    reactive_y = reactive_y_f, reactive_status = reactive(status()$eset_fdata_fig)
  )

  # ============ sample space - scatter plot ============
  # r_sample_fig <- reactiveVal()# DONOT REMOVE
  s_sample_fig <- meta_scatter_module(
    "sample_space", reactive_meta = pdata, reactive_expr = expr, combine = "pheno", source = "scatter_meta_sample",
    reactive_x = reactive_x_s, reactive_y = reactive_y_s, reactive_status = reactive(status()$eset_pdata_fig)
  )
  
  # ============ level 1 selection - forward to dynamic heatmap ============
  
  tab_rows_fdata <- reactiveVal(TRUE)
  tab_rows_pdata <- reactiveVal(TRUE)
  notNullAndPosLength <- function(x) !is.null(x) && length(x) > 0
  
  ## ===== selection from correlation heatmap - only sample =====

  observeEvent(s_cor_heatmap(), {
    if (notNullAndPosLength(s_cor_heatmap()$brushed$col)) {
      tab_rows_pdata(s_cor_heatmap()$brushed$col)
    } else if (notNullAndPosLength(s_cor_heatmap()$clicked)) {
      tab_rows_pdata(s_cor_heatmap()$clicked[["col"]])
    } else
      tab_rows_pdata(TRUE)    
  })
  
  ## ============== selection from heatmap - sample sand feature ========

  observeEvent(s_heatmap(), {
    # fdata
    if (notNullAndPosLength(s_heatmap()$brushed$row)) {
      tab_rows_fdata(s_heatmap()$brushed$row)
    } else if (notNullAndPosLength(s_heatmap()$clicked)) {
      tab_rows_fdata(s_heatmap()$clicked[["row"]])
    } else {
      tab_rows_fdata(TRUE)
    }
    # pdata
    if (notNullAndPosLength(s_heatmap()$brushed$col)) {
      tab_rows_pdata(s_heatmap()$brushed$col)
    } else if (notNullAndPosLength(s_heatmap()$clicked)) {
      tab_rows_pdata(s_heatmap()$clicked[["col"]])
    } else {
      tab_rows_pdata(TRUE)
    }
  })
  
  ## ============== selection from scatter plots - feature only ========

  observeEvent(c(s_feature_fig()), {
    if (notNullAndPosLength(s_feature_fig()$selected)) {
      tab_rows_fdata( s_feature_fig()$selected )
    } else if ( notNullAndPosLength(s_feature_fig()$clicked) ) {
      tab_rows_fdata( s_feature_fig()$clicked )
    } else
      tab_rows_fdata(TRUE)
  })
  
  ## ============= selection from scatter plots - sample only ========

  observeEvent(s_sample_fig(), {
    if (notNullAndPosLength(s_sample_fig()$selected)) {
      tab_rows_pdata( s_sample_fig()$selected )
    } else if ( notNullAndPosLength(s_sample_fig()$clicked) ) {
      tab_rows_pdata( s_sample_fig()$clicked )
    } else
      tab_rows_pdata(TRUE)
  })
  
  ## ============ selection from feature table - feature only ============
  ## need one more level of nesting to avoid circular dependency
  tab_rows_fdata2 <- reactiveVal(TRUE)
  
  observeEvent(tab_rows_fdata(), {
    tab_rows_fdata2( tab_rows_fdata() )
  })
  
  observeEvent(tab_fd(), {
    if (notNullAndPosLength(tab_fd())) {
      if ( length(tab_fd()) > 2 )
        tab_rows_fdata2( tab_fd() ) else
          tab_rows_fdata2(TRUE)
    } else
      tab_rows_fdata2(TRUE)
  })
  
  ## ============ selection from sample table - sample only ============
  ## need one more level of nesting to avoid circular dependency
  
  tab_rows_pdata2 <- reactiveVal(TRUE)
  
  observeEvent(tab_rows_pdata(), {
    tab_rows_pdata2( tab_rows_pdata() )
  })
  
  observeEvent(tab_pd(), {
    if (notNullAndPosLength(tab_pd())) {
      if ( length(tab_pd()) > 2 )
        tab_rows_pdata2( tab_pd() ) else
          tab_rows_pdata2(TRUE)
    } else
      tab_rows_pdata2(TRUE)
  })

  ## tables
  tab_pd <- dataTable_module(
    "tab_pheno", reactive_data = pdata,
    tab_status = reactive(status()$eset_pdata_tab), tab_rows = tab_rows_pdata
  )
  tab_fd <- dataTable_module(
    "tab_feature", reactive_data = fdata,
    tab_status = reactive(status()$eset_fdata_tab), tab_rows = tab_rows_fdata
  )
  tab_expr <- dataTable_module(
    "tab_expr", reactive_data = reactive({
      req(expr())
      cbind(data.frame(feature = rownames(expr()), expr()))
    }), tab_status = reactive(status()$eset_exprs_tab),
    tab_rows = tab_rows_fdata, selector = FALSE)

  ### return selected feature and samples
  selectedFeatures <- reactiveVal()
  selectedSamples <- reactiveVal()

  observeEvent(s_cor_heatmap(), {
    if (!is.null(s_cor_heatmap()$brushed$col)) {
      selectedSamples(s_cor_heatmap()$brushed$col)
    } else if (!is.null(s_cor_heatmap()$clicked)) {
      selectedSamples(s_cor_heatmap()$clicked["col"])  # ?? else set to character(0)
    } else 
      selectedSamples(character(0))
  })

  observeEvent(s_heatmap(), {

    if (!is.null(s_heatmap()$brushed$row)) {
      selectedFeatures(s_heatmap()$brushed$row)
    } else if (!is.null(s_heatmap()$clicked)) {
      selectedFeatures(s_heatmap()$clicked["row"]) # ?? else set to character(0)
    } else
      selectedFeatures(character(0))

    if (!is.null(s_heatmap()$brushed$col)) {
      selectedSamples(s_heatmap()$brushed$col)
    } else if (!is.null(s_heatmap()$clicked)) {
      selectedSamples(s_heatmap()$clicked["col"]) # ?? else set to character(0)
    } else 
      selectedSamples(character(0))
  })

  observeEvent(s_feature_fig(), {
    if (!is.null(s_feature_fig()$selected) && length(s_feature_fig()$selected) > 0) {
      selectedFeatures( s_feature_fig()$selected )
    } else if ( !is.null(s_feature_fig()$clicked) )
      selectedFeatures( s_feature_fig()$clicked )
  })

  observeEvent(s_sample_fig(), {
    if (!is.null(s_sample_fig()$selected) && length(s_sample_fig()$selected) > 0) {
      selectedSamples( s_sample_fig()$selected )
    } else if ( !is.null(s_sample_fig()$clicked) )
      selectedSamples( s_sample_fig()$clicked )
  })

  observeEvent(tab_pd(), {
    # sso <- selectedSamples()
    # if (length(tab_pd()) < length(sso))
      selectedSamples(tab_pd()) 
    } )
  
  singleTrue <- function(x) !is.null(x) && is.logical(x) && length(x) == 1 && any(x)
  observeEvent(tab_fd(), {
    i1 <- length(tab_fd()) < length(tab_rows_fdata())
    i2 <- singleTrue(tab_rows_fdata())
    i3 <- !singleTrue(tab_fd())
    if (  (i1 || i2) && i3  )
      selectedFeatures( tab_fd() )
    })
  observeEvent(tab_expr(), {
    i1 <- length(tab_expr()) < length(tab_rows_fdata())
    i2 <- singleTrue(tab_rows_fdata())
    i3 <- !singleTrue(tab_expr())
    if (  (i1 || i2) && i3  )
      selectedFeatures(tab_expr())
    })

  # GS List
  tab_gslist <- gslist_module(
    "gsList", reactive_i = tab_rows_fdata, reactive_featureData = fdata
  )

  observeEvent(tab_gslist(), {
    req(tab_gslist())
    selectedFeatures( tab_gslist() )
  })

  # ============= status for snapshot ============
  observe({
    if (!is.null(tb <- status()$eset_active_tab))
      updateNavbarPage(session = session, inputId = "eset", selected = tb)
  })

  na2null <- function(x) {
    if (is.null(x) || is.na(x) || length(x) == 0)
      return(NULL)
    x
  }
  observe({
    tab_rows_fdata( status()$eset_fdata_tabrows )
    tab_rows_pdata( status()$eset_pdata_tabrows )
    selectedSamples( na2null( status()$eset_selected_samples ) )
    selectedFeatures( status()$eset_selected_features )
  })

############## dynamic heatmap function start ##################

  hdmat <- reactive({

    # only run when tab is active
    req(input$eset == "Dynamic heatmap")

    req(e0 <- expr())
    req(fd <- fdata())
    req(pd <- pdata())
    if (!isTRUE(tab_rows_fdata2()))
      req(length(tab_rows_fdata2()) > 2) # at least 3 features
    if (!isTRUE(tab_rows_pdata()))
      req(length(tab_rows_pdata()) > 2) # at least 3 samples

    if (length(tab_rows_fdata2()) > 2) {
      fd <- fd[tab_rows_fdata2(), ]
      e0 <- e0[tab_rows_fdata2(), ]
    }

    if (length(tab_rows_pdata2()) > 2) {
      pd <- pd[tab_rows_pdata2(), ]
      e0 <- e0[, tab_rows_pdata2()]
    }

    list(expr = e0, fd = fd, pd = pd)
    })

  s_dyn_heatmap <- iheatmapModule(
    'dynheatmapViewer',
    mat = reactive(hdmat()$expr), pd = reactive(hdmat()$pd), fd = reactive(hdmat()$fd),
    status = reactive(status()$eset_dyn_heatmap), fill.NA = FALSE
  )

  observeEvent(s_dyn_heatmap(), {
    
    if (!is.null(s_dyn_heatmap()$brushed$row)) {
      selectedFeatures(s_dyn_heatmap()$brushed$row)
    } else if (!is.null(s_dyn_heatmap()$clicked)) {
      selectedFeatures(s_dyn_heatmap()$clicked["row"]) 
    } # ?? else set to character(0)
    #else
      # selectedFeatures(character(0))

    if (!is.null(s_dyn_heatmap()$brushed$col)) {
      selectedSamples(s_dyn_heatmap()$brushed$col)
    } else if (!is.null(s_dyn_heatmap()$clicked)) {
      selectedSamples(s_dyn_heatmap()$clicked["col"]) 
    } 
    # ?? else set to character(0)
    # else 
      # selectedSamples(character(0))
  })

  ############## dynamic heatmap function end ##################

  reactive({
    
    l <- list(
      feature = selectedFeatures(),
      sample = selectedSamples()
    )

    sta <- list(
      eset_active_tab = input$eset,
      eset_pdata_tab = attr(tab_pd(), "status"), # -> pdata_tab
      eset_fdata_tab = attr(tab_fd(), "status"), # -> fdata_tab
      eset_exprs_tab = attr(tab_expr(), "status"), # -> exprs_tab
      eset_fdata_fig = attr(s_feature_fig(), "status"), # -> fdata_fig
      eset_pdata_fig = attr(s_sample_fig(), "status"), # -> pdata_fig
      eset_heatmap = attr(s_heatmap(), "status"),
      eset_fdata_tabrows = tab_rows_fdata(),
      eset_pdata_tabrows = tab_rows_pdata(),
      eset_selected_samples = c(selectedSamples()),
      eset_selected_features = c(selectedFeatures())
    )
    if (sta$eset_active_tab != "Feature table" )# fig tab
      sta$eset_fdata_tab$rows_selected <- NULL
    if (sta$eset_active_tab != "Sample table" )
      sta$eset_pdata_tab$rows_selected <- NULL
    attr(l, "status") <- sta
    l
  })

  }) # end moduleServer
}

