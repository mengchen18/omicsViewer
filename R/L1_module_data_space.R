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
    tabPanel('Feature', meta_scatter_ui(ns('feature_space'))),
    tabPanel("Feature table", dataTable_ui(ns("tab_feature")))    ,
    tabPanel("Sample", meta_scatter_ui(ns("sample_space"))),
    tabPanel("Sample table", dataTable_ui(ns("tab_pheno"))),
    tabPanel(
      "Cor",
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
        column(
          12, iheatmapOutput(id = ns("corheatmapViewer"))
        )
      )
    ),
    tabPanel(
      "Heatmap",
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
        column(
          12, iheatmapOutput(id = ns("heatmapViewer"))
        )
      )
    ),
    tabPanel(
      "Dynamic heatmap",
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
        column(
          12, iheatmapOutput(id = ns("dynheatmapViewer"))
        )
      )
    ),
    tabPanel("Expression", dataTable_ui(ns("tab_expr"))),
    tabPanel("GSList", gslist_ui(ns("gsList")))
  )
}

#' @description Application level 1 module - data space
#' @param input input
#' @param output output
#' @param session session
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
  input, output, session, expr, pdata, fdata, 
  reactive_x_s = reactive(NULL), reactive_y_s = reactive(NULL),
  reactive_x_f = reactive(NULL), reactive_y_f = reactive(NULL),
  cormat = reactive(NULL), status = reactive(NULL)
) {
  
  ns <- session$ns

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

  s_cor_heatmap <- callModule(
    iheatmapModule, 'corheatmapViewer', mat = cmat, pd = pdata, fd = pdata,
    status = reactive(status()$eset_cor_heatmap), fill.NA = FALSE
  )
  
  # # heatmap
  s_heatmap <- callModule(
    iheatmapModule, 'heatmapViewer', mat = expr, pd = pdata, fd = fdata, 
    status = reactive(status()$eset_heatmap)
  )

  # ### feature space
  # r_feature_fig <- reactiveVal()
  s_feature_fig <- callModule(
    meta_scatter_module, id = "feature_space", reactive_meta = fdata, reactive_expr = expr,
    combine = "feature", source = "scatter_meta_feature", reactive_x = reactive_x_f,
    reactive_y = reactive_y_f, reactive_status = reactive(status()$eset_fdata_fig)
  )

  # sample space
  # r_sample_fig <- reactiveVal()# DONOT REMOVE
  s_sample_fig <- callModule(
    meta_scatter_module, id = "sample_space", reactive_meta = pdata, reactive_expr = expr, combine = "pheno", source = "scatter_meta_sample",
    reactive_x = reactive_x_s, reactive_y = reactive_y_s, reactive_status = reactive(status()$eset_pdata_fig)
  )

  tab_rows_fdata <- reactiveVal(TRUE)
  tab_rows_pdata <- reactiveVal(TRUE)
  notNullAndPosLength <- function(x) !is.null(x) && length(x) > 0

  observeEvent(s_cor_heatmap(), {
    if (notNullAndPosLength(s_cor_heatmap()$brushed$col)) {
      tab_rows_pdata(s_cor_heatmap()$brushed$col)
    } else if (notNullAndPosLength(s_cor_heatmap()$clicked)) {
      tab_rows_pdata(s_cor_heatmap()$clicked[["col"]])
    } else
      tab_rows_pdata(TRUE)    
  })

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

  observeEvent(c(s_feature_fig()), {
    if (notNullAndPosLength(s_feature_fig()$selected)) {
      tab_rows_fdata( s_feature_fig()$selected )
    } else if ( notNullAndPosLength(s_feature_fig()$clicked) ) {
      tab_rows_fdata( s_feature_fig()$clicked )
    } else
      tab_rows_fdata(TRUE)
  })

  observeEvent(s_sample_fig(), {
    if (notNullAndPosLength(s_sample_fig()$selected)) {
      tab_rows_pdata( s_sample_fig()$selected )
    } else if ( notNullAndPosLength(s_sample_fig()$clicked) ) {
      tab_rows_pdata( s_sample_fig()$clicked )
    } else
      tab_rows_pdata(TRUE)
  })

  ## tables
  tab_pd <- callModule(
    dataTable_module, id = "tab_pheno",  reactive_data = pdata,
    tab_status = reactive(status()$eset_pdata_tab), tab_rows = tab_rows_pdata
  )
  tab_fd <- callModule(
    dataTable_module, id = "tab_feature",  reactive_data = fdata,
    tab_status = reactive(status()$eset_fdata_tab), tab_rows = tab_rows_fdata
  )
  tab_expr <- callModule(
    dataTable_module, id = "tab_expr",  reactive_data = reactive({
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
  tab_gslist <- callModule(
    gslist_module, id = "gsList", reactive_i = tab_rows_fdata, reactive_featureData = fdata
  )

  observeEvent(tab_gslist(), {
    req(tab_gslist())
    selectedFeatures( tab_gslist() )
  })

  #### status for snapshot #####
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

    if (length(tab_rows_fdata()) > 2) {
      fd <- fd[tab_rows_fdata(), ]
      e0 <- e0[tab_rows_fdata(), ]
    }

    if (length(tab_rows_pdata()) > 2) {
      pd <- pd[tab_rows_pdata(), ]
      e0 <- e0[, tab_rows_pdata()]
    }

    list(expr = e0, fd = fd, pd = pd)
    })

  s_dyn_heatmap <- callModule(
    iheatmapModule, 'dynheatmapViewer', 
    mat = reactive(hdmat()$expr), pd = reactive(hdmat()$pd), fd = reactive(hdmat()$fd),
    status = reactive(status()$eset_dyn_heatmap), fill.NA = FALSE
  )

  observeEvent(s_dyn_heatmap(), {

    if (!is.null(s_dyn_heatmap()$brushed$row)) {
      selectedFeatures(s_dyn_heatmap()$brushed$row)
    } else if (!is.null(s_dyn_heatmap()$clicked)) {
      selectedFeatures(s_dyn_heatmap()$clicked["row"]) # ?? else set to character(0)
    } else
      selectedFeatures(character(0))

    if (!is.null(s_dyn_heatmap()$brushed$col)) {
      selectedSamples(s_dyn_heatmap()$brushed$col)
    } else if (!is.null(s_dyn_heatmap()$clicked)) {
      selectedSamples(s_dyn_heatmap()$clicked["col"]) # ?? else set to character(0)
    } else 
      selectedSamples(character(0))
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
}

