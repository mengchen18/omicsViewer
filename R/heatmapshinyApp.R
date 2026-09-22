#' @description interactive heatmap
#' @param x a matrix object or \code{ExpressionSet} or \code{SummarizedExperiment}
#' @param fData feature data, ignored if x is an ExpressionSet
#' @param pData phenotype data, ignored if x is an ExpressionSet
#' @param impute whether impute the expression matrix
#' @importFrom matrixStats rowSums2 rowVars
#' @importFrom Biobase fData pData exprs
#' 
iheatmap <- function(x, fData = NULL, pData = NULL, impute = FALSE) {
  
  if (inherits(x, "ExpressionSet") || inherits(x, "xcmsFeatureSet")) {
    fData <- fData(x)
    pData <- pData(x)
    x <- exprs(x)
  }
  
  ir <- unique(
    c(which(rowSums2(!is.na(x)) == 0), 
      which(rowVars(x) == 0))
  )
  if (length(ir) > 0) {
    fData <- fData[-ir, ]
    x <- x[-ir, ]
  }
  
  if (impute) {
    x <- apply(x, 1, function(xx) {
      xx[is.na(xx)] <- min(xx, na.rm = TRUE)*0.9
      xx
    })
    x <- t(x)
  }
  
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel = sidebarPanel(
        tabsetPanel(
          tabPanel("Parameters", iheatmapInput(id = "test")),
          tabPanel("Legend", iheatmapLegend(id = "test"))
        )
      ),
      mainPanel = mainPanel(
        iheatmapOutput(id = "test")
      )
    )
  )
  
  server <- function(input, output) {
    iheatmapModule(
      'test',
      mat = reactive(x),
      pd = reactive(pData),
      fd = reactive(fData))
  }
  
  shinyApp(ui, server)
}



#' @description iheatmap input
#' @param id id
#' @param scaleOn the default scale, row, column or none
#' @describeIn iheatmap input for heatmap
#' @importFrom shinycssloaders withSpinner
iheatmapInput <- function(id, scaleOn = "row") {
  
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        6, 
        # clustering column/distance/linkage
        selectInput(ns("colSortBy"), "Sorting columns by", choices = NULL, selectize = TRUE),
        conditionalPanel(
          sprintf("input['%s'] == 'hierarchical cluster'", ns("colSortBy")), 
          selectInput(ns("clusterColDist"), "Distance", 
                      choices = c("Pearson correlation", "Euclidean", "Maximum", "Manhattan", "Canberra", "Binary", 
                                  "Minkowski", "Spearman correlation"), 
                      selectize = TRUE),
          selectInput(ns("clusterColLink"), "Linkage", 
                      choices = c("ward.D", "ward.D2", "single", "complete", "average", 
                                  "mcquitty", "median", "centroid"), 
                      selectize = TRUE)
        ),
        hr(),
        # clustering row/distance/linkage
        selectInput(ns("rowSortBy"), "Sorting rows by", choices = NULL, selectize = TRUE, multiple = FALSE),
        conditionalPanel(
          sprintf("input['%s'] == 'hierarchical cluster'", ns("rowSortBy")), 
          selectInput(ns("clusterRowDist"), "Distance", 
                      choices = c("Pearson correlation", "Euclidean", "Maximum", "Manhattan", "Canberra", "Binary", 
                                  "Minkowski", "Spearman correlation"), 
                      selectize = TRUE),
          selectInput(ns("clusterRowLink"), "Linkage", 
                      choices = c("ward.D", "ward.D2", "single", "complete", "average", 
                                  "mcquitty", "median", "centroid"), 
                      selectize = TRUE)
        ),
        hr(),        
        # scale none/row/col
        selectInput(ns("scale"), "Scale on", 
                    choices = c("row", "none", "column"), selected = scaleOn,
                    selectize = TRUE)
      ),
      column(
        6, 
        # row annotations (1 vs 1+; numeric vs categorith)
        selectInput(ns("annotCol"), label = "Column annotations", choices = NULL, multiple = TRUE),
        
        # column annotations
        selectInput(ns("annotRow"), label = "Row annotations", choices = NULL, multiple = TRUE),
        
        # column annotations
        selectInput(ns("tooltipInfo"), label = "Tooltips", choices = NULL, multiple = TRUE),
        hr(),
        
        # margin
        sliderInput(ns("marginBottom"), "Bottom margin", min = 1, max = 20, value = 4),
        sliderInput(ns("marginRight"), "Right margin", min = 1, max = 20, value = 4),
        
        # color pallete
        selectInput(ns("heatmapColors"), label = "Heatmap color panel", choices = c(
          "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn"
        ), selected = "RdYlBu", selectize = TRUE, multiple = FALSE)   
      )
    )
  )
}

#' @description iheatmap Output
#' @describeIn iheatmap
#' 
iheatmapOutput <- function(id) {
  ns <- NS(id)   
  tagList(    
    fluidRow(
      uiOutput(ns("topleft")),
      uiOutput(ns("topmiddle")),
      uiOutput(ns("column_dend_ui")),
      uiOutput(ns("middleleft")),
      uiOutput(ns("middlemiddle")),
      uiOutput(ns("column_sideColor_ui")),
      uiOutput(ns("row_dend_ui")),
      uiOutput(ns("row_sideColor_ui")),
      shinycssloaders::withSpinner(
        uiOutput(ns("heatmap_ui")), type = 8, color = "green"
      )
    ),
    shinyPlotTooltipsUI(ns("tlp_heatmap"))
  )
}

#' @description iheatmap legend
#' @describeIn iheatmap
#' 
iheatmapLegend <- function(id) {
  ns <- NS(id) 
  tagList(
    # verbatimTextOutput(ns("info")),
    uiOutput(ns("key_heatmap_ui")),
    uiOutput(ns("key_colSideCor_ui")),
    uiOutput(ns("key_rowSideCor_ui"))
  )
}

#' @description iheatmap clear function
#' @describeIn iheatmap iheatmap clear
#'
iheatmapClear <- function(id) {
  ns <- NS(id)
  # actionButton(ns("clear"), "Clear selection")
  actionBttn(ns("clear"), "Clear figure selection", style = "minimal", color = "primary", size = "xs")
}


#' @description iheatmap module
#' @param id module id
#' @param mat expression matrix
#' @param pd phenotype data
#' @param fd feature data
#' @param fill.NA fill NA? TRUE or FALSE
#' @param status heatmap states
#' @param store Optional child view of the canonical widget store
#'   (\code{\link{widget_store_child}}) for this heatmap instance. When
#'   given, every parameter-panel widget (palette, scaling, margins,
#'   sorting, clustering, annotations, tooltips) is registered as an
#'   agent-controllable binding (control plane, plan section 6, S3+S4).
#' @importFrom RColorBrewer brewer.pal
#' @name iheatmap
#'
iheatmapModule <- function(
  id, mat, pd, fd, status = reactive(NULL), fill.NA = TRUE,
  store = NULL
  ) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns
  
  matr <- reactive({
    req(mat())
    if (any(is.na(mat())) && fill.NA )
      r <- fillNA(mat()) else
        r <- mat()
    r
    })

  rowDendrogram <- reactive({
    req(mat())
    attr(mat(), "rowDendrogram")
    })

  colDendrogram <- reactive({
    req(mat())
    attr(mat(), "colDendrogram")
    })

  clsRow <- reactive({
    if (nrow(matr()) > 750 || !is.null(rowDendrogram()))
      return("none")
    "hierarchical cluster"
  })

  rdg <- reactive({
    if (is.null(rowDendrogram()) || is.null(names(rowDendrogram())))
      return(NULL)
    x <- rowDendrogram()
    names(x) <- paste("HCL", names(x))
    x
    })

  cdg <- reactive({
    if (is.null(colDendrogram()) || is.null(names(colDendrogram())))
      return(NULL)
    x <- colDendrogram()
    names(x) <- paste("HCL", names(x))
    x
    })
  
  gsdf <- reactive({
    attr(fd(), "GS")
  })
  
  fdColWithGS <- reactive({
    req(fd())
    c(colnames(fd()), paste0("GS|", levels(gsdf()$gsId)))
  })
  # ######## update selectize input ########
  observe({ 
    req(pd())
    updateSelectInput(session, "annotCol", choices = colnames(pd())) 
    })
  observe({
    req(pd())
    updateSelectizeInput(session, "annotRow", choices = fdColWithGS(), server = TRUE)
    }) 
  observe({
    req(pd())
    updateSelectInput(
      session, "colSortBy", choices = c(names(cdg()), "hierarchical cluster", "none", colnames(pd()))) 
    })
  observe({
    cs <- c(names(rdg()), "none", "hierarchical cluster", fdColWithGS())
    ss <- clsRow()
    if (ss == "none" && cs[[1]] != "none")
      ss <- cs[[1]]
    updateSelectizeInput( session, "rowSortBy", choices = cs, selected = ss, server = TRUE )
    })
  observe({
    req(pd())
    updateSelectInput(session, "tooltipInfo", choices = c(colnames(fd()), colnames(pd())) ) 
    })
  
  # ######## prepare heatmap data ########
  mm <- reactive({
    req(input$scale)
    req(matr())
    if (input$scale == "row") {
      mm <- t(scale(t(matr()))) 
      brk <- c(min(mm, na.rm = TRUE), seq(-2, 2, length.out = 99), max(mm, na.rm = TRUE))
    } else if (input$scale == "column") {
      mm <- scale(matr())
      brk <- c(min(mm, na.rm = TRUE), seq(-2, 2, length.out = 99), max(mm, na.rm = TRUE))
    } else {
      mm <- matr()
      brk <- seq(min(mm, na.rm = TRUE), max(mm, na.rm = TRUE), length.out = 101)
    }    
    list(mat = mm, breaks = brk)
  })
  
  pre_hcl <- reactiveVal()
  pre_ord <- reactiveVal()
  pre_hcl_col <- reactiveVal()
  pre_ord_col <- reactiveVal()

  # ------------------------------------------------------------------
  # Canonical widget-store bindings (control plane, plan section 6, S3+S4).
  # The generic agent tier drives these keys through store_apply; user
  # edits sync back below. S4 completes the parameter panel: sorting,
  # clustering, and annotation widgets join the S3 palette/scale/margins.
  # The status observer below still restores legacy snapshot state
  # directly (kept for pre-S3 snapshots; both paths are idempotent).
  # ------------------------------------------------------------------
  if (!is.null(store)) {
    .heatmap_palettes <- c(
      "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn")
    .heatmap_distances <- c(
      "Pearson correlation", "Euclidean", "Maximum", "Manhattan",
      "Canberra", "Binary", "Minkowski", "Spearman correlation")
    .heatmap_linkages <- c(
      "ward.D", "ward.D2", "single", "complete", "average",
      "mcquitty", "median", "centroid")
    # Defensive choice sources for the store validators: they read module
    # reactives WITHOUT req(), because a shiny.validation condition raised
    # inside a choices_provider would abort the whole store_apply
    # transaction (tryCatch(error=) cannot catch it). They mirror exactly
    # what the choices-update observers offer in the UI.
    .heatmap_hcl_names <- function(which) {
      if (is.null(mat())) return(character(0))
      dl <- attr(mat(), which)
      if (is.null(dl) || is.null(names(dl))) return(character(0))
      paste("HCL", names(dl))
    }
    .heatmap_fd_cols <- function() {
      if (is.null(fd())) return(character(0))
      gs <- attr(fd(), "GS")
      gsn <- if (is.null(gs)) character(0) else paste0("GS|", levels(gs$gsId))
      c(colnames(fd()), gsn)
    }
    .heatmap_col_choices <- function()
      c(.heatmap_hcl_names("colDendrogram"), "hierarchical cluster", "none",
        colnames(pd()))
    .heatmap_row_choices <- function()
      c(.heatmap_hcl_names("rowDendrogram"), "none", "hierarchical cluster",
        .heatmap_fd_cols())
    store_register(
      store,
      widget_binding("heatmap_colors", "select", label = "Heatmap color panel",
        help = paste("Diverging RColorBrewer palette for the heatmap cells;",
                     "one of", paste(.heatmap_palettes, collapse = ", ")),
        values = .heatmap_palettes),
      widget_binding("scale", "select", label = "Scale on",
        help = "Center and scale values by row, by column, or not at all",
        values = c("row", "none", "column")),
      widget_binding("margin_bottom", "integer", label = "Bottom margin",
        help = "Bottom plot margin in lines, from 1 through 20",
        min = 1L, max = 20L),
      widget_binding("margin_right", "integer", label = "Right margin",
        help = "Right plot margin in lines, from 1 through 20",
        min = 1L, max = 20L),
      widget_binding("col_sort_by", "select", label = "Sort columns by",
        help = paste("Column ordering: a precomputed dendrogram (HCL entries),",
                     "'hierarchical cluster', 'none', or a phenotype",
                     "annotation column"),
        choices_provider = function(v) .heatmap_col_choices()),
      widget_binding("cluster_col_dist", "select", label = "Column distance",
        help = "Distance measure for column hierarchical clustering",
        values = .heatmap_distances),
      widget_binding("cluster_col_link", "select", label = "Column linkage",
        help = "Linkage method for column hierarchical clustering",
        values = .heatmap_linkages),
      widget_binding("row_sort_by", "select", label = "Sort rows by",
        help = paste("Row ordering: a precomputed dendrogram (HCL entries),",
                     "'none', 'hierarchical cluster', a feature annotation",
                     "column, or a gene set (GS| prefixed)"),
        choices_provider = function(v) .heatmap_row_choices()),
      widget_binding("cluster_row_dist", "select", label = "Row distance",
        help = "Distance measure for row hierarchical clustering",
        values = .heatmap_distances),
      widget_binding("cluster_row_link", "select", label = "Row linkage",
        help = "Linkage method for row hierarchical clustering",
        values = .heatmap_linkages),
      widget_binding("annot_col", "multi_select", label = "Column annotations",
        help = paste("Phenotype annotation columns drawn as color bars",
                     "above the heatmap; empty clears the bars"),
        choices_provider = function(v) colnames(pd())),
      widget_binding("annot_row", "multi_select", label = "Row annotations",
        help = paste("Feature annotation columns (or GS| gene sets) drawn as",
                     "color bars beside the heatmap; empty clears the bars"),
        choices_provider = function(v) .heatmap_fd_cols()),
      widget_binding("tooltip_info", "multi_select", label = "Tooltips",
        help = paste("Feature and phenotype columns shown in the hover",
                     "tooltip; empty clears the tooltip"),
        choices_provider = function(v) c(colnames(fd()), colnames(pd())))
    )
    .heatmap_store_root <- if (is.null(store$parent)) store else store$parent
    .heatmap_store_keys <- c(
      heatmap_colors = "heatmapColors", scale = "scale",
      margin_bottom = "marginBottom", margin_right = "marginRight",
      col_sort_by = "colSortBy", row_sort_by = "rowSortBy",
      cluster_col_dist = "clusterColDist", cluster_col_link = "clusterColLink",
      cluster_row_dist = "clusterRowDist", cluster_row_link = "clusterRowLink",
      annot_col = "annotCol", annot_row = "annotRow",
      tooltip_info = "tooltipInfo")
    # Create the epoch reactive ONCE and keep strong references to every
    # store-glue observer: observers whose dependencies are only weakly
    # held by the reactive graph are otherwise garbage collected between
    # flushes, silently killing later store -> UI pushes.
    .heatmap_epoch <- store_epoch(store)
    .heatmap_store_observers <- list()
    .heatmap_keep <- function(obs) {
      .heatmap_store_observers[[length(.heatmap_store_observers) + 1L]] <<- obs
      invisible(obs)
    }

    # UI -> store: mirror user edits (store_sync_from_ui is acknowledgement
    # aware, so widget confirmations of external writes are not overrides).
    .heatmap_keep(observeEvent(input$heatmapColors, {
      store_sync_from_ui(store, "heatmap_colors", input$heatmapColors)
    }, ignoreInit = TRUE))
    .heatmap_keep(observeEvent(input$scale, {
      store_sync_from_ui(store, "scale", input$scale)
    }, ignoreInit = TRUE))
    .heatmap_keep(observeEvent(input$marginBottom, {
      store_sync_from_ui(store, "margin_bottom", input$marginBottom)
    }, ignoreInit = TRUE))
    .heatmap_keep(observeEvent(input$marginRight, {
      store_sync_from_ui(store, "margin_right", input$marginRight)
    }, ignoreInit = TRUE))
    .heatmap_keep(observeEvent(input$colSortBy, {
      store_sync_from_ui(store, "col_sort_by", input$colSortBy)
    }, ignoreInit = TRUE))
    .heatmap_keep(observeEvent(input$rowSortBy, {
      store_sync_from_ui(store, "row_sort_by", input$rowSortBy)
    }, ignoreInit = TRUE))
    .heatmap_keep(observeEvent(input$clusterColDist, {
      store_sync_from_ui(store, "cluster_col_dist", input$clusterColDist)
    }, ignoreInit = TRUE))
    .heatmap_keep(observeEvent(input$clusterColLink, {
      store_sync_from_ui(store, "cluster_col_link", input$clusterColLink)
    }, ignoreInit = TRUE))
    .heatmap_keep(observeEvent(input$clusterRowDist, {
      store_sync_from_ui(store, "cluster_row_dist", input$clusterRowDist)
    }, ignoreInit = TRUE))
    .heatmap_keep(observeEvent(input$clusterRowLink, {
      store_sync_from_ui(store, "cluster_row_link", input$clusterRowLink)
    }, ignoreInit = TRUE))
    # multi-select inputs report character(0) when cleared; observeEvent
    # treats that as a change (only NULL is ignored), so clearing syncs too
    .heatmap_keep(observeEvent(input$annotCol, {
      store_sync_from_ui(store, "annot_col", input$annotCol)
    }, ignoreInit = TRUE))
    .heatmap_keep(observeEvent(input$annotRow, {
      store_sync_from_ui(store, "annot_row", input$annotRow)
    }, ignoreInit = TRUE))
    .heatmap_keep(observeEvent(input$tooltipInfo, {
      store_sync_from_ui(store, "tooltip_info", input$tooltipInfo)
    }, ignoreInit = TRUE))

    # Seed unset keys with the widget defaults once the inputs exist, so
    # discovery tools report real current values from the start. Restores
    # and agent applies that land first win (seeding skips held keys).
    .heatmap_seeded <- FALSE
    .heatmap_keep(observe({
      if (.heatmap_seeded) return(NULL)
      vals <- list(
        heatmap_colors = input$heatmapColors, scale = input$scale,
        margin_bottom = input$marginBottom, margin_right = input$marginRight,
        col_sort_by = input$colSortBy, row_sort_by = input$rowSortBy,
        cluster_col_dist = input$clusterColDist,
        cluster_col_link = input$clusterColLink,
        cluster_row_dist = input$clusterRowDist,
        cluster_row_link = input$clusterRowLink,
        annot_col = input$annotCol %||% character(0),
        annot_row = input$annotRow %||% character(0),
        tooltip_info = input$tooltipInfo %||% character(0))
      if (any(vapply(vals, is.null, logical(1))))
        return(NULL)
      .heatmap_seeded <<- TRUE
      # store_read on a child view keys results by FULL canonical ids
      held <- store_read(store, names(vals))
      patch <- vals[vapply(names(vals), function(k)
        is.null(held[[paste0(store$prefix, ".", k)]]), logical(1))]
      if (length(patch))
        tryCatch(store_apply(store, patch, origin = "system", strict = FALSE),
                 error = function(e) NULL)
    }))

    # Store -> UI push for external writes only (pending entries mark
    # them); user clicks sync back through the observers above.
    # Shape mirrors meta_scatter's proven axis-mode push: epoch watch +
    # return()-based guards (no for/next control flow inside observe()).
    # rowSortBy/annotRow are server-side selectize inputs; their update
    # goes through updateSelectizeInput exactly like the status path.
    .heatmap_push_input <- function(key, input_id, value) {
      if (startsWith(key, "margin_")) {
        updateSliderInput(session, input_id, value = value)
      } else if (input_id %in% c("rowSortBy", "annotRow")) {
        updateSelectizeInput(session, input_id, selected = value)
      } else {
        updateSelectInput(session, input_id, selected = value)
      }
    }
    .heatmap_keep(observe({
      .heatmap_epoch()
      vals <- store_read(store, names(.heatmap_store_keys))
      pending <- .heatmap_store_root$pending
      invisible(lapply(names(.heatmap_store_keys), function(key) {
        full <- paste0(store$prefix, ".", key)
        value <- vals[[full]]
        if (is.null(value) || is.null(pending[[full]]))
          return(NULL)
        .heatmap_push_input(key, .heatmap_store_keys[[key]], value)
      }))
    }))
  }

  # Status restore: translated into one store transaction when the module
  # is store-backed (values identical to the widget_store snapshot apply
  # idempotently; diff-only writes leave untouched keys alone). The direct
  # update*Input path is kept for store-less callers (standalone heatmap
  # app) and covers values the transactional validators reject per key
  # (strict = FALSE), e.g. annotation columns from a changed dataset.
  .heatmap_status_patch <- function(s) {
    patch <- list(
      heatmap_colors = s$heatmapColors,
      scale = s$scale,
      margin_bottom = s$marginBottom,
      margin_right = s$marginRight,
      col_sort_by = s$colSortBy,
      row_sort_by = s$rowSortBy,
      cluster_col_dist = s$clusterColDist,
      cluster_col_link = s$clusterColLink,
      cluster_row_dist = s$clusterRowDist,
      cluster_row_link = s$clusterRowLink,
      annot_col = null2empty(s$annotCol),
      annot_row = null2empty(s$annotRow),
      tooltip_info = null2empty(s$tooltipInfo))
    # drop absent scalar keys; multi_select keys may legitimately be empty
    # (cleared selections restore as character(0))
    keep <- vapply(names(patch), function(k) {
      if (k %in% c("annot_col", "annot_row", "tooltip_info")) TRUE
      else !is.null(patch[[k]]) && length(patch[[k]]) > 0L
    }, logical(1))
    patch[keep]
  }
  observeEvent(status(), {
    if (is.null(status()))
      return(NULL)
    if (!is.null(store)) {
      tryCatch(
        store_apply(store, .heatmap_status_patch(status()),
                    origin = "restore", strict = FALSE),
        error = function(e) NULL)
    }
    updateSelectInput(session, "annotCol", selected = null2empty(status()$annotCol) )
    updateSelectizeInput(session, "annotRow", selected  = null2empty(status()$annotRow) )
    updateSelectInput(session, "colSortBy", selected  = status()$colSortBy)
    updateSelectizeInput(session, "rowSortBy", selected  = status()$rowSortBy)
    updateSelectInput(session, "tooltipInfo", selected  = null2empty(status()$tooltipInfo) )
    updateSelectInput(session, "heatmapColors", selected  = status()$heatmapColors)
    updateSelectInput(session, "scale", selected  = status()$scale)
    updateSelectInput(session, "clusterColDist", selected  = status()$clusterColDist)
    updateSelectInput(session, "clusterColLink", selected  = status()$clusterColLink)
    updateSelectInput(session, "clusterRowDist", selected  = status()$clusterRowDist)
    updateSelectInput(session, "clusterRowLink", selected  = status()$clusterRowLink)
    updateSliderInput(session, "marginRight", value = status()$marginRight)
    updateSliderInput(session, "marginBottom", value = status()$marginBottom)
    pre_hcl(status()$rowDendrogram)
    pre_ord(status()$rowOrder)
    pre_hcl_col(status()$colDendrogram)
    pre_ord_col(status()$colOrder)
    })

  rowSB <- eventReactive(list(
    input$rowSortBy, mm()$mat, input$clusterRowDist, input$clusterRowLink, status()
    ), {
    if (!is.null(pre_hcl()) & !is.null(pre_ord())) {
      return(list(
        ord = pre_ord(), hcl = pre_hcl()
        ))
    }
    req(input$rowSortBy)    
    hcl_r <- NULL
    ord_r <- seq_len(nrow(mm()$mat))
    if (input$rowSortBy %in% names(rdg())) {
      return(rdg()[[input$rowSortBy]])
    } else if (!is.null(gsdf()) && grepl("GS\\|", input$rowSortBy)[1]) {
      d <- gsdf()$featureId[gsdf()$gsId == sub("^GS\\|", "", input$rowSortBy)]
      d <- as.numeric( rownames(fd()) %fin% as.character(d) )
      ord_r <- order(d)
    } else if (!input$rowSortBy %in% c("", "none", "hierarchical cluster")) {
      req(input$rowSortBy %in% colnames(fd()))
      ord_r <- order(fd()[, input$rowSortBy])
    } else if (input$rowSortBy == "hierarchical cluster" && nrow(mm()$mat) > 2) { #clusterRow
      dd <- tolower(strsplit(input$clusterRowDist, " ")[[1]][1])
      hcl_r <- hclust(adist(mm()$mat, method = dd), method = input$clusterRowLink)
      ord_r <- hcl_r$order 
      hcl_r <- as.dendrogram(hcl_r)
    }
    list(ord = ord_r, hcl = hcl_r)
  })
  
  colSB <- eventReactive(list(    
    input$colSortBy, mm()$mat, input$clusterColDist, input$clusterColLink, status()
    ), {        
    if (!is.null(pre_hcl_col()) & !is.null(pre_ord_col())) {
      return(list(
        ord = pre_ord_col(), hcl = pre_hcl_col()
        ))
    }
    req(input$colSortBy)
    hcl_c <- NULL
    ord_c <- seq_len(ncol(mm()$mat))
    if (input$colSortBy %in% names(cdg())) {
      return(cdg()[[input$colSortBy]])
    } else if (!input$colSortBy %in% c("", "none", "hierarchical cluster")) {
      req(input$colSortBy %in% colnames(pd()))
      ord_c <- order(pd()[, input$colSortBy])
    } else if (input$colSortBy == "hierarchical cluster" && ncol(mm()$mat) > 2) {
      dd <- tolower(strsplit(input$clusterColDist, " ")[[1]][1])
      hcl_c <- hclust(adist(t(mm()$mat), method = dd), method = input$clusterColLink)
      ord_c <- hcl_c$order
      hcl_c <- as.dendrogram(hcl_c)
    }
    list(ord = ord_c, hcl = hcl_c)
  })
  
  hm <- reactive({
    mmo <- mm()$mat[rowSB()$ord, colSB()$ord]
    list(
      dend_c = colSB()$hcl,
      dend_r = rowSB()$hcl,
      mat = t(mmo),
      ord_c = colSB()$ord,
      ord_r = rowSB()$ord,
      brk =  mm()$breaks
    )
  })
  ######## render plot ########
  dat_colSideCol <- reactive({
    req(input$annotCol)    
    addHeatmapAnnotation(pd()[hm()$ord_c, input$annotCol],  var.name = input$annotCol)
  })
  output$colSideCol <- renderPlot({
    par(mar = c(0, 0, 0, input$marginRight))
    addHeatmapAnnotation_plot( dat_colSideCol(), xlim = ranges$x-0.5)
  })
  
  dat_rowSideCol <- reactive({
    req(input$annotRow)
    ic <- grepl("GS\\|", input$annotRow)
    am <- NULL
    if (any(ic)) {
      am <- cbind(am, vapply(input$annotRow[ic], function(i) {
        d <- gsdf()$featureId[gsdf()$gsId == sub("^GS\\|", "", i)]
        rownames(fd()) %fin% as.character(d)
      }, FUN.VALUE = logical(nrow(fd()))))
      colnames(am) <- input$annotRow[ic]
    }
    if (any(!ic)) {
      am2 <- fd()[, input$annotRow[!ic], drop = FALSE]
      if (is.null(am))
        am <- am2 else
          am <- cbind(am, am2)
    }
    am <- as.data.frame(am[, input$annotRow, drop = FALSE])
    addHeatmapAnnotation(am[hm()$ord_r, , drop = FALSE], column = FALSE, var.name = input$annotRow)
  })
  output$rowSideCol <- renderPlot({
    par(mar= c(input$marginBottom, 0, 0, 0))
    addHeatmapAnnotation_plot(dat_rowSideCol(), ylim = ranges$y-0.5)
  })
  
  output$heatmap <- renderPlot({
    par(mar = c(input$marginBottom, 0, 0, input$marginRight))
    req(hm()$mat)
    image(hm()$mat, x = seq_len(nrow(hm()$mat)), y = seq_len(ncol(hm()$mat)), xlim = ranges$x, ylim = ranges$y,
          col = heatColor(), axes = FALSE, xlab = "", ylab = "", breaks = hm()$brk)
    
    if (ranges$x[2] - ranges$x[1] <= 60) {
      irn <- .rg(ranges$x, rownames(hm()$mat))
      mtext(side = 1, at = irn$at, text = irn$lab, las = 2, line = 0.5)
      abline(v = seq(ranges$x[1], ranges$x[2], by = 1), col = "white")
    }
    
    if (ranges$y[2] - ranges$y[1] <= 100) {
      rrn <- .rg(ranges$y, colnames(hm()$mat))
      abline(h = seq(ranges$y[1], ranges$y[2], by = 1), col = "white")
      mtext(side = 4, at = rrn$at, text = rrn$lab, las = 2, line = 0.5)
    }    
  })
  
  output$dendCol <- renderPlot({
    req(hm()$dend_c)
    par(mar = c(0, 0, 1, input$marginRight))
    plot(hm()$dend_c, xaxs="i", yaxs = "i", axes = FALSE, xlim = ranges$x, center = TRUE)
    axis(side = 4)
  })
  
  output$dendRow <- renderPlot({
    req(hm()$dend_r)
    par(mar = c(input$marginBottom, 1, 0, 0))
    plot(hm()$dend_r, horiz = TRUE, yaxs="i", xaxs = "i", axes = FALSE, ylim = ranges$y)
    axis(side = 1)
  })
  
  ######## update range - heatmap ########
  ranges <- reactiveValues(x = NULL, y = NULL)
  observe({
    req(hm()$mat)
    if (!is.null(status()$ranges_x))
      ranges$x <- status()$ranges_x else 
        ranges$x <- c(0, nrow(hm()$mat))+0.5
    if (!is.null(status()$ranges_y))
      ranges$y <- status()$ranges_y else 
        ranges$y <- c(0, ncol(hm()$mat))+0.5
  })
  
  .rg <- function(x, tx) {
    v1 <- ceiling(x[1])
    v2 <- floor(x[2])
    list(at = v1:v2, lab = tx[v1:v2])
  }
  heatColor <- reactive({
    colorRampPalette(
      rev(brewer.pal(n = 7, name = input$heatmapColors))
    )(100)
  })

  observeEvent(input$heatmap_dblclick, {
    brush <- input$heatmap_brush
    if (!is.null(brush)) {
      ranges$x <- c(max(floor(brush$xmin)-0.5, 0.5), min(round(brush$xmax)+0.5, nrow(hm()$mat) + 0.5))
      ranges$y <- c(max(floor(brush$ymin)-0.5, 0.5), min(round(brush$ymax)+0.5, ncol(hm()$mat) + 0.5))
    } else {
      ranges$x <- c(0, nrow(hm()$mat))+0.5
      ranges$y <- c(0, ncol(hm()$mat))+0.5
    }
  })

  ######## update range - dendrogram ########
  observeEvent(input$dendRow_dblclick, {
    brush <- input$dendRow_brush
    if (!is.null(brush)) {
      ranges$y <- c(max(floor(brush$ymin)-0.5, 0.5), min(ceiling(brush$ymax)+0.5, ncol(hm()$mat) + 0.5))
    } else {
      ranges$y <- c(0, ncol(hm()$mat))+0.5
    }
  })
  observeEvent(input$dendCol_dblclick, {
    brush <- input$dendCol_brush
    if (!is.null(brush)) {
      ranges$x <- c(max(floor(brush$xmin)-0.5, 0.5), min(ceiling(brush$xmax)+0.5, nrow(hm()$mat) + 0.5))
    } else {
      ranges$x <- c(0, nrow(hm()$mat))+0.5
    }
  })

  ######## update range - sidebar ########
  observeEvent(input$rowSideCol_dblclick, {
    brush <- input$rowSideCol_brush
    if (!is.null(brush)) {
      ranges$y <- c(max(floor(brush$ymin)-0.5, 0.5), min(ceiling(brush$ymax)+0.5, ncol(hm()$mat) + 0.5))
    } else {
      ranges$y <- c(0, ncol(hm()$mat))+0.5
    }
  })
  observeEvent(input$colSideCol_dblclick, {
    brush <- input$colSideCol_brush
    if (!is.null(brush)) {
      ranges$x <- c(max(floor(brush$xmin)-0.5, 0.5), min(ceiling(brush$xmax)+0.5, nrow(hm()$mat) + 0.5))
    } else {
      ranges$x <- c(0, nrow(hm()$mat))+0.5
    }
  })

  ######### render keys ##########
  output$key_heatmap <- renderPlot(
    heatmapKey(range(hm()$mat, na.rm = TRUE), heatColor())
  )
  output$key_heatmap_ui <- renderUI({
    plotOutput(ns("key_heatmap"), height = "45px")
  })

  output$key_colSideCor <- renderPlot({
    lg <- dat_colSideCol()
    req(lg$key)
    if (lg$type == 2) {
      sideCorKey(x = lg$key, label = lg$var.name)
    } else {
      graphics::layout(matrix(seq_along(lg$key), 1))
      for (i in seq_along(lg$key))
        sideCorKey(x = lg$key[[i]], label = lg$var.name[i])
    }
  })
  output$key_colSideCor_ui <- renderUI({
    if (is.null(dat_colSideCol()$key))
      return()
    tagList(
      h4("Column bars"),
      plotOutput(ns("key_colSideCor"), height = 266)
    )
  })

  output$key_rowSideCor <- renderPlot({
    lg <- dat_rowSideCol()
    req(lg$key)
    if (lg$type == 2) {
      sideCorKey(x = lg$key, label = lg$var.name)
    } else {
      graphics::layout(matrix(seq_along(lg$key), 1))
      for (i in seq_along(lg$key))
        sideCorKey(x = lg$key[[i]], label = lg$var.name[i])
    }
  })
  output$key_rowSideCor_ui <- renderUI({
    if (is.null(dat_rowSideCol()$key))
      return()
    tagList(
      h4("Row bars"),
      plotOutput(ns("key_rowSideCor"), height = 266)
    )
  })

  ######## place holder plots ##########
  output$empty1 <- renderPlot({
    par(mar = c(0, 0, 0, 0))
    plot(0, axes = FALSE, col = NA)
  })
  output$empty2 <- renderPlot({
    par(mar = c(0, 0, 0, 0))
    plot(0, axes = FALSE, col = NA)
  })
  output$empty3 <- renderPlot({
    par(mar = c(0, 0, 0, 0))
    plot(0, axes = FALSE, col = NA)
  })
  output$empty4 <- renderPlot({
    par(mar = c(0, 0, 0, 0))
    plot(0, axes = FALSE, col = NA)
  })

  ######### dynamic UI render; dynamic layout ##########
  show_col_dend <- reactiveVal(TRUE)
  observe( show_col_dend(
    input$colSortBy == "hierarchical cluster" || grepl("^HCL ", input$colSortBy)
    ))

  show_col_sideColor <- reactiveVal(FALSE)
  observe( show_col_sideColor(length(input$annotCol) != 0) )

  show_row_dend <- reactiveVal(TRUE)
  observe( show_row_dend(
      input$rowSortBy == "hierarchical cluster" || grepl("^HCL ", input$rowSortBy)
    ))

  show_row_sideColor <- reactiveVal(FALSE)
  observe( show_row_sideColor(length(input$annotRow) != 0) )

  width_left <- reactive( 2 )
  width_mid <- reactive ({
    if (length(input$annotRow) <= 3)
      r <- 1 else if (length(input$annotRow) < 6)
        r <- 2 else
          r <- 3
  })
  width_right <- reactive({
    if (show_row_sideColor() & show_row_dend()) {
      wid <- 12 - width_mid() - width_left()
    } else if (!show_row_sideColor() & show_row_dend()) {
      wid <- 12 - width_left()
    } else if (show_row_sideColor() & !show_row_dend()) {
      wid <- 12 - width_mid()
    } else
      wid <- 12
    wid
  })

  # top left
  height_colSideCol <- reactive({
    if ( dat_colSideCol()$type == 1 )
      return("50px")
    paste0(20*length(input$annotCol), "px")
  })

  output$topleft <- renderUI({
    if (!show_col_dend() || !show_row_dend())
      return(NULL)
    column(width_left(), plotOutput(ns("empty1"), height = "100px" ), offset = 0, style='padding:0px;')
  })
  # top middle
  output$topmiddle <- renderUI({
    if (!show_col_dend() || !show_row_sideColor())
      return(NULL)
    column(width_mid(), plotOutput(ns("empty2"), height = "100px" ), offset = 0, style='padding:0px;')
  })
  # top right
  output$column_dend_ui <- renderUI({
    if (!show_col_dend())
      return(NULL)
    column(width_right(), plotOutput(ns("dendCol"),
                                     dblclick = ns("dendCol_dblclick"),
                                     brush = brushOpts(id = ns("dendCol_brush"), resetOnNew = TRUE),
                                     height = "100px" ), offset = 0, style='padding:0px;')
  })

  # middle left
  output$middleleft <- renderUI({
    if (!show_col_sideColor() || !show_row_dend() )
      return(NULL)
    column(width_left(), plotOutput(ns("empty3"), height = height_colSideCol() ), offset = 0, style='padding:0px;')
  })
  # middle middle
  output$middlemiddle <- renderUI({
    if (!show_col_sideColor() || !show_row_sideColor() )
      return(NULL)
    column(width_mid(), plotOutput(ns("empty4"), height = height_colSideCol() ), offset = 0, style='padding:0px;')
  })
  # middle right
  output$column_sideColor_ui <- renderUI({
    if ( !show_col_sideColor() )
      return(NULL)
    column(width_right(), plotOutput(ns("colSideCol"),
                                     click = ns("colSideCol_click"),
                                     dblclick = ns("colSideCol_dblclick"),
                                     hover = hoverOpts(id = ns("colSideCol_hover"), delay = 150),
                                     brush = brushOpts(id = ns("colSideCol_brush"), resetOnNew = TRUE),
                                     height = height_colSideCol() ), offset = 0, style='padding:0px;')
  })
  # bottom left
  output$row_dend_ui <- renderUI({
    if ( ! show_row_dend() )
      return(NULL)
    column(width_left(), plotOutput(ns("dendRow"),
                                    dblclick = ns("dendRow_dblclick"),
                                    brush = brushOpts(id = ns("dendRow_brush"), resetOnNew = TRUE),
                                    height = "800px" ), offset = 0, style='padding:0px;')
  })
  # bottom middle
  output$row_sideColor_ui <- renderUI({
    if ( !show_row_sideColor() )
      return(NULL)
    column(width_mid(), plotOutput(ns("rowSideCol"),
                                   click = ns("rowSideCol_click"),
                                   dblclick = ns("rowSideCol_dblclick"),
                                   hover = hoverOpts(id = ns("rowSideCol_hover"), delay = 150),
                                   brush = brushOpts(id = ns("rowSideCol_brush"), resetOnNew = TRUE),
                                   height = "800px" ), offset = 0, style='padding:0px;')
  })
  # bottom right
  output$heatmap_ui <- renderUI({
    column(width_right(),
           plotOutput(ns("heatmap"), click = ns("heatmap_click"),
                      dblclick = ns("heatmap_dblclick"),
                      hover = hoverOpts(id = ns("heatmap_hover"), delay = 150),
                      brush = brushOpts(id = ns("heatmap_brush"), resetOnNew = TRUE),
                      height = "800px"),
           offset = 0, style='padding:0px;')

  })

  ############## tooltips ################

  dat_tooltip <- reactive({
    req(pd())
    req(fd())
    list(
      pd = pd()[hm()$ord_c, colnames(pd()) %in% input$tooltipInfo, drop = FALSE],
      fd = fd()[hm()$ord_r, colnames(fd()) %in% input$tooltipInfo, drop = FALSE]
    )
  })
  tooltip_mouse <- reactiveVal(NULL)
  # col side color tooltips
  observeEvent(list(input$annotCol, input$annotRow), {
    tooltip_mouse(NULL)
  })
  observe({
    res <- NULL
    req(input$colSideCol_hover)
    req(input$colSideCol_click)
    hover_x <- ceiling(input$colSideCol_hover$x)
    hover_y <- ceiling(input$colSideCol_hover$y)
    x <- ceiling(input$colSideCol_click$x)
    y <- ceiling(input$colSideCol_click$y)
    if (hover_x == x && hover_y == y && dat_colSideCol()$type %in% 2:3) {
      if (dat_colSideCol()$type == 2) {
        varname <- dat_colSideCol()$var.name
        key <- dat_colSideCol()$key
        leg <- names(dat_colSideCol()$key$color)[x]
      } else {
        varname <- rownames(dat_colSideCol()$cmat)[y]
        key <- dat_colSideCol()$key[[match(varname, dat_colSideCol()$var.name)]]
        color <- match(dat_colSideCol()$cmat[y, x], dat_colSideCol()$cmat[y, ])
        leg <- names(key$color)[color]
      }
      res <- sprintf("<b>%s:</b> %s", varname, leg)
    }
    tooltip_mouse(res)
  })

  # row side color tooltips
  observe({
    res <- NULL
    req(input$rowSideCol_hover)
    req(input$rowSideCol_click)
    hover_x <- ceiling(input$rowSideCol_hover$x)
    hover_y <- ceiling(input$rowSideCol_hover$y)
    x <- ceiling(input$rowSideCol_click$x)
    y <- ceiling(input$rowSideCol_click$y)
    if (hover_x == x && hover_y == y && dat_rowSideCol()$type %in% 2:3) {
      if (dat_rowSideCol()$type == 2) {
        varname <- dat_rowSideCol()$var.name
        key <- dat_rowSideCol()$key
        leg <- names(dat_rowSideCol()$key$color)[y]
      } else {
        varname <- rownames(dat_rowSideCol()$cmat)[x]
        key <- dat_rowSideCol()$key[[match(varname, dat_rowSideCol()$var.name)]]
        color <- match(dat_rowSideCol()$cmat[x, y], dat_rowSideCol()$cmat[x, ])
        leg <- names(key$color)[color]
      }
      res <- sprintf("<b>%s:</b> %s", varname, leg)
    }
    tooltip_mouse(res)
  })

  clickedName <- reactiveVal(NULL)
  # heatmap tooltips
  observe({
    res <- NULL
    req(input$heatmap_hover)
    req(input$heatmap_click)
    hover_x <- round(input$heatmap_hover$x)
    hover_y <- round(input$heatmap_hover$y)
    x <- round(input$heatmap_click$x)
    y <- round(input$heatmap_click$y)

    if (hover_x == x && hover_y == y ) {

      l <- c(dat_tooltip()$pd[x, , drop = FALSE], dat_tooltip()$fd[y, , drop = FALSE])
      tl <- vapply(names(l), function(x) {
        if (is.numeric(cv <- l[[x]]))
          cv <- round(cv, digits = 2)
        sprintf("<b>%s:</b> %s", x, cv)
      }, character(1))
      tl <- paste(tl, collapse = "<br>")
      res <- sprintf(
        "<b>Sample:</b> %s <br> <b>Feature:</b> %s",
        rownames(hm()$mat)[x], colnames(hm()$mat)[y]
      )
      res <- paste(res, tl, sep = "<br>")
    }

    clickedName( c(col = rownames(hm()$mat)[x], row = colnames(hm()$mat)[y] ) )
    tooltip_mouse(res)
  })

  shinyPlotTooltips("tlp_heatmap", points = tooltip_mouse)

  ############### return ##############

  brushedValues <- reactive({
    brush <- input$heatmap_brush
    l <- list(col = NULL, row = NULL)
    if (!is.null(brush)) {
      x1 <- ceiling(max(round(brush$xmin), 0.5))
      x2 <- floor(min(round(brush$xmax), nrow(hm()$mat)))
      y1 <- ceiling(max(round(brush$ymin), 0.5))
      y2 <- floor(min(round(brush$ymax), ncol(hm()$mat)))
      l <- list(
        row = colnames(hm()$mat)[y1:y2],
        col = rownames(hm()$mat)[x1:x2]
      )
    }
    l
  })

  selVal <- reactiveValues(
    clicked = NULL,
    selected = list(col = NULL, row = NULL)
  )
  observeEvent(input$clear, {
    selVal$clicked <- NULL#
    selVal$selected <- list(col = NULL, row = NULL)
  })
  observeEvent(list(clickedName(), brushedValues()), {
    selVal$clicked <- clickedName()# l[v_scatter()$clicked]
    selVal$selected <- brushedValues()#l[v_scatter()$selected]
  })

  reactive({
    r <- list(
      clicked = selVal$clicked, #clickedName(),
      brushed = selVal$selected) # brushedValues())
    attr(r, "status") <- list(
      annotCol = input$annotCol,
      annotRow = input$annotRow,
      colSortBy = input$colSortBy,
      rowSortBy = input$rowSortBy,
      tooltipInfo = input$tooltipInfo,
      marginRight = input$marginRight,
      marginBottom = input$marginBottom,
      heatmapColors = input$heatmapColors,
      scale = input$scale,
      clusterColDist = input$clusterColDist,
      clusterColLink = input$clusterColLink,
      clusterRowDist = input$clusterRowDist,
      clusterRowLink = input$clusterRowLink,
      rowDendrogram = rowSB()$hcl,
      rowOrder = rowSB()$ord,
      colDendrogram = colSB()$hcl,
      colOrder =  colSB()$ord,
      ranges_x = ranges$x,
      ranges_y = ranges$y
      )
    r
  })

  }) # end moduleServer
}

