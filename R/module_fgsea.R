#' @description Utility - fgsea shiny ui
#' @param id id
enrichment_fgsea_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Module description for AI browsers and screen readers
    div(class = "sr-only", id = ns("module-help"),
      tags$h4("About Fast Gene Set Enrichment Analysis (fGSEA)"),
      tags$p("fGSEA is a computational method that determines whether predefined gene sets show statistically significant, concordant differences between two biological states. Unlike ORA which only uses a cutoff, fGSEA uses the entire ranked list of genes, making it more sensitive to subtle but coordinated changes in gene expression across a pathway."),
      tags$h4("When to use fGSEA"),
      tags$p("Use fGSEA when you have a ranked list of all genes/proteins (e.g., ranked by fold change or t-statistic) and want to identify which gene sets are enriched at the top or bottom of your ranking. This is especially useful when there are no clear cutoffs or when you want to detect pathway-level changes that might not be apparent from individual gene changes."),
      tags$h4("How to interpret results"),
      tags$p("The Normalized Enrichment Score (NES) indicates the degree of enrichment: positive NES means the gene set is enriched at the top of your ranked list, negative NES means enrichment at the bottom. The p-value and FDR indicate statistical significance. The leading edge genes are the core subset of genes driving the enrichment signal. The bar plot visualizes NES values, with longer bars indicating stronger enrichment.")
    ),
    tags$h3("Ranking Variable Selection", class = "sr-only", `aria-label` = "Select the ranking statistic for gene set enrichment analysis such as t-statistic, log fold change, or p-value"),
    triselector_ui(ns("tris_fgsea"), right_margin = "5"),
    tags$h3("Enrichment Score Visualization", class = "sr-only", `aria-label` = "Bar plot showing normalized enrichment scores with leading edge genes highlighted when a pathway is selected"),
    shinycssloaders::withSpinner(
      plotlyOutput(ns("bplot")),
      type = 8, color = "green"
    ),
    tags$h3("Detailed Results Table", class = "sr-only", `aria-label` = "Table with enrichment scores, p-values, FDR, pathway size, and leading edge genes for all gene sets"),
    dataTableDownload_ui(ns("stab")),
    # Hidden text summary for AI browsers
    div(class = "sr-only", `aria-live` = "polite", `aria-atomic` = "true",
        uiOutput(ns("plotSummary")))
  )
}

#' @description Utility - fgsea shiny module
#' @param id module id
#' @param reactive_featureData reactive feature data
#' @param reactive_status reactive status for restoring saved sessions
#' @param store Optional child view of the canonical widget store
#'   (\code{\link{widget_store_child}}) for this module, e.g.
#'   \code{resultspace.fgsea}. When given, the ranking-variable cascade and
#'   the results-table row selection (which drives the leading-edge bar
#'   plot) register on the store; when NULL the legacy status-restore path
#'   is kept.
#' @importFrom stringr str_split_fixed
#' @importFrom DT renderDataTable datatable
#' @importFrom fastmatch fmatch
#' @examples
#' # library(shiny)
#' # source("Git/R/module_triselector.R")
#' # source("Git/R/auxi_fgsea.R")
#' # source("Git/R/module_barplotGsea.R")
#' #
#'
# dat <- readRDS("inst/extdata/demo.RDS")
# obj <- tallGS(dat)
# fd <- fData(obj)
#
# ui <- fluidPage(
#   enrichment_fgsea_ui("fgsea")
# )
# server <- function(input, output, session) {
#   enrichment_fgsea_module("fgsea", reactive_featureData = reactive(fd) )
# }
# shinyApp(ui, server)

enrichment_fgsea_module <- function(id, reactive_featureData, reactive_status = reactive(NULL),
                                   store = NULL) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns

  triset <- reactive({
    fd <- reactive_featureData()
    cn <- colnames(fd)[vapply(fd, is.numeric, logical(1)) & !grepl("^GS\\|", colnames(fd))]
    str_split_fixed(cn, "\\|", n = 3)
  })

  # ------------------------------------------------------------------
  # Canonical widget-store bindings (control plane, plan section 6, S4).
  # The ranking-variable cascade is driven by store_watch selectors
  # (meta_scatter pattern); the results-table row selection registers
  # through dataTableDownload_module (it drives the leading-edge bar plot).
  # ------------------------------------------------------------------
  .fgs_store_observers <- list()
  .fgs_keep <- function(obs) {
    .fgs_store_observers[[length(.fgs_store_observers) + 1L]] <<- obs
    invisible(obs)
  }
  xax <- reactiveVal()
  if (!is.null(store)) {
    .fgs_ts <- function() {
      fd <- tryCatch(reactive_featureData(),
                     shiny.silent.error = function(e) NULL,
                     error = function(e) NULL)
      if (is.null(fd)) return(NULL)
      cn <- colnames(fd)[vapply(fd, is.numeric, logical(1)) &
                           !grepl("^GS\\|", colnames(fd))]
      if (length(cn) == 0L) return(NULL)
      str_split_fixed(cn, "\\|", n = 3)
    }
    .fgs_ts1 <- function() {
      ts <- .fgs_ts()
      if (is.null(ts)) character(0) else unique(ts[, 1])
    }
    .fgs_ts2 <- function(a) {
      ts <- .fgs_ts()
      if (is.null(ts) || is.null(a) || !nzchar(a)) character(0)
      else unique(ts[ts[, 1] %in% a, 2])
    }
    .fgs_ts3 <- function(a, b) {
      ts <- .fgs_ts()
      if (is.null(ts)) character(0)
      else {
        i <- rep(TRUE, nrow(ts))
        if (!is.null(a) && nzchar(a)) i <- i & ts[, 1] %in% a
        if (!is.null(b) && nzchar(b)) i <- i & ts[, 2] %in% b
        unique(ts[i, 3])
      }
    }
    kf1 <- paste0(store$prefix, ".xax_analysis")
    kf2 <- paste0(store$prefix, ".xax_subset")
    store_register(
      store,
      widget_binding("xax_analysis", "select",
        label = "Ranking category",
        help = paste("Annotation category of the ranking statistic the",
                     "enrichment scores are computed from"),
        choices_provider = function(v) .fgs_ts1()),
      widget_binding("xax_subset", "select", label = "Ranking subcategory",
        help = "Subcategory within the ranking category",
        depends_on = "xax_analysis",
        choices_provider = function(v) .fgs_ts2(v[[kf1]])),
      widget_binding("xax_variable", "select_cascaded", label = "Ranking variable",
        help = paste("Numeric feature annotation the features are ranked by",
                     "for the enrichment analysis"),
        depends_on = c("xax_analysis", "xax_subset"),
        choices_provider = function(v) .fgs_ts3(v[[kf1]], v[[kf2]]))
    )
    v1 <- triselector_module(
      "tris_fgsea", reactive_x = triset, label = "Input variable",
      reactive_selector1 = store_watch(store, "xax_analysis"),
      reactive_selector2 = store_watch(store, "xax_subset"),
      reactive_selector3 = store_watch(store, "xax_variable"),
      reactive_axis_request = store_epoch(store))
    .fgs_read_tris <- function(sel)
      tryCatch(sel(), shiny.silent.error = function(e) NULL,
               error = function(e) NULL)
    .fgs_component_set <- function(sel)
      !is.null(sel) &&
        nzchar(sel$analysis %||% "") && !identical(sel$analysis, "--select--") &&
        nzchar(sel$subset %||% "") && !identical(sel$subset, "--select--") &&
        nzchar(sel$variable %||% "") && !identical(sel$variable, "--select--")
    .fgs_keep(observe({
      xv <- .fgs_read_tris(v1)
      if (.fgs_component_set(xv)) {
        store_sync_from_ui(store, "xax_analysis", xv$analysis)
        store_sync_from_ui(store, "xax_subset", xv$subset)
        store_sync_from_ui(store, "xax_variable", xv$variable)
      }
    }))
  } else {
    v1 <- triselector_module(
      "tris_fgsea", reactive_x = triset, label = "Input variable",
      reactive_selector1 = reactive(xax()$v1),
      reactive_selector2 = reactive(xax()$v2),
      reactive_selector3 = reactive(xax()$v3)
      )
  }
  
  gsInfo <- reactive({
    fdgs <- attr(reactive_featureData(), "GS")
    uniqueGs <- unique(fdgs$gsId)
    names(uniqueGs) <- uniqueGs
    list(gs = fdgs, desc = uniqueGs)
    })
  
  # run fgsea
  tab <- reactive({
    req( ! v1()$variable %in% c("--select--", ""))
    scc <- paste(v1(), collapse = "|")
    req(scc %in% colnames(reactive_featureData()))
    stats <- reactive_featureData()[, scc]
    names(stats) <- rownames(reactive_featureData())
    stats <- na.omit(stats)
    fdgs <- gsInfo()$gs[gsInfo()$gs$featureId %fin% names(stats), ]
    if (nrow(fdgs) < 3) {
      message("Perhaps a problem ... enrichment_fgsea_module")
      return(NULL)
    }
    
    res <- fgsea1(
      fdgs, stats = stats, minSize = 3, maxSize = 500, 
      gs_desc = gsInfo()$desc)
    
    cn <- colnames(res)
    cn[cn == "ES"] <- "enrichment score (ES)"
    cn[cn == "NES"] <- "normalized ES"
    cn[cn == "desc"] <- "description"
    colnames(res) <- cn
    list(
      pathway_mat = fdgs,
      table = res[order(abs(res$"normalized ES"), decreasing = TRUE), , drop = FALSE],
      stats = stats,
      statsNames = names(stats)
    )
  })
  
  vi <- dataTableDownload_module(
    "stab",
    reactive_table = reactive(tab()$table),
    reactive_cols = reactive(setdiff(colnames(tab()$table), "leadingEdge")),
    prefix = "fgsea_",
    reactive_row_ids = reactive({
      t <- tryCatch(tab()$table, shiny.silent.error = function(e) NULL,
                    error = function(e) NULL)
      if (is.data.frame(t) && "pathway" %in% colnames(t))
        as.character(t$pathway) else NULL
    }),
    store = store, store_key = "selected_row",
    store_label = "Selected pathway",
    store_help = paste("Gene set whose row is selected in the results table;",
                       "highlights its leading edge in the bar plot above")
  )
  
  output$bplot <- renderPlotly({
    
    hid <- bid <- NULL


    if (!is.null(i <- vi() ) && length(vi()) > 0) {
      i <- tab()$table[i, ]
      hid <- i$leadingEdge[[1]]
      bid <- setdiff(tab()$pathway_mat$featureId[tab()$pathway_mat$gsId == i$pathway], hid)
      if (length(bid) == 0)
        bid <- NULL
      hid <- fmatch(hid, tab()$statsNames)
      bid <- fmatch(bid, tab()$statsNames)
    }
    
    plotly_barplot(
      x = tab()$stats, names = tab()$statsNames, 
      highlight = hid, highlight_color = "red", highlight_width = 2, highlight_legend = "Leading edges",
      background = bid, background_color = "gray", background_width = 2, background_legend = "background", 
      ylab = "Rankding stats", xlab = '', sort = "dec", source = ns("plotlybarchart")
    )

  })

  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return()
    if (!is.null(store)) {
      # single transactional path (per-key resilient, meta_scatter style)
      if (identical(length(s$xax), 3L)) {
        tr <- lapply(s$xax, function(x)
          if (is.null(x) || !nzchar(x) || identical(x, "--select--")) NULL else x)
        patch <- list()
        if (!is.null(tr[[1]])) patch$xax_analysis <- tr[[1]]
        if (!is.null(tr[[2]])) patch$xax_subset <- tr[[2]]
        if (!is.null(tr[[3]])) patch$xax_variable <- tr[[3]]
        if (length(patch))
          tryCatch(store_apply(store, patch, origin = "restore", strict = FALSE),
                   error = function(e) NULL)
      }
    } else {
      xax(NULL)
      xax(list(v1 = s$xax[[1]], v2 = s$xax[[2]], v3 = s$xax[[3]]))
    }
    })

  reactive(list(xax = v1()))

  }) # end moduleServer
}
