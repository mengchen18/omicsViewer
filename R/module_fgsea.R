#' @description Utility - fgsea shiny ui
#' @param id id
enrichment_fgsea_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Module description for AI browsers and screen readers
    div(class = "sr-only", id = ns("module-help"),
      tags$h4("About Fast Gene Set Enrichment Analysis (fGSEA)"),
      tags$p("fGSEA is a computational method that determines whether predefined gene sets show statistically significant, concordant differences between two biological states. Unlike ORA which only uses a cutoff, fGSEA uses the entire ranked list of genes, making it more sensitive to subtle but coordinated changes in gene expression across a pathway. After a run, the top rows of the results table are readable by the AI assistant through the panel_results state section; column sorting is click-only."),
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
#' @param context Optional panel context bundle (\code{\link{panel_context}};
#'   WP-10, plan D13 - built at the L1 call site). The module consumes
#'   the child store (\code{context$store}) and the per-instance
#'   retention helper (\code{context$keep}) from it; the feature data
#'   and the sliced status keep arriving through their formals. NULL
#'   (default) keeps the direct-instantiation path.
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
                                   store = NULL, context = NULL) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns

  # ------------------------------------------------------------------
  # WP-10 (panel debt, D13): the panel context, built at the L1 call
  # site (WP-11 collapses the per-panel builds into the registry
  # loop). Context wins when present; the legacy formals remain the
  # direct-instantiation path (the WP-5..WP-9 precedent). This module
  # consumes store + keep from the context; reactive_featureData /
  # reactive_status keep arriving through their formals (nothing
  # moved in this WP).
  # ------------------------------------------------------------------
  if (!is.null(context)) {
    store <- context$store %||% store
  }

  triset <- reactive({
    fd <- reactive_featureData()
    cn <- colnames(fd)[vapply(fd, is.numeric, logical(1)) & !grepl("^GS\\|", colnames(fd))]
    str_split_fixed(cn, "\\|", n = 3)
  })

  # ------------------------------------------------------------------
  # Canonical widget-store bindings (control plane, plan section 6, S4).
  # The ranking-variable cascade is driven by store_watch selectors
  # (meta_scatter pattern); the results-table row selection registers
  # through dataTableDownload_module (it drives the leading-edge bar
  # plot).
  #
  # WP-10 (panel debt): retention is the shared per-instance
  # panel_keep() factory (D3 - the context's closure when
  # context-built, else a fresh local one for direct instantiation;
  # the hand-rolled .fgs_keep list is gone), and the cascade
  # boilerplate (guarded fd triset source, ts1/ts2/ts3 choices
  # providers, the three xax_* registrations, triselector wiring
  # through store_watch/store_epoch and the UI->store sync) is the
  # shared panel_triset_bind helper (D4). Labels/help stay
  # panel-owned; the registered keys/kinds/labels/help/depends_on are
  # byte-identical to the hand-rolled block (frozen pre-migration
  # golden, tests/test_panelRegistry.R). The triset source keeps this
  # panel's OWN numeric-column derivation (numeric, non-GS feature
  # annotations - a ranking statistic must be numeric; NOT ORA's
  # trisetter(combine="none") whole-column inventory). The store-path
  # triselector reads the helper's guarded reactive(.ts()) instead of
  # the module's unguarded triset reactive (the PD-6 sanctioned
  # difference, identical for every real data flow); the legacy
  # no-store branch keeps the triset reactive above.
  # ------------------------------------------------------------------
  keep <- if (is.null(context)) panel_keep() else context$keep
  xax <- reactiveVal()
  if (!is.null(store)) {
    v1 <- panel_triset_bind(
      store, keep,
      tri_id = "tris_fgsea", tri_label = "Input variable",
      triset_source = function() {
        fd <- tryCatch(reactive_featureData(),
                       shiny.silent.error = function(e) NULL,
                       error = function(e) NULL)
        if (is.null(fd)) return(NULL)
        cn <- colnames(fd)[vapply(fd, is.numeric, logical(1)) &
                             !grepl("^GS\\|", colnames(fd))]
        if (length(cn) == 0L) return(NULL)
        str_split_fixed(cn, "\\|", n = 3)
      },
      bindings = list(
        analysis = list(
          label = "Ranking category",
          help = paste("Annotation category of the ranking statistic the",
                       "enrichment scores are computed from")),
        subset = list(
          label = "Ranking subcategory",
          help = "Subcategory within the ranking category"),
        variable = list(
          label = "Ranking variable",
          help = paste("Numeric feature annotation the features are ranked by",
                       "for the enrichment analysis"))))
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

  # WP-7 (F4): mirror the computed result table into a reactiveVal behind
  # the same visibility gate the outputs use (R-H3 posture: never compute
  # fgsea while the tab is hidden). The module return then carries the
  # table as an ATTRIBUTE (the data-space quickViews pattern) - attributes
  # are dropped by .sanitize_widget_state (lapply) and as.list(), so
  # snapshots and the panels state section never see the raw table; only
  # the result-space reader (get_omics_viewer_state
  # sections='panel_results') reads it. The reactiveVal read is cheap -
  # the module return gains no tab() dependency, so the fgsea computation
  # stays visibility-gated. tab() req()s a settled ranking variable; while
  # unset the observer simply suspends and fgs_tab stays NULL.
  fgs_tab <- reactiveVal(NULL)
  keep(observe({
    if (!output_visible(session, ns("stab-table")))
      return(NULL)
    t <- tab()
    fgs_tab(if (is.list(t) && is.data.frame(t$table)) t$table else NULL)
  }))

  reactive({
    l <- list(xax = v1())
    attr(l, "table") <- tryCatch(
      fgs_tab(),
      shiny.silent.error = function(e) NULL,
      error = function(e) NULL)
    l
  })

  }) # end moduleServer
}
