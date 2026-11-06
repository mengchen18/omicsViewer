#' @description Utility enrichment analysis shiny ui
#' @param id id
#' @importFrom DT dataTableOutput
enrichment_analysis_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Module description for AI browsers and screen readers
    div(class = "sr-only", id = ns("module-help"),
      tags$h4("About Over-Representation Analysis (ORA)"),
      tags$p("Over-Representation Analysis identifies gene sets, pathways, and functional categories that are statistically enriched in your selected features. It uses the hypergeometric test to determine if your selected genes/proteins overlap with known biological gene sets more often than would be expected by chance. After a run, the top rows of the results table are readable by the AI assistant through the panel_results state section; column sorting is click-only."),
      tags$h4("When to use ORA"),
      tags$p("Use this analysis when you have a list of interesting features (e.g., significantly changed genes, top-ranked proteins) and want to discover which biological processes, molecular functions, cellular components, or pathways are associated with them. This helps interpret your results in a biological context."),
      tags$h4("How to interpret results"),
      tags$p("The results table shows enriched gene sets ranked by p-value. Lower p-values indicate stronger enrichment. The FDR (False Discovery Rate) column provides multiple-testing corrected p-values - typically FDR < 0.05 is considered significant. The 'overlap' column shows how many of your selected features belong to each gene set, and the 'odds ratio' indicates the strength of enrichment.")
    ),
    uiOutput(ns("error")) %>%
      tagAppendAttributes(`aria-live` = "assertive", `aria-atomic` = "true"),
    tags$h3("Gene Set Collection Selection", class = "sr-only", `aria-label` = "Select gene set database for over-representation analysis such as GO, KEGG, or custom gene sets"),
    triselector_ui(ns("tris_ora"), right_margin = "5"),
    tags$h3("Enrichment Results", class = "sr-only", `aria-label` = "Table showing enriched gene sets with p-values, FDR, odds ratios, and overlap statistics from hypergeometric test"),
    dataTableDownload_ui(ns("stab")),
    tags$h3("Gene Overlap Details", class = "sr-only", `aria-label` = "Detailed list of overlapping genes between selected features and enriched gene set with annotations"),
    dataTableDownload_ui(ns("overlapTab"))
  )
}

#' @description Utility enrichment analysis shiny module
#' @param id module id
#' @param reactive_featureData reactive feature data
#' @param reactive_i reactive index of rows to be selected (for ORA)
#' @param store Optional child view of the canonical widget store
#'   (\code{\link{widget_store_child}}) for this module, e.g.
#'   \code{resultspace.ora}. When given, the collapse-features cascade and
#'   the results-table row selection (which drives the overlap-genes
#'   table) register on the store; when NULL the legacy status-restore
#'   path is kept.
#' @param reactive_visible Optional reactive logical overriding the
#'   browser-reported tab visibility in the R-H3 gates, for headless
#'   regression tests (MockShinySession cannot report the
#'   \code{output_*_hidden} clientData flag). NULL (default) keeps the
#'   real \code{\link{output_visible}} read.
#' @param context Optional panel context bundle (\code{\link{panel_context}};
#'   WP-9, plan D13 - built at the L1 call site). The module consumes
#'   the child store (\code{context$store}) and the per-instance
#'   retention helper (\code{context$keep}) from it; feature data, the
#'   selection and the sliced status keep arriving through their
#'   formals. NULL (default) keeps the direct-instantiation path.
#' @importFrom fastmatch fmatch
#' @importFrom stats cutree
#' @examples
#' #' # source("Git/R/auxi_fgsea.R")
#' # source("Git/R/auxi_vectORA.R")
#' # source("Git/R/module_barplotGsea.R")
# dat <- readRDS("inst/extdata/demo.RDS")
# obj <- tallGS(dat)
# fd <- Biobase::fData(obj)
# fdgs <- attr(fd, "GS")
# selected_ids <- rownames(fd)[fd$`PCA|All|PC1(10.1%)` > 0.02]
# ui <- fluidPage(
#   enrichment_analysis_ui("ea")
# )
# server <- function(input, output, session) {
#   enrichment_analysis_module("ea",
#     reactive_featureData = reactive(fd), reactive_i = reactive(selected_ids)
#   )
# }
# shinyApp(ui, server)


enrichment_analysis_module <- function(
  id, reactive_featureData, reactive_i, reactive_status = reactive(NULL),
  store = NULL, reactive_visible = NULL, context = NULL
) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns
  
  # ------------------------------------------------------------------
  # WP-9 (panel debt, D13): the panel context, built at the L1 call
  # site (WP-11 collapses the per-panel builds into the registry
  # loop). Context wins when present; the legacy formals remain the
  # direct-instantiation path (the WP-5..WP-8 precedent). This module
  # consumes store + keep from the context; reactive_featureData /
  # reactive_i / reactive_status keep arriving through their formals
  # (nothing moved in this WP).
  # ------------------------------------------------------------------
  if (!is.null(context)) {
    store <- context$store %||% store
  }
  
  reactive_pathway <- reactive({
    attr(reactive_featureData(), "GS")
  })
  
  triset <- reactive({
    trisetter(meta = reactive_featureData(), combine = "none")
  })

  # ------------------------------------------------------------------
  # Canonical widget-store bindings (control plane, plan section 6, S4).
  # The collapse-features cascade is driven by store_watch selectors
  # (meta_scatter pattern); the results-table row selection registers
  # through dataTableDownload_module (it drives the overlap-genes table).
  #
  # WP-9 (panel debt): retention is the shared per-instance
  # panel_keep() factory (D3 - the context's closure when
  # context-built, else a fresh local one for direct instantiation;
  # the hand-rolled .ora_keep list is gone), and the cascade
  # boilerplate (guarded fd triset source, ts1/ts2/ts3 choices
  # providers, the three xax_* registrations, triselector wiring
  # through store_watch/store_epoch and the UI->store sync) is the
  # shared panel_triset_bind helper (D4). Labels/help stay
  # panel-owned; the registered keys/kinds/labels/help/depends_on are
  # byte-identical to the hand-rolled block (frozen pre-migration
  # golden, tests/test_panelRegistry.R). The store-path triselector
  # reads the helper's guarded reactive(.ts()) instead of the module's
  # unguarded triset reactive (the PD-6 sanctioned difference,
  # identical for every real data flow); the legacy no-store branch
  # keeps the triset reactive above.
  # ------------------------------------------------------------------
  keep <- if (is.null(context)) panel_keep() else context$keep
  xax <- reactiveVal()
  if (!is.null(store)) {
    v1 <- panel_triset_bind(
      store, keep,
      tri_id = "tris_ora", tri_label = "Collapse features on",
      triset_source = function() {
        fd <- tryCatch(reactive_featureData(),
                       shiny.silent.error = function(e) NULL,
                       error = function(e) NULL)
        if (is.null(fd)) return(NULL)
        trisetter(meta = fd, combine = "none")
      },
      bindings = list(
        analysis = list(
          label = "Collapse category",
          help = paste("Annotation category whose values the input features",
                       "are collapsed on before testing over-representation")),
        subset = list(
          label = "Collapse subcategory",
          help = "Subcategory within the collapse category"),
        variable = list(
          label = "Collapse variable",
          help = paste("Variable the input features are collapsed on",
                       "(feature annotation column)"))))
  } else {
    v1 <- triselector_module(
      "tris_ora", reactive_x = triset, label = "Collapse features on",
      reactive_selector1 = reactive(xax()$v1),
      reactive_selector2 = reactive(xax()$v2),
      reactive_selector3 = reactive(xax()$v3)
      )
  }

  size_bg <- reactiveVal()
  rii <- reactiveVal()
  reactive_pathway_collapsed <- reactiveVal( NULL )
  col_key <- reactiveVal( NULL )

  # Tab-visibility source for the R-H3 gates: the browser's
  # output_<id>_hidden clientData flag, overridable through the
  # reactive_visible argument for headless regression tests
  # (MockShinySession cannot report the flag).
  .ora_visible <- function() {
    if (!is.null(reactive_visible)) isTRUE(reactive_visible())
    else output_visible(session, ns("stab-table"))
  }

  # R-H3 (2026-10-09 fix): a plain observe(), NOT observeEvent. While
  # the tab is hidden the visibility gate below is the observer's ONLY
  # reactive dependency, so the flip back to visible re-runs it and the
  # computation catches up with the current selection before the table
  # renders. The previous observeEvent form ran the handler inside
  # isolate() (bindEvent.Observer), so the clientData read could never
  # re-trigger it: selection changes arriving while the tab was hidden
  # were dropped and the tab kept showing the previous selection's
  # results (volcano rect selection -> stale ORA). The gate must stay
  # FIRST - reads below it only register as dependencies when they
  # execute. Headless sessions report no flag and stay visible
  # (output_visible).
  observe({

    if (!.ora_visible())
      return(NULL)

    req(reactive_i())
    req(reactive_featureData())
    req(rp <- reactive_pathway())

    # Historical contract: an UNSET collapse triple means "no collapse" -
    # run ORA on the raw feature ids. The unified triselector returns NULL
    # while the cascade has never settled a variable (cold start: only the
    # analysis defaults; the variable sits at "--select--" and never
    # commits without allow_unset), so a req(v1()$variable) here suspended
    # this observer permanently until a manual pick - rii() was never
    # written and the whole tab rendered blank with zero diagnostics
    # (regression vs the pre-unification triselector, which preselected
    # "--select--" as a real value and drove the no-collapse branch on
    # every selection change). v1() is read in this observer, so a later
    # manual pick re-fires it into the collapse branch.
    v <- v1()
    if (is.null(v) || v$variable %in% c("", "--select--")) {
      size_bg( nrow(reactive_featureData()) )
      reactive_pathway_collapsed(NULL)
      col_key( NULL )
      # plain observe(): never read back a reactiveVal this observer
      # writes (the old observeEvent form isolated these reads) - the
      # self-dependency would loop the observer forever
      ri <- reactive_i()
      if ( length(ri) <= 1 )
        rii(NULL) else if ( length(ri) <= 3 )
          rii("notest") else
            rii(ri)
      return()
    }

    cs <- do.call(paste, list(v, collapse = "|"))
    if (!cs %in% colnames(reactive_featureData()))
      return(NULL)
    val <- reactive_featureData()[, cs]
    names(val) <- rownames(reactive_featureData())
    # R-H2: NA/"" collapse values must not count as pseudo-genes (they used
    # to inflate size_bg and leak NA into rii(), which crashed the eager
    # oraTab observer with `NA == "notest"` and closed the session).
    val <- val[!is.na(val) & val != ""]
    ck <- val[reactive_i()]
    ck <- ck[!is.na(ck)]
    col_key( ck )

    rk <- unique(ck)
    if ( length(rk) <= 1)
      rii(NULL) else if (length(rk) <=  3)
        rii("notest") else
          rii(rk)

    rp$featureId <- as.factor( val[ as.character( rp$featureId ) ] )
    rp <- rp[!is.na(rp$featureId), ]
    reactive_pathway_collapsed( unique(rp) )

    size_bg( length(unique(val)) )
    })
  
  OT <- reactive({
  
    req(size_bg())
    req(rii())
  
    notest <- "No geneset has been tested, please try to include more input feature IDs!" 

    if (identical(rii()[1], "notest"))
      return(notest)
  
    if (is.null(reactive_pathway_collapsed()))
      rp <- reactive_pathway() else
        rp <- reactive_pathway_collapsed()
  
    tab <- vectORATall(rp, i = rii(), background = size_bg())
  
    if (is.null(tab))
      return(notest)
  
    ic <- which(vapply(tab, function(x) is.numeric(x) & !is.integer(x), logical(1)))  
    tab[, ic] <- lapply(tab[, ic], signif, digits = 3)  
    tab <- tab[which(tab$p.adjusted < 0.1 | tab$p.value < 0.05 | tab$OR >= 3), ]    
  
    if (nrow(tab) > 3) {
      # R-H3: stable order + cap the clustered rows (jaccard/hclust on
      # 1000+ enriched sets was the second hot spot)
      tab <- tab[order(tab$p.value), , drop = FALSE]
      ncl <- min(nrow(tab), ORA_CLUSTER_MAX_ROWS)
      hcl <- hclust(jaccardList(tab$overlap_ids[seq_len(ncl)]))
      cls <- cutree(hcl, h = 0.5)
      tab$desc <- paste("cluster",
                        c(cls, rep(0L, nrow(tab) - ncl)), tab$desc, sep = "_")
    }
  
    tab
  })
  
  oraTab <- reactiveVal( NULL )
  observe({
    # R-H3: no ORA while the tab is hidden (vectORATall + jaccard
    # clustering on every selection change was ~20 s at 1000 sets).
    # Plain observe(): the .ora_visible() read is a real dependency, so
    # the visibility flip re-renders whatever the collapse observer
    # above produced.
    if (!.ora_visible())
      return(NULL)
    if (is.null(rii())) {
      # R-M6: the selection shrinking below the testable minimum must not
      # leave the previous results table on screen - show an explicit
      # no-test message (only once something was shown; not on cold start).
      if (!is.null(isolate(oraTab())))
        oraTab("Too few feature IDs are selected to test over-representation, please try to include more input feature IDs!")
      return(NULL)
    }
    oraTab( OT() )
  })
  
  
  output$errorMsg <- renderText({
    req(is.character(oraTab()))
    oraTab()
  })
  
  output$error <- renderUI(
    verbatimTextOutput(ns("errorMsg"))
  )
  
  vi <- dataTableDownload_module(
    "stab",
    reactive_table = reactive({
      req(is.data.frame(oraTab()))
      oraTab()
    }),
    reactive_cols = reactive( setdiff(colnames(oraTab()), "overlap_ids") ),
    prefix = "ORA_", sortBy = "p.value", decreasing = FALSE, pageLength = ENRICHMENT_TABLE_PAGE_LENGTH,
    reactive_row_ids = reactive({
      t <- tryCatch(oraTab(), shiny.silent.error = function(e) NULL,
                    error = function(e) NULL)
      if (is.data.frame(t) && "pathway" %in% colnames(t))
        as.character(t$pathway) else NULL
    }),
    store = store, store_key = "selected_row",
    store_label = "Selected pathway",
    store_help = paste("Gene set whose row is selected in the results table;",
                       "drives the overlap-genes table below")
  )

  hd <- reactive({
    req(is.data.frame(oraTab()))
    req( i <- vi() )

    ii <- grep("^General", colnames(reactive_featureData()), ignore.case = TRUE)
    if (length(ii) == 0)
      ii <- seq_len( min(3, ncol(reactive_featureData())) )
    i <- oraTab()[i, ]
    hid <- i$overlap_ids[[1]]
    req(hid)
    if (!is.null(ck <- col_key()))
      hid <- names(ck)[ck %in% hid]
    df1 <- reactive_featureData()[hid, ii, drop = FALSE]
    df1 <- cbind(Overlap = "+", df1)
    apath <- reactive_pathway()[reactive_pathway()$gsId == i$pathway, ]
    aid <- setdiff(apath$featureId, hid)
    if (length(aid) > 0) {
      df2 <- reactive_featureData()[aid, ii, drop = FALSE]
      df2 <- cbind(Overlap = "", df2)
      df1 <- rbind(df1, df2)
    }    
    df1
    })

  vi2 <- dataTableDownload_module(
    "overlapTab",
    reactive_table = hd,
    prefix = "ORA_overlapGenes_", pageLength = ENRICHMENT_TABLE_PAGE_LENGTH
  )

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

  # WP-7 (F4): the module return carries the current ORA result table as
  # an ATTRIBUTE (the data-space quickViews pattern): attributes are
  # dropped by .sanitize_widget_state (lapply) and as.list(), so snapshots
  # and the panels state section never see the raw table - only the
  # result-space reader (get_omics_viewer_state sections='panel_results')
  # reads it. oraTab() is a cheap reactiveVal read (visibility-gated
  # computation, R-H3); a character value is the panel's own no-test
  # message and is passed through for the excerpt to convey.
  reactive({
    l <- list(xax = v1())
    attr(l, "table") <- tryCatch(
      oraTab(),
      shiny.silent.error = function(e) NULL,
      error = function(e) NULL)
    l
  })

  }) # end moduleServer
}
