#' @description Utility enrichment analysis shiny ui
#' @param id id
#' @importFrom DT dataTableOutput
enrichment_analysis_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Module description for AI browsers and screen readers
    div(class = "sr-only", id = ns("module-help"),
      tags$h4("About Over-Representation Analysis (ORA)"),
      tags$p("Over-Representation Analysis identifies gene sets, pathways, and functional categories that are statistically enriched in your selected features. It uses the hypergeometric test to determine if your selected genes/proteins overlap with known biological gene sets more than would be expected by chance."),
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
  store = NULL
) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns
  
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
  # ------------------------------------------------------------------
  .ora_store_observers <- list()
  .ora_keep <- function(obs) {
    .ora_store_observers[[length(.ora_store_observers) + 1L]] <<- obs
    invisible(obs)
  }
  xax <- reactiveVal()
  if (!is.null(store)) {
    .ora_ts <- function() {
      fd <- tryCatch(reactive_featureData(),
                     shiny.silent.error = function(e) NULL,
                     error = function(e) NULL)
      if (is.null(fd)) return(NULL)
      ts <- tryCatch(trisetter(meta = fd, combine = "none"),
                     shiny.silent.error = function(e) NULL,
                     error = function(e) NULL)
      if (is.null(ts) || !is.matrix(ts) || nrow(ts) == 0L) NULL else ts
    }
    .ora_ts1 <- function() {
      ts <- .ora_ts()
      if (is.null(ts)) character(0) else unique(ts[, 1])
    }
    .ora_ts2 <- function(a) {
      ts <- .ora_ts()
      if (is.null(ts) || is.null(a) || !nzchar(a)) character(0)
      else unique(ts[ts[, 1] %in% a, 2])
    }
    .ora_ts3 <- function(a, b) {
      ts <- .ora_ts()
      if (is.null(ts)) character(0)
      else {
        i <- rep(TRUE, nrow(ts))
        if (!is.null(a) && nzchar(a)) i <- i & ts[, 1] %in% a
        if (!is.null(b) && nzchar(b)) i <- i & ts[, 2] %in% b
        unique(ts[i, 3])
      }
    }
    ko1 <- paste0(store$prefix, ".xax_analysis")
    ko2 <- paste0(store$prefix, ".xax_subset")
    store_register(
      store,
      widget_binding("xax_analysis", "select",
        label = "Collapse category",
        help = paste("Annotation category whose values the input features",
                     "are collapsed on before testing over-representation"),
        choices_provider = function(v) .ora_ts1()),
      widget_binding("xax_subset", "select", label = "Collapse subcategory",
        help = "Subcategory within the collapse category",
        depends_on = "xax_analysis",
        choices_provider = function(v) .ora_ts2(v[[ko1]])),
      widget_binding("xax_variable", "select_cascaded", label = "Collapse variable",
        help = paste("Variable the input features are collapsed on",
                     "(feature annotation column)"),
        depends_on = c("xax_analysis", "xax_subset"),
        choices_provider = function(v) .ora_ts3(v[[ko1]], v[[ko2]]))
    )
    v1 <- triselector_module(
      "tris_ora", reactive_x = triset, label = "Collapse features on",
      reactive_selector1 = store_watch(store, "xax_analysis"),
      reactive_selector2 = store_watch(store, "xax_subset"),
      reactive_selector3 = store_watch(store, "xax_variable"),
      reactive_axis_request = store_epoch(store))
    # UI -> store: settled triples only (WP2 store_bind_triselector)
    store_bind_triselector(store,
      keys = c(analysis = "xax_analysis", subset = "xax_subset",
               variable = "xax_variable"),
      sel = v1, keep = .ora_keep)
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

  observeEvent(list(
    reactive_i(),
    reactive_featureData(),
    reactive_pathway(),
    v1()
    ), {

    req(reactive_i())
    req(reactive_featureData())
    req(rp <- reactive_pathway())
    req(v1()$variable)

    if (v1()$variable %in% c("", "--select--")) {
      size_bg( nrow(reactive_featureData()) )
      reactive_pathway_collapsed(NULL)
      rii(reactive_i())
      col_key( NULL )
      if ( length(rii()) <= 1 )
        rii(NULL) else if ( length(rii()) <= 3 )
          rii("notest")
      return()
    }
    
    cs <- do.call(paste, list(v1(), collapse = "|"))
    if (!cs %in% colnames(reactive_featureData()))
      return(NULL)
    val <- reactive_featureData()[, cs]
    names(val) <- rownames(reactive_featureData())
    ck <- val[reactive_i()]
    col_key( ck )

    rii(unique(ck))
    if ( length(rii()) <= 1)
      rii(NULL) else if (length(rii()) <=  3)
        rii("notest")

    rp$featureId <- as.factor( val[ as.character( rp$featureId ) ] )
    reactive_pathway_collapsed( unique(rp) )

    size_bg( length(unique(val)) )
    })
  
  OT <- reactive({
  
    req(size_bg())
    req(rii())
  
    notest <- "No geneset has been tested, please try to include more input feature IDs!" 

    if (rii()[1] == "notest")
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
      hcl <- hclust(jaccardList(tab$overlap_ids))
      cls <- cutree(hcl, h = 0.5)
      tab$desc <- paste("cluster", cls, tab$desc, sep = "_")    
    }    
  
    tab
  })
  
  oraTab <- reactiveVal( NULL )
  observe({
    if (!is.null(rii()))
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

  reactive(list(xax = v1()))

  }) # end moduleServer
}
