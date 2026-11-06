#' @description Utility sample general ui
#' @param id id
#' 
sample_general_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Module description for AI browsers and screen readers
    div(class = "sr-only", id = ns("module-help"),
      tags$h4("About Sample Group Analysis"),
      tags$p("This module provides comprehensive statistical analysis to compare and visualize sample groups or patient cohorts. It offers multiple analysis types including beeswarm plots for continuous variables, contingency table analysis for categorical variables, and survival analysis for time-to-event data. This is the sample-centric complement to the feature analysis module."),
      tags$h4("When to use sample analysis"),
      tags$p("Use this module when you want to analyze relationships between sample metadata variables, compare patient groups, test associations between clinical variables, or perform survival analysis stratified by sample characteristics. This is essential for understanding sample clustering, identifying clinical correlations, and validating patient groupings."),
      tags$h4("How to interpret results"),
      tags$p("The visualization type depends on your selected variables: beeswarm plots show distributions of continuous variables across groups with statistical tests, contingency tables show relationships between categorical variables with chi-square or Fisher's exact test, and survival curves show time-to-event outcomes stratified by groups with log-rank test. The results table provides detailed statistical metrics appropriate for each analysis type.")
    ),
    fluidRow(
      shinydashboard::box(
        batch_comparison_ui(ns("batch_comp")),
        height = "490px",
        width = 12
      ),
      column(12, style = "margin-top: 0px;", triselector_ui(ns("tris_sample_general"), right_margin = "5")),
      shinydashboard::box(
        dataTableDownload_ui(ns("msatab")) ,
        height = "500px",
        width = 5
      ),
      shinydashboard::box(
        # WP6: the attr4 panel is STATIC (it was re-created by every
        # renderUI run, rebuilding its widgets on unchanged view types);
        # the plot itself swaps through the uiOutput below
        tagList(
          column(11, uiOutput(ns("sample_general_plot"))),
          column(1, attr4selector_ui(ns("a4_gp"), circle = FALSE, right = TRUE))
        ),
        height = "500px",
        width = 7
      )
    )
  )
}

#' @description Utility sample general module
#' @param id module id
#' @param reactive_phenoData reactive phenotype data
#' @param reactive_expr reactive expression data
#' @param reactive_j index for which row in phenotype data should be highlighted/selected
#' @param reactive_status saved status to restore
#' @param store Optional child view of the canonical widget store
#'   (\code{\link{widget_store_child}}) for this module, e.g.
#'   \code{resultspace.sample_general}. When given, the module's
#'   user-editable surface (the link-variable cascade, the embedded
#'   attribute-4 panel, the survival censor slider and the batch
#'   comparison switches/rows) registers on the store; when NULL the
#'   legacy status-restore path is kept.
#' @param context Optional panel context bundle (\code{\link{panel_context}};
#'   plan D13). When given, the module consumes the child store
#'   (\code{context$store}) and the per-instance retention helper
#'   (\code{context$keep}) from it; data/selection/status keep arriving
#'   through their formals (nothing moved into the module in WP-8). All
#'   legacy formals keep working for direct instantiation.
#' @examples
#' #' # library(shiny)
#' # #
#' # source("Git/R/module_triselector.R")
#' # source("Git/R/module_scatter.R")
#' # source("Git/R/module_contTableStats.R")
#' # source("Git/R/module_survival.R")
#' # source("Git/R/module_figureAttr4.R")
#' # source("Git/R/auxi_figureAttr4.R")
#' #
#' # dat <- readRDS("Dat/exampleEset.RDS")
#' # pd <- pData(dat)
#' #
#' # ui <- fluidPage(
#' #   sample_general_ui("sample_general")
#' # )
#' #
#' # server <- function(input, output, session) {
#   sample_general_module("sample_general", reactive_phenoData = reactive(pd), #' reactive_j = reactive(sample(rownames(pd), size = 20)) )
#' # }
#' #
#' # shinyApp(ui, server)
#'
sample_general_module <- function(id, reactive_phenoData, reactive_expr,
  reactive_j = reactive(NULL), reactive_status = reactive(NULL),
  store = NULL, context = NULL) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns

  # ------------------------------------------------------------------
  # WP-8 (panel debt, D13): the panel context, built at the L1 call
  # site (WP-11 collapses the per-panel builds into the registry
  # loop). Context wins when present; the legacy formals remain the
  # direct-instantiation path (the feature/ptm/geneshot precedent).
  # This module consumes store + keep from the context; data/selection/
  # reactive_status keep arriving through their formals (nothing moved
  # in this WP).
  # ------------------------------------------------------------------
  if (!is.null(context)) {
    store <- context$store %||% store
  }
  # WP-8 (panel debt, D3): ONE retention closure for the whole module.
  # The organic-growth split (a store-glue .sg_keep list inside the
  # store branch plus a separate .sg_view_keep list for the view-type
  # carrier observer) is unified: the context's panel_keep() closure
  # when context-built, else a fresh local one for direct
  # instantiation. Retention is retention (the WP-7 pattern).
  keep <- if (is.null(context)) panel_keep() else context$keep

  triset <- reactive({
    trisetter(meta = reactive_phenoData(), expr = reactive_expr(), combine = "pheno")
  })

  xax <- reactiveVal()

  # ------------------------------------------------------------------
  # Canonical widget-store bindings (control plane, plan section 6, S4):
  # the link-variable cascade (beeswarm / contingency / survival switch).
  # The attr4 panel registers under <prefix>.attr4.* inside
  # attr4selector_module; the survival submodule's censor slider under
  # <prefix>.survival_censor; the batch comparison switches/rows under
  # <prefix>.batch_* - all through this panel's child store, all
  # untouched by WP-8. Choices providers mirror triset() WITHOUT req()
  # (a shiny.validation condition inside a provider aborts the whole
  # store_apply transaction) and include the Surv category, which decides
  # the survival view (the deliberate opposite of feature-general's
  # minus-Surv cascade).
  #
  # WP-8 (panel debt): retention is the shared per-instance
  # panel_keep() factory (D3), and the cascade boilerplate (guarded
  # pheno triset source KEEPING the Surv rows - Kaplan-Meier -,
  # ts1/ts2/ts3 choices providers, the three xax_* registrations, the
  # triselector wiring through store_watch/store_epoch and the
  # UI->store sync) is the shared panel_triset_bind helper (D4). The
  # .sg_apply_triple dual-path helper (store / legacy xax) stays
  # verbatim (the module's internal legacy fallback, not L1 glue);
  # labels/help stay panel-owned; the registered keys/kinds/labels/
  # help/depends_on are byte-identical to the hand-rolled block (frozen
  # pre-migration golden, tests/test_panelRegistry.R). The trio
  # registers FIRST so the binding registration order stays
  # xax_analysis, xax_subset, xax_variable, then attr4/survival/batch
  # exactly as the pre-WP-8 module produced (store_registry_view/
  # store_read iterate in registration order).
  # ------------------------------------------------------------------
  if (!is.null(store)) {
    v1 <- panel_triset_bind(
      store, keep,
      tri_id = "tris_sample_general", tri_label = "Link selection to",
      triset_source = function() {
        if (is.null(reactive_phenoData())) return(NULL)
        ts <- tryCatch(
          trisetter(meta = reactive_phenoData(), expr = reactive_expr(),
                    combine = "pheno"),
          shiny.silent.error = function(e) NULL, error = function(e) NULL)
        if (is.null(ts) || !is.matrix(ts) || nrow(ts) == 0L) NULL else ts
      },
      bindings = list(
        analysis = list(
          label = "Link variable category",
          help = paste("Annotation category of the variable the selected",
                       "samples are analyzed against")),
        subset = list(
          label = "Link variable subcategory",
          help = "Subcategory within the link-variable category"),
        variable = list(
          label = "Link variable",
          help = paste("Variable driving the analysis view: numeric gives a",
                       "beeswarm, categorical a contingency table, Surv a",
                       "Kaplan-Meier curve"))),
      allow_unset = TRUE)
    # programmatic triple writes (status restore, batch link) share one
    # helper; user picks sync back through the bind observer inside
    # panel_triset_bind (store_bind_triselector, WP2 semantics)
    .sg_apply_triple <- function(p1, p2, p3, origin = "system") {
      patch <- stats::setNames(list(p1, p2, p3),
                               c("xax_analysis", "xax_subset", "xax_variable"))
      tryCatch(store_apply(store, patch, origin = origin, strict = FALSE),
               error = function(e) NULL)
    }
  } else {
    .sg_apply_triple <- function(p1, p2, p3, origin = "system") {
      xax(NULL)
      xax(list(v1 = p1, v2 = p2, v3 = p3))
    }
    v1 <- triselector_module(
      "tris_sample_general", reactive_x = triset, label = "Link selection to",
      reactive_selector1 = reactive(xax()$v1),
      reactive_selector2 = reactive(xax()$v2),
      reactive_selector3 = reactive(xax()$v3), allow_unset = TRUE)
  }

  attr4select_status <- reactiveVal()
  # R-M11 (todo 4.3): the attr4 panel and the link variable resolve
  # through varSelector (Feature|Auto via the expression matrix); the old
  # reactive_input built cbind(pData, t(expr)) - one column PER FEATURE -
  # on every selection change to extract single columns.
  attr4select <- attr4selector_module(
    "a4_gp", reactive_meta = reactive_phenoData,
    reactive_expr = reactive_expr,
    reactive_triset = triset, reactive_status = attr4select_status,
    store = store
  )

  pheno <- reactive({
    tv <- tryCatch(v1(), shiny.silent.error = function(e) NULL,
                   error = function(e) NULL)
    if (is.null(tv) || tv$variable %in% c("", "--select--"))
      return(NULL)
    val <- tryCatch(
      varSelector(tv, expr = reactive_expr(), meta = reactive_phenoData()),
      error = function(e) NULL)
    if (is.null(val))
      return(NULL)
    
    if (v1()$analysis == "Surv") {
      type <- "surv"
    } else if (is.numeric(val)) {
      type <- "beeswarm"  
    } else if (is.character(val) || is.factor(val)) {
      type <- "table"
    } else {
      warnings("Unknown type of val: sample_general_module, return NULL!")
      return(NULL)
    }
    list(value = val, type = type)
  })
  
  select <- reactive({
    req(reactive_j())
    req(reactive_phenoData())
    select <- rep("Unselected", nrow(reactive_phenoData()))
    if (!is.null(reactive_j()))
      select[rownames(reactive_phenoData()) %in% reactive_j()] <- "selected"
    select
  })
  
  # WP6: derive the view type once through a value-deduped reactiveVal
  # carrier (one kept observer) - pheno() recomputes on every link change
  # but only its $type decides the view, so the container must not be
  # rebuilt while the type is unchanged (remounts reset plot state and
  # rebuild the attr4 panel). Observer-GC rule: kept referenced - on
  # the module's ONE panel_keep closure since WP-8 (retention is
  # retention; the split keep lists were organic growth).
  .sg_view_type <- reactiveVal(NULL)
  keep(observe({
    vt <- tryCatch(pheno()$type, shiny.silent.error = function(e) NULL,
                   error = function(e) NULL)
    if (!identical(vt, .sg_view_type()))
      .sg_view_type(vt)
  }))

  output$sample_general_plot <- renderUI({
    vt <- .sg_view_type()
    if (identical(vt, "beeswarm"))
      return(plotly_scatter_ui(ns("sample_general_beeswarm")))
    if (identical(vt, "table"))
      return(factorIndependency_ui(ns("sample_general_contab")))
    if (identical(vt, "surv"))
      return(survival_ui(ns("sample_general_surv")))
    NULL
  })
  
  ## beeswarm
  # showRegLine <- reactiveVal(FALSE)
  htestV1 <- reactiveVal()
  htestV2 <- reactiveVal()
  .sg_keep_observers <- list()
  vs_scatter <- plotly_scatter_module(
    "sample_general_beeswarm",
    reactive_param_plotly_scatter = reactive({
      req(reactive_j())
      req(pheno()$value)
      tooltips <- attr4select$tooltips
      if (is.null(tooltips))
        tooltips <- rownames(reactive_phenoData())
      l <- list(
        x = select(),
        y = pheno()$value,
        xlab = "",
        ylab = do.call(paste, list(v1(), collapse = "|")),
        tooltips = tooltips
      )
      l$color <- attr4select$color
      l$shape <- attr4select$shape
      l$size <- attr4select$size
      l$highlight <- attr4select$highlight
      l$highlightName <- attr4select$highlightName
      l
    }),
    reactive_regLine = reactive(FALSE), # showRegLine,
    reactive_checkpoint = reactive(pheno()$type == "beeswarm"),
    htest_var1 = htestV1, htest_var2 = htestV2)
  # L7: keep the user's t-test group picks across re-renders (observer
  # retained per the observer-GC rule)
  .sg_keep_observers$htest <- observe({
    v <- vs_scatter()
    if (!is.null(v$htest_V1)) htestV1(v$htest_V1)
    if (!is.null(v$htest_V2)) htestV2(v$htest_V2)
  })

  # cont table stats
  factorIndependency_module("sample_general_contab",
             x = select, y = reactive(pheno()$value),
             reactive_checkpoint = reactive(pheno()$type == "table")
  )

  ## survival
  survival_module('sample_general_surv',
             reactive_resp = reactive(pheno()$value), reactive_strata = select,
             reactive_checkpoint = reactive(pheno()$type == "surv"),
             store = store
  )

  ## table
  metatab <- reactive({
    req(reactive_j())
    tab <- reactive_phenoData()
    kp <- strsplit(colnames(tab), split = "\\|")[[1]][[1]]
    kp <- sprintf("^%s\\|", kp)
    tab <- tab[, grep(kp, colnames(tab)), drop = FALSE]
    tab <- tab[reactive_j(), , drop = FALSE]
    ic <- vapply(tab, is.numeric, logical(1)) & vapply(tab, is.integer, logical(1))
    tab[ic] <- lapply(tab[ic], signif, digits = 2)
    colnames(tab) <- sub('General\\|All\\|', "", colnames(tab))
    tab
  })

  dataTableDownload_module(
    "msatab", reactive_table = metatab, prefix = "SampleTable_"
  )

  ## batch comparison module
  # Convert reactive_j (row names) to indices for batch_comparison_module
  reactive_i_samples <- reactive({
    req(reactive_j())
    req(reactive_phenoData())
    which(rownames(reactive_phenoData()) %in% reactive_j())
  })

  # Call batch comparison module
  batch_comp_selected <- batch_comparison_module(
    "batch_comp",
    reactive_expr = reactive_expr,
    reactive_phenoData = reactive_phenoData,
    reactive_featureData = reactive(NULL),  # No feature data available in this context
    reactive_i_samples = reactive_i_samples,
    store = store
  )

  # Update triselector when a row is selected from batch comparison
  observeEvent(batch_comp_selected(), {
    selected <- batch_comp_selected()
    req(!is.null(selected))

    if (selected$source == "phenotype") {
      # Parse variable_name to extract category|subcategory|variable
      var_name <- selected$data$variable_name
      parts <- strsplit(var_name, "\\|")[[1]]

      if (length(parts) >= 3) {
        .sg_apply_triple(parts[1], parts[2], parts[3])
      }
    } else if (selected$source == "features") {
      # For features: set to "Feature" -> "Auto" -> [feature_name]
      feature_name <- selected$data$feature_name
      .sg_apply_triple("Feature", "Auto", feature_name)
    }
  })

  # WP-8 (panel debt): the cascade UI->store sync (store_bind_triselector,
  # WP2 semantics: settled triples only, mid-cascade echoes never enter
  # the store, store pushes acknowledged through the same path) moved
  # into panel_triset_bind above. The cascade has no meaningful default,
  # so it is not seeded and fills on the first pick.

  ## save and restore status
  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return()
    if (identical(length(s$xax), 3L)) {
      tr <- lapply(s$xax, function(x)
        if (is.null(x) || !nzchar(x) || identical(x, "--select--")) NULL else x)
      if (!any(vapply(tr, is.null, logical(1)))) {
        .sg_apply_triple(tr[[1]], tr[[2]], tr[[3]], origin = "restore")
      } else if (is.null(tr[[3]])) {
        # saved unset: remove a link variable set after the snapshot
        # (todo 2.3); the store path pushes the canonical unset
        if (!is.null(store)) {
          patch <- list()
          if (!is.null(tr[[1]])) patch$xax_analysis <- tr[[1]]
          if (!is.null(tr[[2]])) patch$xax_subset <- tr[[2]]
          patch$xax_variable <- "--select--"
          tryCatch(store_apply(store, patch, origin = "restore", strict = FALSE),
                   error = function(e) NULL)
        } else {
          xax(NULL)
        }
      }
    }
    })

  observeEvent(reactive_status(), {
    if (is.null(s <- reactive_status()))
      return()
    attr4select_status(NULL)
    attr4select_status(s$attr4)    
    })

  ## return status ##
  rv <- reactiveValues()
  observe( rv$xax <- v1() )
  observe( rv$attr4 <- attr4select$status )
  reactive(
    reactiveValuesToList(rv)
    )

  }) # end moduleServer
}
