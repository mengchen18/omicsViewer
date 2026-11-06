#' @title Utility L1 result space ui
#' @param id id
#' @return A Shiny UI definition (\code{tagList}).
#' @importFrom shinythemes shinytheme
L1_result_space_ui <- function(id) {
  ns <- NS(id)

  tagList(
    uiOutput(ns("optTabs"))
  )
}

#' @title Application level 1 module - result space
#' @param id module id
#' @param reactive_expr expression matrix
#' @param reactive_phenoData phentype data
#' @param reactive_featureData feature data
#' @param reactive_i row ID/name of rows selected
#' @param reactive_highlight col ID/name of columns selected
#' @param additionalTabs additional tabs added to "Analyst" panel
#' @param object originally loaded object, mostly an \code{ExpressionSet} or \code{SummarizedExperiment} object
#' @param status intial status
#' @param store Root of the canonical widget store
#'   (\code{\link{widget_store_new}}). Result-space child stores are
#'   derived from it (\code{resultspace.*}) and the analyst tab navbar is
#'   registered for agent control.
#' @return A Shiny module server function, called for its side effects
#'   (renders the result-space panels and registers their widget-store
#'   bindings).
#' @keywords internal
L1_result_space_module <- function(
  id,
  reactive_expr, reactive_phenoData, reactive_featureData,
  reactive_i = reactive(NULL),
  reactive_highlight = reactive(NULL),
  additionalTabs = NULL,
  object = NULL, status = reactive(NULL),
  store = NULL
) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ------------------------------------------------------------------
    # Canonical widget-store bindings (control plane, plan section 6, S4):
    # the analyst tab navbar itself, so the agent can switch analysis
    # tabs exactly like a user. This block is L1's OWN state (plan
    # appendix "L1-owned"); the per-panel child stores are created by
    # the WP-11 registry loop below, in registry (server-composition)
    # order - the pre-WP-11 hardcoded child-store lines dissolved there.
    # ------------------------------------------------------------------
    store_rs <- NULL  # referenced by the status-restore observer even when store is NULL
    if (!is.null(store)) {
      store_rs <- widget_store_child(store, "resultspace")
      # tab titles in renderUI order; providers must never req()
      # WP-1 (P2): thin wrapper over the ONE availability function
      # (rs_tab_set); fd acquisition stays here, the availability
      # predicates/ordering/uniquing are the registry's
      .rs_tab_choices <- function() {
        fd <- tryCatch(reactive_featureData(),
                       shiny.silent.error = function(e) NULL,
                       error = function(e) NULL)
        rs_tab_set(fd, additional = if (length(additionalTabs) > 0)
          vapply(additionalTabs, function(lo) lo$tabName, character(1)))
      }
      store_register(
        store_rs,
        widget_binding("analyst_tab", "navbar", label = "Analysis tab",
          help = paste("Active tab of the analysis (result-space) navbar:",
                       "Feature, ORA, fGSEA, Response, StringDB, SeqLogo,",
                       "Geneshot, or Sample (dataset-dependent)"),
          choices_provider = function(v) .rs_tab_choices())
      )
      .rs_store_observers <- list()
      .rs_keep <- function(obs) {
        .rs_store_observers[[length(.rs_store_observers) + 1L]] <<- obs
        invisible(obs)
      }
      .rs_root_store <- if (is.null(store$parent)) store else store$parent
      .rs_keep(observeEvent(input$analyst, {
        if (!is.null(input$analyst))
          store_sync_from_ui(store_rs, "analyst_tab", input$analyst)
      }, ignoreInit = TRUE))
      .rs_seeded <- FALSE
      .rs_keep(observe({
        if (.rs_seeded) return(NULL)
        if (is.null(input$analyst)) return(NULL)
        .rs_seeded <<- TRUE
        store_seed(store_rs, list(analyst_tab = input$analyst))
      }))
      .rs_epoch <- store_epoch(store_rs)
      .rs_keep(observe({
        .rs_epoch()
        tb <- store_read(store_rs, "analyst_tab")[[1]]
        if (!is.null(tb) &&
            !is.null(.rs_root_store$pending[[paste0(store_rs$prefix, ".analyst_tab")]]))
          updateNavbarPage(session = session, inputId = "analyst", selected = tb)
      }))
    }

    # ------------------------------------------------------------------
    # WP-11 (plan section 4): L1 SERVER-side composition from the panel
    # registry. The eight per-panel call sites with their bespoke formals
    # and per-panel context builds (WP-3..WP-10, D13) collapse into ONE
    # loop over panel_records: per record, create the child store view
    # (dose keeps its resultspace.dose_response namespace, D2/D18), build
    # the panel_context (the established pattern: data/selection/status/
    # session/active), and invoke the record's server with its module id
    # plus the legacy formals (the modules still consume data through
    # those; context data fields stay dormant). Module returns are
    # collected into a map keyed by status_key for the status builder
    # and the enrichment reader below.
    #
    # ORDER INVARIANT (call order = snapshot/enumeration order): the
    # registry list order IS the module invocation order, which is the
    # binding registration order - the store_snapshot bytes (values+
    # epochs+actions, names AND order) and the capability-record order
    # depend on it. The pre-WP-11 literal call order was feature_general,
    # fgsea, ora, stringdb, sample_general, geneshotTab, ptm, rescurve;
    # the registry records were reordered to exactly that (data fix in
    # auxi_panelRegistry.R), never the frozen goldens (D12).
    # ------------------------------------------------------------------
    .rs_active_tab <- reactive(status()$analyst_active_tab)
    .rs_sources <- list(
      expr = reactive_expr, pdata = reactive_phenoData,
      fdata = reactive_featureData, object = object,
      feature = reactive_i, sample = reactive_highlight)
    # which context data/selection fields each panel receives (the exact
    # per-call-site shapes; Filter(is.function) keeps a NULL `object` out
    # of the dose context, PD-17)
    .rs_ctx_data <- list(
      feature_general = c("expr", "pdata", "fdata"),
      fgsea = "fdata",
      ora = "fdata",
      stringdb = "fdata",
      sample_general = c("expr", "pdata"),
      geneshotTab = "fdata",
      ptm = "fdata",
      rescurve = c("expr", "pdata", "fdata", "object"))
    .rs_ctx_selection <- list(
      feature_general = c("feature", "sample"),
      fgsea = "feature",
      ora = "feature",
      stringdb = "feature",
      sample_general = "sample",
      geneshotTab = "feature",
      ptm = "feature",
      rescurve = c("feature", "sample"))
    # the legacy formals per module, verbatim from the dissolved call
    # sites (D13 argument coalescing: context wins for store/keep, the
    # explicit formals remain the direct-instantiation path). status
    # slices come from the record's status_key (the special-arg glue is
    # gone: string ids moved in WP-4, dose attrs in WP-3)
    .rs_server_args <- list(
      feature_general = function(slice) list(
        reactive_expr = reactive_expr,
        reactive_i = reactive_i,
        reactive_highlight = reactive_highlight,
        reactive_phenoData = reactive_phenoData,
        reactive_featureData = reactive_featureData,
        reactive_status = slice),
      fgsea = function(slice) list(
        reactive_featureData = reactive_featureData,
        reactive_status = slice),
      ora = function(slice) list(
        reactive_i = reactive_i,
        reactive_featureData = reactive_featureData,
        reactive_status = slice),
      # string consumes everything through its context since WP-4 (the
      # ids derivation moved into the module)
      stringdb = function(slice) list(),
      sample_general = function(slice) list(
        reactive_phenoData = reactive_phenoData,
        reactive_expr = reactive_expr,
        reactive_j = reactive_highlight,
        reactive_status = slice),
      geneshotTab = function(slice) list(
        fdata = reactive_featureData,
        feature_selected = reactive_i,
        reactive_status = slice),
      ptm = function(slice) list(
        fdata = reactive_featureData,
        feature_selected = reactive_i,
        reactive_status = slice
        # background = reactive(attr(object(), "ptm.seq.window"))
        # background = reactive({
        #   i <- grep("^PTMSeq\\|", colnames(reactive_featureData()))
        #   reactive_featureData()[[i]]
        # })
        ),
      # dose: the S6.6_drc derivation moved into the module in WP-3
      # (over context$data$object); no reactive_status formal is passed
      rescurve = function(slice) list(
        reactive_expr = reactive_expr,
        reactive_i = reactive_i,
        reactive_phenoData = reactive_phenoData,
        reactive_featureData = reactive_featureData))

    .rs_panel_results <- list()
    # static guard: the loop tables must cover exactly the registry's
    # module ids (actionable boot failure when a future record is added
    # without its loop entry - the failures below would otherwise be
    # silent empty contexts / "attempt to apply non-function")
    .rs_loop_ids <- vapply(panel_records, function(r) r$module_id, "")
    for (.rs_tbl in list(.rs_server_args, .rs_ctx_data, .rs_ctx_selection)) {
      if (!setequal(.rs_loop_ids, names(.rs_tbl)))
        stop("L1 result space: the WP-11 registry loop tables must cover ",
             "exactly the registry module ids (",
             paste(.rs_loop_ids, collapse = ", "), "); got: ",
             paste(names(.rs_tbl), collapse = ", "), call. = FALSE)
    }
    for (r in panel_records) {
      .rs_child <- if (is.null(store)) NULL else
        widget_store_child(store, r$store_id)
      .rs_slice <- if (is.null(r$status_key)) reactive(NULL) else
        reactive(status()[[r$status_key]])
      .rs_ctx <- panel_context(
        id = r$module_id, tab_value = r$tab_value, store = .rs_child,
        data = Filter(is.function,
                      .rs_sources[.rs_ctx_data[[r$module_id]]]),
        selection = Filter(is.function,
                           .rs_sources[.rs_ctx_selection[[r$module_id]]]),
        status = .rs_slice,
        active_tab = .rs_active_tab,
        session = session)
      .rs_panel_results[[r$status_key %||% r$module_id]] <- do.call(
        r$server,
        c(list(r$module_id), .rs_server_args[[r$module_id]](.rs_slice),
          list(store = .rs_child, context = .rs_ctx)))
    }

    # ------------------------------------------------------------------
    # WP-7 (F4): bounded read access to the enrichment result tables.
    # The reader runs from tool contexts (no reactive context), so every
    # read is isolate()d + tryCatch-guarded; the modules expose their
    # current tables as attributes on their returns (never as status
    # elements - snapshots must not grow result tables). Returned as an
    # attribute on this module's status value, mirroring the data-space
    # quickViews pattern; as.list()/.sanitize_widget_state drop it.
    # L1-owned (plan appendix); since WP-11 the ora/fgsea module returns
    # come from the registry loop's result map (identical reactive
    # values - the loop replaced the v3/v2 call sites 1:1).
    # ------------------------------------------------------------------
    .read_one_enrichment_table <- function(module_return, method) {
      tab <- tryCatch(
        attr(shiny::isolate(module_return()), "table"),
        shiny.silent.error = function(e) NULL,
        error = function(e) NULL)
      agent_enrichment_table_excerpt(tab, method)
    }
    read_enrichment_results <- function(method = NULL, n = 15L) {
      m <- tolower(trimws(if (is.null(method)) "" else as.character(method)[1]))
      if (!m %in% c("", "ora", "fgsea"))
        stop("Unknown enrichment method: ", m,
             ". Use 'ora', 'fgsea', or omit it for both.")
      out <- list()
      if (m %in% c("", "ora"))
        out$ora <- .read_one_enrichment_table(
          .rs_panel_results$analyst_ora, "ORA")
      if (m %in% c("", "fgsea"))
        out$fgsea <- .read_one_enrichment_table(
          .rs_panel_results$analyst_fgsea, "fGSEA")
      out
    }

    #
    additional_module_results <- list()
    if (length(additionalTabs) > 0) {
      for (lo in additionalTabs) {
        module_args <- list(
          pdata = reactive_phenoData, fdata = reactive_featureData, expr = reactive_expr,
          feature_selected = reactive_i, sample_selected = reactive_highlight, object = object
        )
        if (isTRUE(lo$stateful)) {
          module_args$reactive_status <- reactive(
            status()$analyst_additional[[lo$moduleName]]
          )
        }
        additional_module_results[[lo$moduleName]] <- do.call(
          lo$moduleServer, c(list(lo$moduleName), module_args)
        )
      }
    }

    #### status for snapshot #####
    safe_module_state <- function(x) {
      tryCatch(.sanitize_widget_state(x()), error = function(e) NULL)
    }

    observe({
      if (!is.null(tb <- status()$analyst_active_tab)) {
        if (!is.null(store_rs)) {
          # single transactional path (per-key resilient, meta_scatter style)
          tryCatch(store_apply(store_rs, list(analyst_tab = tb),
                               origin = "restore", strict = FALSE),
                   error = function(e) NULL)
        } else {
          updateNavbarPage(session = session, inputId = "analyst", selected = tb)
        }
      }
    })
    ####

    output$optTabs <- renderUI({
      titleTabs <- list(
        title = tagList(
          tags$span(class = "omicsviewer-nav-title", "Analysis"),
          tags$span(class = "omicsviewer-nav-sub", " \u2014 analyse the current selection on the fly")
        ),
        windowTitle = "Analysis",
        id = ns("analyst"),
        theme = shinytheme("spacelab"),
        tabPanel(
          tags$span(icon("chart-column"), "Feature", title = "Selected-feature analysis.\nCompares the abundance of the currently selected feature(s) across sample groups:\nbeeswarm/box plots with group-comparison statistics (t-test or Wilcoxon for two groups,\nANOVA/Kruskal-Wallis for more), or ROC / precision-recall curves when a binary\noutcome is chosen.\nInterpretation: the box spans the interquartile range with the median as center line.\nA small p-value (e.g. < 0.05) indicates a statistically significant group difference -\nalways judge effect size and overlap on the plot too. For ROC/PR curves: AUC near 1.0\nmeans the feature separates the two classes well; 0.5 is no better than chance."),
          value = "Feature",
          tags$h2("Feature Analysis", class = "sr-only", `aria-label` = "Statistical analysis and visualization of selected features including boxplots with group comparisons and ROC curves for binary outcomes"),
          feature_general_ui(ns("feature_general"))
        )
      )
      sampleAnalyst <- list(
        tabPanel(
          tags$span(icon("arrows-down-to-people"), "Sample", title = "Selected-sample analysis.\nAnalyses the currently selected samples against the rest of the cohort:\nbeeswarm plots for continuous variables, contingency tables of categorical\nannotations (chi-square / Fisher's exact test), and Kaplan-Meier survival\ncurves with the log-rank test when time-to-event data are available.\nInterpretation: a small contingency p-value means the selected samples differ in\nthat annotation more than expected by chance. Kaplan-Meier curves show the\nevent-free proportion over time per group; the log-rank p-value tests whether the\ncurves differ. Check group sizes (n) - small groups make tests unreliable."),
          value = "Sample",
          tags$h2("Sample Analysis", class = "sr-only", `aria-label` = "Statistical analysis of selected samples including group comparisons, contingency tables, and Kaplan-Meier survival curves"),
          sample_general_ui(ns("sample_general"))
        )
      )

      ### geneshot
      geneshot <- list(
        tabPanel(
          tags$span(icon("book-open"), "Geneshot", title = "Literature-based gene discovery.\nSearches PubMed-associated literature (AutoRIF) for genes co-mentioned with your\nsearch term(s); separate multiple terms with ';'. Useful for hypothesis generation.\nInterpretation: the scatter plot shows publication count (x) versus the share of a\ngene's publications mentioning your term (y). Genes in the upper right are both\nwell studied and specifically associated with your topic. A '+' marks genes that\nare present in your current feature selection. Ranking reflects literature bias,\nnot statistical evidence from your data."),
          value = "Geneshot",
          tags$h2("Geneshot Literature Search", class = "sr-only", `aria-label` = "Literature-based gene discovery using PubMed co-mention analysis to find genes associated with search terms"),
          geneshot_ui(ns("geneshotTab"))
        )
      )
      ### end

      optionalTabs <- list()

      if (!is.null(attr(reactive_featureData(), "GS"))) {
        optionalTabs <- c(optionalTabs, list(tabPanel(
          tags$span(icon("object-group"), "ORA", title = "Over-representation analysis (hypergeometric test).\nTests whether your currently selected features belong to a gene set (pathway, GO\ncategory, custom set) more often than expected by chance. Needs a well-defined\nselection, e.g. significant hits from the Feature tab.\nInterpretation: the p-value is the chance probability of the observed overlap; FDR\ncorrects for testing many gene sets (FDR < 0.05 is a common significance cutoff);\nthe odds ratio / gene ratio measures enrichment strength. The second table lists\nthe overlapping features per gene set. Results depend on selection size and the\nbackground gene set collection."),
          value = "ORA",
          tags$h2("Over-Representation Analysis", class = "sr-only", `aria-label` = "Gene set over-representation analysis using hypergeometric test to identify enriched pathways and functional categories"),
          enrichment_analysis_ui(ns("ora"))
        )))
        optionalTabs <- c(optionalTabs, list(tabPanel(
          tags$span(icon("arrow-down-short-wide"), "fGSEA", title = "Fast gene set enrichment analysis (GSEA on a ranked list).\nUnlike ORA, fGSEA uses ALL features, ranked by a statistic you choose (e.g. fold\nchange or signed p-value), and asks whether a gene set's members cluster at the\ntop or bottom of that ranking - no hard significance cutoff needed.\nInterpretation: NES (normalized enrichment score) > 0 means the set is enriched\namong up-ranked features, NES < 0 depleted/enriched at the other end; padj < 0.05\nis a common significance threshold; the leading edge lists the members driving\nthe enrichment signal."),
          value = "fGSEA",
          tags$h2("Fast Gene Set Enrichment Analysis", class = "sr-only", `aria-label` = "Ranked gene set enrichment analysis computing normalized enrichment scores and identifying leading edge genes"),
          enrichment_fgsea_ui(ns("fgsea"))
        )))
      }

      if (any(grepl("^ResponseCurve\\|", colnames(reactive_featureData())))) {
        optionalTabs <- c(optionalTabs, list(tabPanel(
          tags$span(icon("capsules"), "Response", title = "Dose-response curves (4-parameter logistic fit).\nFits abundance vs concentration for the selected features and estimates EC50/IC50\n(concentration giving half-maximal effect), Hill slope (curve steepness) and the\ntop/bottom response plateaus.\nInterpretation: a lower EC50/IC50 means higher potency; the Hill slope shows how\nabruptly the response switches; the parameter table gives fitted values with\nconfidence intervals. Points are observed data, the line the model fit."),
          value = "Response",
          tags$h2("Dose-Response Curves", class = "sr-only", `aria-label` = "Dose-response curve fitting with EC50 and IC50 estimation using 4-parameter logistic regression model"),
          dose_response_ui(ns("rescurve"))
        )))
      }

      if (any(grepl("^StringDB\\|", colnames(reactive_featureData())))) {
        optionalTabs <- c(optionalTabs, list(tabPanel(
          tags$span(icon("share-nodes"), "StringDB", title = "Protein-protein interaction network (STRING database).\nQueries STRING for known and predicted interactions among the currently selected\nfeatures (max ~300; set the taxonomy code, e.g. 9606 = human) and draws the\nnetwork - nodes are proteins, edges are interactions. Requires internet access.\nInterpretation: edge confidence reflects combined evidence (experiments,\ndatabases, text mining); dense clusters often correspond to functional modules or\ncomplexes; highly connected hub proteins may be key regulators."),
          value = "StringDB",
          tags$h2("STRING Protein Interaction Network", class = "sr-only", `aria-label` = "Protein-protein interaction network from STRING database showing physical and functional associations"),
          string_ui(ns("stringdb"))
        )))
      }

      if (any(grepl("^SeqLogo\\|", colnames(reactive_featureData())))) {
        optionalTabs <- c(optionalTabs, list(tabPanel(
          tags$span(icon("font"), "SeqLogo", title = "PTM motif enrichment (sequence logo).\nTakes the amino-acid sequence windows around the selected modification sites and\ntests which residues are over-represented at each position, shown as a sequence\nlogo around the modification site (position 0).\nInterpretation: taller letters mean the amino acid (or property class) is more\nconserved/enriched at that position; significant enrichment versus all sites in\nthe dataset hints at the kinase or recognition machinery involved."),
          value = "SeqLogo",
          tags$h2("Post-Translational Modification Motifs", class = "sr-only", `aria-label` = "PTM motif enrichment analysis identifying overrepresented amino acid sequence patterns around modification sites"),
          ptmotif_ui(ns("ptm"))
        )))
      }

      ######
      if (length(additionalTabs) > 0) {
        for (lo in additionalTabs) {
          optionalTabs <- c(optionalTabs, list(tabPanel(lo$tabName, lo$moduleUi(ns(lo$moduleName)))))
        }
      }
      do.call(navbarPage, c(titleTabs, optionalTabs, geneshot, sampleAnalyst))
    })

    reactive({
      additional_state <- lapply(additional_module_results, function(x) {
        if (is.null(x))
          return(NULL)
        tryCatch({
          value <- if (is.function(x)) x() else x
          if (is.list(value) && !is.null(value$state))
            value$state
          else
            attr(value, "status")
        }, error = function(e) NULL)
      })
      if (length(additionalTabs) > 0) {
        for (lo in additionalTabs) {
          if (isTRUE(lo$stateful) && !is.null(additional_state[[lo$moduleName]]))
            next
          additional_state[[lo$moduleName]] <- list(state_support = FALSE)
        }
      }

      # WP-11: the status builder pulls every stateful panel's value from
      # the registry loop's result map (keyed by status_key). The key
      # EMISSION order below is pinned to the historical literal order -
      # analyst_active_tab first, analyst_additional last, and the
      # panels in between in today's order (feature_general,
      # sample_general, gene_shot, fgsea, ora, ptm) - which is NOT the
      # registry order; the snapshot panels$result_space name/order
      # parity is asserted by the WP-11 goldens (test_panelRegistry.R).
      # analyst_dose_response keeps its state_support marker (WP-3; the
      # plot is a computed payload) and analyst_stringdb stays unemitted
      # (WP-4; widget state rides resultspace.stringdb.*).
      .rs_status_emission <- c(
        "analyst_feature_general", "analyst_sample_general",
        "analyst_gene_shot", "analyst_fgsea", "analyst_ora",
        "analyst_ptm")
      out <- list(analyst_active_tab = input$analyst)
      for (.rs_key in .rs_status_emission)
        out[[.rs_key]] <- safe_module_state(.rs_panel_results[[.rs_key]])
      out$analyst_dose_response <- list(state_support = FALSE)
      out$analyst_additional <- additional_state
      # WP-7 (F4): bounded enrichment-result reader for the assistant's
      # panel_results state section (attribute - never serialized into
      # snapshots; the same pattern as the data-space quickViews attr).
      attr(out, "read_enrichment_results") <- read_enrichment_results
      out
    })
  }) # end moduleServer
}
