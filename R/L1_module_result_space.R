#' @title Utility L1 result space ui
#' @param id id
#' @importFrom shinythemes shinytheme
L1_result_space_ui <- function(id) {
  ns <- NS(id)

  tagList(
    uiOutput(ns("optTabs"))
  )
}

#' @title Utility L1 result space ui
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
#' @export
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
    # result-space child stores plus the analyst tab navbar itself, so the
    # agent can switch analysis tabs exactly like a user.
    # ------------------------------------------------------------------
    store_feature_general <- NULL
    store_sample_general <- NULL
    store_ora <- NULL
    store_fgsea <- NULL
    store_stringdb <- NULL
    store_ptm <- NULL
    store_geneshot <- NULL
    if (!is.null(store)) {
      store_feature_general <- widget_store_child(store, "resultspace.feature_general")
      store_sample_general <- widget_store_child(store, "resultspace.sample_general")
      store_ora <- widget_store_child(store, "resultspace.ora")
      store_fgsea <- widget_store_child(store, "resultspace.fgsea")
      store_stringdb <- widget_store_child(store, "resultspace.stringdb")
      store_ptm <- widget_store_child(store, "resultspace.ptm")
      store_geneshot <- widget_store_child(store, "resultspace.geneshot")
      store_rs <- widget_store_child(store, "resultspace")
      # tab titles in renderUI order; providers must never req()
      .rs_tab_choices <- function() {
        fd <- tryCatch(reactive_featureData(),
                       shiny.silent.error = function(e) NULL,
                       error = function(e) NULL)
        tabs <- "Feature"
        if (!is.null(fd)) {
          if (!is.null(attr(fd, "GS")))
            tabs <- c(tabs, "ORA", "fGSEA")
          if (any(grepl("^ResponseCurve\\|", colnames(fd))))
            tabs <- c(tabs, "Response")
          if (any(grepl("^StringDB\\|", colnames(fd))))
            tabs <- c(tabs, "StringDB")
          if (any(grepl("^SeqLogo\\|", colnames(fd))))
            tabs <- c(tabs, "SeqLogo")
        }
        if (length(additionalTabs) > 0)
          tabs <- c(tabs, vapply(additionalTabs,
                                 function(lo) lo$tabName, character(1)))
        c(tabs, "Geneshot", "Sample")
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

    # session restore finished
    v <- feature_general_module("feature_general",
      reactive_expr = reactive_expr,
      reactive_i = reactive_i,
      reactive_highlight = reactive_highlight,
      reactive_phenoData = reactive_phenoData,
      reactive_featureData = reactive_featureData,
      reactive_status = reactive(status()$analyst_feature_general),
      store = store_feature_general
    )

    # session restore finished
    v2 <- enrichment_fgsea_module("fgsea",
      reactive_featureData = reactive_featureData,
      reactive_status = reactive(status()$analyst_fgsea),
      store = store_fgsea
    )

    # session restore finished
    v3 <- enrichment_analysis_module("ora",
      reactive_i = reactive_i, reactive_featureData = reactive_featureData,
      reactive_status = reactive(status()$analyst_ora),
      store = store_ora
    )

    # session restore finished
    v4 <- string_module(
      "stringdb",
      reactive_ids = reactive({
        i <- grep("^StringDB\\|", colnames(reactive_featureData()))
        reactive_featureData()[reactive_i(), i[1]]
      }), reactive_status = reactive(NULL),
      active = reactive(status()$analyst_active_tab == "StringDB"),
      store = store_stringdb
    )

    # session restore finished
    v5 <- sample_general_module(
      "sample_general",
      reactive_phenoData = reactive_phenoData, reactive_expr = reactive_expr,
      reactive_j = reactive_highlight,
      reactive_status = reactive(status()$analyst_sample_general),
      store = store_sample_general
    )

    # session restore finished
    v6 <- geneshot_module(
      "geneshotTab",
      fdata = reactive_featureData,
      feature_selected = reactive_i,
      reactive_status = reactive(status()$analyst_gene_shot),
      store = store_geneshot
    )

    # session restore finished
    v7 <- ptmotif_module(
      "ptm",
      fdata = reactive_featureData,
      feature_selected = reactive_i,
      reactive_status = reactive(status()$analyst_ptm),
      store = store_ptm
      # ,
      # background = reactive( attr(object(), "ptm.seq.window") )
      # background = reactive({
      #   i <- grep("^PTMSeq\\|", colnames(reactive_featureData()))
      #   reactive_featureData()[[i]]
      #   })
    )

    # dose response
    v8 <- dose_response_module(
      "rescurve",
      reactive_expr = reactive_expr,
      reactive_i = reactive_i,
      reactive_phenoData = reactive_phenoData,
      reactive_featureData = reactive_featureData,
      reactive_attr_drc = reactive({
        req(object())
        attr(object(), "S6.6_drc")
      })
    )

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
        title = "Analysis", id = ns("analyst"),
        theme = shinytheme("spacelab"),
        tabPanel(
          "Feature",
          tags$h2("Feature Analysis", class = "sr-only", `aria-label` = "Statistical analysis and visualization of selected features including boxplots with group comparisons and ROC curves for binary outcomes"),
          feature_general_ui(ns("feature_general"))
        )
      )
      sampleAnalyst <- list(
        tabPanel(
          "Sample",
          tags$h2("Sample Analysis", class = "sr-only", `aria-label` = "Statistical analysis of selected samples including group comparisons, contingency tables, and Kaplan-Meier survival curves"),
          sample_general_ui(ns("sample_general"))
        )
      )

      ### geneshot
      geneshot <- list(
        tabPanel(
          "Geneshot",
          tags$h2("Geneshot Literature Search", class = "sr-only", `aria-label` = "Literature-based gene discovery using PubMed co-mention analysis to find genes associated with search terms"),
          geneshot_ui(ns("geneshotTab"))
        )
      )
      ### end

      optionalTabs <- list()

      if (!is.null(attr(reactive_featureData(), "GS"))) {
        optionalTabs <- c(optionalTabs, list(tabPanel(
          "ORA",
          tags$h2("Over-Representation Analysis", class = "sr-only", `aria-label` = "Gene set over-representation analysis using hypergeometric test to identify enriched pathways and functional categories"),
          enrichment_analysis_ui(ns("ora"))
        )))
        optionalTabs <- c(optionalTabs, list(tabPanel(
          "fGSEA",
          tags$h2("Fast Gene Set Enrichment Analysis", class = "sr-only", `aria-label` = "Ranked gene set enrichment analysis computing normalized enrichment scores and identifying leading edge genes"),
          enrichment_fgsea_ui(ns("fgsea"))
        )))
      }

      if (any(grepl("^ResponseCurve\\|", colnames(reactive_featureData())))) {
        optionalTabs <- c(optionalTabs, list(tabPanel(
          "Response",
          tags$h2("Dose-Response Curves", class = "sr-only", `aria-label` = "Dose-response curve fitting with EC50 and IC50 estimation using 4-parameter logistic regression model"),
          dose_response_ui(ns("rescurve"))
        )))
      }

      if (any(grepl("^StringDB\\|", colnames(reactive_featureData())))) {
        optionalTabs <- c(optionalTabs, list(tabPanel(
          "StringDB",
          tags$h2("STRING Protein Interaction Network", class = "sr-only", `aria-label` = "Protein-protein interaction network from STRING database showing physical and functional associations"),
          string_ui(ns("stringdb"))
        )))
      }

      if (any(grepl("^SeqLogo\\|", colnames(reactive_featureData())))) {
        optionalTabs <- c(optionalTabs, list(tabPanel(
          "SeqLogo",
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

      list(
        analyst_active_tab = input$analyst,
        analyst_feature_general = safe_module_state(v),
        analyst_sample_general = safe_module_state(v5),
        analyst_gene_shot = safe_module_state(v6),
        analyst_fgsea = safe_module_state(v2),
        analyst_stringdb = list(
          state_support = FALSE,
          gap = APP_STATE_GAPS$stringdb
        ),
        analyst_ora = safe_module_state(v3),
        analyst_ptm = safe_module_state(v7),
        analyst_dose_response = list(state_support = FALSE),
        analyst_additional = additional_state
      )
    })
  }) # end moduleServer
}
