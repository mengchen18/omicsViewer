#' @description Utility L1 result space ui
#' @param id id
#' @importFrom shinythemes shinytheme
L1_result_space_ui <- function(id) {
  ns <- NS(id)

  tagList(
    uiOutput(ns("optTabs"))
  )
}

#' @description Utility L1 result space ui
#' @param id module id
#' @param reactive_expr expression matrix
#' @param reactive_phenoData phentype data
#' @param reactive_featureData feature data
#' @param reactive_i row ID/name of rows selected
#' @param reactive_highlight col ID/name of columns selected
#' @param additionalTabs additional tabs added to "Analyst" panel
#' @param object originally loaded object, mostly an \code{ExpressionSet} or \code{SummarizedExperiment} object
#' @param status intial status
#' @export
L1_result_space_module <- function(
  id,
  reactive_expr, reactive_phenoData, reactive_featureData,
  reactive_i = reactive(NULL),
  reactive_highlight = reactive(NULL),
  additionalTabs = NULL,
  object = NULL, status = reactive(NULL)
) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # session restore finished
    v <- feature_general_module("feature_general",
      reactive_expr = reactive_expr,
      reactive_i = reactive_i,
      reactive_highlight = reactive_highlight,
      reactive_phenoData = reactive_phenoData,
      reactive_featureData = reactive_featureData,
      reactive_status = reactive(status()$analyst_feature_general)
    )

    # session restore finished
    v2 <- enrichment_fgsea_module("fgsea",
      reactive_featureData = reactive_featureData,
      reactive_status = reactive(status()$analyst_fgsea)
    )

    v3 <- enrichment_analysis_module("ora",
      reactive_i = reactive_i, reactive_featureData = reactive_featureData
    )

    # session restore finished
    v4 <- string_module(
      "stringdb",
      reactive_ids = reactive({
        i <- grep("^StringDB\\|", colnames(reactive_featureData()))
        reactive_featureData()[reactive_i(), i[1]]
      }), reactive_status = reactive(status()$analyst_stringdb),
      active = reactive(status()$analyst_active_tab == "StringDB")
    )

    # session restore finished
    v5 <- sample_general_module(
      "sample_general",
      reactive_phenoData = reactive_phenoData, reactive_expr = reactive_expr,
      reactive_j = reactive_highlight,
      reactive_status = reactive(status()$analyst_sample_general)
    )

    # session restore finished
    v6 <- geneshot_module(
      "geneshotTab",
      fdata = reactive_featureData,
      feature_selected = reactive_i,
      reactive_status = reactive(status()$analyst_gene_shot)
    )

    v7 <- ptmotif_module(
      "ptm",
      fdata = reactive_featureData,
      feature_selected = reactive_i
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
    if (length(additionalTabs) > 0) {
      for (lo in additionalTabs) {
        lo$moduleServer(
          lo$moduleName,
          pdata = reactive_phenoData, fdata = reactive_featureData, expr = reactive_expr,
          feature_selected = reactive_i, sample_selected = reactive_highlight, object = object
        )
      }
    }

    #### status for snapshot #####
    observe({
      if (!is.null(tb <- status()$analyst_active_tab)) {
        updateNavbarPage(session = session, inputId = "analyst", selected = tb)
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
      list(
        analyst_active_tab = input$analyst,
        analyst_feature_general = v(),
        analyst_sample_general = v5(),
        analyst_gene_shot = v6(),
        analyst_fgsea = v2(),
        analyst_stringdb = v4()
      )
    })
  }) # end moduleServer
}
