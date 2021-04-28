#' Utility L1 result space ui
#' @param id id
#' @importFrom shinythemes shinytheme
L1_result_space_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    uiOutput(ns("optTabs"))
  )
}

#' Utility L1 result space ui
#' @param input input
#' @param output output
#' @param session session
#' @param reactive_expr expression matrix
#' @param reactive_phenoData phentype data
#' @param reactive_featureData feature data
#' @param reactive_i row ID/name of rows selected
#' @param reactive_highlight col ID/name of columns selected
#' @param additionalTabs additional tabs added to "Analyst" panel
#' @param object originally loaded object, mostly an ExpressionSet object


L1_result_space_module <- function(
  input, output, session, 
  reactive_expr, reactive_phenoData, reactive_featureData,
  reactive_i = reactive(NULL),  
  reactive_highlight = reactive(NULL),
  additionalTabs = NULL,
  object  = NULL
) {
  ns <- session$ns
  
  v <- callModule(feature_general_module, id = "feature_general", 
                  reactive_expr = reactive_expr, 
                  reactive_i = reactive_i,
                  reactive_highlight = reactive_highlight,
                  reactive_phenoData = reactive_phenoData,
                  reactive_featureData = reactive_featureData)
  
  v2 <- callModule(enrichment_fgsea_module, id = "fgsea", 
                   reactive_featureData = reactive_featureData)
  
  v3 <- callModule(enrichment_analysis_module, id = "ora", 
                   reactive_pathway_mat = reactive({
                     fd <- reactive_featureData()
                     fd[, grep("^GS\\|", colnames(fd))]
                   }), reactive_i = reactive_i,
                   reactive_featureData = reactive_featureData
                   )
  
  v4 <- callModule(
    string_module, id = "stringdb", reactive_ids = reactive({
      i <- grep("^StringDB\\|", colnames(reactive_featureData()))
      reactive_featureData()[reactive_i(), i[1]]
    }))
  
  v5 <- callModule(
    sample_general_module, id = "sample_general", 
    reactive_phenoData = reactive_phenoData, 
    reactive_j = reactive_highlight
  )

  v6 <- callModule(
    geneshot_module, id = "geneshotTab", 
    fdata = reactive_featureData, 
    feature_selected = reactive_i
    )
  
  v7 <- callModule(
    ptmotif_module, id = "ptm", 
    fdata = reactive_featureData, 
    feature_selected = reactive_i
  )
  
  # 
  if (length(additionalTabs) > 0) {
    for (lo in additionalTabs){
      callModule(
        lo$moduleServer, id = lo$moduleName,
        pdata = reactive_phenoData, fdata = reactive_featureData, expr = reactive_expr, 
        feature_selected = reactive_i, sample_selected = reactive_highlight, object = object
      )
    }
  }
  
  output$optTabs <- renderUI({
    
    titleTabs <- list(
      title = "Analyst",
      theme = shinytheme("spacelab"), 
      tabPanel("Feature general", feature_general_ui(ns("feature_general")))     
    )
    sampleAnalyst <- list(
      tabPanel("Sample general", sample_general_ui(ns("sample_general")))
      )
    
    ### geneshot
    geneshot <- list(
      tabPanel("Geneshot", geneshot_ui(ns("geneshotTab")))
      )
    ### end
    
    optionalTabs <- list()
    
    if (any(grepl("^GS\\|", colnames(reactive_featureData())))) {
      optionalTabs <- c(optionalTabs, list(tabPanel('ORA', enrichment_analysis_ui(ns("ora")))))
      optionalTabs <- c(optionalTabs, list(tabPanel("fGSEA", enrichment_fgsea_ui(ns("fgsea")))))
    }
    
    if (any(grepl("^StringDB\\|", colnames(reactive_featureData()))))
      optionalTabs <- c(optionalTabs, list(tabPanel("StringDB", string_ui(ns("stringdb")))))
    
    if (any(grepl("^PTMSeq\\|", colnames(reactive_featureData()))))
      optionalTabs <- c(optionalTabs, list(tabPanel("PTM motif", ptmotif_ui(ns("ptm")))))
    
    ######
    if (length(additionalTabs) > 0) {
      for (lo in additionalTabs) {
        optionalTabs <- c( optionalTabs, list(tabPanel(lo$tabName, lo$moduleUi(ns(lo$moduleName)))) )
      }
    }
    ######
    
    do.call(navbarPage, c(titleTabs, optionalTabs, geneshot, sampleAnalyst))
  })
}