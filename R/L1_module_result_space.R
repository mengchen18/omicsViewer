#' Utility L1 result space ui
#' @param id id
#' @importFrom shinythemes shinytheme
L1_result_space_ui <- function(id) {
  ns <- NS(id)
  tagList(
    navbarPage(
      "Analyst",
      theme = shinytheme("spacelab"), 
      tabPanel("Feature general", feature_general_ui(ns("feature_general"))),
      tabPanel("Sample general", sample_general_ui(ns("sample_general"))),      
      tabPanel('ORA', enrichment_analysis_ui(ns("ora"))),
      tabPanel("StringDB", string_ui(ns("stringdb"))),
      tabPanel("fGSEA", enrichment_fgsea_ui(ns("fgsea")))
    )
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


L1_result_space_module <- function(
  input, output, session, 
  reactive_expr, reactive_phenoData, reactive_featureData,
  reactive_i = reactive(NULL),  
  reactive_highlight = reactive(NULL)
) {
  
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
                   }), reactive_i = reactive_i)
  
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
}