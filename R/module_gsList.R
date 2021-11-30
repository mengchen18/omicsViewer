#' @description Utility enrichment analysis shiny ui
#' @param id id
#' @importFrom DT dataTableOutput
gslist_ui <- function(id) {
  ns <- NS(id)
  dataTableDownload_ui(ns("stab"))
}

#' @description Utility enrichment analysis shiny module
#' @param input input
#' @param output output
#' @param session session
#' @param reactive_i reactive index of rows to be selected (for ORA)
#' @param reactive_featureData reactive feature data
#' @importFrom fastmatch fmatch

gslist_module <- function(
  input, output, session, reactive_featureData, reactive_i
) {
  
  ns <- session$ns
  
  reactive_pathway <- reactive({    
    req(f1 <- reactive_featureData())
    gss <- attr(f1, "GS")
    s <- cbind(f1[fmatch(gss$featureId, rownames(f1)), grep("^General", colnames(f1))], gss)
    colnames(s)[colnames(s) == "gsId"] <- "Gene-set"
    s
  })
  
  tab <- reactive({
    req(reactive_pathway())
    ic <- setdiff(colnames(reactive_pathway()), "featureId")
    if (length(reactive_i()) == 0 || is.na(reactive_i()))
      return(reactive_pathway()[, ic])
    df <- reactive_pathway()[reactive_pathway()$featureId %fin% reactive_i(), ic]
    req(is.data.frame(df))
    df
    })
  
  callModule(
    dataTableDownload_module, id = "stab", reactive_table = tab, prefix = "gslist_", pageLength = 25
  )
}
