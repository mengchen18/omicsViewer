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
    req(gss)
    s <- cbind(gss, f1[fmatch(gss$featureId, rownames(f1)), grep("^General", colnames(f1)), drop = FALSE])
    colnames(s)[colnames(s) == "gsId"] <- "Gene-set"
    s
  })
  
  tab <- reactive({
    req(reactive_pathway())
    if (length(reactive_i()) == 1 && is.logical(reactive_i()) && reactive_i())
      return(reactive_pathway())    
    if (length(reactive_i()) == 0 || is.na(reactive_i()))
      return(reactive_pathway())
    df <- reactive_pathway()[reactive_pathway()$featureId %fin% reactive_i(), ]
    req(is.data.frame(df))
    df
    })
  
  ii <- callModule(
    dataTableDownload_module, id = "stab", reactive_table = reactive({
      tab()[, setdiff(colnames(tab()), "featureId")]
      }), prefix = "gslist_", pageLength = 25
  )

  reactive({    
    req(ii())
    as.character( tab()$featureId[ii()] )
    })
}
