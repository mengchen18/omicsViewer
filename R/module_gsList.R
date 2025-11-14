#' @description Utility enrichment analysis shiny ui
#' @param id id
#' @importFrom DT dataTableOutput
gslist_ui <- function(id) {
  ns <- NS(id)
  dataTableDownload_ui(ns("stab"))
}

#' @description Utility gene set list shiny module
#' @param id module id
#' @param reactive_featureData reactive feature data
#' @param reactive_i reactive index of rows to be selected
#' @importFrom fastmatch fmatch

gslist_module <- function(
  id, reactive_featureData, reactive_i
) {

  moduleServer(id, function(input, output, session) {

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
    if (length(reactive_i()) == 0 || all(is.na(reactive_i())))
      return(reactive_pathway())
    df <- reactive_pathway()[reactive_pathway()$featureId %fin% reactive_i(), ]
    req(is.data.frame(df))
    df
    })
  
  ii <- dataTableDownload_module(
    "stab", reactive_table = reactive({
      tab()[, setdiff(colnames(tab()), "featureId")]
      }), prefix = "gslist_", pageLength = DEFAULT_TABLE_PAGE_LENGTH_LARGE
  )

  reactive({
    req(ii())
    as.character( tab()$featureId[ii()] )
    })

  }) # end moduleServer
}
