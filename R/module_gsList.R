#' Gene Set List UI Function
#'
#' @description
#' Creates the user interface for the gene set membership display module.
#' Shows which gene sets contain the selected features with downloadable results.
#'
#' @param id Character. Namespace ID for the Shiny module. Must match the ID
#'   used in \code{\link{gslist_module}}.
#'
#' @return
#' A \code{tagList} containing a data table with download functionality showing
#' gene set memberships for selected features.
#'
#' @family enrichment modules
#' @seealso
#' \code{\link{gslist_module}} for the corresponding server logic.
#'
#' @keywords internal
#' @importFrom DT dataTableOutput
#'
gslist_ui <- function(id) {
  ns <- NS(id)
  dataTableDownload_ui(ns("stab"))
}

#' Gene Set List Server Function
#'
#' @description
#' Server logic for the gene set membership display module. Retrieves and
#' displays gene set annotations for selected features from the feature data
#' attributes.
#'
#' @param id Character. Namespace ID for the Shiny module. Must match the ID
#'   used in \code{\link{gslist_ui}}.
#'
#' @param reactive_featureData Reactive expression. Returns a data.frame of
#'   feature metadata with a "GS" attribute containing gene set membership
#'   information. The "GS" attribute should be a data.frame with columns:
#'   \itemize{
#'     \item featureId: Feature identifiers
#'     \item gsId: Gene set identifiers
#'     \item Additional annotation columns (optional)
#'   }
#'
#' @param reactive_i Reactive expression. Returns feature IDs or indices to
#'   display. Can be:
#'   \itemize{
#'     \item Character vector of feature IDs
#'     \item Integer vector of feature indices
#'     \item Logical scalar TRUE (show all features)
#'     \item NULL or NA (show all features)
#'   }
#'
#' @details
#' The module extracts gene set annotations from the "GS" attribute of
#' feature data and filters to show only selected features. Gene sets are
#' displayed with associated feature annotations from columns starting with
#' "General|".
#'
#' @return
#' A reactive expression returning a character vector of feature IDs for
#' the selected table row, or NULL if no row is selected.
#'
#' @family enrichment modules
#' @seealso
#' \code{\link{gslist_ui}} for the corresponding UI function.
#'
#' @keywords internal
#' @importFrom fastmatch fmatch
#'
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
