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
#' @param reactive_status Reactive gene-set table state used for snapshot
#'   restoration.
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
  id, reactive_featureData, reactive_i, reactive_status = reactive(NULL)
) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns
  
  reactive_pathway <- reactive({
    # R-M9: the gene-set long table (up to ~1M rows) must only be built
    # when the panel is actually visible; the visibility flip re-runs this
    # reactive (clientData dependency) before the table paints. Headless
    # sessions report no flag and stay visible.
    if (!output_visible(session, ns("stab-table")))
      return(NULL)
    req(f1 <- reactive_featureData())
    gss <- attr(f1, "GS")
    req(gss)
    kp <- strsplit(colnames(f1), split = "\\|")[[1]][[1]]
    kp <- sprintf("^%s\\|", kp)
    s <- cbind(gss, f1[fmatch(gss$featureId, rownames(f1)), grep(kp, colnames(f1)), drop = FALSE])
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
      }), tab_status = reactive_status, prefix = "gslist_",
    pageLength = DEFAULT_TABLE_PAGE_LENGTH_LARGE,
    reactive_row_ids = reactive({
      if (is.null(tab()) || nrow(tab()) == 0) return(NULL)
      # R-M8: gsId was renamed to "Gene-set" in reactive_pathway(); the
      # old tab()$gsId was NULL and produced garbage ids
      paste(tab()[["Gene-set"]], tab()$featureId, sep = "|")
    })
  )

  reactive({
    v <- ii()
    req(v)
    out <- as.character( tab()$featureId[v] )
    # R-M8: carry the DT status so the snapshot save path
    # (eset_gslist_tab) is not dead
    attr(out, "status") <- attr(v, "status")
    out
    })

  }) # end moduleServer
}
