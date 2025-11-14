#' Dose Response Curve UI Function
#'
#' @description
#' Creates the user interface for the dose response curve visualization module.
#' Displays dose-response curves with fitted parameters and feature information
#' in a multi-panel layout.
#'
#' @param id Character. Namespace ID for the Shiny module. Must match the ID
#'   used in \code{\link{dose_response_module}}.
#'
#' @return
#' A \code{fluidRow} containing:
#' \itemize{
#'   \item Plot panel: Dose-response curve visualization
#'   \item Parameter table: Fitted curve parameters (hill slope, EC50, etc.)
#'   \item Feature table: Selected feature information
#' }
#'
#' @family analysis modules
#' @seealso
#' \code{\link{dose_response_module}} for the corresponding server logic.
#' \code{\link{plotDCMat}} for the underlying plotting function.
#'
#' @keywords internal
#' @importFrom DT dataTableOutput
dose_response_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(
      7, plotOutput(ns("plot"), height = DOSE_RESPONSE_PLOT_HEIGHT)
    ),
    column(
      5, dataTableDownload_ui(ns("param")), style = paste0("margin-top: ", DOSE_RESPONSE_TABLE_MARGIN_TOP, ";")
    ),
    column(
      12, dataTableDownload_ui(ns("feature"))
    )
  )
}

#' Dose Response Curve Server Function
#'
#' @description
#' Server logic for the dose response curve visualization module. Fits
#' dose-response curves using the drc package and displays fitted parameters
#' and plots for selected features.
#'
#' @param id Character. Namespace ID for the Shiny module. Must match the ID
#'   used in \code{\link{dose_response_ui}}.
#'
#' @param reactive_expr Reactive expression. Returns a numeric matrix with
#'   features as rows and samples as columns containing expression/response values.
#'
#' @param reactive_phenoData Reactive expression. Returns a data.frame of
#'   sample metadata including dose and curve ID columns.
#'
#' @param reactive_featureData Reactive expression. Returns a data.frame of
#'   feature metadata with features as rows.
#'
#' @param reactive_i Reactive expression. Returns a single integer index
#'   indicating which feature row to display. Must have length 1.
#'
#' @param reactive_attr_drc Reactive expression. Returns a named character vector
#'   with elements "dose_col" and "curveid_col" specifying the column names
#'   in phenoData for dose values and curve identifiers.
#'
#' @details
#' The module fits 4-parameter or 5-parameter log-logistic curves to the
#' dose-response data and extracts parameters including EC50, hill slope,
#' minimum and maximum responses. Only single feature selection is supported.
#'
#' @return
#' NULL (invisibly). The module displays plots and tables but does not return
#' a reactive value.
#'
#' @family analysis modules
#' @seealso
#' \code{\link{dose_response_ui}} for the corresponding UI function.
#' \code{\link{plotDCMat}} for the plotting function.
#' \code{\link{drmMat}} for curve fitting.
#'
#' @keywords internal
#' @importFrom fastmatch fmatch

dose_response_module <- function(
    id,
    reactive_expr,
    reactive_phenoData,
    reactive_featureData,
    reactive_i,
    reactive_attr_drc
) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns
  
  dr1 <- reactive({
    req(length(reactive_i()) == 1)
    req(reactive_featureData())
    req(reactive_attr_drc())
    v <- reactive_attr_drc()
    list(
      col_dose = paste("General", "All", v["dose_col"], sep = "|"),
      col_curveid = paste("General", "All", v["curveid_col"], sep = "|")
      )
  })
  
  tabs <- reactive({
    req(dr1())
    plotDCMat(
      expr = reactive_expr(), 
      pd = reactive_phenoData(), 
      fd = reactive_featureData(), 
      featid = reactive_i(),
      dose.var = dr1()$col_dose, 
      curve.var = dr1()$col_curveid,
      only.par = TRUE
    )
  })
  
  output$plot <- renderPlot({
    req(dr1())
    plotDCMat(
      expr = reactive_expr(), 
      pd = reactive_phenoData(), 
      fd = reactive_featureData(), 
      featid = reactive_i(),
      dose.var = dr1()$col_dose, 
      curve.var = dr1()$col_curveid
    )
  })
  
  dataTableDownload_module(
    "feature", reactive_table = reactive(tabs()$featInfo), prefix = "Response_featureInfo_"
  )
  dataTableDownload_module(
    "param", reactive_table = reactive(tabs()$par), prefix = "Response_par_"
  )

  }) # end moduleServer
}
