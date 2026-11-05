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
  tagList(
    # Module description for AI browsers and screen readers
    div(class = "sr-only", id = ns("module-help"),
      tags$h4("About Dose-Response Curve Analysis"),
      tags$p("Dose-response analysis models the relationship between drug/compound concentration and biological response using a 4-parameter logistic curve. This analysis fits sigmoid curves to your data and estimates key pharmacological parameters including EC50 (half-maximal effective concentration), IC50 (half-maximal inhibitory concentration), Hill slope, and maximum/minimum response values."),
      tags$h4("When to use dose-response analysis"),
      tags$p("Use this analysis when you have measured responses across multiple concentrations of a compound or treatment. This is essential for drug discovery, toxicology studies, and understanding concentration-dependent biological effects. It helps determine the potency of compounds and compare their effectiveness."),
      tags$h4("How to interpret results"),
      tags$p("The curve plot shows the fitted sigmoid relationship between concentration (x-axis, usually log-scale) and response (y-axis). EC50/IC50 represents the concentration producing 50% of maximum effect - lower values indicate higher potency. The Hill slope indicates the steepness of the curve. The parameter table shows fitted values with their confidence intervals. Top and bottom parameters represent the maximum and minimum response plateaus.")
    ),
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
#'   in phenoData for dose values and curve identifiers. Defaults to an empty
#'   reactive (the module mounts idle); since WP-3 the L1 call site passes no
#'   explicit derivation - the attribute is derived inside the module from
#'   \code{context$data$object} instead.
#' @param store Optional child view of the canonical widget store
#'   (\code{\link{widget_store_child}}) for this module, e.g.
#'   \code{widget_store_child(store, "resultspace.dose_response")}. When
#'   given (directly or through \code{context}), the parameter and feature
#'   tables register their selected-row bindings on it (\code{param_table} /
#'   \code{feature_table}; the string \code{selected_row} pattern); when NULL
#'   the tables keep the legacy store-less behavior (WP-3, P4).
#' @param context Optional panel context bundle (\code{\link{panel_context}};
#'   plan D13). When given, \code{context$store} takes over the store wiring
#'   and the \code{S6.6_drc} attribute derivation is computed from
#'   \code{context$data$object} (moved in verbatim from the L1 call site in
#'   WP-3). All legacy formals keep working for direct instantiation.
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
    reactive_attr_drc = shiny::reactive(NULL),
    store = NULL,
    context = NULL
) {

  moduleServer(id, function(input, output, session) {

  ns <- session$ns

  # WP-3 (P4): the panel joins the widget store. Resolution per plan
  # D13: the context (when given) owns the store; the explicit formal is
  # the legacy / direct-instantiation path.
  if (!is.null(context))
    store <- context$store %||% store

  # WP-3 (P4): the S6.6_drc attribute derivation moved in VERBATIM from
  # the L1 call site (plan WP-3 risk-table rule: identical expressions
  # moved across the module boundary). The source is the context's
  # object reactive; without a context the legacy reactive_attr_drc
  # formal is used (direct instantiation).
  .dose_object <- if (is.null(context)) NULL else context$data$object
  attr_drc <- if (is.null(.dose_object)) reactive_attr_drc else reactive({
    req(.dose_object())
    # sqlite-revival WP-3: S6.6_drc rides through getAx (the default
    # seam) so SQLite connections serve it from the attrs table;
    # containers without a getAx method fall back to plain attr()
    # (the pre-2.1.81 behavior)
    tryCatch(getAx(.dose_object(), "S6.6_drc"),
             error = function(e) attr(.dose_object(), "S6.6_drc"))
  })
  
  dr1 <- reactive({
    req(length(reactive_i()) == 1)
    req(reactive_featureData())
    req(attr_drc())
    v <- attr_drc()
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
    "feature", reactive_table = reactive(tabs()$featInfo), prefix = "Response_featureInfo_",
    reactive_row_ids = reactive({
      fi <- tryCatch(tabs()$featInfo,
                     shiny.silent.error = function(e) NULL,
                     error = function(e) NULL)
      if (is.data.frame(fi) && nrow(fi) > 0L) rownames(fi) else NULL
    }),
    store = store, store_key = "feature_table",
    store_label = "Selected feature row",
    store_help = paste("Feature whose row is selected in the dose-response",
                       "feature-information table; the table follows the",
                       "current single-feature selection. PRECONDITION:",
                       "requires exactly one feature selected and",
                       "dose-response metadata (S6.6_drc) on the dataset")
  )
  dataTableDownload_module(
    "param", reactive_table = reactive(tabs()$par), prefix = "Response_par_",
    reactive_row_ids = reactive({
      pt <- tryCatch(tabs()$par,
                     shiny.silent.error = function(e) NULL,
                     error = function(e) NULL)
      if (is.data.frame(pt) && "param" %in% colnames(pt))
        as.character(pt$param) else NULL
    }),
    store = store, store_key = "param_table",
    store_label = "Selected parameter",
    store_help = paste("Parameter row currently selected in the",
                       "dose-response parameter table (hill slope, response",
                       "plateaus, inflection point, fit statistics).",
                       "PRECONDITION: requires exactly one feature selected",
                       "and dose-response metadata (S6.6_drc) on the dataset")
  )

  }) # end moduleServer
}
