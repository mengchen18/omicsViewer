
#' Launch the omicsViewer Shiny Application
#'
#' @description
#' Starts an interactive Shiny application for exploring omics data, including visualization
#' of expression matrices, feature and sample metadata, statistical analyses, and functional
#' enrichment results. The viewer supports both \code{ExpressionSet} and
#' \code{SummarizedExperiment} objects.
#'
#' @param dir Character. Path to directory containing the \code{ExpressionSet} or
#'   \code{SummarizedExperiment} object saved as .RDS file. Provide only the directory path,
#'   not the full file path. The viewer will list all compatible files in this directory.
#' @param additionalTabs List. Optional custom Shiny modules to add as tabs in the "Analyst" panel.
#'   Each element should be a list with: \code{tabName} (character), \code{moduleName} (character),
#'   \code{moduleUi} (UI function), and \code{moduleServer} (server function).
#' @param filePattern Character. Regular expression pattern to filter files displayed in the
#'   directory. Default: \code{".(RDS|DB|SQLITE|SQLITE3)$"} (case-insensitive).
#' @param ESVObj ExpressionSet or SummarizedExperiment. Optional pre-loaded object to view
#'   directly without file selection. If provided, the file dropdown will show "ESVObj.RDS".
#' @param esetLoader Function. Custom loader for reading saved objects. Default: \code{readESVObj}.
#'   Should accept a file path and return an ExpressionSet or SummarizedExperiment object.
#' @param exprsGetter Function. Extracts expression matrix from loaded object.
#'   Default: \code{getExprs}. Should return a numeric matrix.
#' @param pDataGetter Function. Extracts phenotype/sample metadata. Default: \code{getPData}.
#'   Should return a data.frame with rownames matching sample names.
#' @param fDataGetter Function. Extracts feature metadata. Default: \code{getFData}.
#'   Should return a data.frame with rownames matching feature names.
#' @param defaultAxisGetter Function. Determines default axes for plots. Takes two arguments:
#'   \code{x} (the loaded object) and \code{what} (one of "sx", "sy", "fx", "fy" for
#'   sample/feature space x/y axes). Should return column name from metadata.
#' @param appName Character. Application title displayed in the UI. Default: "omicsViewer".
#' @param appVersion Character or package_version. Version number displayed in UI.
#'   Default: current package version.
#'
#' @export
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @examples
#' if (interactive()) {
#'   # Basic usage with example data
#'   omicsViewer(system.file("extdata", package = "omicsViewer"))
#'
#'   # With pre-loaded object
#'   packdir <- system.file("extdata", package = "omicsViewer")
#'   eset <- readRDS(file.path(packdir, "exampleEset.RDS"))
#'   omicsViewer(packdir, ESVObj = eset)
#' }
#'
#' @return NULL (invisibly). Launches the Shiny application. The app runs until stopped by the user.
#'
#' @seealso
#' \code{\link{prepOmicsViewer}} for preparing data objects for visualization.
#' \code{\link{app_module}} for the main application module (developers only).
#' 
omicsViewer <- function(
  dir, additionalTabs = NULL, filePattern = ".(RDS|DB|SQLITE|SQLITE3)$", ESVObj = NULL,
  esetLoader = readESVObj, 
  exprsGetter = getExprs, pDataGetter = getPData, fDataGetter = getFData, 
  defaultAxisGetter = getAx,
  appName = "omicsViewer", appVersion = packageVersion("omicsViewer")
  ) {
  
  app <- list(
    ui = fluidPage(
      app_ui("app")
    ),
    server = function(input, output, session, aTabs = additionalTabs,
                      f_eset = esetLoader, f_exprs = exprsGetter, f_pd = pDataGetter, f_fd = fDataGetter,
                      axg = defaultAxisGetter) {
      app_module(
        "app", .dir = reactive(dir), additionalTabs = aTabs, filePattern = filePattern,
        esetLoader = f_eset, exprsGetter = f_exprs, pDataGetter = f_pd, fDataGetter = f_fd,
        defaultAxisGetter = axg, appName = appName, appVersion = appVersion, ESVObj = reactive(ESVObj)
        )
    }
  )
  runApp(app)
}


