
#' Start omicsViewer
#' @param dir directory to the \code{ExpressionSet} or \code{SummarizedExperiment} object. Only give the directory
#'  in this argument, not the .rds file.
#' @param additionalTabs additional tabs added to "Analyst" panel
#' @param filePattern file pattern to be displayed.
#' @param esetLoader function to load the eset object, if an RDS file, should be "readRDS"
#' @param exprsGetter function to get the expression matrix from eset
#' @param pDataGetter function to get the phenotype data from eset
#' @param fDataGetter function to get the feature data from eset
#' @param defaultAxisGetter function to get the default axes to be visualized. It should be a function with two 
#'   arguments: x - the object loaded to the viewer; what - one of "sx", "sy", "fx" and "fy", representing the 
#'   sample space x-axis, sample space y-axis, feature space x-axis and feature space y-axis respectively.
#' @param appName name of the application
#' @param appVersion version of the application
#' @param ESVObj the ESV object
#' @export
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#' 1
#' ## To start the shiny app: 
#' # omicsViewer(
#' #  system.file("extdata", package = "omicsViewer")
#' # )
#' @return do not return values
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
      callModule(
        app_module, id = "app", .dir = reactive(dir), additionalTabs = aTabs, filePattern = filePattern,
        esetLoader = f_eset, exprsGetter = f_exprs, pDataGetter = f_pd, fDataGetter = f_fd,
        defaultAxisGetter = axg, appName = appName, appVersion = appVersion, ESVObj = reactive(ESVObj)
        )
    }
  )
  runApp(app)
}


