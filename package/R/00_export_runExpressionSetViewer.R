
#' Start ExpressionSetViewer
#' @param dir directory to the ExpressionSet object. Only give the directoy
#'  in this argument, not the .rds file.
#' @param additionalTabs additional tabs added to "Analyst" panel
#' @param esetLoader function to load the eset object, if an RDS file, should be "readRDS"
#' @param exprsGetter function to get the expression matrix from eset
#' @param pDataGetter function to get the phenotype data from eset
#' @param fDataGetter function to get the feature data from eset
#' @export
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @examples
#' 1
#' ## To start the shiny app: 
#' # ExpressionSetViewer(
#' # 	ExpressionSetViewer(system.file("extdata", package = "ExpressionSetViewer"))
#' # )
#' 
ExpressionSetViewer <- function(
  dir, additionalTabs = NULL, 
  esetLoader = readRDS, exprsGetter = exprs, pDataGetter = pData, fDataGetter = fData
  ) {
  
  app <- list(
    ui = fluidPage(
      app_ui("app")
    ),
    server = function(input, output, session, aTabs = additionalTabs,
                      f_eset = esetLoader, f_exprs = exprsGetter, f_pd = pDataGetter, f_fd = fDataGetter) {
      callModule(
        app_module, id = "app", dir = reactive(dir), additionalTabs = aTabs,
        esetLoader = f_eset, exprsGetter = f_exprs, pDataGetter = f_pd, fDataGetter = f_fd
        )
    }
  )
  runApp(app)
}


