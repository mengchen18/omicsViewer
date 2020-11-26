
#' Start ExpressionSetViewer
#' @param dir directory to the ExpressionSet object. Only give the directoy
#'  in this argument, not the .rds file.
#' @export
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @examples
#' 1
#' ## To start the shiny app: 
#' # ExpressionSetViewer(
#' # 	"/media/share_baybioms/Projects/008_Bioinformatics/B032_ExpressionSetViewer/Dat/"
#' # )
#' 
ExpressionSetViewer <- function(dir) {
  app <- list(
    ui = fluidPage(
      app_ui("app")
    ),
    server = function(input, output, session) {
      callModule(app_module, id = "app", dir = reactive(dir))
    }
  )  
  runApp(app)
}


