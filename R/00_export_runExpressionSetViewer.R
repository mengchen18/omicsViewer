
#' Start ExpressionSetViewer
#' @param dir directory to the ExpressionSet object
#' @export
#' @examples
#' 1
#' # v <- testf("/media/share_baybioms/Projects/008_Bioinformatics/B032_ExpressionSetViewer/Dat/")
#' 
runExpressionSetViewer <- function(dir) {
  app <- list(
    ui = fluidPage(
      app_ui("app")
    ),
    server = function(input, output, session) {
      callModule(app_module, id = "app", dir = reactive(dir))
    }
  )
  
  shiny::runApp(app)
}


