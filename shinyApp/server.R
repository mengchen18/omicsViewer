library(ExpressionSetViewer)
v <- "/media/ExpressionSetViewerData"
server <- function(input, output, session) {
	observe({
		print(
			list.files(v)
			)
		})
  callModule(ExpressionSetViewer:::app_module, id = "app", dir = reactive(v))
}


