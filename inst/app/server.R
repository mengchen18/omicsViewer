library(omicsViewer)
v <- "/media/ExpressionSetViewerData"
server <- function(input, output, session) {
	observe({
		print(
			list.files(v)
			)
		})
  callModule(omicsViewer::app_module, id = "app", .dir = reactive(v), esetLoader = readESVObj)
}


