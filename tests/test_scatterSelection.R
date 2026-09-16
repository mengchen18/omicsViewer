library(shiny)
library(omicsViewer)
library(unittest, quietly = TRUE)

# Plotly owns the visible box/lasso. If a selection event invalidates the plot
# parameters, renderPlotly recreates the figure and erases that selection box.
redraw_count <- 0L
scatter_result <- NULL
trace(
  omicsViewer:::plotly_scatter,
  tracer = quote(redraw_count <<- redraw_count + 1L),
  print = FALSE,
  where = asNamespace("omicsViewer")
)

expr <- matrix(rnorm(40), nrow = 4, dimnames = list(paste0("f", 1:4), paste0("s", 1:10)))
pdata <- data.frame(
  `General|All|x` = 1:10,
  `General|All|y` = (1:10) * 2,
  check.names = FALSE,
  row.names = colnames(expr)
)

scatterApp <- function(input, output, session) {
  scatter_result <<- omicsViewer:::meta_scatter_module(
    "scatter",
    reactive_meta = reactive(pdata),
    reactive_expr = reactive(expr),
    combine = "pheno",
    source = "selectiontest"
  )
}

testServer(scatterApp, {
  session$setInputs(
    `scatter-tris_main_scatter1-analysis` = "General",
    `scatter-tris_main_scatter1-subset` = "All",
    `scatter-tris_main_scatter1-variable` = "x",
    `scatter-tris_main_scatter2-analysis` = "General",
    `scatter-tris_main_scatter2-subset` = "All",
    `scatter-tris_main_scatter2-variable` = "y",
    `scatter-main_scatterOutput-showRegLine` = FALSE
  )
  session$flushReact()
  session$flushOutput()

  before <- redraw_count
  ok(before > 0L, "scatter test draws an initial Plotly figure")

  session$setInputs(
    `plotly_selected-selectiontest` = '[{"x":2,"y":4},{"x":3,"y":6}]'
  )
  session$flushReact()
  session$flushOutput()

  ok(
    ut_cmp_identical(redraw_count, before),
    "a Plotly selection does not redraw and erase the selection box"
  )
  ok(
    ut_cmp_identical(scatter_result()$selected, c("s2", "s3")),
    "the semantic Plotly selection is still forwarded to downstream modules"
  )
})

untrace(omicsViewer:::plotly_scatter, where = asNamespace("omicsViewer"))
