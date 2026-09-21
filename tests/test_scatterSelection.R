library(shiny)
library(omicsViewer)
library(unittest, quietly = TRUE)

# Plotly owns the visible box/lasso. If a selection event invalidates the plot
# parameters, renderPlotly recreates the figure and erases that selection box.
redraw_count <- 0L
last_in_selection <- NULL
last_y_axis <- NULL
scatter_result <- NULL
trace(
  omicsViewer:::plotly_scatter,
  tracer = quote({
    redraw_count <<- redraw_count + 1L
    last_in_selection <<- inSelection
    last_y_axis <<- ylab
  }),
  print = FALSE,
  where = asNamespace("omicsViewer")
)

expr <- matrix(rnorm(40), nrow = 4, dimnames = list(paste0("f", 1:4), paste0("s", 1:10)))
pdata <- data.frame(
  `General|All|x` = 1:10,
  `General|All|y` = (1:10) * 2,
  `General|All|z` = (1:10) * 4,
  check.names = FALSE,
  row.names = colnames(expr)
)

scatterApp <- function(input, output, session) {
  scatter_result <<- omicsViewer:::meta_scatter_module(
    "scatter",
    reactive_meta = reactive(pdata),
    reactive_expr = reactive(expr),
    combine = "pheno",
    source = "selectiontest",
    store = omicsViewer:::widget_store_child(
      omicsViewer:::widget_store_new(), "test.scatter")
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

  after_axis_change <- redraw_count
  session$setInputs(
    `scatter-tris_main_scatter2-variable` = "z"
  )
  session$flushReact()
  session$flushOutput()
  ok(
    ut_cmp_identical(redraw_count > after_axis_change, TRUE),
    "changing a scatter axis redraws the new figure"
  )
  ok(
    ut_cmp_identical(last_in_selection, NA),
    "selection opacity is not inherited by a different scatter figure"
  )
})

# Snapshot restoration may emphasize a semantic selection, but only after the
# restored triselectors reach the axes saved with that selection.
restoration_status <- NULL
restorationApp <- function(input, output, session) {
  restoration_status <<- reactiveVal(NULL)
  omicsViewer:::meta_scatter_module(
    "scatter",
    reactive_meta = reactive(pdata),
    reactive_expr = reactive(expr),
    combine = "pheno",
    source = "restoretest",
    reactive_status = restoration_status,
    store = omicsViewer:::widget_store_child(
      omicsViewer:::widget_store_new(), "test.restore")
  )
}

testServer(restorationApp, {
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

  restoration_status(list(
    axisMode = "custom",
    xax = list(analysis = "General", subset = "All", variable = "x"),
    yax = list(analysis = "General", subset = "All", variable = "z"),
    showRegLine = FALSE,
    attr4 = NULL,
    selection_clicked = NULL,
    selection_selected = c("s2", "s3"),
    selectByCorner = FALSE
  ))
  # testServer does not always flush programmatic selectize updates made by
  # module restoration; drive the restored value explicitly to exercise the
  # final axis/selection state in the same reactive graph.
  session$setInputs(
    `scatter-tris_main_scatter2-variable` = "z"
  )
  session$flushReact()
  session$flushOutput()
  session$flushReact()
  session$flushOutput()

  ok(
    ut_cmp_identical(last_in_selection, c(2L, 3L)),
    "a restored selection is emphasized on its restored axes"
  )
})

untrace(omicsViewer:::plotly_scatter, where = asNamespace("omicsViewer"))
