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
      omicsViewer:::widget_store_new(), "test.scatter"),
    selection = omicsViewer:::selection_port(
      omicsViewer:::selection_store_new(), "sample")
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
      omicsViewer:::widget_store_new(), "test.restore"),
    selection = omicsViewer:::selection_port(
      omicsViewer:::selection_store_new(), "sample")
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

# ---------------------------------------------------------------------------
# Volcano corner selections through the unified selection bus. The corner
# auto-selection is server-owned: its emphasis is resolved IN-RENDER from
# the live rectangles (same paint as an axis or cutoff change), and its
# propagation rides the selection bus like every other source. These rows
# regress the reported defect: switching the y axis log.fdr -> log.pvalue
# grew the corner selection in the tables while the figure kept NA
# emphasis, and the clear button could not clear a corner selection (the
# corner observer re-fired on the clear_counter invalidation and
# resurrected it).
fd2 <- data.frame(
  row.names = paste0("f", 1:4),
  check.names = FALSE,
  `ttest|KO_vs_WT|mean.diff` = c(2, -2, 0.5, -0.4),
  `ttest|KO_vs_WT|log.fdr` = c(3, 3, 0.5, 0.4),
  `ttest|KO_vs_WT|log.pvalue` = c(6, 6, 1.5, 1.4)
)
expr2 <- matrix(rnorm(4 * 10), nrow = 4,
                dimnames = list(paste0("f", 1:4), paste0("s", 1:10)))

corner_bus <- NULL
cornerApp <- function(input, output, session) {
  corner_bus <<- omicsViewer:::selection_port(
    omicsViewer:::selection_store_new(), "feature")
  omicsViewer:::meta_scatter_module(
    "scatter",
    reactive_meta = reactive(fd2),
    reactive_expr = reactive(expr2),
    combine = "feature",
    source = "cornerbus",
    store = omicsViewer:::widget_store_child(
      omicsViewer:::widget_store_new(), "test.corner"),
    selection = corner_bus
  )
}

testServer(cornerApp, {
  session$setInputs(
    `scatter-tris_main_scatter1-analysis` = "ttest",
    `scatter-tris_main_scatter1-subset` = "KO_vs_WT",
    `scatter-tris_main_scatter1-variable` = "mean.diff",
    `scatter-tris_main_scatter2-analysis` = "ttest",
    `scatter-tris_main_scatter2-subset` = "KO_vs_WT",
    `scatter-tris_main_scatter2-variable` = "log.fdr",
    `scatter-a4selector-xcut` = "log10(2)",
    `scatter-a4selector-ycut` = "-log10(0.05)",
    `scatter-a4selector-scorner` = "volcano",
    `scatter-main_scatterOutput-showRegLine` = FALSE
  )
  session$flushReact(); session$flushOutput()
  session$flushReact(); session$flushOutput()

  ok(
    ut_cmp_identical(corner_bus$read()$ids, c("f1", "f2")),
    "the volcano corner selection is reported to the selection bus"
  )
  ok(
    identical(corner_bus$read()$origin, "corner"),
    "the corner is the reporting origin"
  )
  ok(
    ut_cmp_identical(last_in_selection, c(1L, 2L)),
    "the load-time corner selection is emphasized in the figure"
  )
  ok(
    ut_cmp_identical(corner_bus$read()$mirror, c("f1", "f2")),
    "the corner selection mirrors into the table row highlight"
  )

  paints_before_switch <- redraw_count
  session$setInputs(`scatter-tris_main_scatter2-variable` = "log.pvalue")
  session$flushReact(); session$flushOutput()
  session$flushReact(); session$flushOutput()

  ok(
    ut_cmp_identical(corner_bus$read()$ids, c("f1", "f2", "f3", "f4")),
    "switching log.fdr -> log.pvalue re-derives the (larger) corner selection"
  )
  ok(
    ut_cmp_identical(last_in_selection, c(1L, 2L, 3L, 4L)),
    "the new figure paints the new corner emphasis (the reported bug)"
  )
  ok(
    ut_cmp_identical(redraw_count - paints_before_switch, 1L),
    "the axis switch plus corner reselection is a single paint"
  )

  # a browser lasso on the switched figure reports through the bus and
  # does not repaint (Plotly owns the visible selection)
  paints_before_lasso <- redraw_count
  session$setInputs(`plotly_selected-cornerbus` = '[{"x":2,"y":6}]')
  session$flushReact(); session$flushOutput()
  ok(
    ut_cmp_identical(corner_bus$read()$ids, "f1"),
    "a browser lasso reports through the bus and overrides the corner"
  )
  ok(
    ut_cmp_identical(redraw_count, paints_before_lasso),
    "a browser lasso does not repaint the figure"
  )

  # clear empties the bus record (tables and result space follow) and
  # disengages the corner; an axis switch must not resurrect it
  session$setInputs(`scatter-clear` = 1L)
  session$flushReact(); session$flushOutput()
  session$flushReact(); session$flushOutput()
  ok(
    length(corner_bus$read()$ids) == 0L && identical(corner_bus$read()$origin, "clear"),
    "clear empties the selection bus record"
  )
  ok(
    ut_cmp_identical(corner_bus$read()$mirror, TRUE),
    "clear un-filters the table mirror"
  )
  ok(
    ut_cmp_identical(last_in_selection, NA),
    "clear drops the figure emphasis"
  )

  session$setInputs(`scatter-tris_main_scatter2-variable` = "log.fdr")
  session$flushReact(); session$flushOutput()
  session$flushReact(); session$flushOutput()
  ok(
    length(corner_bus$read()$ids) == 0L,
    "an axis switch after clear does not resurrect the corner selection"
  )

  # a genuine corner edit re-engages the corner
  session$setInputs(`scatter-a4selector-scorner` = "None")
  session$flushReact(); session$flushOutput()
  session$setInputs(`scatter-a4selector-scorner` = "volcano")
  session$flushReact(); session$flushOutput()
  session$flushReact(); session$flushOutput()
  ok(
    ut_cmp_identical(corner_bus$read()$ids, c("f1", "f2")),
    "a corner edit re-engages the corner auto-selection"
  )
})

untrace(omicsViewer:::plotly_scatter, where = asNamespace("omicsViewer"))
