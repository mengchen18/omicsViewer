library(omicsViewer)
library(unittest, quietly = TRUE)
library(shiny)

# R-H5: after the table DATA changes, the sync observer used to re-derive
# the selected id from the OLD positional input$table_rows_selected
# against the NEW table and write it to the store as a user edit -
# silently switching the linked view (e.g. the selected ORA pathway).
# The epoch guard must hold the report until the browser re-reports under
# the new table (or DT clears it).

dtd_host <- function() {
  tab <- reactiveVal(data.frame(pathway = paste0("gs", 1:5), p = 1:5))
  ids <- reactiveVal(paste0("gs", 1:5))
  function(input, output, session) {
    root <- omicsViewer:::widget_store_new()
    omicsViewer:::dataTableDownload_module(
      "dtd", reactive_table = tab, reactive_row_ids = ids,
      prefix = "t_", store = omicsViewer:::widget_store_child(root, "tab"),
      store_key = "sel")
    session$userData$root <- root
    session$userData$tab <- tab
    session$userData$ids <- ids
  }
}

val <- NULL
shiny::testServer(dtd_host(), {
  for (i in 1:6) session$flushReact()
  # user picks row 3 -> gs3
  session$setInputs(`dtd-table_rows_selected` = 3)
  for (i in 1:4) session$flushReact()
  v1 <- omicsViewer:::store_read(session$userData$root, "tab.sel")[["tab.sel"]]
  # table data changes underneath (new ranking, ids reshuffled: gs3 now row 5)
  session$userData$tab(data.frame(pathway = paste0("gs", c(1, 2, 4, 5, 3)),
                                  p = c(5, 4, 3, 2, 1)))
  session$userData$ids(paste0("gs", c(1, 2, 4, 5, 3)))
  for (i in 1:4) session$flushReact()
  # the OLD positional input (3) still reports; it must NOT be re-mapped
  # against the new table (which would select gs4 and switch the linked
  # view). The store keeps gs3 until the browser re-reports.
  v2 <- omicsViewer:::store_read(session$userData$root, "tab.sel")[["tab.sel"]]
  val <<- c(v1 = v1, v2 = v2)
})
ok(identical(val[["v1"]], "gs3"),
   "dataTableDownload: user click reports the row id (R-H5 setup)")
ok(identical(val[["v2"]], "gs3"),
   "dataTableDownload: stale positional report after a table data change does not switch the selected id (R-H5)")
