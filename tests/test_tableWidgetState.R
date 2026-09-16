library(shiny)
library(omicsViewer)
library(unittest, quietly = TRUE)

pd <- data.frame(
  `General|All|group` = rep(c("A", "B"), each = 15),
  check.names = FALSE,
  row.names = paste0("feature_", sprintf("%03d", 1:30))
)

table_state <- list(
  start = 25L,
  length = 25L,
  order = list(c(1L, "asc")),
  columns = list(list(search = list(search = ""))),
  time = 12345
)

saved_status <- NULL
tableApp <- function(input, output, session) {
  table_result <<- omicsViewer:::dataTable_module(
    "dt",
    reactive_data = reactive(pd),
    tab_status = reactive(NULL),
    tab_rows = reactive(TRUE)
  )
}

testServer(tableApp, {
  session$setInputs(
    `dt-multisel` = TRUE,
    `dt-table_rows_selected` = 25L,
    `dt-table_state` = table_state
  )
  session$flushReact()
  session$flushOutput()
  saved_status <<- attr(table_result(), "status")
})

ok(
  ut_cmp_identical(saved_status$start, 25L),
  "data table captures the selected page"
)
ok(
  ut_cmp_identical(saved_status$length, 25L),
  "data table captures the selected page length"
)
ok(
  ut_cmp_identical(saved_status$selected_rows, "feature_025"),
  "data table captures a stable selected row identifier"
)
ok(
  ut_cmp_identical(saved_status$order, list(c(1L, "asc"))),
  "data table captures column ordering"
)
ok(
  ut_cmp_identical(saved_status$time, NULL),
  "data table excludes non-portable DataTable state"
)
