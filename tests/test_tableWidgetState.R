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

# ---- WP8: dataTable page + column_filters store bindings ----------------
widget_store_new <- omicsViewer:::widget_store_new
widget_store_child <- omicsViewer:::widget_store_child
store_read <- omicsViewer:::store_read
store_apply <- omicsViewer:::store_apply

store_saved <- NULL
tableStoreApp <- function(input, output, session) {
  store <- widget_store_new()
  .dt_store <- widget_store_child(store, "dataspace.tab_pheno")
  table_result <<- omicsViewer:::dataTable_module(
    "dt",
    reactive_data = reactive(pd),
    tab_status = reactive(NULL),
    tab_rows = reactive(TRUE),
    store = .dt_store
  )
  store_saved <<- store
}

testServer(tableStoreApp, {
  session$setInputs(
    `dt-multisel` = TRUE,
    `dt-table_rows_selected` = 25L,
    `dt-table_state` = table_state
  )
  session$flushReact()
  session$flushOutput()
})

reg_keys <- names(store_saved$bindings)
ok(
  setequal(intersect(reg_keys, c("dataspace.tab_pheno.page",
                                 "dataspace.tab_pheno.column_filters")),
           c("dataspace.tab_pheno.page", "dataspace.tab_pheno.column_filters")),
  "dataTable module registers page and column_filters bindings"
)
vals <- store_read(store_saved)
ok(
  identical(vals$dataspace.tab_pheno.page, 2L),
  "browser-reported table state syncs the page into the store"
)
ok(
  identical(vals$dataspace.tab_pheno.column_filters,
            setNames(character(0), character(0))),
  "empty column filters sync as an empty mapping"
)

# agent apply lands in the store (browser push needs a real session)
applied_keys <- character()
testServer(tableStoreApp, {
  session$setInputs(`dt-multisel` = FALSE)
  session$flushReact()
  session$flushOutput()
  receipt <- store_apply(
    store_saved,
    list(`dataspace.tab_pheno.column_filters` = list(`General|All|group` = "A"),
         `dataspace.tab_pheno.page` = 3L),
    origin = "agent"
  )
  applied_keys <<- unlist(receipt$applied)
  session$flushReact()
  session$flushOutput()
})
vals2 <- store_read(store_saved)
ok(
  setequal(applied_keys, c("dataspace.tab_pheno.column_filters",
                           "dataspace.tab_pheno.page")),
  "agent apply receipt reports the pushed keys"
)
ok(
  identical(vals2$dataspace.tab_pheno.column_filters, c(`General|All|group` = "A")) &&
    identical(vals2$dataspace.tab_pheno.page, 3L),
  "agent applies of filters/page land in the store"
)
