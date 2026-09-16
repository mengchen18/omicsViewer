library(shiny)
library(shinytest2)
library(Biobase)

test_that("widget-only snapshots preserve active tabs and DataTable state", {
  skip_on_cran()

  snapshot_dir <- file.path(tempdir(), "omicsviewer-widget-state")
  dir.create(snapshot_dir, recursive = TRUE, showWarnings = FALSE)

  set.seed(1234)
  expr <- matrix(
    rnorm(60 * 8),
    nrow = 60,
    dimnames = list(paste0("feature_", sprintf("%03d", 1:60)), paste0("sample_", 1:8))
  )
  pdata <- data.frame(
    `General|All|group` = rep(c("A", "B"), each = 4),
    check.names = FALSE,
    row.names = colnames(expr)
  )
  fdata <- data.frame(
    `General|All|score` = seq_len(60),
    check.names = FALSE,
    row.names = rownames(expr)
  )
  eset <- ExpressionSet(
    assayData = expr,
    phenoData = AnnotatedDataFrame(pdata),
    featureData = AnnotatedDataFrame(fdata)
  )

  widget_app <- shinyApp(
    ui = fluidPage(omicsViewer:::app_ui("app")),
    server = function(input, output, session) {
      omicsViewer:::app_module(
        "app",
        .dir = reactive(snapshot_dir),
        ESVObj = reactive(eset),
        appName = "omicsViewer widget-state test"
      )
    }
  )

  app <- AppDriver$new(
    widget_app,
    name = "widget-state",
    seed = 1234L,
    timeout = 15000L,
    height = 900L,
    width = 1400L
  )
  on.exit(app$stop(), add = TRUE)

  app$wait_for_idle()

  # Move to the feature table and advance to page three. Active-page and table
  # pagination state must round-trip; stable selected-row IDs are covered by
  # tests/test_tableWidgetState.R.
  app$set_inputs(`app-dataspace-eset` = "Feature table")
  app$wait_for_idle()
  app$wait_for_value(input = "app-dataspace-tab_feature-table_state", ignore = list(NULL))

  # Advance with the real pagination control. This browser-owned widget value must
  # survive both snapshot creation and restoration.
  app$click(selector = "#app-dataspace-tab_feature-table .next")
  app$wait_for_idle()
  app$click(selector = "#app-dataspace-tab_feature-table .next")
  app$wait_for_idle()

  saved_state <- app$get_value(input = "app-dataspace-tab_feature-table_state")
  testthat::expect_equal(saved_state$start, 50L)
  testthat::expect_equal(saved_state$length, 25L)

  app$click("app-snapshot")
  app$set_inputs(`app-snapshot_name` = "phase4")
  app$click("app-snapshot_save")
  app$wait_for_idle()

  snapshot_path <- file.path(
    snapshot_dir,
    "ESVSnapshot_ESVObj.RDS_phase4.ESS"
  )
  stopifnot(file.exists(snapshot_path))
  snapshot <- readRDS(snapshot_path)

  stopifnot(
    identical(snapshot$policy, "widget-only"),
    identical(snapshot$app$data_active_tab, "Feature table"),
    identical(snapshot$panels$data_space$eset_active_tab, "Feature table"),
    identical(snapshot$panels$data_space$eset_fdata_tab$start, 50L),
    identical(snapshot$panels$data_space$eset_fdata_tab$length, 25L),
    is.null(snapshot$panels$result_space$analyst_feature_general$htestV1),
    is.null(snapshot$panels$result_space$analyst_feature_general$htestV2),
    is.null(snapshot$panels$result_space$analyst_gene_shot$rif),
    identical(snapshot$gaps$stringdb, omicsViewer:::APP_STATE_GAPS$stringdb)
  )

  # Restore through the same snapshot table and verify the active tab and saved
  # DataTable page return without storing any result payload.
  app$set_inputs(`app-dataspace-eset` = "Feature")
  app$wait_for_idle()
  app$click("app-snapshot")
  app$wait_for_idle()
  app$click(selector = "#app-tab_saveSS table tbody tr:first-child td:first-child")
  app$wait_for_idle(timeout = 15000L)
  stopifnot(
    identical(app$get_value(input = "app-dataspace-eset"), "Feature table"),
    identical(app$get_value(input = "app-dataspace-tab_feature-table_state")$start, 50L)
  )
})
