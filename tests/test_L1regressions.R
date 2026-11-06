# todo 4.6 module-level regression probes for the L1 data/result-space
# shells. These run the real modules headless through testServer.
#
# WP-11 (panel debt): extended with the registry-composition probes -
# the L1 server body composes from panel_records (one loop), so this
# suite pins that the composed L1 boots headless WITH a store and that
# the status builder emits the frozen analyst_* key set in order. The
# suite loads the working tree (pkgload; the test_panelRegistry
# convention) because the WP-11 code under test is in R/ - fall back to
# the installed package when pkgload is unavailable.

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkg_root <- if (file.exists("DESCRIPTION")) "."
  else if (file.exists(file.path("..", "DESCRIPTION"))) ".." else "."
  pkgload::load_all(pkg_root, quiet = TRUE, attach_testthat = FALSE)
} else {
  library(omicsViewer)
}
library(unittest, quietly = TRUE)
library(shiny)

# ------------------------------------------------------------------
# L1 (result space): store_rs was referenced unconditionally by the
# status-restore observer although it is only assigned when a store is
# given - exported use with store = NULL plus a status carrying
# analyst_active_tab errored the observer ("object 'store_rs' not
# found"), closing the session.
# ------------------------------------------------------------------
dat <- readRDS(system.file("extdata", "demo.RDS", package = "omicsViewer"))
l1_err <- NULL
l1_alive <- NULL
shiny::testServer(function(input, output, session) {
  omicsViewer:::L1_result_space_module(
    "rs",
    reactive_expr = reactive(Biobase::exprs(dat)),
    reactive_phenoData = reactive(Biobase::pData(dat)),
    reactive_featureData = reactive(Biobase::fData(dat)),
    status = reactive(list(analyst_active_tab = "Sample")),
    store = NULL)
}, {
  l1_err <<- tryCatch({
    for (i in 1:10) session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  l1_alive <<- !session$isClosed()
})
ok(is.null(l1_err),
   "L1 result space: store=NULL + status restore raises no error (4.6/L1)")
ok(isTRUE(l1_alive),
   "L1 result space: session stays alive with store=NULL (4.6/L1)")

# ------------------------------------------------------------------
# WP-11: the registry-composed L1 boots headless WITH a store and the
# status builder emits the frozen key set in exact order (the snapshot
# panels$result_space schema). The full byte-level goldens (store
# snapshot, registration order, payload content) live in
# tests/test_panelRegistry.R's WP-11 section; this is the L1 shell's own
# quick probe.
# ------------------------------------------------------------------
L1_W11_STATUS_NAMES <- c(
  "analyst_active_tab", "analyst_feature_general",
  "analyst_sample_general", "analyst_gene_shot", "analyst_fgsea",
  "analyst_ora", "analyst_ptm", "analyst_dose_response",
  "analyst_additional")
l1_store <- omicsViewer:::widget_store_new()
l1_ret <- NULL
l1_w11_err <- NULL
l1_w11_names <- NULL
shiny::testServer(function(input, output, session) {
  l1_ret <<- omicsViewer:::L1_result_space_module(
    "rs",
    reactive_expr = reactive(Biobase::exprs(dat)),
    reactive_phenoData = reactive(Biobase::pData(dat)),
    reactive_featureData = reactive(Biobase::fData(dat)),
    status = reactive(NULL),
    store = l1_store)
}, {
  l1_w11_err <<- tryCatch({
    for (i in 1:3) session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  l1_w11_names <<- tryCatch(
    names(shiny::isolate(l1_ret())), error = function(e) NULL)
})
ok(is.null(l1_w11_err),
   "L1 result space: the WP-11 registry-composed L1 boots headless with a store")
ok(identical(l1_w11_names, L1_W11_STATUS_NAMES),
   "L1 result space: the status builder emits the frozen analyst_* key set in order (WP-11)")
ok(identical(length(omicsViewer:::store_registry_view(l1_store)), 72L),
   "L1 result space: the composed store carries the frozen 72-record enumeration (WP-11/D12; the full frozen vector lives in test_panelRegistry.R)")

# ------------------------------------------------------------------
# WP-11: with additional tabs the legacy 5-field loop stays its own
# loop (Q1/D15 - never absorbed); the composed status payload keeps
# appending the plugin state after the built-in keys.
# ------------------------------------------------------------------
.at_rec <- list(
  tabName = "Custom", moduleName = "customMod",
  moduleUi = function(id) shiny::div(paste0("custom-", id)),
  moduleServer = function(id, pdata, fdata, expr, feature_selected,
                          sample_selected, object, ...)
    shiny::reactive(list(state = list(k = 1L))))
l1_at_names <- NULL
l1_at_err <- NULL
shiny::testServer(function(input, output, session) {
  v <- omicsViewer:::L1_result_space_module(
    "rs",
    reactive_expr = reactive(Biobase::exprs(dat)),
    reactive_phenoData = reactive(Biobase::pData(dat)),
    reactive_featureData = reactive(Biobase::fData(dat)),
    additionalTabs = list(.at_rec),
    status = reactive(NULL),
    store = omicsViewer:::widget_store_new())
  l1_at_err <<- tryCatch({
    session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  l1_at_names <<- tryCatch(names(shiny::isolate(v())),
                           error = function(e) NULL)
}, { })
ok(is.null(l1_at_err) && identical(l1_at_names, L1_W11_STATUS_NAMES),
   "L1 result space: additionalTabs compose unchanged through the WP-11 loop (Q1: the legacy loop stays)")
