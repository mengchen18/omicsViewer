library(omicsViewer)
library(unittest, quietly = TRUE)
library(shiny)

# todo 4.6 module-level regression probes for the L1 data/result-space
# shells. These run the real modules headless through testServer.

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
