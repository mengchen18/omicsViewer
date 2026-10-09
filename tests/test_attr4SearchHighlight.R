## attr4 Search cascade -> searchon select -> gene highlight (open circle)
##
## Regression board for the 2026-10 fix: the searchValue wiring from the
## searchon select had lived in a `debounce(foo, 1000)`-create-and-discard
## observe that was removed as dead code (ddde9d4) - but shiny's debounce()
## itself installs a tracker observer that eagerly runs its argument, so the
## removal silently killed gene highlighting in every attr4 panel (the
## volcano's "circle the gene with an open ring" feature). Re-installed as a
## plain observeEvent; a select needs no debounce.
##
##   1. module: picking a gene sets params$highlight (+ highlightName)
##   2. module: deselecting clears params$highlight
##   3. module: unsetting the Search cascade clears params$highlight
##   4. module: params$status carries searchValue (snapshot restore path)
##   5. meta_scatter: the volcano params carry highlight/highlightName
##   6. plotly_scatter: highlight draws the circle-open trace at the gene

local({
  library(unittest, quietly = TRUE)
  pkgload::load_all(".", quiet = TRUE, attach_testthat = FALSE)
})
dat <- readRDS("inst/extdata/demo.RDS")
fd <- fData(dat)
ts_demo <- do.call(rbind, strsplit(colnames(fd), "|", fixed = TRUE))

## browser model (updateSelectInput ack loop; test_attr4TooltipDefault pattern)
Q <- list()
restore_fns <- local({
  imp <- parent.env(asNamespace("omicsViewer"))
  fns <- lapply(c("updateSelectInput", "updateSelectizeInput"), function(f) local({
    fn <- f; orig <- get(fn, envir = imp); force(orig)
    if (bindingIsLocked(fn, imp)) unlockBinding(fn, imp)
    assign(fn, function(session, inputId, ..., selected = NULL) {
      if (is.character(selected) && length(selected) == 1L && nzchar(selected))
        Q[[length(Q) + 1L]] <<- list(id = session$ns(inputId), value = selected)
      orig(session = session, inputId = inputId, ..., selected = selected)
    }, envir = imp)
    function() {
      if (bindingIsLocked(fn, imp)) unlockBinding(fn, imp)
      assign(fn, orig, envir = imp)
      lockBinding(fn, imp)
    }
  }))
  fns
})
ack <- function(session, rounds = 40) {
  for (r in seq_len(rounds)) {
    q <- Q; Q <<- list()
    final <- list()
    for (m in q) final[[m$id]] <- m$value
    final <- final[!vapply(names(final), function(id)
      identical(session$input[[id]], final[[id]]), logical(1))]
    if (!length(final)) break
    do.call(session$setInputs, final)
    session$flushReact()
  }
  for (i in 1:5) session$flushReact()
}
warm_connect <- function(session) {
  v <- list()
  for (grp in c("selectColorUI", "selectShapeUI", "selectSizeUI",
                "selectTooltipUI", "selectSearchCol"))
    for (w in c("analysis", "subset", "variable"))
      v[[paste0("a4-", grp, "-", w)]] <- ""
  v[["a4-xcut"]] <- "log10(2)"; v[["a4-ycut"]] <- "-log10(0.05)"
  v[["a4-scorner"]] <- "None"
  do.call(session$setInputs, v)
  session$flushReact()
  ack(session)
}
# pick the Search cascade column like a user (level by level)
pick_search_col <- function(session) {
  session$setInputs("a4-selectSearchCol-analysis" = "General")
  session$flushReact(); ack(session)
  session$setInputs("a4-selectSearchCol-subset" = "All")
  session$flushReact(); ack(session)
  session$setInputs("a4-selectSearchCol-variable" = "Gene.name")
  session$flushReact(); ack(session); ack(session)
}

## ---- 1-4. module level ---------------------------------------------------
root <- omicsViewer:::widget_store_new()
st <- omicsViewer:::widget_store_child(root, "test.a4")
res <- new.env(parent = emptyenv())
shiny::testServer(function(input, output, session) {
  res$a4 <<- omicsViewer:::attr4selector_module(
    "a4", reactive_meta = reactive(fd), reactive_expr = reactive(NULL),
    reactive_triset = reactive(ts_demo), store = st, default_tooltip = TRUE)
}, {
  warm_connect(session)
  pick_search_col(session)
  gene <- as.character(fd[3, "General|All|Gene.name"])
  session$setInputs("a4-searchon" = gene)
  session$flushReact(); ack(session)
  res$picked <<- gene
  res$highlight <<- shiny::isolate(res$a4$highlight)
  res$highlightName <<- shiny::isolate(res$a4$highlightName)
  res$expected <<- which(fd[, "General|All|Gene.name"] %in% gene)
  # 2. deselect everything
  session$setInputs("a4-searchon" = character(0))
  session$flushReact(); ack(session)
  res$highlight_deselect <<- shiny::isolate(res$a4$highlight)
  # reselect, then 3. unset the Search cascade
  session$setInputs("a4-searchon" = gene)
  session$flushReact(); ack(session)
  session$setInputs("a4-selectSearchCol-variable" = "--select--")
  session$flushReact(); ack(session); ack(session)
  res$highlight_unset <<- shiny::isolate(res$a4$highlight)
  # 4. status mirror (reselect first)
  session$setInputs("a4-selectSearchCol-variable" = "Gene.name")
  session$flushReact(); ack(session)
  session$setInputs("a4-searchon" = gene)
  session$flushReact(); ack(session)
  res$status <<- shiny::isolate(res$a4$status)
})
ok(ut_cmp_identical(res$highlight, res$expected),
   "search - picking a gene sets params$highlight to the feature index")
ok(ut_cmp_equal(res$highlightName, "Gene.name"),
   "search - highlightName carries the Search column variable")
ok(ut_cmp_equal(is.null(res$highlight_deselect), TRUE),
   "search - deselecting all genes clears params$highlight")
ok(ut_cmp_equal(is.null(res$highlight_unset), TRUE),
   "search - unsetting the Search cascade clears params$highlight")
ok(ut_cmp_equal(res$status$searchValue, res$picked),
   "search - params$status carries searchValue (snapshot restore path)")

## ---- 5. meta_scatter (volcano) integration -------------------------------
CAP <- new.env(); CAP$last <- NULL
orig_psm <- omicsViewer:::plotly_scatter_module
assignInNamespace("plotly_scatter_module", function(id, reactive_param_plotly_scatter, ...) {
  wrapped <- reactive({
    p <- reactive_param_plotly_scatter()
    CAP$last <<- p
    p
  })
  orig_psm(id, wrapped, ...)
}, ns = "omicsViewer")
mres <- new.env(parent = emptyenv())
mstore <- omicsViewer:::widget_store_new()
mst <- omicsViewer:::widget_store_child(mstore, "dataspace.feature_space")
shiny::testServer(function(input, output, session)
  omicsViewer:::meta_scatter_module(
    "fs", reactive_meta = reactive(fd), reactive_expr = reactive(exprs(dat)),
    combine = "feature", source = "fs",
    reactive_x = reactive("ttest|RE_vs_ME|mean.diff"),
    reactive_y = reactive("ttest|RE_vs_ME|log.fdr"), store = mst,
    selection = omicsViewer:::selection_port(
      omicsViewer:::selection_store_new(), "feature")),
  {
    init <- list("fs-axisMode" = "quick",
                 "fs-a4selector-xcut" = "log10(2)",
                 "fs-a4selector-ycut" = "-log10(0.05)",
                 "fs-a4selector-scorner" = "None")
    for (t in c("fs-tris_main_scatter1", "fs-tris_main_scatter2"))
      for (w in c("analysis", "subset", "variable"))
        init[[paste0(t, "-", w)]] <- ""
    for (g in c("selectColorUI", "selectShapeUI", "selectSizeUI",
                "selectTooltipUI", "selectSearchCol"))
      for (w in c("analysis", "subset", "variable"))
        init[[paste0("fs-a4selector-", g, "-", w)]] <- ""
    do.call(session$setInputs, init)
    for (i in 1:5) session$flushReact()
    ack(session)
    mres$highlight_before <<- CAP$last$highlight
    session$setInputs("fs-a4selector-selectSearchCol-analysis" = "General")
    session$flushReact(); ack(session)
    session$setInputs("fs-a4selector-selectSearchCol-subset" = "All")
    session$flushReact(); ack(session)
    session$setInputs("fs-a4selector-selectSearchCol-variable" = "Gene.name")
    session$flushReact(); ack(session)
    gene <- as.character(fd[3, "General|All|Gene.name"])
    session$setInputs("fs-a4selector-searchon" = gene)
    session$flushReact(); ack(session)
    mres$highlight <<- CAP$last$highlight
    mres$highlightName <<- CAP$last$highlightName
    mres$expected <<- which(fd[, "General|All|Gene.name"] %in% gene)
  })
assignInNamespace("plotly_scatter_module", orig_psm, ns = "omicsViewer")
ok(ut_cmp_equal(is.null(mres$highlight_before), TRUE),
   "volcano - no highlight before a gene is searched")
ok(ut_cmp_identical(mres$highlight, mres$expected),
   "volcano - the searched gene reaches the scatter params as highlight")
ok(ut_cmp_equal(mres$highlightName, "Gene.name"),
   "volcano - the scatter params carry the highlightName")

## ---- 6. plotly_scatter open-circle trace ---------------------------------
pp <- omicsViewer:::plotly_scatter(
  x = fd[, "ttest|RE_vs_ME|mean.diff"], y = fd[, "ttest|RE_vs_ME|log.fdr"],
  xlab = "mean.diff", ylab = "log.fdr",
  highlight = which(fd[, "General|All|Gene.name"] %in% "MAP1LC3B"),
  highlightName = "Gene.name")
hl <- Filter(function(tr)
  identical(tr$type, "scatter") && grepl("circle-open", tr$marker$symbol %||% "", fixed = TRUE),
  pp$fig$x$attrs)
ok(ut_cmp_equal(length(hl), 1L),
   "plotly - exactly one circle-open highlight trace is added")
ok(ut_cmp_equal(
  c(hl[[1]]$x, hl[[1]]$y),
  c(fd[3, "ttest|RE_vs_ME|mean.diff"], fd[3, "ttest|RE_vs_ME|log.fdr"])),
  "plotly - the open circle sits on the highlighted gene's coordinates")

invisible(lapply(restore_fns, function(f) f()))
