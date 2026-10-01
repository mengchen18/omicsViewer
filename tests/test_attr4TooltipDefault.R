## Tooltips attribute cascade default (attr4selector_module default_tooltip)
##
## The data-space scatters pass default_tooltip = TRUE: when the annotation
## first arrives and the Tooltips cascade has never been set (no store value,
## no user selection, no restore), the FIRST annotation column whose variable
## matches ATTR4_TOOLTIP_DEFAULT_PATTERN ("symbol|name", case-insensitive) is
## preselected via a one-shot store_seed. Everything that lands first - a
## restore, an agent apply, a user pick - wins, and an explicit unset is
## never re-seeded. Dataset switches re-arm the seed.
##
## Regression board:
##   1. pattern matching + column-order priority ("whichever comes first")
##   2. no-match annotations keep the cascade unset
##   3. default_tooltip = FALSE (right-panel instances) change nothing
##   4. other attr4 groups (colour/shape/size/search) are never preselected
##   5. explicit user pick and explicit unset both beat the default;
##      the unset is not re-seeded
##   6. store writes (agent / snapshot restores) that land first beat the
##      default (restore-first-wins)
##   7. dataset switch re-seeds the new dataset's default
##   8. demo integration: meta_scatter boots with General|All|Gene.name

local({
  library(unittest, quietly = TRUE)
  pkgload::load_all(".", quiet = TRUE, attach_testthat = FALSE)
})
dat <- readRDS("inst/extdata/demo.RDS")
fd <- fData(dat)
pd <- pData(dat)
ex <- exprs(dat)

## triset fixtures --------------------------------------------------------
ts_demo <- do.call(rbind, strsplit(colnames(fd), "|", fixed = TRUE))
# "name" column BEFORE the "symbol" column: whichever comes first wins
ts_name_first <- rbind(
  c("General", "All", "Protein.ID"),
  c("General", "All", "Gene.name"),
  c("General", "Annotation", "SYMBOL"))
# symbol first
ts_symbol_first <- rbind(
  c("General", "All", "Protein.ID"),
  c("General", "Annotation", "SYMBOL"),
  c("General", "All", "Gene.name"))
# no match at all
ts_none <- rbind(
  c("General", "All", "Protein.ID"),
  c("ttest", "A_vs_B", "log.fdr"),
  c("PCA", "All", "PC1(10%)"))
# second dataset for the switch probe
ts_switch <- rbind(
  c("General", "All", "Protein.ID"),
  c("General", "Annotation", "Symbol"))

## browser model ---------------------------------------------------------
## updateSelectInput / updateSelectizeInput pushes are queued with their
## FULL id; the ack loop applies each batch like a browser round trip
## (last write per widget wins). Adapted from test_renderStability.R.
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
# a real browser reports every input once on connect
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

## module-level boot ------------------------------------------------------
# Boots attr4selector_module on a fresh child store and returns an env with
# $store (root), $vals (leaf attr4 values), $status (params$status probe)
# and $a4 (params). `script` runs INSIDE the testServer body after the
# connect/ack warm-up and can drive further interactions.
boot_a4 <- function(triset, default_tooltip = TRUE, script = NULL,
                    triset_rv = NULL) {
  root <- omicsViewer:::widget_store_new()
  st <- omicsViewer:::widget_store_child(root, "test.a4")
  res <- new.env(parent = emptyenv())
  tsr <- if (is.null(triset_rv)) shiny::reactive(triset) else triset_rv
  shiny::testServer(function(input, output, session) {
    res$a4 <<- omicsViewer:::attr4selector_module(
      "a4", reactive_meta = reactive(fd), reactive_expr = reactive(NULL),
      reactive_triset = tsr, store = st, default_tooltip = default_tooltip)
  }, {
    warm_connect(session)
    if (!is.null(script)) script(res, session, root, st)
    leaf <- omicsViewer:::widget_store_child(root, "test.a4.attr4")
    keys <- c("tooltip_analysis", "tooltip_subset", "tooltip_variable",
              "color_variable", "shape_variable", "size_variable",
              "search_variable")
    v <- omicsViewer:::store_read(leaf, keys)
    # store_read returns FULL canonical keys - normalise to leaf ids
    names(v) <- keys
    res$vals <<- v
    res$inputs <<- list(
      analysis = input[["a4-selectTooltipUI-analysis"]],
      subset = input[["a4-selectTooltipUI-subset"]],
      variable = input[["a4-selectTooltipUI-variable"]])
    res$status <<- res$a4$status
    res$tooltips <<- res$a4$tooltips
    res$root <<- root
  })
  res
}
tip <- function(res)
  paste(unlist(res$vals[
    c("tooltip_analysis", "tooltip_subset", "tooltip_variable")]), collapse = "|")

## ---- 1. pattern matching + column order ---------------------------------
r <- boot_a4(ts_demo)
ok(ut_cmp_equal(tip(r), "General|All|Gene.name"),
   "default - demo annotation preselects the first symbol/name column")
ok(ut_cmp_equal(r$inputs$variable, "Gene.name"),
   "default - the tooltip variable widget shows the default")

r <- boot_a4(ts_name_first)
ok(ut_cmp_equal(tip(r), "General|All|Gene.name"),
   "default - column order wins (Gene.name before SYMBOL)")

r <- boot_a4(ts_symbol_first)
ok(ut_cmp_equal(tip(r), "General|Annotation|SYMBOL"),
   "default - symbol column wins when it comes first")

r <- boot_a4(ts_none)
ok(ut_cmp_equal(is.null(r$vals$tooltip_variable), TRUE),
   "default - no matching column leaves the cascade unset")
ok(ut_cmp_equal(is.null(r$status$selectTooltip), TRUE),
   "default - no matching column yields no tooltips status")

## ---- 2. params follow the default ----------------------------------------
r <- boot_a4(ts_demo)
ok(ut_cmp_identical(c(r$tooltips),
                    as.character(fd[, "General|All|Gene.name"])),
   "default - params$tooltips carries the Gene.name column")
ok(ut_cmp_equal(identical(attr(r$tooltips, "label"), "General|All|Gene.name"),
                TRUE),
   "default - params$tooltips keeps its plot label")

## ---- 3. other groups are never preselected ------------------------------
r <- boot_a4(ts_demo)
ok(ut_cmp_equal(is.null(r$vals$color_variable), TRUE),
   "scope - colour cascade stays unset")
ok(ut_cmp_equal(is.null(r$vals$shape_variable), TRUE),
   "scope - shape cascade stays unset")
ok(ut_cmp_equal(is.null(r$vals$size_variable), TRUE),
   "scope - size cascade stays unset")
ok(ut_cmp_equal(is.null(r$vals$search_variable), TRUE),
   "scope - search cascade stays unset")

## ---- 4. default_tooltip = FALSE keeps the legacy behaviour ---------------
r <- boot_a4(ts_demo, default_tooltip = FALSE)
ok(ut_cmp_equal(is.null(r$vals$tooltip_variable), TRUE),
   "legacy - default_tooltip=FALSE leaves tooltips unset")

## ---- 5. user pick / unset beat the default -------------------------------
r <- boot_a4(ts_demo, script = function(res, session, root, st) {
  session$setInputs("a4-selectTooltipUI-variable" = "Protein.ID")
  ack(session)
})
ok(ut_cmp_equal(r$vals$tooltip_variable, "Protein.ID"),
   "user pick - a manual variable pick replaces the default")
ok(ut_cmp_equal(length(r$tooltips), nrow(fd)),
   "user pick - params$tooltips follows the picked column")

r <- boot_a4(ts_demo, script = function(res, session, root, st) {
  session$setInputs("a4-selectTooltipUI-variable" = "--select--")
  session$flushReact()
  session$setInputs("a4-scorner" = "None")  # unrelated observer re-run
  ack(session); ack(session)
})
ok(ut_cmp_equal(is.null(r$vals$tooltip_variable), TRUE) ||
      identical(r$vals$tooltip_variable, "--select--"),
   "unset - --select-- clears the store key (no re-seed)")
ok(ut_cmp_equal(r$status$selectTooltip$variable, "--select--"),
   "unset - the committed status carries the explicit unset marker")

## ---- 6. store writes (agent / snapshot restores) win over the default ----
# a restore after the default settled replaces it (store_apply semantics;
# the load-time seed can never clobber a restore because seeds skip held
# keys - store_seed line 707)
r <- boot_a4(ts_demo, script = function(res, session, root, st) {
  leaf <- omicsViewer:::widget_store_child(root, "test.a4.attr4")
  omicsViewer:::store_apply(leaf, list(
    tooltip_analysis = "General", tooltip_subset = "All",
    tooltip_variable = "Protein.ID"), origin = "restore", strict = FALSE)
  session$flushReact()
  ack(session)
})
ok(ut_cmp_equal(r$vals$tooltip_variable, "Protein.ID"),
   "restore - a restore apply replaces the seeded default")
ok(ut_cmp_equal(r$inputs$variable, "Protein.ID"),
   "restore - the widget follows the restored tooltip")

## ---- 7. dataset switch re-seeds the new dataset's default ----------------
r <- local({
  root <- omicsViewer:::widget_store_new()
  st <- omicsViewer:::widget_store_child(root, "test.a4")
  tsr <- shiny::reactiveVal(ts_demo)
  res <- new.env(parent = emptyenv())
  shiny::testServer(function(input, output, session) {
    res$a4 <<- omicsViewer:::attr4selector_module(
      "a4", reactive_meta = reactive(NULL), reactive_expr = reactive(NULL),
      reactive_triset = tsr, store = st, default_tooltip = TRUE)
  }, {
    warm_connect(session)
    leaf <- omicsViewer:::widget_store_child(root, "test.a4.attr4")
    res$before <<- omicsViewer:::store_read(leaf, "tooltip_variable")[[1]]
    # dataset switch: reset + new triset (the app's L0 order)
    omicsViewer:::store_reset(root)
    shiny::isolate(tsr(ts_switch))
    session$flushReact()
    warm_connect(session)
    res$after <<- omicsViewer:::store_read(leaf, "tooltip_variable")[[1]]
  })
  res
})
ok(ut_cmp_equal(r$before, "Gene.name"),
   "dataset switch - dataset A settled on its default")
ok(ut_cmp_equal(r$after, "Symbol"),
   "dataset switch - dataset B re-seeds its own default")

## ---- 8. meta_scatter integration (demo wiring) --------------------------
sc <- omicsViewer:::widget_store_new()
scs <- omicsViewer:::widget_store_child(sc, "dataspace.feature_space")
sres <- new.env(parent = emptyenv())
shiny::testServer(function(input, output, session)
  omicsViewer:::meta_scatter_module(
    "fs", reactive_meta = reactive(fd), reactive_expr = reactive(ex),
    combine = "feature", source = "fs",
    reactive_x = reactive("ttest|RE_vs_ME|mean.diff"),
    reactive_y = reactive("ttest|RE_vs_ME|log.fdr"),
    store = scs,
    selection = omicsViewer:::selection_port(
      omicsViewer:::selection_store_new(), "feature")),
  {
    session$setInputs("fs-axisMode" = "quick")
    session$setInputs("fs-a4selector-xcut" = "log10(2)")
    session$setInputs("fs-a4selector-ycut" = "-log10(0.05)")
    session$setInputs("fs-a4selector-scorner" = "None")
    session$flushReact()
    ack(session); ack(session)
    leaf <- omicsViewer:::widget_store_child(sc, "dataspace.feature_space.attr4")
    sres$vals <<- omicsViewer:::store_read(leaf, c(
      "tooltip_analysis", "tooltip_subset", "tooltip_variable",
      "color_variable"))
  })
ok(ut_cmp_equal(
  paste(unlist(sres$vals[1:3]), collapse = "|"), "General|All|Gene.name"),
  "meta_scatter - demo boots with the Gene.name tooltip default")
ok(ut_cmp_equal(is.null(sres$vals[[4]]), TRUE),
  "meta_scatter - colour cascade still starts unset")

invisible(lapply(restore_fns, function(f) f()))
