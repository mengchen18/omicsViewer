# Render/echo stability harness for the store-driven scatter stack
# (promoted from HANDOVER.md Appendix A).
#
# A "paint" is one evaluation of the params reactive that a VISIBLE
# renderPlotly would consume:
#   - before the WP1 render barrier: every evaluation of the raw params
#     reactive (checkpoint-gated) - the historic baseline semantics;
#   - after the WP1 render barrier (plotly_scatter_module gains a
#     reactive_ready argument and commits params into a reactiveVal): the
#     counter applies the same gates the module's commit observer applies
#     (checkpoint + reactive_ready) and skips evaluations whose params are
#     identical() to the last counted one, mirroring the reactiveVal
#     dedupe a visible renderPlotly sits behind.
# The mode is chosen by inspecting the module's formals, so this file
# reproduces the documented 11/35 baseline on pre-WP1 code unchanged.
#
# Rows owned by WP4 (store hygiene), WP5 (epoch scoping) and WP6
# (renderUI rebuilds) are printed but not asserted until those work
# packages land; flip their assert = FALSE to TRUE with the matching WP.

suppressMessages({
  library(shiny)
  library(Biobase)
  library(unittest, quietly = TRUE)
  pkgload::load_all(".", quiet = TRUE, attach_testthat = FALSE)
})
dat <- readRDS("inst/extdata/demo.RDS")
fd <- fData(dat); pd <- pData(dat); ex <- exprs(dat)
`%||%` <- function(a, b) if (is.null(a)) b else a

PH <- new.env()
PH$paints <- list(); PH$q <- list(); PH$keep <- list(); PH$ui <- 0L; PH$rows <- list()

## — instrumentation 1: eager paint counter on every plotly_scatter_module instance.
##   The observer MUST stay referenced or it is garbage collected between flushes.
.scatter_has_barrier <- "reactive_ready" %in% names(formals(
  get("plotly_scatter_module", envir = asNamespace("omicsViewer"))))
local({
  orig <- omicsViewer:::plotly_scatter_module
  assignInNamespace("plotly_scatter_module", function(id, reactive_param_plotly_scatter,
      reactive_regLine = reactive(FALSE), reactive_checkpoint = reactive(TRUE),
      htest_var1 = reactive(NULL), htest_var2 = reactive(NULL),
      reactive_ready = reactive(TRUE)) {
    full_id <- getDefaultReactiveDomain()$ns(id)
    last_p <- NULL
    PH$keep[[length(PH$keep) + 1L]] <- observe({
      if (!isTRUE(tryCatch(reactive_checkpoint(), error = function(e) FALSE))) return()
      if (.scatter_has_barrier &&
          !isTRUE(tryCatch(reactive_ready(), error = function(e) FALSE))) return()
      p <- tryCatch(reactive_param_plotly_scatter(),
                    shiny.silent.error = function(e) NULL, error = function(e) NULL)
      if (is.null(p)) return()
      if (.scatter_has_barrier && identical(p, last_p)) return()
      last_p <<- p
      PH$paints[[length(PH$paints) + 1L]] <- list(id = full_id,
        x = attr(p$x, "label") %||% "", y = attr(p$y, "label") %||% "",
        rects = length(p$rect))
    })
    if (.scatter_has_barrier)
      orig(id, reactive_param_plotly_scatter, reactive_regLine, reactive_checkpoint,
           htest_var1, htest_var2, reactive_ready)
    else
      orig(id, reactive_param_plotly_scatter, reactive_regLine, reactive_checkpoint,
           htest_var1, htest_var2)
  }, ns = "omicsViewer")
})

## — instrumentation 2: queue every server->client select update with its FULL id.
##   Patch the package's imports env (not shiny's namespace), and force(orig): a lazy
##   promise would resolve to the wrapper itself and recurse.
local({
  imp <- parent.env(asNamespace("omicsViewer"))
  for (f in c("updateSelectInput", "updateSelectizeInput")) local({
    fn <- f; orig <- get(fn, envir = imp); force(orig)
    if (bindingIsLocked(fn, imp)) unlockBinding(fn, imp)
    assign(fn, function(session, inputId, ..., selected = NULL) {
      if (is.character(selected) && length(selected) == 1L)
        PH$q[[length(PH$q) + 1L]] <- list(id = session$ns(inputId), value = selected)
      orig(session = session, inputId = inputId, ..., selected = selected)
    }, envir = imp)
  })
})

## — instrumentation 3: count plot-container rebuilds (each renderUI run calls one *_ui)
for (f in c("plotly_scatter_ui", "plotly_boxplot_ui", "plot_roc_pr_ui")) local({
  fn <- f; orig <- get(fn, envir = asNamespace("omicsViewer"))
  assignInNamespace(fn, function(id, ...) { PH$ui <- PH$ui + 1L; orig(id, ...) },
                    ns = "omicsViewer")
})
PH$a4 <- 0L
local({
  orig <- get("attr4selector_ui", envir = asNamespace("omicsViewer"))
  assignInNamespace("attr4selector_ui", function(id, ...) { PH$a4 <- PH$a4 + 1L; orig(id, ...) },
                    ns = "omicsViewer")
})

## — browser model: acknowledge queued widget updates in round trips.
##    split : plain selectInputs ack together; server-side selectize (`variable') next trip
##    serial: one widget per round trip (worst case)
ph_ack <- function(session, model = "split", rounds = 30) {
  for (r in seq_len(rounds)) {
    q <- PH$q; PH$q <- list()
    # the browser applies the WHOLE batch to the DOM first (last write per widget wins),
    # then reports the resulting values — never an older message after a newer reply
    final <- list(); for (m in q) final[[m$id]] <- m$value
    final <- final[!vapply(names(final), function(id)
      identical(session$input[[id]], final[[id]]), logical(1))]
    if (!length(final)) break
    ids <- names(final)
    groups <- if (model == "serial") as.list(ids) else
      Filter(length, list(ids[!grepl("-variable$", ids)], ids[grepl("-variable$", ids)]))
    for (g in groups) {
      do.call(session$setInputs, final[g]); session$flushReact()
    }
  }
  for (i in 1:5) session$flushReact()
}
## A real browser reports every input once on connect; triselectors wait for that.
ph_init_inputs <- function(session, prefixes, tris) {
  v <- list()
  for (p in prefixes) {
    v[[paste0(p, "-axisMode")]] <- "quick"
    for (a in c("a4selector", "a4_gf", "a4_gp")) {
      v[[paste0(p, "-", a, "-xcut")]] <- "log10(2)"
      v[[paste0(p, "-", a, "-ycut")]] <- "-log10(0.05)"
      v[[paste0(p, "-", a, "-scorner")]] <- "None"
    }
  }
  for (t in tris) for (w in c("analysis", "subset", "variable")) v[[paste0(t, "-", w)]] <- ""
  do.call(session$setInputs, v)
  for (i in 1:5) session$flushReact()
}
A4 <- c("selectColorUI", "selectShapeUI", "selectSizeUI", "selectTooltipUI", "selectSearchCol")
ph_row <- function(scenario, model, measured, target, detail = "", assert = TRUE) {
  ok <- if (is.function(target)) target(measured) else identical(as.integer(measured), as.integer(target))
  if (isTRUE(assert))
    unittest::ok(ok, sprintf("%s [%s]", scenario, model))
  PH$rows[[length(PH$rows) + 1L]] <- data.frame(scenario = scenario, model = model,
    measured = as.character(measured),
    target = if (is.function(target)) attr(target, "label") else as.character(target),
    pass = ok, assert = assert, detail = detail, stringsAsFactors = FALSE)
}
ph_paints <- function(id) Filter(function(p) identical(p$id, id), PH$paints)
ph_desc <- function(pp) paste(vapply(pp, function(p) sprintf("[%s | %s | r%d]",
  sub(".*\\|", "", p$x), sub(".*\\|", "", p$y), p$rects), ""), collapse = " ")
triple <- function(st, k) paste(unlist(store_read(st, paste0(k, c("_analysis", "_subset", "_variable")))), collapse = "|")
## ----------------------------------------------------------------- A. quick-view switches
scenario_quick <- function(model) {
  store <- widget_store_new(); st <- widget_store_child(store, "dataspace.feature_space")
  testServer(function(input, output, session)
  meta_scatter_module("df", reactive_meta = reactive(fd), reactive_expr = reactive(ex),
    combine = "feature", source = "df",
    reactive_x = reactive("ttest|RE_vs_ME|mean.diff"),
    reactive_y = reactive("ttest|RE_vs_ME|log.fdr"), store = st), {
  ph_init_inputs(session, "df", c("df-tris_main_scatter1", "df-tris_main_scatter2",
                                  paste0("df-a4selector-", A4)))
  ph_ack(session, model)
  go <- function(label, x, y) {
    PH$paints <- list(); PH$q <- list()
    xs <- strsplit(x, "|", fixed = TRUE)[[1]]; ys <- strsplit(y, "|", fixed = TRUE)[[1]]
    store_apply(st, list(x_analysis = xs[1], x_subset = xs[2], x_variable = xs[3],
      y_analysis = ys[1], y_subset = ys[2], y_variable = ys[3], axis_mode = "quick"),
      origin = "system")
    session$flushReact(); ph_ack(session, model)
    pp <- ph_paints("df-main_scatterOutput")
    ph_row(paste("quick:", label), model, length(pp), 1L, ph_desc(pp))
    shown <- if (length(pp)) pp[[length(pp)]] else list(x = "", y = "")
    ph_row(paste("quick:", label, "-> final figure is the target"), model,
           identical(c(shown$x, shown$y), c(x, y)), TRUE,
           sprintf("shown x=%s y=%s", shown$x, shown$y))
  }
  go("volcano RE_vs_ME -> volcano MT_vs_WT", "ttest|MT_vs_WT|mean.diff", "ttest|MT_vs_WT|log.fdr")
  go("volcano -> Cor|MDR", "Cor|MDR|R", "Cor|MDR|logP")
  go("Cor|MDR -> Cor|Doubleing.Time", "Cor|Doubleing.Time|R", "Cor|Doubleing.Time|logP")
  go("Cor -> volcano MT_vs_WT", "ttest|MT_vs_WT|mean.diff", "ttest|MT_vs_WT|log.fdr")
  go("volcano -> PCA", "PCA|All|PC1(10.5%)", "PCA|All|PC2(7.2%)")
})
}

## ----------------------------------------------------------------- B. custom-mode edits
scenario_custom <- function(model) {
  store <- widget_store_new(); st <- widget_store_child(store, "dataspace.feature_space")
  testServer(function(input, output, session)
  meta_scatter_module("df", reactive_meta = reactive(fd), reactive_expr = reactive(ex),
    combine = "feature", source = "df",
    reactive_x = reactive("ttest|RE_vs_ME|mean.diff"),
    reactive_y = reactive("ttest|RE_vs_ME|log.fdr"), store = st), {
  ph_init_inputs(session, "df", c("df-tris_main_scatter1", "df-tris_main_scatter2",
                                  paste0("df-a4selector-", A4)))
  ph_ack(session, model)
  edit <- function(label, w, value, target) {
    PH$paints <- list(); PH$q <- list()
    do.call(session$setInputs, stats::setNames(list(value), paste0("df-tris_main_scatter1-", w)))
    ph_ack(session, model)
    pp <- ph_paints("df-main_scatterOutput")
    ph_row(paste("custom:", label), model, length(pp), target,
           paste("store x =", triple(st, "x"), ph_desc(pp)))
  }
  edit("x variable -> log.pvalue", "variable", "log.pvalue", 1L)
  edit("x subset -> MT_vs_WT", "subset", "MT_vs_WT", 1L)
  tgt <- function(n) n >= 1L; attr(tgt, "label") <- ">=1 and store follows"
  edit("x analysis -> Cor", "analysis", "Cor", tgt)
  ok_store <- startsWith(triple(st, "x"), "Cor|")
  ph_row("custom: store follows analysis edit", model, ok_store, TRUE, triple(st, "x"))
})
}

## ----------------------------------------------------------------- C. right panel rebuilds
scenario_feature_general <- function() {
  store <- widget_store_new(); st <- widget_store_child(store, "resultspace.feature_general")
  ri <- NULL
  cat_col <- strsplit("General|All|Cell.line", "|", fixed = TRUE)[[1]]
  testServer(function(input, output, session) {
    ri <<- reactiveVal(c(1L, 2L))
    feature_general_module("fg", reactive_expr = reactive(ex), reactive_i = ri,
      reactive_phenoData = reactive(pd), reactive_featureData = reactive(fd), store = st)
  }, {
    session$setInputs(`fg-internal_radio` = "Bees")
    ph_init_inputs(session, character(0), c("fg-tris_feature_general", paste0("fg-a4_gf-", A4)))
    for (i in 1:5) { session$flushReact(); try(session$output$`fg-feature_general_plot`, silent = TRUE) }
    # no link variable picked: the settled-triple module output is NULL, and
    # the panel must still render its fallback view (relative-abundance
    # boxplot), exactly as the historical "--select--" placeholder did -
    # a NULL v1() that aborts the pheno chain leaves a BLANK right panel
    fallback <- try(session$output$`fg-feature_general_plot`, silent = TRUE)
    ph_row("feature_general: no link variable -> fallback view renders", "-",
           !is.null(fallback), TRUE)
    do.call(session$setInputs, stats::setNames(as.list(cat_col),
      paste0("fg-tris_feature_general-", c("analysis", "subset", "variable"))))
    for (i in 1:5) { session$flushReact(); try(session$output$`fg-feature_general_plot`, silent = TRUE) }
    PH$ui <- 0L; PH$paints <- list()
    ri(c(3L, 4L))
    for (i in 1:5) { session$flushReact(); try(session$output$`fg-feature_general_plot`, silent = TRUE) }
    # WP6 target: renderUI must not rebuild on an unchanged view type
    ph_row("feature_general: selection change, same view type -> container rebuilds",
           "-", PH$ui, 0L, assert = FALSE)
    ph_row("feature_general: selection change, same view type -> paints", "-",
           length(ph_paints("fg-feature_general_beeswarm")), 1L)
    # a restore of a different analysis, acknowledged one widget per round trip
    PH$q <- list()
    store_apply(st, list(xax_analysis = "PCA", xax_subset = "All",
                         xax_variable = "PC1(10.5%)"), origin = "restore")
    session$flushReact(); ph_ack(session, "serial")
    ol <- store$override_log
    ph_row("component-set module: restore logs no false user override", "serial",
           length(ol), 0L,
           if (length(ol)) sprintf("%s pending=%s user=%s", ol[[1]]$id, ol[[1]]$pending, ol[[1]]$user) else "")
    # WP4 target: seeds must not leave un-acknowledgeable pending entries
    left <- names(Filter(Negate(is.null), store$pending))
    ph_row("seed leaves no un-acknowledgeable pending", "-", length(left), 0L,
           paste(left, collapse = ", "))
  })
}

## ----------------------------------------------------------------- C2. sample_general rebuilds
scenario_sample_general <- function() {
  store <- widget_store_new(); st <- widget_store_child(store, "resultspace.sample_general")
  testServer(function(input, output, session)
  sample_general_module("sg", reactive_phenoData = reactive(pd), reactive_expr = reactive(ex),
    reactive_j = reactive(colnames(ex)[1:5]), store = st), {
    set <- function(col) do.call(session$setInputs, stats::setNames(
      as.list(strsplit(col, "|", fixed = TRUE)[[1]]),
      paste0("sg-tris_sample_general-", c("analysis", "subset", "variable"))))
    pump <- function() for (i in 1:5) {
      session$flushReact(); try(session$output$`sg-sample_general_plot`, silent = TRUE) }
    set("General|All|MDR"); pump()
    PH$ui <- 0L; PH$a4 <- 0L
    set("General|All|Doubleing.Time"); pump()  # numeric -> numeric: same view type
    # WP6 targets
    ph_row("sample_general: link change, same view type -> container rebuilds", "-", PH$ui, 0L,
           assert = FALSE)
    ph_row("sample_general: link change, same view type -> attr4 panel rebuilds", "-", PH$a4, 0L,
           assert = FALSE)
  })
}

## ----------------------------------------------------------------- D. store fan-out
scenario_fanout <- function() {
  store <- widget_store_new()
  s <- lapply(c(df = "dataspace.feature_space", ds = "dataspace.sample_space",
                fg = "resultspace.feature_general", sg = "resultspace.sample_general"),
              function(p) widget_store_child(store, p))
  testServer(function(input, output, session) {
    meta_scatter_module("df", reactive_meta = reactive(fd), reactive_expr = reactive(ex),
      combine = "feature", source = "df", reactive_x = reactive("ttest|RE_vs_ME|mean.diff"),
      reactive_y = reactive("ttest|RE_vs_ME|log.fdr"), store = s$df)
    meta_scatter_module("ds", reactive_meta = reactive(pd), reactive_expr = reactive(ex),
      combine = "pheno", source = "ds", reactive_x = reactive("PCA|All|PC1(10.5%)"),
      reactive_y = reactive("PCA|All|PC2(7.2%)"), store = s$ds)
    feature_general_module("fg", reactive_expr = reactive(ex), reactive_i = reactive(5L),
      reactive_phenoData = reactive(pd), reactive_featureData = reactive(fd), store = s$fg)
    sample_general_module("sg", reactive_phenoData = reactive(pd), reactive_expr = reactive(ex),
      reactive_j = reactive(colnames(ex)[1:5]), store = s$sg)
  }, {
    session$setInputs(`fg-internal_radio` = "Bees")
    tris <- c(paste0(rep(c("df", "ds"), each = 2), "-tris_main_scatter", 1:2),
              "fg-tris_feature_general", "sg-tris_sample_general",
              as.vector(outer(c("df-a4selector", "ds-a4selector", "fg-a4_gf", "sg-a4_gp"), A4, paste, sep = "-")))
    ph_init_inputs(session, c("df", "ds", "fg", "sg"), tris)
    ph_ack(session, "split")
    PH$q <- list()
    store_apply(s$df, list(x_subset = "MT_vs_WT", y_subset = "MT_vs_WT"), origin = "system")
    session$flushReact()
    owner <- sub("-x", "", vapply(PH$q, `[[`, "", "id"))
    # WP5 target: epoch scoping must keep uninvolved modules quiet
    ph_row("fan-out: updates sent to modules NOT involved in the switch", "-",
           sum(owner != "df"), 0L,
           paste(names(table(owner)), table(owner), sep = "=", collapse = ", "),
           assert = FALSE)
  })
}

## ------------------------------------------------- A2. selection persistence
# The volcano corner auto-selection must follow quick-view switches: the
# new view's corner genes replace the old selection, leaving the volcano
# (e.g. to a correlation view) clears it, and re-entering a volcano
# re-engages it. These invariants regressed repeatedly (stale-echo
# adoptions, GC-swallowed clears, ignoreNULL-dead clear branches) while
# every paint-count row stayed green - selection state needs its own rows.
scenario_selection <- function(model) {
  store <- widget_store_new(); st <- widget_store_child(store, "dataspace.feature_space")
  corner_genes <- function(sub) {
    md <- fd[[sprintf("ttest|%s|mean.diff", sub)]]
    lf <- fd[[sprintf("ttest|%s|log.fdr", sub)]]
    ok <- !is.na(md) & !is.na(lf)
    rn <- rownames(fd)
    sort(unique(c(rn[ok & md > log10(2) & lf > -log10(0.05)],
                 rn[ok & md < -log10(2) & lf > -log10(0.05)])))
  }
  testServer(function(input, output, session) {
    s_fig <- meta_scatter_module("df", reactive_meta = reactive(fd), reactive_expr = reactive(ex),
      combine = "feature", source = "df",
      reactive_x = reactive("ttest|RE_vs_ME|mean.diff"),
      reactive_y = reactive("ttest|RE_vs_ME|log.fdr"), store = st)
    exported <<- list(fig = s_fig)
  }, {
    ph_init_inputs(session, "df", c("df-tris_main_scatter1", "df-tris_main_scatter2",
                                    paste0("df-a4selector-", A4)))
    ph_ack(session, model)
    settled_sel <- function() tryCatch(exported$fig()$selected,
      shiny.silent.error = function(e) character(0), error = function(e) character(0))
    switch_to <- function(x, y) {
      PH$paints <- list(); PH$q <- list()
      xs <- strsplit(x, "|", fixed = TRUE)[[1]]; ys <- strsplit(y, "|", fixed = TRUE)[[1]]
      store_apply(st, list(x_analysis = xs[1], x_subset = xs[2], x_variable = xs[3],
        y_analysis = ys[1], y_subset = ys[2], y_variable = ys[3], axis_mode = "quick"),
        origin = "system")
      session$flushReact(); ph_ack(session, model)
    }
    sel_row <- function(label, expected) {
      got <- settled_sel()
      ok <- if (length(expected)) identical(sort(got), expected) else length(got) == 0L
      ph_row(paste("selection:", label), model, ok, TRUE,
        paste0("n=", length(got), " expected=", length(expected)))
    }
    # load view settles first
    ph_ack(session, model)
    sel_row("load (volcano RE_vs_ME)", corner_genes("RE_vs_ME"))
    switch_to("ttest|RE_vs_LE|mean.diff", "ttest|RE_vs_LE|log.fdr")
    sel_row("volcano -> volcano RE_vs_LE", corner_genes("RE_vs_LE"))
    switch_to("ttest|RE_vs_ME|mean.diff", "ttest|RE_vs_ME|log.fdr")
    sel_row("back to volcano RE_vs_ME", corner_genes("RE_vs_ME"))
    switch_to("Cor|MDR|R", "Cor|MDR|logP")
    sel_row("leaving volcano clears", character(0))
    switch_to("ttest|RE_vs_LE|mean.diff", "ttest|RE_vs_LE|log.fdr")
    sel_row("re-entry re-engages (RE_vs_LE)", corner_genes("RE_vs_LE"))
  })
}

for (m in c("split", "serial")) { scenario_quick(m); scenario_custom(m); scenario_selection(m) }
scenario_feature_general()
scenario_sample_general()
scenario_fanout()

res <- do.call(rbind, PH$rows)
op <- options(width = 250)
print(res[, c("pass", "measured", "target", "model", "scenario")], right = FALSE, row.names = FALSE)
cat("\nDetails:\n"); for (i in seq_len(nrow(res))) if (nzchar(res$detail[i]))
  cat(sprintf("  %-62s %-6s %s\n", res$scenario[i], res$model[i], res$detail[i]))
cat(sprintf("\n%d / %d targets met (%d asserted)\n", sum(res$pass), nrow(res), sum(res$assert)))
