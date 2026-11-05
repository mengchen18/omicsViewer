# WP-0 unit suite for the result-space panel registry socket
# (R/auxi_panelRegistry.R; PANEL_DEBT_PLAN.md sections 3-4, WP-0).
#
# Covers: panel_record validation (malformed records rejected with
# actionable errors), the seed registry content (verbatim names from the
# plan appendix), rs_tab_set availability truth tables against the real
# demo-data conditions, panel_context field presence/types/defaults,
# panel_keep retention (the deterministic gc() repro from AGENTS.md /
# plan section 5.4), panel_triset_bind equivalence against a verbatim
# copy of today's module_ora.R cascade block plus functional UI->store
# coverage (incl. allow_unset), and the D12 enumeration freeze: the
# agent-visible store-key enumeration of the REAL result-space
# composition booted from demo.RDS. Every later panel WP asserts the
# frozen vector UNCHANGED; only WP-3 may grow it (by exactly the dose
# keys, in its own reviewed commit).
#
# The package is NOT reinstalled for this suite; load the working tree
# (the agentChannel / scatterDefaults convention).

library(unittest, quietly = TRUE)

if (!requireNamespace("pkgload", quietly = TRUE)) {
  ok(TRUE, "panel registry tests skipped because pkgload is unavailable")
  quit(save = "no", status = 0)
}

pkg_root <- if (file.exists("DESCRIPTION")) "." else
  if (file.exists(file.path("..", "DESCRIPTION"))) ".." else "."
pkgload::load_all(pkg_root, quiet = TRUE, attach_testthat = FALSE)

suppressMessages({
  library(shiny)
  library(Biobase)
})

panel_record <- omicsViewer:::panel_record
panel_records <- omicsViewer:::panel_records  # access forces the delayed binding
panel_keep <- omicsViewer:::panel_keep
panel_triset_bind <- omicsViewer:::panel_triset_bind
panel_context <- omicsViewer:::panel_context
rs_tab_set <- omicsViewer:::rs_tab_set
widget_store_new <- omicsViewer:::widget_store_new
widget_store_child <- omicsViewer:::widget_store_child
widget_binding <- omicsViewer:::widget_binding
store_register <- omicsViewer:::store_register
store_watch <- omicsViewer:::store_watch
store_epoch <- omicsViewer:::store_epoch
store_read <- omicsViewer:::store_read
store_apply <- omicsViewer:::store_apply
store_reset <- omicsViewer:::store_reset
store_restore <- omicsViewer:::store_restore
store_bind_triselector <- omicsViewer:::store_bind_triselector
store_registry_view <- omicsViewer:::store_registry_view
trisetter <- omicsViewer:::trisetter

.err <- function(expr) tryCatch(force(expr), error = function(e) conditionMessage(e))

dat <- readRDS(system.file("extdata", "demo.RDS", package = "omicsViewer"))
fd0 <- Biobase::fData(dat)
fd0_no_opt <- fd0[, !grepl("^ResponseCurve\\||^StringDB\\||^SeqLogo\\|",
                           colnames(fd0)), drop = FALSE]
fd0_gs <- fd0
attr(fd0_gs, "GS") <- data.frame(featureId = "g1", gsId = "gs1")

## ----------------------------------------------- panel_record validation ----
ok(grepl("module_id", .err(panel_record(
     "", "x.y", "T", "L", "i", 1, function(fd, pd) TRUE, identity, identity))),
   "panel_record: empty module_id rejected, error names the field")
ok(grepl("module_id", .err(panel_record(
     "a.b", "x.y", "T", "L", "i", 1, function(fd, pd) TRUE, identity, identity))),
   "panel_record: dotted module_id rejected (single-segment namespace ids)")
ok(grepl("store_id", .err(panel_record(
     "a", "xy", "T", "L", "i", 1, function(fd, pd) TRUE, identity, identity))),
   "panel_record: non-dotted store_id rejected")
ok(grepl("store_id", .err(panel_record(
     "a", ".leading.dot", "T", "L", "i", 1, function(fd, pd) TRUE, identity, identity))),
   "panel_record: malformed store_id rejected")
ok(grepl("tab_value", .err(panel_record(
     "a", "x.y", NA, "L", "i", 1, function(fd, pd) TRUE, identity, identity))),
   "panel_record: NA tab_value rejected")
ok(grepl("label", .err(panel_record(
     "a", "x.y", "T", 1, "i", 1, function(fd, pd) TRUE, identity, identity))),
   "panel_record: non-character label rejected")
ok(grepl("icon", .err(panel_record(
     "a", "x.y", "T", "L", "", 1, function(fd, pd) TRUE, identity, identity))),
   "panel_record: empty icon rejected")
ok(grepl("order", .err(panel_record(
     "a", "x.y", "T", "L", "i", 1.5, function(fd, pd) TRUE, identity, identity))),
   "panel_record: fractional order rejected")
ok(grepl("order", .err(panel_record(
     "a", "x.y", "T", "L", "i", c(1L, 2L), function(fd, pd) TRUE, identity, identity))),
   "panel_record: vector order rejected")
ok(grepl("requires", .err(panel_record(
     "a", "x.y", "T", "L", "i", 1, "not-a-fun", identity, identity))),
   "panel_record: non-function requires rejected")
ok(grepl("ui", .err(panel_record(
     "a", "x.y", "T", "L", "i", 1, function(fd, pd) TRUE, NULL, identity))),
   "panel_record: non-function ui rejected")
ok(grepl("server", .err(panel_record(
     "a", "x.y", "T", "L", "i", 1, function(fd, pd) TRUE, identity, NULL))),
   "panel_record: non-function server rejected")
ok(grepl("stateful", .err(panel_record(
     "a", "x.y", "T", "L", "i", 1, function(fd, pd) TRUE, identity, identity,
     stateful = NA))),
   "panel_record: NA stateful rejected")
ok(grepl("status_key", .err(panel_record(
     "a", "x.y", "T", "L", "i", 1, function(fd, pd) TRUE, identity, identity,
     status_key = ""))),
   "panel_record: empty status_key rejected")
r_ok <- panel_record(
  "a", "x.y", "T", "L", "i", 2, function(fd, pd) TRUE, identity, identity,
  stateful = FALSE)
ok(is.list(r_ok) && r_ok$order == 2L && !r_ok$stateful && is.null(r_ok$status_key),
   "panel_record: valid record returned as a plain validated list")

## -------------------------------------------------- seed registry content ----
ok(identical(names(panel_records),
             c("feature_general", "ora", "fgsea", "rescurve", "stringdb",
               "ptm", "geneshotTab", "sample_general")),
   "registry: 8 records keyed by module id in composition order")
ok(identical(unname(vapply(panel_records, function(r) r$tab_value, "")), c(
  "Feature", "ORA", "fGSEA", "Response", "StringDB", "SeqLogo",
  "Geneshot", "Sample")),
   "registry: tab values verbatim, registry order == today's composition order")
ok(identical(unname(vapply(panel_records, function(r) r$store_id, "")), c(
  "resultspace.feature_general", "resultspace.ora", "resultspace.fgsea",
  "resultspace.dose_response", "resultspace.stringdb", "resultspace.ptm",
  "resultspace.geneshot", "resultspace.sample_general")),
   "registry: store prefixes verbatim incl. the NEW resultspace.dose_response (D2/D18)")
ok(identical(unname(vapply(panel_records, function(r) r$status_key, "")), c(
  "analyst_feature_general", "analyst_ora", "analyst_fgsea",
  "analyst_dose_response", "analyst_stringdb", "analyst_ptm",
  "analyst_gene_shot", "analyst_sample_general")),
   "registry: status keys verbatim (analyst_gene_shot != module id, historical)")
ok(identical(unname(vapply(panel_records, function(r) r$order, 1L)), 1:8),
   "registry: order field 1..8")
ok(identical(
  unname(vapply(panel_records, function(r) isTRUE(r$stateful), TRUE)),
  c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE)),
   "registry: stateful flags - dose/STRING FALSE (current semantics; fixes are WP-3/4)")
ok(identical(
  unname(vapply(panel_records, function(r) r$icon, "")),
  c("chart-column", "object-group", "arrow-down-short-wide", "capsules",
    "share-nodes", "font", "book-open", "arrows-down-to-people")),
   "registry: icons match today's L1 navbar composition verbatim")
ok(identical(panel_records$ora$ui, omicsViewer:::enrichment_analysis_ui) &&
   identical(panel_records$ora$server, omicsViewer:::enrichment_analysis_module) &&
   identical(panel_records$geneshotTab$server, omicsViewer:::geneshot_module) &&
   identical(panel_records$rescurve$server, omicsViewer:::dose_response_module) &&
   identical(panel_records$sample_general$ui, omicsViewer:::sample_general_ui),
   "registry: ui/server are the module functions themselves")

## predicates replicate today's conditions exactly (NULL-safe)
ok(isTRUE(panel_records$feature_general$requires(NULL, NULL)) &&
   isTRUE(panel_records$geneshotTab$requires(NULL, NULL)) &&
   isTRUE(panel_records$sample_general$requires(NULL, NULL)),
   "registry: pinned/always-on panels available even with fd NULL")
ok(isFALSE(panel_records$ora$requires(NULL, NULL)) &&
   isFALSE(panel_records$ora$requires(fd0, NULL)) &&
   isTRUE(panel_records$ora$requires(fd0_gs, NULL)) &&
   isTRUE(panel_records$fgsea$requires(fd0_gs, NULL)),
   "registry: ORA/fGSEA require a non-NULL fd carrying the GS attr")
ok(isFALSE(panel_records$rescurve$requires(fd0_no_opt, NULL)) &&
   isTRUE(panel_records$rescurve$requires(fd0, NULL)) &&
   isFALSE(panel_records$stringdb$requires(fd0_no_opt, NULL)) &&
   isTRUE(panel_records$stringdb$requires(fd0, NULL)) &&
   isFALSE(panel_records$ptm$requires(fd0_no_opt, NULL)) &&
   isTRUE(panel_records$ptm$requires(fd0, NULL)),
   "registry: Response/StringDB/SeqLogo require their column families")

## ------------------------------------------------- rs_tab_set truth table ----
# demo data: no GS attr; 32 ResponseCurve, 1 StringDB, 1 SeqLogo columns
# (verified 2026-11-05, plan section 3.2 parity obligation)
ok(ut_cmp_identical(
  rs_tab_set(fd0),
  c("Feature", "Response", "StringDB", "SeqLogo", "Geneshot", "Sample")),
  "rs_tab_set: demo fd truth table")
ok(ut_cmp_identical(
  rs_tab_set(fd0_gs),
  c("Feature", "ORA", "fGSEA", "Response", "StringDB", "SeqLogo",
    "Geneshot", "Sample")),
  "rs_tab_set: synthetic GS fd adds ORA/fGSEA in positions 2-3")
ok(ut_cmp_identical(
  rs_tab_set(NULL), c("Feature", "Geneshot", "Sample")),
  "rs_tab_set: NULL fd degrades to Feature/Geneshot/Sample")
ok(ut_cmp_identical(
  rs_tab_set(fd0_no_opt), c("Feature", "Geneshot", "Sample")),
  "rs_tab_set: fd without optional column families degrades")
ok(ut_cmp_identical(
  rs_tab_set(fd0, c("C1", "C2")),
  c("Feature", "Response", "StringDB", "SeqLogo", "C1", "C2",
    "Geneshot", "Sample")),
  "rs_tab_set: additional tabs compose between core panels and the pinned pair")
ok(ut_cmp_identical(
  rs_tab_set(fd0, c("Sample", "Geneshot", "C1")),
  c("Feature", "Response", "StringDB", "SeqLogo", "Sample", "Geneshot", "C1")),
  "rs_tab_set: unique() semantics for duplicated additional names")
ok(ut_cmp_identical(
  rs_tab_set(fd0, "Feature"),
  c("Feature", "Response", "StringDB", "SeqLogo", "Geneshot", "Sample")),
  "rs_tab_set: additional duplicating a core tab is uniqued away")
ok(ut_cmp_identical(
  rs_tab_set(fd0, NULL), rs_tab_set(fd0)),
  "rs_tab_set: NULL additional handled like the default")
ok(ut_cmp_identical(
  rs_tab_set(fd0),
  unique(c("Feature", "Response", "StringDB", "SeqLogo", "Geneshot", "Sample"))),
  "rs_tab_set: byte-equal to today's L0 agent_analysis_tabs vector for demo data")
# a failed fd read (hidden tab / erroring reactive) must degrade, not throw
.errfd <- structure(list(1), class = c("panel_registry_errfd", "list"))
dimnames.panel_registry_errfd <- function(x) stop("boom")
ok(ut_cmp_identical(
  rs_tab_set(.errfd), c("Feature", "Geneshot", "Sample")),
  "rs_tab_set: erroring fd read degrades to Feature/Geneshot/Sample (guarded predicates)")

## ------------------------------------------------------ panel_context ----
ctx <- panel_context("ora")
ok(setequal(
  names(ctx),
  c("id", "tab_value", "data", "selection", "store", "status", "keep",
    "triset", "table", "logger", "session", "active")),
  "panel_context: every section-3.3 field present (and nothing more)")
ok(identical(ctx$id, "ora") && is.null(ctx$tab_value) && is.null(ctx$store) &&
   is.null(ctx$session) && identical(ctx$data, list()) &&
   identical(ctx$selection, list()),
   "panel_context: identity/store/session/data/selection defaults")
ok(is.null(shiny::isolate(ctx$status())),
   "panel_context: status defaults to reactive(NULL)")
ok(is.null(shiny::isolate(ctx$active())),
   "panel_context: active defaults to reactive(NULL) without active_tab")
ok(is.function(ctx$logger) && is.null(ctx$logger("any_event", list(a = 1))) &&
   is.null(ctx$logger()),
   "panel_context: logger defaults to the dormant no-op (D14) - accepts anything")
ok(is.function(ctx$keep) && !withVisible(ctx$keep(1))$visible,
   "panel_context: keep is a closure returning its observer invisibly")
ok(!identical(ctx$keep, panel_context("ora")$keep),
   "panel_context: keep is per-instance (fresh panel_keep() per context, D3)")
ok(grepl("data", .err(panel_context("x", data = list(bad = reactive(1))))),
   "panel_context: unknown data element rejected with the field named")
ok(grepl("selection$feature", .err(panel_context("x", selection = list(feature = 1))),
   fixed = TRUE),
   "panel_context: non-reactive selection element rejected")
ok(grepl("status", .err(panel_context("x", status = NULL))),
   "panel_context: non-reactive status rejected")
ok(grepl("active_tab", .err(panel_context("x", active_tab = "nope"))),
   "panel_context: non-reactive active_tab rejected")
ok(grepl("tab_value", .err(panel_context("x", active_tab = reactive("ORA")))),
   "panel_context: active_tab without tab_value rejected (never a silent logical(0) from active; WP-0 review NIT-a)")
ok(grepl("no store", .err(ctx$triset(
     "t", "l", function() NULL,
     list(analysis = list(label = "a", help = "b"),
          subset = list(label = "a", help = "b"),
          variable = list(label = "a", help = "b"))))),
   "panel_context: triset helper refuses store-less contexts with an actionable error")

## active derivation (reactive flip needs a session domain)
act <- new.env()
shiny::testServer(function(input, output, session) {
  rv <- reactiveVal("StringDB")
  act$ctx <- panel_context("stringdb", tab_value = "StringDB",
                           active_tab = reactive(rv()))
  act$flip <- function() rv("ORA")
}, {
  session$flushReact()
  ok(isTRUE(shiny::isolate(act$ctx$active())),
     "panel_context: active TRUE when the session tab equals tab_value")
  act$flip()
  session$flushReact()
  ok(isFALSE(shiny::isolate(act$ctx$active())),
     "panel_context: active flips FALSE when the session tab changes")
})

## ------------------------------------- panel_keep retention (gc repro) ----
# AGENTS.md observer-GC lore / plan section 5.4: Shiny holds observers
# weakly; glue-style observers whose ONLY reference is the keep closure
# must survive a deterministic gc() and keep firing. (The collection of
# UNkept observers is weak-ref timing dependent and survives in short
# mock sessions, so this asserts the retention guarantee - the direction
# the keep idiom owns.)
#
# Structural assertion (WP-0 review MINOR-2): the gc case above passes
# even for a NON-storing keep (PD-8 - an unkept control observer also
# survived), so the deterministic guarantee is asserted here directly:
# every kept observer must land in the closure's list. The gc case
# below stays as integration documentation of the retention idiom.
keep_struct <- new.env()
shiny::testServer(function(input, output, session) {
  k <- panel_keep()
  v <- shiny::observe({ NULL })
  k(v)
  keep_struct$k <- k
  keep_struct$v <- v
}, {
  session$flushReact()
})
ok(
  identical(get("obs", envir = environment(keep_struct$k)),
            list(keep_struct$v)),
  "panel_keep: structural retention - kept observers land in the closure's list (deterministic; fails for a non-storing keep)"
)

gcenv <- new.env()
gcenv$counts <- rep(0L, 5L)
app_gc <- function(input, output, session) {
  root <- widget_store_new()
  holder <- new.env(parent = emptyenv())
  holder$root <- root
  holder$keeps <- lapply(seq_len(5L), function(i) panel_keep())
  for (i in seq_len(5L)) {
    kp <- holder$keeps[[i]]
    local({
      ii <- i
      ch <- widget_store_child(root, paste0("gctest", ii))
      store_register(ch, widget_binding("x", "string", label = "x"))
      ep <- store_epoch(ch)  # module-level pattern: created once, referenced
      kp(shiny::observe({
        ep()                   # epoch-driven store glue dependency
        gcenv$counts[ii] <<- gcenv$counts[ii] + 1L
      }))
    })
  }
  gcenv$holder <- holder  # observers referenced ONLY through the keeps
}
shiny::testServer(app_gc, {
  session$flushReact()
  gcenv$base <- gcenv$counts
  # drop every reference the test frame holds, then collect twice
  rm(list = ls(all.names = TRUE, envir = environment()), envir = environment())
  gc(); gc()
  for (i in seq_len(5L)) {
    ch <- widget_store_child(gcenv$holder$root, paste0("gctest", i))
    store_apply(ch, list(x = "written"), origin = "system")
  }
  session$flushReact()
})
ok(ut_cmp_identical(gcenv$base, rep(1L, 5L)),
   "panel_keep: glue observers ran once before collection")
ok(ut_cmp_identical(gcenv$counts - gcenv$base, rep(1L, 5L)),
   "panel_keep: kept observers still fire after rm + gc() x2 (retention, P1)")

## ------------------------------- panel_triset_bind vs verbatim ORA block ----
# Verbatim copy of today's module_ora.R store branch (the cascade
# boilerplate) driven over the demo fdata, against the helper over the
# same data: identical registered keys/kinds/labels/help/depends_on and
# identical choices-provider outputs (plan WP-0 test contract).
#
# FROZEN GOLDEN FIXTURE (WP-0 review NIT-c): when WP-9 replaces the real
# module_ora.R cascade with panel_triset_bind, this copy STAYS HERE
# verbatim - the equivalence claim must not go circular against the very
# code it once mirrored. Do not "simplify" it to call the helper.
tsenv <- new.env()
app_ts <- function(input, output, session) {
  reactive_featureData <- reactive(fd0)

  ## ---- VERBATIM module_ora.R (cascade block, store branch) ----
  triset <- reactive({
    trisetter(meta = reactive_featureData(), combine = "none")
  })
  .ora_store_observers <- list()
  .ora_keep <- function(obs) {
    .ora_store_observers[[length(.ora_store_observers) + 1L]] <<- obs
    invisible(obs)
  }
  store <- widget_store_child(widget_store_new(), "resultspace.ora")
  .ora_ts <- function() {
    fd <- tryCatch(reactive_featureData(),
                   shiny.silent.error = function(e) NULL,
                   error = function(e) NULL)
    if (is.null(fd)) return(NULL)
    ts <- tryCatch(trisetter(meta = fd, combine = "none"),
                   shiny.silent.error = function(e) NULL,
                   error = function(e) NULL)
    if (is.null(ts) || !is.matrix(ts) || nrow(ts) == 0L) NULL else ts
  }
  .ora_ts1 <- function() {
    ts <- .ora_ts()
    if (is.null(ts)) character(0) else unique(ts[, 1])
  }
  .ora_ts2 <- function(a) {
    ts <- .ora_ts()
    if (is.null(ts) || is.null(a) || !nzchar(a)) character(0)
    else unique(ts[ts[, 1] %in% a, 2])
  }
  .ora_ts3 <- function(a, b) {
    ts <- .ora_ts()
    if (is.null(ts)) character(0)
    else {
      i <- rep(TRUE, nrow(ts))
      if (!is.null(a) && nzchar(a)) i <- i & ts[, 1] %in% a
      if (!is.null(b) && nzchar(b)) i <- i & ts[, 2] %in% b
      unique(ts[i, 3])
    }
  }
  ko1 <- paste0(store$prefix, ".xax_analysis")
  ko2 <- paste0(store$prefix, ".xax_subset")
  store_register(
    store,
    widget_binding("xax_analysis", "select",
      label = "Collapse category",
      help = paste("Annotation category whose values the input features",
                   "are collapsed on before testing over-representation"),
      choices_provider = function(v) .ora_ts1()),
    widget_binding("xax_subset", "select", label = "Collapse subcategory",
      help = "Subcategory within the collapse category",
      depends_on = "xax_analysis",
      choices_provider = function(v) .ora_ts2(v[[ko1]])),
    widget_binding("xax_variable", "select_cascaded", label = "Collapse variable",
      help = paste("Variable the input features are collapsed on",
                   "(feature annotation column)"),
      depends_on = c("xax_analysis", "xax_subset"),
      choices_provider = function(v) .ora_ts3(v[[ko1]], v[[ko2]]))
  )
  v1 <- triselector_module(
    "tris_ora", reactive_x = triset, label = "Collapse features on",
    reactive_selector1 = store_watch(store, "xax_analysis"),
    reactive_selector2 = store_watch(store, "xax_subset"),
    reactive_selector3 = store_watch(store, "xax_variable"),
    reactive_axis_request = store_epoch(store))
  store_bind_triselector(store,
    keys = c(analysis = "xax_analysis", subset = "xax_subset",
             variable = "xax_variable"),
    sel = v1, keep = .ora_keep)
  tsenv$storeA <- store
  tsenv$vA <- v1

  ## ---- the helper, same panel, same strings ----
  storeB <- widget_store_child(widget_store_new(), "resultspace.ora")
  .keepB <- panel_keep()
  vB <- panel_triset_bind(
    storeB, .keepB, "trisB", "Collapse features on",
    triset_source = function() {
      fd <- tryCatch(reactive_featureData(),
                     shiny.silent.error = function(e) NULL,
                     error = function(e) NULL)
      if (is.null(fd)) return(NULL)
      trisetter(meta = fd, combine = "none")
    },
    bindings = list(
      analysis = list(
        label = "Collapse category",
        help = paste("Annotation category whose values the input features",
                     "are collapsed on before testing over-representation")),
      subset = list(
        label = "Collapse subcategory",
        help = "Subcategory within the collapse category"),
      variable = list(
        label = "Collapse variable",
        help = paste("Variable the input features are collapsed on",
                     "(feature annotation column)"))))
  tsenv$storeB <- storeB
  tsenv$vB <- vB

  ## ---- allow_unset variant (feature/sample-general shape) ----
  storeC <- widget_store_child(widget_store_new(), "resultspace.fgsea")
  vC <- panel_triset_bind(
    storeC, panel_keep(), "trisC", "Input variable",
    triset_source = function() {
      fd <- tryCatch(reactive_featureData(),
                     shiny.silent.error = function(e) NULL,
                     error = function(e) NULL)
      if (is.null(fd)) return(NULL)
      trisetter(meta = fd, combine = "none")
    },
    bindings = list(
      analysis = list(label = "Ranking category", help = "a"),
      subset = list(label = "Ranking subcategory", help = "b"),
      variable = list(label = "Ranking variable", help = "c")),
    allow_unset = TRUE)
  tsenv$storeC <- storeC
  tsenv$vC <- vC
}

ts_demo <- trisetter(meta = fd0, combine = "none")
a1 <- unique(ts_demo[, 1])[1]
a2 <- unique(ts_demo[ts_demo[, 1] == a1, 2])[1]
v3 <- unique(ts_demo[ts_demo[, 1] == a1 & ts_demo[, 2] == a2, 3])[1]

shiny::testServer(app_ts, {
  session$flushReact()
  rootA <- tsenv$storeA$parent
  rootB <- tsenv$storeB$parent
  ok(ut_cmp_identical(names(rootA$bindings), names(rootB$bindings)),
     "panel_triset_bind: registers the same canonical keys in the same order")
  for (k in names(rootA$bindings)) {
    ba <- rootA$bindings[[k]]; bb <- rootB$bindings[[k]]
    ok(identical(ba$kind, bb$kind) && identical(ba$label, bb$label) &&
       identical(ba$help, bb$help) && identical(ba$depends_on, bb$depends_on),
       paste0("panel_triset_bind: binding ", k,
              " has identical kind/label/help/depends_on"))
  }
  k1 <- "resultspace.ora.xax_analysis"
  k2 <- "resultspace.ora.xax_subset"
  k3 <- "resultspace.ora.xax_variable"
  states <- list(
    list(),
    stats::setNames(list(a1), k1),
    stats::setNames(list(a1, a2), c(k1, k2)),
    stats::setNames(list(NULL, a2), c(k1, k2)),
    stats::setNames(list(a1, NULL), c(k1, k2)),
    stats::setNames(list("zzz", "zzz"), c(k1, k2)))
  all_eq <- TRUE; any_nonempty <- FALSE
  for (st in states) {
    for (kk in c(k1, k2, k3)) {
      oa <- shiny::isolate(rootA$bindings[[kk]]$choices_provider(st))
      ob <- shiny::isolate(rootB$bindings[[kk]]$choices_provider(st))
      if (!identical(oa, ob)) all_eq <- FALSE
      if (length(oa) > 0L) any_nonempty <- TRUE
    }
  }
  ok(all_eq,
     "panel_triset_bind: choices providers identical to the verbatim block over all cascade states")
  ok(any_nonempty && length(shiny::isolate(rootB$bindings[[k1]]$choices_provider(list()))) == length(unique(ts_demo[, 1])),
     "panel_triset_bind: providers return the real demo cascade (not the guarded empty fallback)")
  ok(is.function(tsenv$vB) && is.function(tsenv$vA),
     "panel_triset_bind: returns the triselector module value (reactive)")

  ## functional: user picks -> committed triple -> store (UI->store glue)
  session$setInputs(`trisB-analysis` = a1)
  session$flushReact()
  session$setInputs(`trisB-subset` = a2)
  session$flushReact()
  session$setInputs(`trisB-variable` = v3)
  session$flushReact()
  ok(ut_cmp_identical(
    unname(unlist(store_read(tsenv$storeB, c("xax_analysis", "xax_subset", "xax_variable")))),
    c(a1, a2, v3)),
     "panel_triset_bind: settled user picks sync into the store")
  ok(ut_cmp_identical(unname(unlist(tsenv$vB())), c(a1, a2, v3)),
     "panel_triset_bind: module value reports the settled triple")

  ## allow_unset threading: settle, then unset the variable
  session$setInputs(`trisC-analysis` = a1)
  session$flushReact()
  session$setInputs(`trisC-subset` = a2)
  session$flushReact()
  session$setInputs(`trisC-variable` = v3)
  session$flushReact()
  ok(ut_cmp_identical(
    unname(unlist(store_read(tsenv$storeC, c("xax_analysis", "xax_subset", "xax_variable")))),
    c(a1, a2, v3)),
     "panel_triset_bind: allow_unset cascade settles like a fixed one")
  session$setInputs(`trisC-variable` = "--select--")
  session$flushReact()
  ok(is.null(store_read(tsenv$storeC, "xax_variable")[[1]]) &&
     identical(store_read(tsenv$storeC, "xax_analysis")[[1]], a1),
     "panel_triset_bind: allow_unset --select-- clears the variable key, keeps upstream")
  ok(identical(tsenv$vC()$variable, "--select--"),
     "panel_triset_bind: allow_unset commits the explicit unset marker")
})

ok(grepl("tri_id", .err(panel_triset_bind(
     widget_store_child(widget_store_new(), "x"), panel_keep(), NULL, "l",
     function() NULL, list(analysis = list(label = "a", help = "b"),
                            subset = list(label = "a", help = "b"),
                            variable = list(label = "a", help = "b"))))),
   "panel_triset_bind: missing tri_id rejected with an actionable error")
ok(grepl("tri_label", .err(panel_triset_bind(
     widget_store_child(widget_store_new(), "x"), panel_keep(), "t", NULL,
     function() NULL, list(analysis = list(label = "a", help = "b"),
                            subset = list(label = "a", help = "b"),
                            variable = list(label = "a", help = "b"))))),
   "panel_triset_bind: NULL tri_label rejected, validated like tri_id (WP-0 review NIT-b)")
ok(grepl("tri_label", .err(panel_triset_bind(
     widget_store_child(widget_store_new(), "x"), panel_keep(), "t", "",
     function() NULL, list(analysis = list(label = "a", help = "b"),
                            subset = list(label = "a", help = "b"),
                            variable = list(label = "a", help = "b"))))),
   "panel_triset_bind: empty tri_label rejected")
ok(grepl("bindings", .err(panel_triset_bind(
     widget_store_child(widget_store_new(), "x"), panel_keep(), "t", "l",
     function() NULL, list(analysis = list(label = "a", help = "b"))))),
   "panel_triset_bind: incomplete bindings rejected with the field named")
ok(grepl("label", .err(panel_triset_bind(
     widget_store_child(widget_store_new(), "x"), panel_keep(), "t", "l",
     function() NULL,
     list(analysis = list(label = "a"),
          subset = list(label = "a", help = "b"),
          variable = list(label = "a", help = "b"))))),
   "panel_triset_bind: binding without help rejected")
ok(grepl("triset_source", .err(panel_triset_bind(
     widget_store_child(widget_store_new(), "x"), panel_keep(), "t", "l",
     NULL, list(analysis = list(label = "a", help = "b"),
                subset = list(label = "a", help = "b"),
                variable = list(label = "a", help = "b"))))),
   "panel_triset_bind: non-function triset_source rejected")
ok(grepl("store", .err(panel_triset_bind(
     NULL, panel_keep(), "t", "l", function() NULL,
     list(analysis = list(label = "a", help = "b"),
          subset = list(label = "a", help = "b"),
          variable = list(label = "a", help = "b"))))),
   "panel_triset_bind: missing store rejected with an actionable error")

## --------------------------- D12 enumeration freeze (plan section 5.1) ----
# The agent-visible store-key enumeration of the REAL result-space
# composition (L1 instantiates every module server unconditionally, so
# ORA/fGSEA register even though the demo ships no GS attr). This vector
# is the permanent regression gate: every later panel WP asserts it
# UNCHANGED. WP-3 (its own reviewed commit, D12) grew it by EXACTLY the
# two dose keys - what the dose tables' dataTableDownload_module calls
# register. The updated vector was verified byte-equal to a dput of the
# live post-WP-3 composition (72 ids; captured 2026-11-05 from the
# freeze harness below), not hand-typed: the equality assertion IS the
# proof.
PANEL_REGISTRY_FREEZE_IDS <- sort(unique(c(
  "resultspace.analyst_tab",
  "resultspace.dose_response.feature_table",
  "resultspace.dose_response.param_table",
  "resultspace.feature_general.xax_analysis",
  "resultspace.feature_general.xax_subset",
  "resultspace.feature_general.xax_variable",
  "resultspace.feature_general.plot_type",
  "resultspace.feature_general.regression_line",
  paste0("resultspace.feature_general.attr4.", c(
    "color_analysis", "color_subset", "color_variable", "scorner",
    "search_analysis", "search_subset", "search_value",
    "search_variable", "shape_analysis", "shape_subset",
    "shape_variable", "size_analysis", "size_subset", "size_variable",
    "tooltip_analysis", "tooltip_subset", "tooltip_variable",
    "xcut", "ycut")),
  "resultspace.ora.xax_analysis",
  "resultspace.ora.xax_subset",
  "resultspace.ora.xax_variable",
  "resultspace.ora.selected_row",
  "resultspace.fgsea.xax_analysis",
  "resultspace.fgsea.xax_subset",
  "resultspace.fgsea.xax_variable",
  "resultspace.fgsea.selected_row",
  "resultspace.stringdb.taxonomy",
  "resultspace.stringdb.show_labels",
  "resultspace.stringdb.selected_row",
  "resultspace.ptm.xax_analysis",
  "resultspace.ptm.xax_subset",
  "resultspace.ptm.xax_variable",
  "resultspace.geneshot.term",
  "resultspace.geneshot.xax_analysis",
  "resultspace.geneshot.xax_subset",
  "resultspace.geneshot.xax_variable",
  "resultspace.sample_general.xax_analysis",
  "resultspace.sample_general.xax_subset",
  "resultspace.sample_general.xax_variable",
  "resultspace.sample_general.survival_censor",
  paste0("resultspace.sample_general.", c(
    "batch_feature_selected_row", "batch_phenotype_selected_row",
    "batch_show_features", "batch_show_phenotype")),
  paste0("resultspace.sample_general.attr4.", c(
    "color_analysis", "color_subset", "color_variable", "scorner",
    "search_analysis", "search_subset", "search_value",
    "search_variable", "shape_analysis", "shape_subset",
    "shape_variable", "size_analysis", "size_subset", "size_variable",
    "tooltip_analysis", "tooltip_subset", "tooltip_variable",
    "xcut", "ycut")))))

esv <- new.env()
esv$store <- widget_store_new()
app_l1 <- function(input, output, session) {
  omicsViewer:::L1_result_space_module(
    "rs",
    reactive_expr = reactive(Biobase::exprs(dat)),
    reactive_phenoData = reactive(Biobase::pData(dat)),
    reactive_featureData = reactive(fd0),
    status = reactive(NULL),
    store = esv$store)
}
l1_err <- NULL; l1_alive <- NULL
provider_identity <- NULL
shiny::testServer(app_l1, {
  l1_err <<- tryCatch({
    for (i in 1:3) session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  l1_alive <<- !session$isClosed()
  # WP-0 review MINOR-1 guard, in its WP-1 form: the live analyst_tab
  # choices provider is invoked with the real store state. Before WP-1
  # this was compared against rs_tab_set(fd0) (a genuine two-
  # implementation check); WP-1 makes .rs_tab_choices delegate to
  # rs_tab_set, so that comparison went circular and the provider is now
  # pinned against the FROZEN hand-typed golden vector below (captured
  # 2026-11-05 from the pre-WP-1 hand-written provider over demo fd,
  # byte-for-byte). The vector stays as the permanent fixture.
  provider_identity <<- tryCatch({
    root <- esv$store
    identical(
      root$bindings[["resultspace.analyst_tab"]]$choices_provider(
        store_read(root, names(root$bindings))),
      c("Feature", "Response", "StringDB", "SeqLogo", "Geneshot", "Sample"))
  }, error = function(e) conditionMessage(e))
})
ok(is.null(l1_err) && isTRUE(l1_alive),
   "freeze harness: real L1 composition boots headless with a fresh store")
ok(isTRUE(provider_identity),
   "provider identity: the live analyst_tab choices provider equals the frozen golden vector (hand-typed; the rs_tab_set comparison went circular at WP-1 and was replaced - recorded in PANEL_DEBT_STATUS)")
live_ids <- sort(unique(vapply(
  store_registry_view(esv$store), function(r) r$id, character(1))))
# WP-3's D12 update folded the dose keys into the literal vector above
# (captured from the live composition); the gate is plain equality again
ok(ut_cmp_identical(live_ids, PANEL_REGISTRY_FREEZE_IDS),
   "D12 enumeration freeze: agent-visible result-space store keys unchanged (72 ids incl. the two WP-3 dose keys)")
live_ns <- unname(vapply(panel_records, function(r) r$store_id, ""))
ok(all(startsWith(live_ids, "resultspace.")) &&
   all(vapply(live_ns, function(p)
     any(startsWith(live_ids, paste0(p, "."))), TRUE)) &&
   any(startsWith(live_ids, "resultspace.dose_response.")) &&
   any(live_ids == "resultspace.analyst_tab"),
   "freeze cross-check: every registry store prefix is represented; dose present since WP-3 (its child registers unconditionally, like ORA/fGSEA on GS-less data)")

## ------------------- WP-3 (P4): dose-response joins the store ----
# The first snapshot-additive adoption: the two dataTableDownload tables
# register selected-row bindings under the NEW resultspace.dose_response
# namespace (D2/D18) - exactly the two keys above, nothing else. The
# S6.6_drc attribute derivation moved into the module from the L1 call
# site (verbatim expressions over context$data$object); these cases pin
# the moved glue end-to-end over the demo data (32 ResponseCurve
# columns; S6.6_drc = c(dose_col = "Dose", curveid_col = "Curve")).
dose_response_module <- omicsViewer:::dose_response_module
agent_capability_records <- omicsViewer:::agent_capability_records
.w3_params <- c("b_hill.slope", "c_min.response", "d_max.response",
                "e_inflection", "pval", "log.pval", "pseudo.rsq", "range")
# the sanctioned WP-3 additions (now folded into the frozen vector above;
# the D12 update commit's own subject)
PANEL_W3_DOSE_IDS <- sort(c(
  "resultspace.dose_response.feature_table",
  "resultspace.dose_response.param_table"))
# shared boot factory: the module wired exactly like the L1 call site
# does since WP-3 (context-built store + object in the context data)
.w3_mkapp <- function(child, object = dat) {
  force(child); force(object)
  function(input, output, session) {
    dose_response_module(
      "rescurve",
      reactive_expr = reactive(Biobase::exprs(dat)),
      reactive_phenoData = reactive(Biobase::pData(dat)),
      reactive_featureData = reactive(fd0),
      reactive_i = reactive(1L),
      context = panel_context(
        id = "rescurve", tab_value = "Response", store = child,
        data = list(object = reactive(object))))
  }
}

# --- context adoption: bindings, live choices, agent round-trip
w3 <- new.env()
w3$store <- widget_store_new()
w3$child <- widget_store_child(w3$store, "resultspace.dose_response")
w3_err <- NULL; w3_view <- NULL; w3_receipt <- NULL; w3_applied <- NULL
w3_bad <- NULL
shiny::testServer(.w3_mkapp(w3$child), {
  w3_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  # registry view (inside the session: the choices providers read live
  # reactives - outside a reactive context they degrade to unavailable)
  w3_view <<- store_registry_view(w3$store, "resultspace.dose_response")
  # agent round-trip exactly like the other panels' tables
  w3_receipt <<- tryCatch(
    store_apply(w3$child, list(param_table = "pval"), origin = "agent"),
    error = function(e) conditionMessage(e))
  w3_applied <<- store_read(w3$child, "param_table")[[1]]
  w3_bad <<- tryCatch(
    store_apply(w3$child, list(feature_table = "GHOST_FEATURE"),
                origin = "agent"),
    error = function(e) conditionMessage(e))
  session$flushReact()
})
ok(is.null(w3_err),
   "WP-3 dose: module boots headless with a context-built store child")
ok(ut_cmp_identical(sort(vapply(w3_view, function(r) r$id, "")),
                    PANEL_W3_DOSE_IDS),
   "WP-3 dose: exactly the two sanctioned selected-row bindings register")
ok(all(vapply(w3_view, function(r) identical(r$kind, "select"), TRUE)),
   "WP-3 dose: both bindings are select-kind (string selected_row pattern)")
w3_par <- Filter(function(r) grepl("param_table$", r$id), w3_view)[[1]]
w3_feat <- Filter(function(r) grepl("feature_table$", r$id), w3_view)[[1]]
ok(ut_cmp_identical(w3_par$allowed_values, .w3_params),
   "WP-3 dose: parameter choices live via the moved-in drc derivation (exact demo rows pinned)")
ok(ut_cmp_identical(w3_feat$allowed_values, rownames(fd0)[1]),
   "WP-3 dose: feature-table choice is the selected feature id (demo feature 1)")
ok(!is.character(w3_receipt) &&
   ut_cmp_identical(w3_receipt$applied, "resultspace.dose_response.param_table") &&
   ut_cmp_identical(w3_receipt$diff$`resultspace.dose_response.param_table`, "pval"),
   "WP-3 dose: agent store_apply round-trips with a receipt (param_table)")
ok(ut_cmp_identical(w3_applied, "pval"),
   "WP-3 dose: the applied value is readable back through the store")
ok(is.character(w3_bad) && grepl("Unknown value for resultspace.dose_response.feature_table", w3_bad),
   "WP-3 dose: a value outside the live choices is rejected (validation over the moved-in derivation)")

# --- user clicks sync into the store (the tables' acknowledgement path)
w3u <- new.env()
w3u$store <- widget_store_new()
w3u$child <- widget_store_child(w3u$store, "resultspace.dose_response")
w3u_err <- NULL; w3u_click <- NULL
shiny::testServer(.w3_mkapp(w3u$child), {
  w3u_err <<- tryCatch({
    # warm the sync observers before the click (AGENTS.md testServer
    # lore; bare observers fire on the first flush regardless)
    session$flushReact()
    session$setInputs(`rescurve-param-table_rows_selected` = 2L)
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w3u_click <<- store_read(w3u$child, "param_table")[[1]]
})
ok(is.null(w3u_err),
   "WP-3 dose: table click path stays alive headless")
ok(ut_cmp_identical(w3u_click, "c_min.response"),
   "WP-3 dose: a user row click syncs the semantic row id into the store")

# --- dataset-switch reset: L0's store_reset prefix "resultspace" now
# actually covers dose (P4's invisible-before reset gap). Runs inside a
# live session (the reset bumps the store's reactive values)
w3r <- new.env()
w3r$store <- widget_store_new()
w3r$child <- widget_store_child(w3r$store, "resultspace.dose_response")
w3r_err <- NULL; w3r_reset <- NULL; w3r_after <- NULL
shiny::testServer(.w3_mkapp(w3r$child), {
  w3r_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    store_apply(w3r$child,
                list(param_table = "pval",
                     feature_table = rownames(fd0)[1]), origin = "agent")
    session$flushReact()
    # L0's exact dataset-switch call (todo 2.5)
    w3r_reset <<- store_reset(w3r$store, c("dataspace", "resultspace"))
    session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w3r_after <<- store_read(w3r$child, c("param_table", "feature_table"))
})
ok(is.null(w3r_err), "WP-3 dose: dataset-switch reset path stays alive")
ok(ut_cmp_identical(sort(intersect(w3r_reset$reset, PANEL_W3_DOSE_IDS)),
                    PANEL_W3_DOSE_IDS),
   "WP-3 dose: dataset-switch reset clears the dose keys (resultspace prefix coverage)")
ok(ut_cmp_identical(w3r_after,
                    list(`resultspace.dose_response.param_table` = NULL,
                         `resultspace.dose_response.feature_table` = NULL)),
   "WP-3 dose: dose keys read NULL after the reset")

# --- legacy paths: store-less boots idle; the explicit reactive_attr_drc
# formal still feeds the tables (direct instantiation, D13)
w3l <- new.env()
w3l$store <- widget_store_new()
w3l$child <- widget_store_child(w3l$store, "resultspace.dose_response")
w3l_err <- NULL; w3l_view <- NULL
shiny::testServer(function(input, output, session) {
  dose_response_module(
    "rescurve",
    reactive_expr = reactive(Biobase::exprs(dat)),
    reactive_phenoData = reactive(Biobase::pData(dat)),
    reactive_featureData = reactive(fd0),
    reactive_i = reactive(1L),
    reactive_attr_drc = reactive(c(dose_col = "Dose", curveid_col = "Curve")),
    store = w3l$child)
}, {
  w3l_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w3l_view <<- store_registry_view(w3l$store, "resultspace.dose_response")
})
ok(is.null(w3l_err),
   "WP-3 dose: legacy direct instantiation (explicit reactive_attr_drc + store) boots")
ok(ut_cmp_identical(
  Filter(function(r) grepl("param_table$", r$id), w3l_view)[[1]]$allowed_values,
  .w3_params),
   "WP-3 dose: the legacy reactive_attr_drc formal still drives the choices (no context)")

w3n <- new.env()
w3n_err <- NULL
shiny::testServer(function(input, output, session) {
  dose_response_module(
    "rescurve",
    reactive_expr = reactive(Biobase::exprs(dat)),
    reactive_phenoData = reactive(Biobase::pData(dat)),
    reactive_featureData = reactive(fd0),
    reactive_i = reactive(1L))
}, {
  w3n_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
})
ok(is.null(w3n_err),
   "WP-3 dose: store-less legacy use boots exactly as before (no store, no context)")

# --- the derivation sources context$data$object: without the S6.6_drc
# attribute the choices stay cold and agent writes are rejected with the
# WP-A/F3 unavailable reason (same class as string pre-Run)
dat_nodrc <- dat
attr(dat_nodrc, "S6.6_drc") <- NULL
w3c <- new.env()
w3c$store <- widget_store_new()
w3c$child <- widget_store_child(w3c$store, "resultspace.dose_response")
w3c_err <- NULL; w3c_bad <- NULL; w3c_allowed <- NULL
shiny::testServer(.w3_mkapp(w3c$child, dat_nodrc), {
  w3c_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w3c_allowed <<- Filter(function(r) grepl("param_table$", r$id),
    store_registry_view(w3c$store, "resultspace.dose_response"))[[1]]$allowed_values
  w3c_bad <<- tryCatch(
    store_apply(w3c$child, list(param_table = "pval"), origin = "agent"),
    error = function(e) conditionMessage(e))
})
ok(is.null(w3c_err),
   "WP-3 dose: an object without S6.6_drc boots idle (no error)")
ok(identical(w3c_allowed, character(0)),
   "WP-3 dose: choices stay cold without the drc attribute (derivation reads context$data$object)")
ok(is.character(w3c_bad) && grepl("Choices are currently unavailable", w3c_bad),
   "WP-3 dose: cold choices reject agent writes (WP-A/F3 unavailable semantics)")

# --- capability records: additive, mapped to the Response tab (D10)
w3_caps <- Filter(function(r) grepl("^resultspace\\.dose_response\\.", r$id),
                  agent_capability_records(esv$store))
ok(ut_cmp_identical(sort(vapply(w3_caps, function(r) r$id, "")),
                    PANEL_W3_DOSE_IDS),
   "WP-3 dose: exactly two capability records appear (additive; the enumeration freeze bounds the rest)")
ok(ut_cmp_identical(unique(vapply(w3_caps, function(r) r$tab, "")), "Response"),
   "WP-3 dose: capability records map to the Response tab (the one AGENT_TAB_LABELS entry)")

# --- rollback direction: a snapshot carrying dose keys restores into a
# store WITHOUT the dose bindings through unknown-id degradation
# (auxi_widgetStore store_restore; pre-WP-3 code reading post-WP-3
# snapshots). The app-level both-direction cases live in
# test_snapshotRoundTrip.R.
w3_old <- widget_store_new()
store_register(widget_store_child(w3_old, "resultspace.ora"),
               widget_binding("xax_analysis", "select", values = c("ttest", "PCA")))
w3_rr <- tryCatch(
  store_restore(w3_old, list(
    values = list(`resultspace.ora.xax_analysis` = "ttest",
                  `resultspace.dose_response.param_table` = "pval",
                  `resultspace.dose_response.feature_table` = "X1"),
    unset = character(0)), replace = TRUE),
  error = function(e) conditionMessage(e))
ok(!is.character(w3_rr) &&
   ut_cmp_identical(sort(w3_rr$unknown_ids), PANEL_W3_DOSE_IDS) &&
   ut_cmp_identical(w3_rr$applied, "resultspace.ora.xax_analysis") &&
   ut_cmp_identical(store_read(w3_old, "resultspace.ora.xax_analysis")[[1]], "ttest"),
   "WP-3 dose: dose keys in a snapshot degrade to unknown_ids on a dose-less store (rollback-safe restore)")

## ------------------- WP-2: additionalTabs entry validation (P8) ----
# panel_validate_additional_tabs checks EXACTLY the 5 fields the two L1
# loops and the L0 availability read (tabName/moduleName non-empty scalar
# chr, moduleUi/moduleServer functions, optional scalar-logical stateful)
# and nothing stricter: extra fields, unnamed top-level lists and
# NA/absent/NULL stateful pass because the loops tolerate them (isTRUE()
# semantics). Malformed input errors with the entry index + field +
# expected type. The guard runs at BOTH entries: omicsViewer() (first
# statement, before runApp) and app_module() (before moduleServer).
panel_validate_additional_tabs <- omicsViewer:::panel_validate_additional_tabs

.at_good <- list(
  tabName = "Custom", moduleName = "customMod",
  moduleUi = function(id) div(paste0("custom-", id)),
  moduleServer = function(id, pdata, fdata, expr, feature_selected,
                          sample_selected, object, ...) paste("srv", id))
.at <- function(...) utils::modifyList(.at_good, list(...))   # merge by name
.at_wo <- function(fld) { r <- .at_good; r[[fld]] <- NULL; r }  # field removed
.at_st_null <- { r <- .at_good; r["stateful"] <- list(NULL); r }  # explicit NULL element

ok(is.null(panel_validate_additional_tabs(NULL)),
   "additionalTabs: NULL passes (the fast path)")
ok(is.null(panel_validate_additional_tabs(list())),
   "additionalTabs: empty list passes (length(x) > 0 guards in the loops)")
ok(is.null(panel_validate_additional_tabs(list(.at_good))),
   "additionalTabs: well-formed record passes")
ok(is.null(panel_validate_additional_tabs(list(.at(stateful = TRUE)))),
   "additionalTabs: optional stateful = TRUE passes")
ok(is.null(panel_validate_additional_tabs(list(.at(stateful = NA)))),
   "additionalTabs: NA stateful passes (isTRUE(NA) is FALSE in the loops - no stricter)")
ok(is.null(panel_validate_additional_tabs(list(.at_st_null))),
   "additionalTabs: explicit NULL stateful passes (indistinguishable from absent via $)")
ok(is.null(panel_validate_additional_tabs(list(.at(extra_field = 1L)))),
   "additionalTabs: extra fields pass (the loops never read them)")
ok(is.null(panel_validate_additional_tabs(list(.at_good, .at(tabName = "Second")))),
   "additionalTabs: multiple well-formed records pass")

## negative paths: every bad shape errors with index + field + expected type
ok(grepl("additionalTabs must be NULL or a list of tab records; got a numeric of length 1",
         .err(panel_validate_additional_tabs(42))),
   "additionalTabs: non-list top level rejected with an actionable message")
ok(grepl("additionalTabs\\[\\[2\\]\\] must be a list \\(a tab record\\); got a character of length 1",
         .err(panel_validate_additional_tabs(list(.at_good, "nope")))),
   "additionalTabs: non-list entry rejected, error names the right index")

for (fld in c("tabName", "moduleName", "moduleUi", "moduleServer")) {
  ok(grepl(paste0("additionalTabs\\[\\[2\\]\\] is missing required field '", fld, "'"),
           .err(panel_validate_additional_tabs(list(.at_good, .at_wo(fld))))),
     paste("additionalTabs: missing", fld, "-> error names index + field"))
}

ok(grepl("additionalTabs\\[\\[2\\]\\]\\$tabName must be a non-empty character string of length 1; got a numeric of length 1",
         .err(panel_validate_additional_tabs(list(.at_good, .at(tabName = 42))))),
   "additionalTabs: numeric tabName rejected, message states the expected type")
ok(grepl("tabName must be a non-empty character string of length 1; got a character of length 2",
         .err(panel_validate_additional_tabs(list(.at_good, .at(tabName = c("a", "b")))))),
   "additionalTabs: length-2 tabName rejected (L0 vapply is character(1))")
ok(grepl("moduleName must be a non-empty character string of length 1; got a character of length 1",
         .err(panel_validate_additional_tabs(list(.at_good, .at(moduleName = ""))))),
   "additionalTabs: empty-string moduleName rejected")
ok(grepl("moduleName must be a non-empty character string of length 1; got a logical of length 1",
         .err(panel_validate_additional_tabs(list(.at_good, .at(moduleName = NA))))),
   "additionalTabs: NA moduleName rejected")
ok(grepl("moduleUi must be a function; got a character of length 2",
         .err(panel_validate_additional_tabs(list(.at_good, .at(moduleUi = c("a", "b")))))),
   "additionalTabs: non-function moduleUi rejected")
ok(grepl("moduleServer must be a function; got a list of length 1",
         .err(panel_validate_additional_tabs(list(.at_good, .at(moduleServer = list(x = 1)))))),
   "additionalTabs: non-function moduleServer rejected")
ok(grepl("stateful must be TRUE or FALSE \\(a logical of length 1\\), or absent/NULL; got a character of length 1",
         .err(panel_validate_additional_tabs(list(.at_good, .at(stateful = "yes"))))),
   "additionalTabs: character stateful rejected")
ok(grepl("stateful must be TRUE or FALSE \\(a logical of length 1\\), or absent/NULL; got a logical of length 2",
         .err(panel_validate_additional_tabs(list(.at_good, .at(stateful = c(TRUE, FALSE)))))),
   "additionalTabs: length-2 stateful rejected")
ok(grepl("additionalTabs\\[\\[3\\]\\] is missing required field 'moduleServer'",
         .err(panel_validate_additional_tabs(list(
           .at_good, .at(tabName = "Second"), .at_wo("moduleServer"))))),
   "additionalTabs: valid entries before an invalid one - the error names the RIGHT index")

## entry wiring: malformed input fails FAST at both entries (before any
## Shiny machinery runs - neither runApp nor moduleServer is reached)
ok(grepl("additionalTabs\\[\\[1\\]\\] is missing required field 'moduleServer'",
         .err(omicsViewer::omicsViewer(dir = tempdir(),
                                       additionalTabs = list(.at_wo("moduleServer"))))),
   "omicsViewer(): malformed additionalTabs fails fast at entry (validator runs before runApp)")
ok(grepl("additionalTabs\\[\\[1\\]\\]\\$moduleUi must be a function",
         .err(omicsViewer:::app_module("app", .dir = reactive(NULL),
                                       additionalTabs = list(.at(moduleUi = "nope"))))),
   "app_module(): entry guard fires before moduleServer (developer embedding path)")

## no behavior change for well-formed input: the real L1 composition
## boots with an additional tab exactly as before (freeze-harness pattern)
esv_at <- new.env()
esv_at$store <- widget_store_new()
esv_at$status <- NULL
app_l1_at <- function(input, output, session) {
  esv_at$status <- omicsViewer:::L1_result_space_module(
    "rs",
    reactive_expr = reactive(Biobase::exprs(dat)),
    reactive_phenoData = reactive(Biobase::pData(dat)),
    reactive_featureData = reactive(fd0),
    additionalTabs = list(.at_good),
    status = reactive(NULL),
    store = esv_at$store)
}
l1_at_err <- NULL; l1_at_status <- NULL
shiny::testServer(app_l1_at, {
  l1_at_err <<- tryCatch({
    session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  l1_at_status <<- tryCatch(shiny::isolate(esv_at$status()),
                            error = function(e) conditionMessage(e))
})
ok(is.null(l1_at_err),
   "boot: well-formed additionalTabs still boots the real L1 unchanged")
ok(ut_cmp_identical(l1_at_status$analyst_additional,
                    list(customMod = list(state_support = FALSE))),
   "boot: the additional module composes into the status payload exactly as before")
