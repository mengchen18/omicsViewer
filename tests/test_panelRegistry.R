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
# WP-11 record-order data fix: the LIST order is the L1 SERVER
# composition order (module invocation = binding registration = snapshot
# byte order; the loop in L1_module_result_space.R iterates it), while
# the TAB composition order lives in the order field (rs_tab_set sorts
# by it - its output is unchanged; the truth table below pins that).
ok(identical(names(panel_records),
             c("feature_general", "fgsea", "ora", "stringdb",
               "sample_general", "geneshotTab", "ptm", "rescurve")),
   "registry: 8 records keyed by module id in L1 server-composition order (WP-11 data fix; the pre-WP-11 literal call order)")
ok(identical(unname(vapply(panel_records, function(r) r$tab_value, "")), c(
  "Feature", "fGSEA", "ORA", "StringDB", "Sample", "Geneshot",
  "SeqLogo", "Response")),
   "registry: tab values verbatim, list order = server-composition order")
ok(identical(sort(unname(vapply(panel_records, function(r) r$tab_value, ""))), sort(c(
  "Feature", "ORA", "fGSEA", "Response", "StringDB", "SeqLogo",
  "Geneshot", "Sample"))),
   "registry: the tab-value SET is the plan's 8 verbatim values")
ok(identical(unname(vapply(panel_records, function(r) r$store_id, "")), c(
  "resultspace.feature_general", "resultspace.fgsea", "resultspace.ora",
  "resultspace.stringdb", "resultspace.sample_general",
  "resultspace.geneshot", "resultspace.ptm",
  "resultspace.dose_response")),
   "registry: store prefixes verbatim incl. the NEW resultspace.dose_response (D2/D18)")
ok(identical(unname(vapply(panel_records, function(r) r$status_key, "")), c(
  "analyst_feature_general", "analyst_fgsea", "analyst_ora",
  "analyst_stringdb", "analyst_sample_general", "analyst_gene_shot",
  "analyst_ptm", "analyst_dose_response")),
   "registry: status keys verbatim (analyst_gene_shot != module id, historical)")
ok(identical(sort(unname(vapply(panel_records, function(r) r$order, 1L))), 1:8),
   "registry: order field is a permutation of 1..8 (one rank per record; the TAB composition rank)")
ok(identical(
  unname(vapply(panel_records[order(vapply(panel_records, function(r) r$order, 1L))],
                function(r) r$tab_value, "")),
  c("Feature", "ORA", "fGSEA", "Response", "StringDB", "SeqLogo",
    "Geneshot", "Sample")),
   "registry: sorting by the order field yields the tab-composition order (what rs_tab_set emits; unchanged by the WP-11 data fix)")
ok(identical(
  unname(vapply(panel_records, function(r) isTRUE(r$stateful), TRUE)),
  c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)),
   "registry: stateful flags - stringdb/dose FALSE (WP-4/WP-3 semantics)")
ok(identical(
  unname(vapply(panel_records, function(r) r$icon, "")),
  c("chart-column", "arrow-down-short-wide", "object-group",
    "share-nodes", "arrows-down-to-people", "book-open", "font",
    "capsules")),
   "registry: icons match today's L1 navbar composition verbatim (one per record)")
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
w3c_err <- NULL; w3c_bad <- NULL; w3c_allowed <- NULL; w3c_rr <- NULL
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
  # NIT-1 (WP-3 review): the no-drc RESTORE direction under the same cold
  # conditions as the apply above — a dose-carrying snapshot degrades
  # per-key through store_restore (a carrier key proves the valid
  # remainder still applies while the cold dose keys are rejected)
  store_register(widget_store_child(w3c$store, "resultspace.ora"),
                 widget_binding("xax_analysis", "select",
                                values = c("ttest", "PCA")))
  w3c_rr <<- tryCatch(
    store_restore(w3c$store, list(
      values = list(`resultspace.ora.xax_analysis` = "ttest",
                    `resultspace.dose_response.param_table` = "pval",
                    `resultspace.dose_response.feature_table` = "X1"),
      unset = character(0)), replace = TRUE),
    error = function(e) conditionMessage(e))
})
ok(is.null(w3c_err),
   "WP-3 dose: an object without S6.6_drc boots idle (no error)")
ok(identical(w3c_allowed, character(0)),
   "WP-3 dose: choices stay cold without the drc attribute (derivation reads context$data$object)")
ok(is.character(w3c_bad) && grepl("Choices are currently unavailable", w3c_bad),
   "WP-3 dose: cold choices reject agent writes (WP-A/F3 unavailable semantics)")
ok(!is.character(w3c_rr) && length(w3c_rr$rejected) == 2L &&
   ut_cmp_identical(sort(vapply(w3c_rr$rejected, `[[`, "", "id")),
                    PANEL_W3_DOSE_IDS) &&
   all(vapply(w3c_rr$rejected,
              function(r) grepl("Choices are currently unavailable", r$reason),
              TRUE)) &&
   ut_cmp_identical(w3c_rr$applied, "resultspace.ora.xax_analysis") &&
   all(vapply(store_read(w3c$store, PANEL_W3_DOSE_IDS), is.null, TRUE)),
   "WP-3 dose: a dose-carrying snapshot restored into a no-drc dataset degrades per-key (dose keys rejected unavailable, valid remainder applies, dose values stay NULL)")

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

## ------------------- WP-4 (P5): STRING snapshot story made honest ----
# The ids derivation moved INTO the module (verbatim from the L1 call
# site, over context$data$fdata + context$selection$feature) and the
# active flag comes from context$active; the dead analyst_stringdb
# status emission is GONE from L1 (app-level round-trip proof lives in
# test_snapshotRoundTrip.R case 13). These cases pin the moved glue at
# module level over the demo data (StringDB|All|ID column; the run
# action's availability is the live consumer of the derived ids).
string_module <- omicsViewer:::string_module
PANEL_W4_STRING_IDS <- sort(c(
  "resultspace.stringdb.taxonomy",
  "resultspace.stringdb.show_labels",
  "resultspace.stringdb.selected_row"))
# shared boot factory: the module wired exactly like the L1 call site
# does since WP-4 (context-built store + the derivation sources in the
# context); the feature selection is a reactiveVal so the case can flip
# it and watch the derived ids go live
.w4_sel <- shiny::reactiveVal(NULL)
.w4_mkapp <- function(child) {
  force(child)
  function(input, output, session) {
    string_module(
      "stringdb",
      context = panel_context(
        id = "stringdb", tab_value = "StringDB", store = child,
        data = list(fdata = reactive(fd0)),
        selection = list(feature = .w4_sel),
        active_tab = reactive("StringDB")))
  }
}

# --- context adoption: bindings, and the moved-in ids derivation LIVE
# through the run action's availability (NULL selection -> no ids;
# feature 1 selected -> the demo's STRING id column yields the ids)
w4 <- new.env()
w4$store <- widget_store_new()
w4$child <- widget_store_child(w4$store, "resultspace.stringdb")
w4_err <- NULL; w4_view <- NULL; w4_act0 <- NULL; w4_act1 <- NULL
shiny::testServer(.w4_mkapp(w4$child), {
  w4_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w4_view <<- store_registry_view(w4$store, "resultspace.stringdb")
  w4_act0 <<- omicsViewer:::store_action_view(
    w4$store, "resultspace.stringdb")[[1]]$available
  .w4_sel(1L)                       # a user picks feature row 1
  session$flushReact(); session$flushReact()
  w4_act1 <<- omicsViewer:::store_action_view(
    w4$store, "resultspace.stringdb")[[1]]$available
})
ok(is.null(w4_err),
   "WP-4 string: module boots headless with a context-built store child")
ok(ut_cmp_identical(sort(vapply(w4_view, function(r) r$id, "")),
                    PANEL_W4_STRING_IDS),
   "WP-4 string: exactly the three frozen bindings register (D12 untouched)")
ok(ut_cmp_identical(
  w4_act0,
  "select features first (set_omics_viewer_state features)"),
   "WP-4 string: without a selection the derived ids are empty (the moved-in derivation reads context$selection$feature)")
ok(ut_cmp_identical(w4_act1, TRUE),
   "WP-4 string: a feature selection makes the derived ids live (the L1-era expression, moved in verbatim)")

# --- legacy paths: the reactive_ids/active/store formals still work
# (direct instantiation, D13); no store and no context boots idle
w4l <- new.env()
w4l$store <- widget_store_new()
w4l$child <- widget_store_child(w4l$store, "resultspace.stringdb")
w4l_err <- NULL; w4l_view <- NULL
shiny::testServer(function(input, output, session) {
  string_module("sdb", reactive_ids = reactive(c("g1", "g2", "g3")),
                store = w4l$child)
}, {
  w4l_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w4l_view <<- store_registry_view(w4l$store, "resultspace.stringdb")
})
ok(is.null(w4l_err),
   "WP-4 string: legacy direct instantiation (explicit reactive_ids + store) boots")
ok(ut_cmp_identical(sort(vapply(w4l_view, function(r) r$id, "")),
                    PANEL_W4_STRING_IDS),
   "WP-4 string: the legacy formals register the same three bindings (signature stability)")

w4n_err <- NULL
shiny::testServer(function(input, output, session) {
  string_module("sdb", reactive_ids = reactive(c("g1", "g2")))
}, {
  w4n_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
})
ok(is.null(w4n_err),
   "WP-4 string: store-less legacy use boots exactly as before (no store, no context)")

# --- guards: neither ids nor context is a hard error; a context without
# the derivation sources is (actionable, per the panel_context school)
ok(grepl("reactive_ids or a context",
         tryCatch(string_module("sdb"), error = function(e) conditionMessage(e))),
   "WP-4 string: neither reactive_ids nor context rejected with an actionable message")
err_w4_ctx <- tryCatch(
  string_module("sdb",
                context = panel_context(id = "stringdb",
                                        data = list(fdata = reactive(fd0)))),
  error = function(e) conditionMessage(e))
ok(grepl("selection[$]feature", err_w4_ctx),
   "WP-4 string: a context missing the derivation sources errors naming them")


## ------------- WP-5: PTM migration to the shared socket (D3/D4/D13) ----
# The FIRST pure migration WP: the panel was already a well-formed store
# citizen, so this is a conversion to the shared machinery (panel_keep
# retention D3, panel_triset_bind cascade D4, context threading D13)
# with ZERO store-key/behavior change. The goldens below were captured
# from the PRE-MIGRATION HEAD (4ff9f8b, 2026-11-06) BEFORE any edit, by
# booting the then-current ptmotif_module exactly like the L1 call site
# (child store resultspace.ptm over demo.RDS fdata, feature 1 selected)
# and dput()ing the outputs - the WP-0 ORA-golden-fixture pattern. They
# are FROZEN: never regenerate them from the migrated code (that would
# make the equivalence claim circular).
ptmotif_module <- omicsViewer:::ptmotif_module
PANEL_W5_PTM_IDS <- sort(c(
  "resultspace.ptm.xax_analysis",
  "resultspace.ptm.xax_subset",
  "resultspace.ptm.xax_variable"))
# frozen pre-migration golden: the registry view over demo fd0 (exactly
# one SeqLogo|All|seq.window column) - id/kind/label/help/depends_on/
# min/max/allowed_values of all three bindings
PANEL_W5_PTM_GOLDEN_VIEW <- list(
  list(id = "resultspace.ptm.xax_analysis", kind = "select",
       label = "Sequence category",
       help = "Annotation category of the sequence-window column",
       depends_on = character(0), min = NULL, max = NULL,
       allowed_values = "SeqLogo"),
  list(id = "resultspace.ptm.xax_subset", kind = "select",
       label = "Sequence subcategory",
       help = "Subcategory within the sequence category",
       depends_on = "resultspace.ptm.xax_analysis", min = NULL, max = NULL,
       allowed_values = "All"),
  list(id = "resultspace.ptm.xax_variable", kind = "select_cascaded",
       label = "Sequence column",
       help = paste("Feature annotation column holding the sequence windows",
                    "around modification sites"),
       depends_on = c("resultspace.ptm.xax_analysis",
                     "resultspace.ptm.xax_subset"),
       min = NULL, max = NULL, allowed_values = "seq.window"))
# frozen pre-migration golden: the auto-seeded store values over demo fd0
PANEL_W5_PTM_GOLDEN_SEED <- list(
  resultspace.ptm.xax_analysis = "SeqLogo",
  resultspace.ptm.xax_subset = "All",
  resultspace.ptm.xax_variable = "seq.window")
# frozen pre-migration golden: the raw binding metadata incl. the
# agent_writable flags (the registry view only shows agent-writable ids;
# this proves all three are in it AND writable)
PANEL_W5_PTM_GOLDEN_META <- list(
  list(id = "resultspace.ptm.xax_analysis", kind = "select",
       label = "Sequence category",
       help = "Annotation category of the sequence-window column",
       depends_on = character(0), agent_writable = TRUE),
  list(id = "resultspace.ptm.xax_subset", kind = "select",
       label = "Sequence subcategory",
       help = "Subcategory within the sequence category",
       depends_on = "resultspace.ptm.xax_analysis", agent_writable = TRUE),
  list(id = "resultspace.ptm.xax_variable", kind = "select_cascaded",
       label = "Sequence column",
       help = paste("Feature annotation column holding the sequence windows",
                    "around modification sites"),
       depends_on = c("resultspace.ptm.xax_analysis",
                     "resultspace.ptm.xax_subset"),
       agent_writable = TRUE))
# frozen pre-migration golden: the real-L1 status payload key order and
# the analyst_ptm entry (module return list(xax = v1()); v1 is NULL
# until the triselector commits through the UI, hence xax = NULL in a
# cold headless boot)
PANEL_W5_PTM_GOLDEN_L1_STATUS <- list(
  names = c("analyst_active_tab", "analyst_feature_general",
            "analyst_sample_general", "analyst_gene_shot", "analyst_fgsea",
            "analyst_ora", "analyst_ptm", "analyst_dose_response",
            "analyst_additional"),
  analyst_ptm = list(xax = NULL))
# synthetic 3-column SeqLogo fdata (the test_agentWidgets shape, plus a
# Kinase subcategory so the cascade has a real branch to pick into)
.w5_fd3 <- data.frame(
  `General|All|Gene.name` = paste0("g", 1:8),
  `SeqLogo|All|win15` = replicate(8, paste(sample(c("A","C","D","E","F","G",
    "H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"), 15, TRUE),
    collapse = "")),
  `SeqLogo|All|win20` = replicate(8, paste(sample(c("A","C","D","E","F","G",
    "H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"), 20, TRUE),
    collapse = "")),
  `SeqLogo|Kinase|win15k` = replicate(8, paste(sample(c("A","C","D","E","F","G",
    "H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"), 15, TRUE),
    collapse = "")),
  row.names = paste0("P", 1:8), check.names = FALSE)
# shared boot factories: context-built (the WP-5 L1 call-site shape;
# holder receives the context so the keep closure is inspectable) and
# legacy direct (formals only - the pre-WP-5 shape, still supported)
.w5_mkctx <- function(child, holder, fd = fd0, feat = 1L, status_payload = NULL) {
  force(child); force(holder)
  function(input, output, session) {
    .fd <- reactive(fd); .ft <- reactive(feat); .st <- reactive(status_payload)
    holder$ctx <- panel_context(
      id = "ptm", tab_value = "SeqLogo", store = child,
      data = list(fdata = .fd), selection = list(feature = .ft),
      status = .st, session = session)
    holder$ret <- ptmotif_module(
      "ptm", fdata = .fd, feature_selected = .ft, reactive_status = .st,
      store = child, context = holder$ctx)
  }
}
.w5_mklegacy <- function(child, fd = fd0, feat = 1L, status_payload = NULL) {
  force(child)
  function(input, output, session) {
    ptmotif_module(
      "ptm", fdata = reactive(fd), feature_selected = reactive(feat),
      reactive_status = reactive(status_payload), store = child)
  }
}
.w5_meta <- function(store) {
  lapply(grep("^resultspace\\.ptm\\.", names(store$bindings), value = TRUE),
         function(k) {
           b <- store$bindings[[k]]
           list(id = k, kind = b$kind, label = b$label, help = b$help,
                depends_on = b$depends_on, agent_writable = b$agent_writable)
         })
}

# --- context boot equals the frozen pre-migration golden, byte-for-byte
w5 <- new.env()
w5$store <- widget_store_new()
w5$child <- widget_store_child(w5$store, "resultspace.ptm")
w5_err <- NULL; w5_view <- NULL; w5_meta <- NULL; w5_seed <- NULL; w5_ret <- NULL
shiny::testServer(.w5_mkctx(w5$child, w5), {
  w5_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  # registry view INSIDE the session (WP-3 lore: choices providers read
  # live reactives; outside a reactive context they degrade to empty)
  w5_view <<- store_registry_view(w5$store, "resultspace.ptm")
  w5_meta <<- .w5_meta(w5$store)
  w5_seed <<- store_read(w5$child)
  w5_ret <<- tryCatch(shiny::isolate(w5$ret()), error = function(e) NULL)
})
ok(is.null(w5_err),
   "WP-5 ptm: module boots headless context-built (the new L1 call-site shape)")
ok(ut_cmp_identical(sort(vapply(w5_view, function(r) r$id, "")), PANEL_W5_PTM_IDS),
   "WP-5 ptm: exactly the three frozen bindings register under resultspace.ptm (D12 untouched)")
ok(ut_cmp_identical(w5_view, PANEL_W5_PTM_GOLDEN_VIEW),
   "WP-5 ptm: registry view byte-identical to the frozen pre-migration golden (id/kind/label/help/depends_on/allowed_values over demo fd)")
ok(ut_cmp_identical(w5_meta, PANEL_W5_PTM_GOLDEN_META),
   "WP-5 ptm: raw binding meta byte-identical incl. agent_writable TRUE x3")
ok(ut_cmp_identical(w5_seed, PANEL_W5_PTM_GOLDEN_SEED),
   "WP-5 ptm: first-choice auto-default seeds the same triple (SeqLogo|All|seq.window)")
ok(ut_cmp_identical(w5_ret, list(xax = NULL)),
   "WP-5 ptm: module return (the analyst_ptm payload source) stays list(xax = NULL) before any UI commit")

# --- legacy direct boot registers identically to the context boot
w5l <- new.env()
w5l$store <- widget_store_new()
w5l$child <- widget_store_child(w5l$store, "resultspace.ptm")
w5l_err <- NULL; w5l_view <- NULL; w5l_meta <- NULL
shiny::testServer(.w5_mklegacy(w5l$child), {
  w5l_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w5l_view <<- store_registry_view(w5l$store, "resultspace.ptm")
  w5l_meta <<- .w5_meta(w5l$store)
})
ok(is.null(w5l_err),
   "WP-5 ptm: legacy direct instantiation (formals only, no context) still boots")
ok(ut_cmp_identical(w5l_view, PANEL_W5_PTM_GOLDEN_VIEW) &&
   ut_cmp_identical(w5l_meta, PANEL_W5_PTM_GOLDEN_META),
   "WP-5 ptm: legacy boot registers byte-identically to the golden (both boot paths agree)")
ok(ut_cmp_identical(w5l_view, w5_view),
   "WP-5 ptm: context boot and legacy boot register identically (D13 argument coalescing)")

# --- behavior pins over the synthetic 3-column fd (pre-migration
# captures (b2)/(c)): auto-seed, restore-first-wins, user-pick-sticks,
# cascade-dependent live choices - all through the CONTEXT boot path
w5b <- new.env()
w5b$store <- widget_store_new()
w5b$child <- widget_store_child(w5b$store, "resultspace.ptm")
w5b_err <- NULL; w5b_seed <- NULL; w5b_pick <- NULL; w5b_ret <- NULL
w5b_ch <- NULL
shiny::testServer(.w5_mkctx(w5b$child, w5b, fd = .w5_fd3), {
  w5b_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w5b_seed <<- store_read(w5b$child)
  # user pick into the Kinase branch sticks over the auto-seed
  session$setInputs(`ptm-tris_seqlogo-analysis` = "SeqLogo")
  session$flushReact()
  session$setInputs(`ptm-tris_seqlogo-subset` = "Kinase")
  session$flushReact()
  session$setInputs(`ptm-tris_seqlogo-variable` = "win15k")
  session$flushReact(); session$flushReact()
  w5b_pick <<- store_read(w5b$child)
  w5b_ret <<- tryCatch(shiny::isolate(w5b$ret()), error = function(e) NULL)
  # state-dependent live choices (store now holds SeqLogo|Kinase|win15k)
  v <- store_registry_view(w5b$store, "resultspace.ptm")
  w5b_ch <<- stats::setNames(lapply(v, `[[`, "allowed_values"),
                            vapply(v, function(r) r$id, ""))
})
ok(is.null(w5b_err), "WP-5 ptm: synthetic-fd context boot stays alive")
ok(ut_cmp_identical(unname(unlist(w5b_seed)), c("SeqLogo", "All", "win15")),
   "WP-5 ptm: auto-default seeds the first cascade row (SeqLogo|All|win15)")
ok(ut_cmp_identical(unname(unlist(w5b_pick)), c("SeqLogo", "Kinase", "win15k")),
   "WP-5 ptm: user pick into the Kinase branch syncs and sticks over the auto-seed")
ok(ut_cmp_identical(unname(unlist(w5b_ret$xax)), c("SeqLogo", "Kinase", "win15k")),
   "WP-5 ptm: module return reports the settled triple after the commit (analyst_ptm payload content)")
ok(ut_cmp_identical(
  w5b_ch,
  list(`resultspace.ptm.xax_analysis` = "SeqLogo",
       `resultspace.ptm.xax_subset` = c("All", "Kinase"),
       `resultspace.ptm.xax_variable` = "win15k")),
   "WP-5 ptm: cascade choices stay state-dependent (subset branches, variable filtered by the Kinase branch)")

# --- legacy status restore path unchanged: an analyst_ptm payload with
# a DIFFERENT triple overrides the auto-seed through
# store_apply(origin = "restore") (pre-migration capture b2)
w5r <- new.env()
w5r$store <- widget_store_new()
w5r$child <- widget_store_child(w5r$store, "resultspace.ptm")
w5r_err <- NULL; w5r_read <- NULL
shiny::testServer(.w5_mkctx(
  w5r$child, w5r, fd = .w5_fd3,
  status_payload = list(xax = c("SeqLogo", "Kinase", "win15k"))), {
  w5r_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w5r_read <<- store_read(w5r$child)
})
ok(is.null(w5r_err), "WP-5 ptm: boot with a legacy analyst_ptm status payload stays alive")
ok(ut_cmp_identical(unname(unlist(w5r_read)), c("SeqLogo", "Kinase", "win15k")),
   "WP-5 ptm: the legacy restore path re-applies the payload triple over the auto-seed (restore-first-wins preserved)")

# --- retention (D3): the shared per-instance panel_keep holds BOTH the
# helper's UI->store sync observer and the module's seed observer
# (structural, deterministic - WP-0 MINOR-2 pattern); after gc() x2 the
# kept sync observer still fires (the WP-0 gc pattern, retention
# direction only per PD-8)
w5k <- new.env()
w5k$store <- widget_store_new()
w5k$child <- widget_store_child(w5k$store, "resultspace.ptm")
w5k_err <- NULL; w5k_nobs <- NULL; w5k_pick <- NULL
shiny::testServer(.w5_mkctx(w5k$child, w5k, fd = .w5_fd3), {
  w5k_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  # structural: everything the module kept lives in the CONTEXT's
  # panel_keep closure (D3 - the retention owner), not a module-local list
  obs <- get("obs", envir = environment(w5k$ctx$keep))
  w5k_nobs <<- length(obs)
  # gc integration (retention direction): collect twice, then a user
  # pick must still reach the store through the kept sync observer
  # (the context - the L1-frame analogue - is the observers' only anchor)
  gc(); gc()
  session$setInputs(`ptm-tris_seqlogo-analysis` = "SeqLogo")
  session$flushReact()
  session$setInputs(`ptm-tris_seqlogo-subset` = "All")
  session$flushReact()
  session$setInputs(`ptm-tris_seqlogo-variable` = "win20")
  session$flushReact(); session$flushReact()
  w5k_pick <<- store_read(w5k$child)
})
ok(is.null(w5k_err), "WP-5 ptm: retention boot stays alive")
ok(ut_cmp_identical(w5k_nobs, 2L),
   "WP-5 ptm: exactly the bind-sync + seed observers land in the context's panel_keep list (D3 unification)")
ok(ut_cmp_identical(unname(unlist(w5k_pick)), c("SeqLogo", "All", "win20")),
   "WP-5 ptm: kept sync observer still fires after gc() x2 (user pick reaches the store)")

# --- app level: the REAL L1 composition's status payload keeps its key
# order and its analyst_ptm entry byte-identical (pre-migration capture
# d); the D12 enumeration freeze earlier in this suite already pins the
# store keys of the same boot (72 ids incl. the three ptm keys)
w5l1 <- new.env()
w5l1$store <- widget_store_new()
app_w5l1 <- function(input, output, session) {
  w5l1$status <- omicsViewer:::L1_result_space_module(
    "rs",
    reactive_expr = reactive(Biobase::exprs(dat)),
    reactive_phenoData = reactive(Biobase::pData(dat)),
    reactive_featureData = reactive(fd0),
    status = reactive(NULL),
    store = w5l1$store)
}
w5l1_err <- NULL; w5l1_payload <- NULL
shiny::testServer(app_w5l1, {
  w5l1_err <<- tryCatch({
    for (i in 1:3) session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w5l1_payload <<- tryCatch(shiny::isolate(w5l1$status()),
                            error = function(e) conditionMessage(e))
})
ok(is.null(w5l1_err),
   "WP-5 ptm: the real L1 composition boots headless post-migration")
ok(ut_cmp_identical(names(w5l1_payload), PANEL_W5_PTM_GOLDEN_L1_STATUS$names),
   "WP-5 ptm: L1 status payload keeps its exact key order (snapshot schema untouched)")
ok(ut_cmp_identical(w5l1_payload$analyst_ptm, PANEL_W5_PTM_GOLDEN_L1_STATUS$analyst_ptm),
   "WP-5 ptm: analyst_ptm payload byte-identical pre/post migration (list(xax = NULL) in a cold boot)")
ok(ut_cmp_identical(
  sort(intersect(vapply(store_registry_view(w5l1$store), function(r) r$id, ""),
                 PANEL_W5_PTM_IDS)),
  PANEL_W5_PTM_IDS),
   "WP-5 ptm: the real L1 boot carries the three ptm keys (the D12 freeze above pins the full 72-id vector)")

## --------------- WP-6: Geneshot migration + P9 (D3/D4/D13) ----
# The SECOND pure migration WP (the WP-5 pattern) PLUS the P9 triset
# dedup: the module's twin cascade derivations (the `triset` reactive
# and the `.gs_ts` choices provider, with divergent empty-data guards)
# collapsed into ONE geneshot_triset(fd). The goldens below were
# captured from the PRE-MIGRATION HEAD (16c6cc6, 2026-11-06) BEFORE any
# edit, by booting the then-current geneshot_module exactly like the L1
# call site (child store resultspace.geneshot over demo.RDS fd0, feature
# 1 selected) and dput()ing the outputs - the WP-0 ORA-golden /
# WP-5 PTM fixture pattern. They are FROZEN: never regenerate them from
# the migrated code (that would make the equivalence claim circular).
geneshot_module <- omicsViewer:::geneshot_module
geneshot_triset <- omicsViewer:::geneshot_triset
PANEL_W6_GS_IDS <- sort(c(
  "resultspace.geneshot.term",
  "resultspace.geneshot.xax_analysis",
  "resultspace.geneshot.xax_subset",
  "resultspace.geneshot.xax_variable"))
# frozen pre-migration golden: BOTH pre-fix P9 derivations return this
# same 4x3 matrix over the demo fd (the derivations only diverged on
# crafted empty data - see the crafted cases below)
PANEL_W6_GS_TRIPLET <- structure(
  c("General", "General", "StringDB", "SeqLogo",
    "All", "All", "All", "All",
    "Gene.name", "Protein.ID", "ID", "seq.window"),
  dim = c(4L, 3L))
# frozen pre-migration golden: the registry view over demo fd0 -
# id/kind/label/help/depends_on/min/max/allowed_values of all four
# bindings (term first: the pre-WP-6 single store_register call produced
# exactly this registration order)
PANEL_W6_GS_GOLDEN_VIEW <- list(
  list(id = "resultspace.geneshot.term", kind = "string",
       label = "Search terms",
       help = paste("Semicolon-separated literature search terms for the",
                    "Geneshot AutoRIF query"),
       depends_on = character(0), min = NULL, max = NULL,
       allowed_values = NULL),
  list(id = "resultspace.geneshot.xax_analysis", kind = "select",
       label = "ID mapper category",
       help = "Annotation category used to map Geneshot genes to features",
       depends_on = character(0), min = NULL, max = NULL,
       allowed_values = c("General", "StringDB", "SeqLogo")),
  list(id = "resultspace.geneshot.xax_subset", kind = "select",
       label = "ID mapper subcategory",
       help = "Subcategory within the ID-mapper category",
       depends_on = "resultspace.geneshot.xax_analysis", min = NULL,
       max = NULL, allowed_values = character(0)),
  list(id = "resultspace.geneshot.xax_variable", kind = "select_cascaded",
       label = "ID mapper column",
       help = paste("Feature annotation column Geneshot genes are matched",
                    "against for the overlap highlight"),
       depends_on = c("resultspace.geneshot.xax_analysis",
                      "resultspace.geneshot.xax_subset"),
       min = NULL, max = NULL,
       allowed_values = c("Gene.name", "Protein.ID", "ID", "seq.window")))
# frozen pre-migration golden: the raw binding metadata incl. the
# agent_writable flags (all four agent-visible AND writable)
PANEL_W6_GS_GOLDEN_META <- list(
  list(id = "resultspace.geneshot.term", kind = "string",
       label = "Search terms",
       help = paste("Semicolon-separated literature search terms for the",
                    "Geneshot AutoRIF query"),
       depends_on = character(0), agent_writable = TRUE),
  list(id = "resultspace.geneshot.xax_analysis", kind = "select",
       label = "ID mapper category",
       help = "Annotation category used to map Geneshot genes to features",
       depends_on = character(0), agent_writable = TRUE),
  list(id = "resultspace.geneshot.xax_subset", kind = "select",
       label = "ID mapper subcategory",
       help = "Subcategory within the ID-mapper category",
       depends_on = "resultspace.geneshot.xax_analysis",
       agent_writable = TRUE),
  list(id = "resultspace.geneshot.xax_variable", kind = "select_cascaded",
       label = "ID mapper column",
       help = paste("Feature annotation column Geneshot genes are matched",
                    "against for the overlap highlight"),
       depends_on = c("resultspace.geneshot.xax_analysis",
                      "resultspace.geneshot.xax_subset"),
       agent_writable = TRUE))
# frozen pre-migration golden: the cold store values (geneshot has NO
# cascade auto-default - the cascade fills on the first pick - and the
# term seeds only once a non-empty term exists)
PANEL_W6_GS_GOLDEN_COLD <- list(
  `resultspace.geneshot.term` = NULL,
  `resultspace.geneshot.xax_analysis` = NULL,
  `resultspace.geneshot.xax_subset` = NULL,
  `resultspace.geneshot.xax_variable` = NULL)
# frozen pre-migration golden: the module return (the analyst_gene_shot
# payload source) in a cold headless boot
PANEL_W6_GS_GOLDEN_RET <- list(xax = NULL, term = NULL)
# frozen pre-migration golden: the real-L1 status payload key order and
# the analyst_gene_shot entry
PANEL_W6_GS_GOLDEN_L1_STATUS <- list(
  names = c("analyst_active_tab", "analyst_feature_general",
            "analyst_sample_general", "analyst_gene_shot", "analyst_fgsea",
            "analyst_ora", "analyst_ptm", "analyst_dose_response",
            "analyst_additional"),
  analyst_gene_shot = list(xax = NULL, term = NULL))
# shared boot factories: context-built (the WP-6 L1 call-site shape;
# holder receives the context so the keep closure is inspectable) and
# legacy direct (formals only - the pre-WP-6 shape, still supported)
.w6_mkctx <- function(child, holder, fd = fd0, feat = 1L,
                      status_payload = NULL) {
  force(child); force(holder)
  function(input, output, session) {
    .fd <- reactive(fd); .ft <- reactive(feat); .st <- reactive(status_payload)
    holder$ctx <- panel_context(
      id = "geneshotTab", tab_value = "Geneshot", store = child,
      data = list(fdata = .fd), selection = list(feature = .ft),
      status = .st, session = session)
    holder$ret <- geneshot_module(
      "geneshotTab", pdata = reactive(NULL), fdata = .fd,
      expr = reactive(NULL), feature_selected = .ft,
      sample_selected = reactive(NULL), object = reactive(NULL),
      reactive_status = .st, store = child, context = holder$ctx)
  }
}
.w6_mklegacy <- function(child, fd = fd0, feat = 1L, status_payload = NULL) {
  force(child)
  function(input, output, session) {
    geneshot_module(
      "geneshotTab", pdata = reactive(NULL), fdata = reactive(fd),
      expr = reactive(NULL), feature_selected = reactive(feat),
      sample_selected = reactive(NULL), object = reactive(NULL),
      reactive_status = reactive(status_payload), store = child)
  }
}
.w6_meta <- function(store) {
  lapply(grep("^resultspace\\.geneshot\\.", names(store$bindings),
              value = TRUE),
         function(k) {
           b <- store$bindings[[k]]
           list(id = k, kind = b$kind, label = b$label, help = b$help,
                depends_on = b$depends_on, agent_writable = b$agent_writable)
         })
}

# --- P9: the ONE derivation equals BOTH pre-fix derivations over the
# demo fd (they were captured separately pre-migration and agreed), and
# the plan-sanctioned unified semantics on the crafted divergent cases
# (the old `triset` reactive returned a 0x3 matrix there; the old .gs_ts
# returned NULL - the shared function adopts the stricter NULL guard)
ok(ut_cmp_identical(geneshot_triset(fd0), PANEL_W6_GS_TRIPLET),
   "WP-6 geneshot P9: geneshot_triset(demo fd) equals BOTH captured pre-fix derivations (they agreed over real data)")
ok(ut_cmp_identical(
  geneshot_triset(data.frame(a = 1:8, b = 9:16)), NULL) &&
    ut_cmp_identical(
      geneshot_triset(data.frame(`GS|x|y` = letters[1:8],
                                 row.names = paste0("F", 1:8),
                                 check.names = FALSE)), NULL) &&
    ut_cmp_identical(geneshot_triset(NULL), NULL),
   "WP-6 geneshot P9: crafted divergent cases (all-numeric / GS-only / NULL fd) pin the UNIFIED stricter guard (NULL; the plan's sanctioned semantic, old triset reactive gave a 0x3 matrix)")
ok(ut_cmp_identical(
  geneshot_triset(data.frame(
    num = 1:8, `GS|x|y` = letters[1:8], `General|All|Gene` = paste0("g", 1:8),
    row.names = paste0("F", 1:8), check.names = FALSE)),
  structure(c("General", "All", "Gene"), dim = c(1L, 3L))),
   "WP-6 geneshot P9: numeric and GS columns are filtered out of the cascade (mixed fd yields the eligible column only)")

# --- context boot equals the frozen pre-migration golden, byte-for-byte
w6 <- new.env()
w6$store <- widget_store_new()
w6$child <- widget_store_child(w6$store, "resultspace.geneshot")
w6_err <- NULL; w6_view <- NULL; w6_meta <- NULL; w6_vals <- NULL
w6_ret <- NULL
shiny::testServer(.w6_mkctx(w6$child, w6), {
  w6_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  # registry view INSIDE the session (WP-3 lore: choices providers read
  # live reactives; outside a reactive context they degrade to empty)
  w6_view <<- store_registry_view(w6$store, "resultspace.geneshot")
  w6_meta <<- .w6_meta(w6$store)
  w6_vals <<- store_read(w6$child)
  w6_ret <<- tryCatch(shiny::isolate(w6$ret()), error = function(e) NULL)
})
ok(is.null(w6_err),
   "WP-6 geneshot: module boots headless context-built (the new L1 call-site shape)")
ok(ut_cmp_identical(sort(vapply(w6_view, function(r) r$id, "")), PANEL_W6_GS_IDS),
   "WP-6 geneshot: exactly the four frozen bindings register under resultspace.geneshot (D12 untouched)")
ok(ut_cmp_identical(w6_view, PANEL_W6_GS_GOLDEN_VIEW),
   "WP-6 geneshot: registry view byte-identical to the frozen pre-migration golden (term first + the cascade trio over demo fd)")
ok(ut_cmp_identical(w6_meta, PANEL_W6_GS_GOLDEN_META),
   "WP-6 geneshot: raw binding meta byte-identical incl. agent_writable TRUE x4")
ok(ut_cmp_identical(w6_vals, PANEL_W6_GS_GOLDEN_COLD),
   "WP-6 geneshot: cold store values byte-identical (no cascade auto-default; term unseeded)")
ok(ut_cmp_identical(w6_ret, PANEL_W6_GS_GOLDEN_RET),
   "WP-6 geneshot: module return (the analyst_gene_shot payload source) stays list(xax = NULL, term = NULL) before any UI commit")

# --- legacy direct boot registers identically to the context boot
w6l <- new.env()
w6l$store <- widget_store_new()
w6l$child <- widget_store_child(w6l$store, "resultspace.geneshot")
w6l_err <- NULL; w6l_view <- NULL; w6l_meta <- NULL
shiny::testServer(.w6_mklegacy(w6l$child), {
  w6l_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w6l_view <<- store_registry_view(w6l$store, "resultspace.geneshot")
  w6l_meta <<- .w6_meta(w6l$store)
})
ok(is.null(w6l_err),
   "WP-6 geneshot: legacy direct instantiation (formals only, no context) still boots")
ok(ut_cmp_identical(w6l_view, PANEL_W6_GS_GOLDEN_VIEW) &&
   ut_cmp_identical(w6l_meta, PANEL_W6_GS_GOLDEN_META),
   "WP-6 geneshot: legacy boot registers byte-identically to the golden (both boot paths agree)")
ok(ut_cmp_identical(w6l_view, w6_view),
   "WP-6 geneshot: context boot and legacy boot register identically (D13 argument coalescing)")

# --- behavior pins (pre-migration captures): user cascade pick + term
# edit sync into the store; the module return reports the settled
# triple; the cascade choices stay state-dependent
w6b <- new.env()
w6b$store <- widget_store_new()
w6b$child <- widget_store_child(w6b$store, "resultspace.geneshot")
w6b_err <- NULL; w6b_vals <- NULL; w6b_ret <- NULL; w6b_ch <- NULL
shiny::testServer(.w6_mkctx(w6b$child, w6b), {
  w6b_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  session$setInputs(`geneshotTab-geneNameCol-analysis` = "General")
  session$flushReact()
  session$setInputs(`geneshotTab-geneNameCol-subset` = "All")
  session$flushReact()
  session$setInputs(`geneshotTab-geneNameCol-variable` = "Gene.name")
  session$flushReact(); session$flushReact()
  session$setInputs(`geneshotTab-term` = "p53; BRCA1")
  session$flushReact(); session$flushReact()
  w6b_vals <<- store_read(w6b$child)
  w6b_ret <<- tryCatch(shiny::isolate(w6b$ret()), error = function(e) NULL)
  v <- store_registry_view(w6b$store, "resultspace.geneshot")
  w6b_ch <<- stats::setNames(lapply(v, `[[`, "allowed_values"),
                            vapply(v, function(r) r$id, ""))
})
ok(is.null(w6b_err), "WP-6 geneshot: user-pick boot stays alive")
ok(ut_cmp_identical(
  unname(unlist(w6b_vals[c("resultspace.geneshot.xax_analysis",
                    "resultspace.geneshot.xax_subset",
                    "resultspace.geneshot.xax_variable")])),
  c("General", "All", "Gene.name")),
   "WP-6 geneshot: the user ID-mapper pick syncs into the store (pre-migration behavior)")
ok(ut_cmp_identical(w6b_vals$resultspace.geneshot.term, "p53; BRCA1"),
   "WP-6 geneshot: the user search-term edit syncs into the store (pre-migration behavior)")
ok(ut_cmp_identical(
  w6b_ret,
  list(xax = list(analysis = "General", subset = "All",
                  variable = "Gene.name"),
       term = "p53; BRCA1")),
   "WP-6 geneshot: module return reports the settled triple + term after the commits (analyst_gene_shot payload content)")
ok(ut_cmp_identical(
  w6b_ch,
  list(`resultspace.geneshot.term` = NULL,
       `resultspace.geneshot.xax_analysis` = c("General", "StringDB", "SeqLogo"),
       `resultspace.geneshot.xax_subset` = "All",
       `resultspace.geneshot.xax_variable` = c("Gene.name", "Protein.ID"))),
   "WP-6 geneshot: cascade choices stay state-dependent (subset narrows, variable filtered by the General branch)")

# --- legacy status restore path unchanged: an analyst_gene_shot payload
# with term + triple re-applies through store_apply(origin = "restore")
# (pre-migration capture; geneshot has no auto-seed, so the restore IS
# the store's initial cascade content)
w6r <- new.env()
w6r$store <- widget_store_new()
w6r$child <- widget_store_child(w6r$store, "resultspace.geneshot")
w6r_err <- NULL; w6r_vals <- NULL; w6r_ret <- NULL
shiny::testServer(.w6_mkctx(
  w6r$child, w6r,
  status_payload = list(term = "restored term",
                        xax = c("StringDB", "All", "ID"))), {
  w6r_err <<- tryCatch({
    for (i in 1:5) session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w6r_vals <<- store_read(w6r$child)
  w6r_ret <<- tryCatch(shiny::isolate(w6r$ret()), error = function(e) NULL)
})
ok(is.null(w6r_err),
   "WP-6 geneshot: boot with a legacy analyst_gene_shot status payload stays alive")
ok(ut_cmp_identical(
  unname(unlist(w6r_vals[c("resultspace.geneshot.xax_analysis",
                    "resultspace.geneshot.xax_subset",
                    "resultspace.geneshot.xax_variable")])),
  c("StringDB", "All", "ID")),
   "WP-6 geneshot: the legacy restore path re-applies the payload triple (single transactional path preserved)")
ok(ut_cmp_identical(w6r_vals$resultspace.geneshot.term, "restored term"),
   "WP-6 geneshot: the legacy restore path re-applies the payload term")
ok(ut_cmp_identical(w6r_ret, list(xax = NULL, term = NULL)),
   "WP-6 geneshot: module return stays cold after a restore-only boot (v1 commits through the UI, not the store)")

# --- retention (D3): the shared per-instance panel_keep holds the
# helper's bind-sync observer AND the module's term sync/seed/push
# observers (structural, deterministic - WP-0 MINOR-2 pattern); after
# gc() x2 a user pick still reaches the store through the kept
# observers (the WP-0 gc pattern, retention direction only per PD-8)
w6k <- new.env()
w6k$store <- widget_store_new()
w6k$child <- widget_store_child(w6k$store, "resultspace.geneshot")
w6k_err <- NULL; w6k_nobs <- NULL; w6k_vals <- NULL
shiny::testServer(.w6_mkctx(w6k$child, w6k), {
  w6k_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  # structural: everything the module kept lives in the CONTEXT's
  # panel_keep closure (D3 - the retention owner), not a module-local list
  obs <- get("obs", envir = environment(w6k$ctx$keep))
  w6k_nobs <<- length(obs)
  # gc integration (retention direction): collect twice, then user
  # picks must still reach the store through the kept observers
  # (the context - the L1-frame analogue - is the observers' only anchor)
  gc(); gc()
  session$setInputs(`geneshotTab-geneNameCol-analysis` = "SeqLogo")
  session$flushReact()
  session$setInputs(`geneshotTab-geneNameCol-subset` = "All")
  session$flushReact()
  session$setInputs(`geneshotTab-geneNameCol-variable` = "seq.window")
  session$flushReact(); session$flushReact()
  session$setInputs(`geneshotTab-term` = "post-gc term")
  session$flushReact(); session$flushReact()
  w6k_vals <<- store_read(w6k$child)
})
ok(is.null(w6k_err), "WP-6 geneshot: retention boot stays alive")
ok(ut_cmp_identical(w6k_nobs, 4L),
   "WP-6 geneshot: exactly the bind-sync + term-sync + term-seed + term-push observers land in the context's panel_keep list (D3 unification)")
ok(ut_cmp_identical(
  unname(unlist(w6k_vals[c("resultspace.geneshot.xax_analysis",
                    "resultspace.geneshot.xax_subset",
                    "resultspace.geneshot.xax_variable")])),
  c("SeqLogo", "All", "seq.window")),
   "WP-6 geneshot: the kept cascade sync observer still fires after gc() x2 (user pick reaches the store)")
ok(ut_cmp_identical(w6k_vals$resultspace.geneshot.term, "post-gc term"),
   "WP-6 geneshot: the kept term sync observer still fires after gc() x2 (term edit reaches the store)")

# --- app level: the REAL L1 composition's status payload keeps its key
# order and its analyst_gene_shot entry byte-identical (pre-migration
# capture); the D12 enumeration freeze earlier in this suite already
# pins the store keys of the same boot (72 ids incl. the four geneshot
# keys)
w6l1 <- new.env()
w6l1$store <- widget_store_new()
app_w6l1 <- function(input, output, session) {
  w6l1$status <- omicsViewer:::L1_result_space_module(
    "rs",
    reactive_expr = reactive(Biobase::exprs(dat)),
    reactive_phenoData = reactive(Biobase::pData(dat)),
    reactive_featureData = reactive(fd0),
    status = reactive(NULL),
    store = w6l1$store)
}
w6l1_err <- NULL; w6l1_payload <- NULL
shiny::testServer(app_w6l1, {
  w6l1_err <<- tryCatch({
    for (i in 1:3) session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w6l1_payload <<- tryCatch(shiny::isolate(w6l1$status()),
                            error = function(e) conditionMessage(e))
})
ok(is.null(w6l1_err),
   "WP-6 geneshot: the real L1 composition boots headless post-migration")
ok(ut_cmp_identical(names(w6l1_payload), PANEL_W6_GS_GOLDEN_L1_STATUS$names),
   "WP-6 geneshot: L1 status payload keeps its exact key order (snapshot schema untouched)")
ok(ut_cmp_identical(w6l1_payload$analyst_gene_shot,
                   PANEL_W6_GS_GOLDEN_L1_STATUS$analyst_gene_shot),
   "WP-6 geneshot: analyst_gene_shot payload byte-identical pre/post migration (list(xax = NULL, term = NULL) in a cold boot)")
ok(ut_cmp_identical(
  sort(intersect(vapply(store_registry_view(w6l1$store), function(r) r$id, ""),
                 PANEL_W6_GS_IDS)),
  PANEL_W6_GS_IDS),
   "WP-6 geneshot: the real L1 boot carries the four geneshot keys (the D12 freeze above pins the full 72-id vector)")

## --------------- WP-7: Feature-general migration (D3/D4/D13) ----
# The LARGEST pure migration (the WP-5/WP-6 pattern): feature-general
# composes attr4selector, boxplot, plotly-scatter and ROC/PR submodules,
# and its store surface is 24 keys (5 own + 19 attr4-registered through
# this panel's child store). SUBMODULES ARE UNTOUCHED (the plan's WP-7
# scope is panel glue only; attr4/plotly/boxplot/ROC internals are
# front-run boundaries of later work, not this WP). The migration
# unifies retention (.fg_keep + .fg_view_keep -> ONE panel_keep
# instance, D3), moves the cascade boilerplate into panel_triset_bind
# (D4: triset_source = trisetter(combine = "pheno") minus Surv rows,
# allow_unset = TRUE) and threads the context (D13).
#
# Golden-fixture provenance (PD-29/PD-33 school): all goldens below
# were captured from PRE-MIGRATION HEAD abc068a BEFORE any edit, by
# booting the then-current feature_general_module exactly like the L1
# call site and dput()ing the outputs. Capture commands (re-runnable;
# the identical script was re-run post-migration and every output was
# byte-identical - the PD-33 pre/post proof):
#   boot: shiny::testServer(function(input, output, session)
#     feature_general_module("feature_general",
#       reactive_expr = reactive(<ex>), reactive_i = reactive(<feat>),
#       reactive_highlight = reactive(NULL),
#       reactive_phenoData = reactive(<pd>),
#       reactive_featureData = reactive(<fd>),
#       reactive_status = reactive(<payload>),
#       store = widget_store_child(widget_store_new(),
#         "resultspace.feature_general")))
#   with <pd>/<ex>/<fd> = demo.RDS pData/exprs/fData (feat 1) OR the
#   synthetic agentWidgets-shaped data below (feat 3); captures:
#   store_registry_view(store, "resultspace.feature_general") (all 24
#   bindings), the raw binding meta incl. agent_writable, store_read,
#   isolate(module return), user-pick/restore triples over both
#   datasets, the real-L1 status payload, and the triset derivation
#   over demo pdata. The fixtures are FROZEN: never regenerate them
#   from the migrated code (that would make the equivalence claim
#   circular).
feature_general_module <- omicsViewer:::feature_general_module

# -------------------------------------------------- frozen goldens ----
# the 24 ids (5 own + 19 attr4, all under the panel's child store)
PANEL_W7_FG_IDS <- sort(c(
  paste0("resultspace.feature_general.", c(
    "xax_analysis", "xax_subset", "xax_variable", "plot_type",
    "regression_line")),
  paste0("resultspace.feature_general.attr4.", c(
    outer(c("color", "shape", "size", "tooltip", "search"),
          c("_analysis", "_subset", "_variable"), FUN = paste0),
    "search_value", "xcut", "ycut", "scorner"))))
# frozen pre-migration golden: the full registry view over the
# synthetic data (the agentWidgets .fg_pd/.fg_ex/.fg_fd shape) - all 24
# bindings, id/kind/label/help/depends_on/min/max/allowed_values, in
# registration order (xax_analysis, xax_subset, xax_variable,
# plot_type, regression_line, then the attr4 registrations flowing
# through this panel's child store)
PANEL_W7_FG_GOLDEN_VIEW <- list(list(id = "resultspace.feature_general.xax_analysis", kind = "select", 
    label = "Link variable category", help = "Annotation category of the variable the selected features are compared across", 
    depends_on = character(0), min = NULL, max = NULL, allowed_values = c("General", 
    "Feature")), list(id = "resultspace.feature_general.xax_subset", 
    kind = "select", label = "Link variable subcategory", help = "Subcategory within the link-variable category", 
    depends_on = "resultspace.feature_general.xax_analysis", 
    min = NULL, max = NULL, allowed_values = character(0)), list(
    id = "resultspace.feature_general.xax_variable", kind = "select_cascaded", 
    label = "Link variable", help = "Variable the selected features are plotted or tested against (phenotype annotation or Feature|Auto| row)", 
    depends_on = c("resultspace.feature_general.xax_analysis", 
    "resultspace.feature_general.xax_subset"), min = NULL, max = NULL, 
    allowed_values = c("group", "score", "F1", "F2", "F3", "F4", 
    "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12", "F13", 
    "F14", "F15", "F16")), list(id = "resultspace.feature_general.plot_type", 
    kind = "enum", label = "Plot type", help = "Beeswarm/box view or ROC/PR classification curves (categorical link variable with few features)", 
    depends_on = character(0), min = NULL, max = NULL, allowed_values = c("Bees", 
    "Curve")), list(id = "resultspace.feature_general.regression_line", 
    kind = "boolean", label = "Regression line", help = "Draw a regression line in the single-feature scatter view", 
    depends_on = character(0), min = NULL, max = NULL, allowed_values = NULL), 
    list(id = "resultspace.feature_general.attr4.color_analysis", 
        kind = "select", label = "Color attribute category", 
        help = "First level of the color attribute cascade (annotation category) to color the figure's points by a group or annotation column", 
        depends_on = character(0), min = NULL, max = NULL, allowed_values = c("General", 
        "Feature")), list(id = "resultspace.feature_general.attr4.color_subset", 
        kind = "select", label = "Color attribute subcategory", 
        help = "Second level of the color attribute cascade to color the figure's points by a group or annotation column", 
        depends_on = "resultspace.feature_general.attr4.color_analysis", 
        min = NULL, max = NULL, allowed_values = character(0)), 
    list(id = "resultspace.feature_general.attr4.color_variable", 
        kind = "select_cascaded", label = "Color attribute variable", 
        help = "Annotation variable to color the figure's points by a group or annotation column", 
        depends_on = c("resultspace.feature_general.attr4.color_analysis", 
        "resultspace.feature_general.attr4.color_subset"), min = NULL, 
        max = NULL, allowed_values = c("group", "score", "F1", 
        "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", 
        "F11", "F12", "F13", "F14", "F15", "F16")), list(id = "resultspace.feature_general.attr4.shape_analysis", 
        kind = "select", label = "Shape attribute category", 
        help = "First level of the shape attribute cascade (annotation category) to give the figure's points a different shape per group", 
        depends_on = character(0), min = NULL, max = NULL, allowed_values = c("General", 
        "Feature")), list(id = "resultspace.feature_general.attr4.shape_subset", 
        kind = "select", label = "Shape attribute subcategory", 
        help = "Second level of the shape attribute cascade to give the figure's points a different shape per group", 
        depends_on = "resultspace.feature_general.attr4.shape_analysis", 
        min = NULL, max = NULL, allowed_values = character(0)), 
    list(id = "resultspace.feature_general.attr4.shape_variable", 
        kind = "select_cascaded", label = "Shape attribute variable", 
        help = "Annotation variable to give the figure's points a different shape per group", 
        depends_on = c("resultspace.feature_general.attr4.shape_analysis", 
        "resultspace.feature_general.attr4.shape_subset"), min = NULL, 
        max = NULL, allowed_values = c("group", "score", "F1", 
        "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", 
        "F11", "F12", "F13", "F14", "F15", "F16")), list(id = "resultspace.feature_general.attr4.size_analysis", 
        kind = "select", label = "Size attribute category", help = "First level of the size attribute cascade (annotation category) to scale the figure's point size by a numeric column", 
        depends_on = character(0), min = NULL, max = NULL, allowed_values = c("General", 
        "Feature")), list(id = "resultspace.feature_general.attr4.size_subset", 
        kind = "select", label = "Size attribute subcategory", 
        help = "Second level of the size attribute cascade to scale the figure's point size by a numeric column", 
        depends_on = "resultspace.feature_general.attr4.size_analysis", 
        min = NULL, max = NULL, allowed_values = character(0)), 
    list(id = "resultspace.feature_general.attr4.size_variable", 
        kind = "select_cascaded", label = "Size attribute variable", 
        help = "Annotation variable to scale the figure's point size by a numeric column", 
        depends_on = c("resultspace.feature_general.attr4.size_analysis", 
        "resultspace.feature_general.attr4.size_subset"), min = NULL, 
        max = NULL, allowed_values = c("group", "score", "F1", 
        "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", 
        "F11", "F12", "F13", "F14", "F15", "F16")), list(id = "resultspace.feature_general.attr4.tooltip_analysis", 
        kind = "select", label = "Tooltips attribute category", 
        help = "First level of the tooltips attribute cascade (annotation category) shown in the figure's hover tooltip", 
        depends_on = character(0), min = NULL, max = NULL, allowed_values = c("General", 
        "Feature")), list(id = "resultspace.feature_general.attr4.tooltip_subset", 
        kind = "select", label = "Tooltips attribute subcategory", 
        help = "Second level of the tooltips attribute cascade shown in the figure's hover tooltip", 
        depends_on = "resultspace.feature_general.attr4.tooltip_analysis", 
        min = NULL, max = NULL, allowed_values = character(0)), 
    list(id = "resultspace.feature_general.attr4.tooltip_variable", 
        kind = "select_cascaded", label = "Tooltips attribute variable", 
        help = "Annotation variable shown in the figure's hover tooltip", 
        depends_on = c("resultspace.feature_general.attr4.tooltip_analysis", 
        "resultspace.feature_general.attr4.tooltip_subset"), 
        min = NULL, max = NULL, allowed_values = c("group", "score", 
        "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", 
        "F10", "F11", "F12", "F13", "F14", "F15", "F16")), list(
        id = "resultspace.feature_general.attr4.search_analysis", 
        kind = "select", label = "Search attribute category", 
        help = "First level of the search attribute cascade (annotation category) whose values the gene highlight picker circles as open-ring markers on the figure", 
        depends_on = character(0), min = NULL, max = NULL, allowed_values = c("General", 
        "Feature")), list(id = "resultspace.feature_general.attr4.search_subset", 
        kind = "select", label = "Search attribute subcategory", 
        help = "Second level of the search attribute cascade whose values the gene highlight picker circles as open-ring markers on the figure", 
        depends_on = "resultspace.feature_general.attr4.search_analysis", 
        min = NULL, max = NULL, allowed_values = character(0)), 
    list(id = "resultspace.feature_general.attr4.search_variable", 
        kind = "select_cascaded", label = "Search attribute variable", 
        help = "Annotation variable whose values the gene highlight picker circles as open-ring markers on the figure", 
        depends_on = c("resultspace.feature_general.attr4.search_analysis", 
        "resultspace.feature_general.attr4.search_subset"), min = NULL, 
        max = NULL, allowed_values = c("group", "score", "F1", 
        "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", 
        "F11", "F12", "F13", "F14", "F15", "F16")), list(id = "resultspace.feature_general.attr4.search_value", 
        kind = "multi_select", label = "Search highlight values", 
        help = "Gene name or other value of the selected Search column to draw as open-ring highlights on the figure; set the Search column first; empty clears the highlight", 
        depends_on = "resultspace.feature_general.attr4.search_variable", 
        min = NULL, max = NULL, allowed_values = character(0)), 
    list(id = "resultspace.feature_general.attr4.xcut", kind = "string", 
        label = "X cutoff", help = "X-axis cutoff for area selection; a number or an R expression such as -log10(0.05)", 
        depends_on = character(0), min = NULL, max = NULL, allowed_values = NULL), 
    list(id = "resultspace.feature_general.attr4.ycut", kind = "string", 
        label = "Y cutoff", help = "Y-axis cutoff for area selection; a number or an R expression such as log10(2)", 
        depends_on = character(0), min = NULL, max = NULL, allowed_values = NULL), 
    list(id = "resultspace.feature_general.attr4.scorner", kind = "select", 
        label = "Selected area", help = "Plot region selected by the cutoffs (corner or volcano mode)", 
        depends_on = character(0), min = NULL, max = NULL, allowed_values = "None"))
# frozen pre-migration golden: the raw binding metadata incl. the
# agent_writable flags (all 24 in the agent surface AND writable)
PANEL_W7_FG_GOLDEN_META <- list(list(id = "resultspace.feature_general.xax_analysis", kind = "select", 
    label = "Link variable category", help = "Annotation category of the variable the selected features are compared across", 
    depends_on = character(0), agent_writable = TRUE), list(id = "resultspace.feature_general.xax_subset", 
    kind = "select", label = "Link variable subcategory", help = "Subcategory within the link-variable category", 
    depends_on = "resultspace.feature_general.xax_analysis", 
    agent_writable = TRUE), list(id = "resultspace.feature_general.xax_variable", 
    kind = "select_cascaded", label = "Link variable", help = "Variable the selected features are plotted or tested against (phenotype annotation or Feature|Auto| row)", 
    depends_on = c("resultspace.feature_general.xax_analysis", 
    "resultspace.feature_general.xax_subset"), agent_writable = TRUE), 
    list(id = "resultspace.feature_general.plot_type", kind = "enum", 
        label = "Plot type", help = "Beeswarm/box view or ROC/PR classification curves (categorical link variable with few features)", 
        depends_on = character(0), agent_writable = TRUE), list(
        id = "resultspace.feature_general.regression_line", kind = "boolean", 
        label = "Regression line", help = "Draw a regression line in the single-feature scatter view", 
        depends_on = character(0), agent_writable = TRUE), list(
        id = "resultspace.feature_general.attr4.color_analysis", 
        kind = "select", label = "Color attribute category", 
        help = "First level of the color attribute cascade (annotation category) to color the figure's points by a group or annotation column", 
        depends_on = character(0), agent_writable = TRUE), list(
        id = "resultspace.feature_general.attr4.color_subset", 
        kind = "select", label = "Color attribute subcategory", 
        help = "Second level of the color attribute cascade to color the figure's points by a group or annotation column", 
        depends_on = "resultspace.feature_general.attr4.color_analysis", 
        agent_writable = TRUE), list(id = "resultspace.feature_general.attr4.color_variable", 
        kind = "select_cascaded", label = "Color attribute variable", 
        help = "Annotation variable to color the figure's points by a group or annotation column", 
        depends_on = c("resultspace.feature_general.attr4.color_analysis", 
        "resultspace.feature_general.attr4.color_subset"), agent_writable = TRUE), 
    list(id = "resultspace.feature_general.attr4.shape_analysis", 
        kind = "select", label = "Shape attribute category", 
        help = "First level of the shape attribute cascade (annotation category) to give the figure's points a different shape per group", 
        depends_on = character(0), agent_writable = TRUE), list(
        id = "resultspace.feature_general.attr4.shape_subset", 
        kind = "select", label = "Shape attribute subcategory", 
        help = "Second level of the shape attribute cascade to give the figure's points a different shape per group", 
        depends_on = "resultspace.feature_general.attr4.shape_analysis", 
        agent_writable = TRUE), list(id = "resultspace.feature_general.attr4.shape_variable", 
        kind = "select_cascaded", label = "Shape attribute variable", 
        help = "Annotation variable to give the figure's points a different shape per group", 
        depends_on = c("resultspace.feature_general.attr4.shape_analysis", 
        "resultspace.feature_general.attr4.shape_subset"), agent_writable = TRUE), 
    list(id = "resultspace.feature_general.attr4.size_analysis", 
        kind = "select", label = "Size attribute category", help = "First level of the size attribute cascade (annotation category) to scale the figure's point size by a numeric column", 
        depends_on = character(0), agent_writable = TRUE), list(
        id = "resultspace.feature_general.attr4.size_subset", 
        kind = "select", label = "Size attribute subcategory", 
        help = "Second level of the size attribute cascade to scale the figure's point size by a numeric column", 
        depends_on = "resultspace.feature_general.attr4.size_analysis", 
        agent_writable = TRUE), list(id = "resultspace.feature_general.attr4.size_variable", 
        kind = "select_cascaded", label = "Size attribute variable", 
        help = "Annotation variable to scale the figure's point size by a numeric column", 
        depends_on = c("resultspace.feature_general.attr4.size_analysis", 
        "resultspace.feature_general.attr4.size_subset"), agent_writable = TRUE), 
    list(id = "resultspace.feature_general.attr4.tooltip_analysis", 
        kind = "select", label = "Tooltips attribute category", 
        help = "First level of the tooltips attribute cascade (annotation category) shown in the figure's hover tooltip", 
        depends_on = character(0), agent_writable = TRUE), list(
        id = "resultspace.feature_general.attr4.tooltip_subset", 
        kind = "select", label = "Tooltips attribute subcategory", 
        help = "Second level of the tooltips attribute cascade shown in the figure's hover tooltip", 
        depends_on = "resultspace.feature_general.attr4.tooltip_analysis", 
        agent_writable = TRUE), list(id = "resultspace.feature_general.attr4.tooltip_variable", 
        kind = "select_cascaded", label = "Tooltips attribute variable", 
        help = "Annotation variable shown in the figure's hover tooltip", 
        depends_on = c("resultspace.feature_general.attr4.tooltip_analysis", 
        "resultspace.feature_general.attr4.tooltip_subset"), 
        agent_writable = TRUE), list(id = "resultspace.feature_general.attr4.search_analysis", 
        kind = "select", label = "Search attribute category", 
        help = "First level of the search attribute cascade (annotation category) whose values the gene highlight picker circles as open-ring markers on the figure", 
        depends_on = character(0), agent_writable = TRUE), list(
        id = "resultspace.feature_general.attr4.search_subset", 
        kind = "select", label = "Search attribute subcategory", 
        help = "Second level of the search attribute cascade whose values the gene highlight picker circles as open-ring markers on the figure", 
        depends_on = "resultspace.feature_general.attr4.search_analysis", 
        agent_writable = TRUE), list(id = "resultspace.feature_general.attr4.search_variable", 
        kind = "select_cascaded", label = "Search attribute variable", 
        help = "Annotation variable whose values the gene highlight picker circles as open-ring markers on the figure", 
        depends_on = c("resultspace.feature_general.attr4.search_analysis", 
        "resultspace.feature_general.attr4.search_subset"), agent_writable = TRUE), 
    list(id = "resultspace.feature_general.attr4.search_value", 
        kind = "multi_select", label = "Search highlight values", 
        help = "Gene name or other value of the selected Search column to draw as open-ring highlights on the figure; set the Search column first; empty clears the highlight", 
        depends_on = "resultspace.feature_general.attr4.search_variable", 
        agent_writable = TRUE), list(id = "resultspace.feature_general.attr4.xcut", 
        kind = "string", label = "X cutoff", help = "X-axis cutoff for area selection; a number or an R expression such as -log10(0.05)", 
        depends_on = character(0), agent_writable = TRUE), list(
        id = "resultspace.feature_general.attr4.ycut", kind = "string", 
        label = "Y cutoff", help = "Y-axis cutoff for area selection; a number or an R expression such as log10(2)", 
        depends_on = character(0), agent_writable = TRUE), list(
        id = "resultspace.feature_general.attr4.scorner", kind = "select", 
        label = "Selected area", help = "Plot region selected by the cutoffs (corner or volcano mode)", 
        depends_on = character(0), agent_writable = TRUE))
# frozen pre-migration golden: the cold store values over the synthetic
# data (xax trio + plot_type NULL; regression_line FALSE is the store's
# boolean-kind default; attr4's search_value character(0) and scorner
# "None" defaults ride along untouched)
PANEL_W7_FG_GOLDEN_COLD <- list(resultspace.feature_general.xax_analysis = NULL, resultspace.feature_general.xax_subset = NULL, 
    resultspace.feature_general.xax_variable = NULL, resultspace.feature_general.plot_type = NULL, 
    resultspace.feature_general.regression_line = FALSE, resultspace.feature_general.attr4.color_analysis = NULL, 
    resultspace.feature_general.attr4.color_subset = NULL, resultspace.feature_general.attr4.color_variable = NULL, 
    resultspace.feature_general.attr4.shape_analysis = NULL, 
    resultspace.feature_general.attr4.shape_subset = NULL, resultspace.feature_general.attr4.shape_variable = NULL, 
    resultspace.feature_general.attr4.size_analysis = NULL, resultspace.feature_general.attr4.size_subset = NULL, 
    resultspace.feature_general.attr4.size_variable = NULL, resultspace.feature_general.attr4.tooltip_analysis = NULL, 
    resultspace.feature_general.attr4.tooltip_subset = NULL, 
    resultspace.feature_general.attr4.tooltip_variable = NULL, 
    resultspace.feature_general.attr4.search_analysis = NULL, 
    resultspace.feature_general.attr4.search_subset = NULL, resultspace.feature_general.attr4.search_variable = NULL, 
    resultspace.feature_general.attr4.search_value = character(0), 
    resultspace.feature_general.attr4.xcut = NULL, resultspace.feature_general.attr4.ycut = NULL, 
    resultspace.feature_general.attr4.scorner = "None")
# frozen pre-migration golden: the module return (the
# analyst_feature_general payload source) in a cold headless boot
PANEL_W7_FG_GOLDEN_RET <- list(xax = NULL, showRegLine = FALSE, attr4 = list(selectColor = NULL, 
    selectShape = NULL, selectSize = NULL, selectTooltip = NULL, 
    searchOnCol = NULL, searchValue = NULL, xcut = NULL, ycut = NULL, 
    acorner = NULL), plotType = NULL)
# frozen pre-migration golden: user picks + agent writes over the
# synthetic data - cascade General|All|group, radio warm->Curve
# (ignoreInit lore), agent regression_line = TRUE - full store values,
# the module return, and the own-binding allowed_values at that state
PANEL_W7_FG_GOLDEN_PICK_VALS <- list(resultspace.feature_general.xax_analysis = "General", resultspace.feature_general.xax_subset = "All", 
    resultspace.feature_general.xax_variable = "group", resultspace.feature_general.plot_type = "Curve", 
    resultspace.feature_general.regression_line = TRUE, resultspace.feature_general.attr4.color_analysis = NULL, 
    resultspace.feature_general.attr4.color_subset = NULL, resultspace.feature_general.attr4.color_variable = NULL, 
    resultspace.feature_general.attr4.shape_analysis = NULL, 
    resultspace.feature_general.attr4.shape_subset = NULL, resultspace.feature_general.attr4.shape_variable = NULL, 
    resultspace.feature_general.attr4.size_analysis = NULL, resultspace.feature_general.attr4.size_subset = NULL, 
    resultspace.feature_general.attr4.size_variable = NULL, resultspace.feature_general.attr4.tooltip_analysis = NULL, 
    resultspace.feature_general.attr4.tooltip_subset = NULL, 
    resultspace.feature_general.attr4.tooltip_variable = NULL, 
    resultspace.feature_general.attr4.search_analysis = NULL, 
    resultspace.feature_general.attr4.search_subset = NULL, resultspace.feature_general.attr4.search_variable = NULL, 
    resultspace.feature_general.attr4.search_value = character(0), 
    resultspace.feature_general.attr4.xcut = NULL, resultspace.feature_general.attr4.ycut = NULL, 
    resultspace.feature_general.attr4.scorner = "None")
PANEL_W7_FG_GOLDEN_PICK_RET <- list(xax = list(analysis = "General", subset = "All", variable = "group"), 
    showRegLine = TRUE, attr4 = list(selectColor = NULL, selectShape = NULL, 
        selectSize = NULL, selectTooltip = NULL, searchOnCol = NULL, 
        searchValue = NULL, xcut = NULL, ycut = NULL, acorner = NULL), 
    plotType = "Curve")
PANEL_W7_FG_GOLDEN_PICK_CH <- list(resultspace.feature_general.xax_analysis = c("General", 
"Feature"), resultspace.feature_general.xax_subset = "All", resultspace.feature_general.xax_variable = c("group", 
"score"), resultspace.feature_general.plot_type = c("Bees", "Curve"
), resultspace.feature_general.regression_line = NULL)
# frozen pre-migration golden: the legacy analyst_feature_general
# restore path (payload with xax triple + showRegLine + plotType
# re-applies through store_apply(origin = "restore")); the module
# return stays cold (v1 commits through the UI, not the store)
PANEL_W7_FG_GOLDEN_RESTORE_VALS <- list(resultspace.feature_general.xax_analysis = "General", resultspace.feature_general.xax_subset = "All", 
    resultspace.feature_general.xax_variable = "score", resultspace.feature_general.plot_type = "Curve", 
    resultspace.feature_general.regression_line = TRUE, resultspace.feature_general.attr4.color_analysis = NULL, 
    resultspace.feature_general.attr4.color_subset = NULL, resultspace.feature_general.attr4.color_variable = NULL, 
    resultspace.feature_general.attr4.shape_analysis = NULL, 
    resultspace.feature_general.attr4.shape_subset = NULL, resultspace.feature_general.attr4.shape_variable = NULL, 
    resultspace.feature_general.attr4.size_analysis = NULL, resultspace.feature_general.attr4.size_subset = NULL, 
    resultspace.feature_general.attr4.size_variable = NULL, resultspace.feature_general.attr4.tooltip_analysis = NULL, 
    resultspace.feature_general.attr4.tooltip_subset = NULL, 
    resultspace.feature_general.attr4.tooltip_variable = NULL, 
    resultspace.feature_general.attr4.search_analysis = NULL, 
    resultspace.feature_general.attr4.search_subset = NULL, resultspace.feature_general.attr4.search_variable = NULL, 
    resultspace.feature_general.attr4.search_value = character(0), 
    resultspace.feature_general.attr4.xcut = NULL, resultspace.feature_general.attr4.ycut = NULL, 
    resultspace.feature_general.attr4.scorner = "None")
PANEL_W7_FG_GOLDEN_RESTORE_RET <- list(xax = NULL, showRegLine = TRUE, attr4 = list(selectColor = NULL, 
    selectShape = NULL, selectSize = NULL, selectTooltip = NULL, 
    searchOnCol = NULL, searchValue = NULL, xcut = NULL, ycut = NULL, 
    acorner = NULL), plotType = NULL)
# frozen pre-migration golden: a saved triple with an UNSET variable
# pushes the canonical unset "--select--" so a restore REMOVES a link
# variable set after the snapshot was taken (todo 2.3 semantics,
# unique to this module among the result-space panels)
PANEL_W7_FG_GOLDEN_UNSET_VALS <- list(resultspace.feature_general.xax_analysis = "General", resultspace.feature_general.xax_subset = "All", 
    resultspace.feature_general.xax_variable = "--select--", 
    resultspace.feature_general.plot_type = "Curve", resultspace.feature_general.regression_line = TRUE)
# frozen pre-migration goldens over demo.RDS: the own-binding registry
# view (xax_analysis General/PCA/Feature; xax_variable head-50 capped
# by store_registry_view), the own-binding cold values, the module
# return, the PCA|All|PC1(10.5%) user pick (values + return + the
# state-dependent choices), the General|All|Cell.line restore triple,
# and two attr4 analysis-level bindings proving the attr4 registrations
# share the demo triset source through this panel's child store
PANEL_W7_FG_DEMO_OWN_VIEW <- list(list(id = "resultspace.feature_general.xax_analysis", kind = "select", 
    label = "Link variable category", help = "Annotation category of the variable the selected features are compared across", 
    depends_on = character(0), min = NULL, max = NULL, allowed_values = c("General", 
    "PCA", "Feature")), list(id = "resultspace.feature_general.xax_subset", 
    kind = "select", label = "Link variable subcategory", help = "Subcategory within the link-variable category", 
    depends_on = "resultspace.feature_general.xax_analysis", 
    min = NULL, max = NULL, allowed_values = character(0)), list(
    id = "resultspace.feature_general.xax_variable", kind = "select_cascaded", 
    label = "Link variable", help = "Variable the selected features are plotted or tested against (phenotype annotation or Feature|Auto| row)", 
    depends_on = c("resultspace.feature_general.xax_analysis", 
    "resultspace.feature_general.xax_subset"), min = NULL, max = NULL, 
    allowed_values = c("Cell.line", "MDR", "Doubleing.Time", 
    "Origin", "TP53.Status", "Gender", "numberOfFeatures", "PC1(10.5%)", 
    "PC2(7.2%)", "PC3(5.4%)", "PC4(4.8%)", "PC5(3.6%)", "PC6(2.9%)", 
    "PC7(2.8%)", "PC8(2.5%)", "Dose", "Curve", "X32.TYW5", "X52.CNOT1", 
    "X55.MAP1LC3B.MAP1LC3B2", "X73.METTL15", "X107.ESPN", "X120.WASH3P.WASH2P.WASH1", 
    "X266.KIF2A", "X267.KIF2A", "X269.DDX39A", "X273.ACOT7", 
    "X274.MYO1C", "X276.MYO1F", "X279.HAX1", "X283.STXBP3", "X288.AP3B1", 
    "X290.APBB1", "X291.LGALS8", "X294.NFKBIE", "X295.PSMD11", 
    "X296.PSMD12", "X297.PSMD9", "X298.ATOX1", "X299.MEN1", "X302.SUPT5H", 
    "X308.HIP1", "X311.WWP2", "X312.CDC7", "X313.PIK3CD", "X318.DCTN6", 
    "X320.WASL", "X321.IPO5", "X322.IPO5", "X327.IGF2BP3")), 
    list(id = "resultspace.feature_general.plot_type", kind = "enum", 
        label = "Plot type", help = "Beeswarm/box view or ROC/PR classification curves (categorical link variable with few features)", 
        depends_on = character(0), min = NULL, max = NULL, allowed_values = c("Bees", 
        "Curve")), list(id = "resultspace.feature_general.regression_line", 
        kind = "boolean", label = "Regression line", help = "Draw a regression line in the single-feature scatter view", 
        depends_on = character(0), min = NULL, max = NULL, allowed_values = NULL))
PANEL_W7_FG_DEMO_OWN_COLD <- list(resultspace.feature_general.xax_analysis = NULL, resultspace.feature_general.xax_subset = NULL, 
    resultspace.feature_general.xax_variable = NULL, resultspace.feature_general.plot_type = NULL, 
    resultspace.feature_general.regression_line = FALSE)
PANEL_W7_FG_DEMO_RET <- list(xax = NULL, showRegLine = FALSE, attr4 = list(selectColor = NULL, 
    selectShape = NULL, selectSize = NULL, selectTooltip = NULL, 
    searchOnCol = NULL, searchValue = NULL, xcut = NULL, ycut = NULL, 
    acorner = NULL), plotType = NULL)
PANEL_W7_FG_DEMO_PICK_VALS <- list(resultspace.feature_general.xax_analysis = "PCA", resultspace.feature_general.xax_subset = "All", 
    resultspace.feature_general.xax_variable = "PC1(10.5%)", 
    resultspace.feature_general.plot_type = NULL, resultspace.feature_general.regression_line = FALSE)
PANEL_W7_FG_DEMO_PICK_RET <- list(xax = list(analysis = "PCA", subset = "All", variable = "PC1(10.5%)"), 
    showRegLine = FALSE, attr4 = list(selectColor = NULL, selectShape = NULL, 
        selectSize = NULL, selectTooltip = NULL, searchOnCol = NULL, 
        searchValue = NULL, xcut = NULL, ycut = NULL, acorner = NULL), 
    plotType = NULL)
PANEL_W7_FG_DEMO_PICK_CH <- list(resultspace.feature_general.xax_analysis = c("General", 
"PCA", "Feature"), resultspace.feature_general.xax_subset = "All", 
    resultspace.feature_general.xax_variable = c("PC1(10.5%)", 
    "PC2(7.2%)", "PC3(5.4%)", "PC4(4.8%)", "PC5(3.6%)", "PC6(2.9%)", 
    "PC7(2.8%)", "PC8(2.5%)"), resultspace.feature_general.plot_type = c("Bees", 
    "Curve"), resultspace.feature_general.regression_line = NULL)
PANEL_W7_FG_DEMO_RESTORE_VALS <- list(resultspace.feature_general.xax_analysis = "General", resultspace.feature_general.xax_subset = "All", 
    resultspace.feature_general.xax_variable = "Cell.line", resultspace.feature_general.plot_type = "Curve", 
    resultspace.feature_general.regression_line = TRUE)
PANEL_W7_FG_DEMO_ATTR4_CH <- list(list(id = "resultspace.feature_general.attr4.color_analysis", 
    kind = "select", label = "Color attribute category", help = "First level of the color attribute cascade (annotation category) to color the figure's points by a group or annotation column", 
    depends_on = character(0), min = NULL, max = NULL, allowed_values = c("General", 
    "PCA", "Feature")), list(id = "resultspace.feature_general.attr4.search_analysis", 
    kind = "select", label = "Search attribute category", help = "First level of the search attribute cascade (annotation category) whose values the gene highlight picker circles as open-ring markers on the figure", 
    depends_on = character(0), min = NULL, max = NULL, allowed_values = c("General", 
    "PCA", "Feature")))
# frozen pre-migration golden: the triset derivation over demo pdata -
# trisetter(combine = "pheno") MINUS the Surv rows (dim 2719x3; levels
# General/PCA/Feature with Surv excluded; the 17 non-Feature rows; the
# 2 Surv rows present in the unfiltered triset; head-50 variable order)
PANEL_W7_FG_TRISET <- list(dim = c(2719L, 3L), analysis_levels = c("General", "PCA", 
"Feature"), non_feature_rows = structure(c("General", "General", 
"General", "General", "General", "General", "General", "PCA", 
"PCA", "PCA", "PCA", "PCA", "PCA", "PCA", "PCA", "General", "General", 
"All", "All", "All", "All", "All", "All", "All", "All", "All", 
"All", "All", "All", "All", "All", "All", "All", "All", "Cell.line", 
"MDR", "Doubleing.Time", "Origin", "TP53.Status", "Gender", "numberOfFeatures", 
"PC1(10.5%)", "PC2(7.2%)", "PC3(5.4%)", "PC4(4.8%)", "PC5(3.6%)", 
"PC6(2.9%)", "PC7(2.8%)", "PC8(2.5%)", "Dose", "Curve"), dim = c(17L, 
3L)), surv_rows_in_full = structure(c("Surv", "Surv", "all", 
"all", "OS", "DFS"), dim = 2:3), variable_head50 = c("Cell.line", 
"MDR", "Doubleing.Time", "Origin", "TP53.Status", "Gender", "numberOfFeatures", 
"PC1(10.5%)", "PC2(7.2%)", "PC3(5.4%)", "PC4(4.8%)", "PC5(3.6%)", 
"PC6(2.9%)", "PC7(2.8%)", "PC8(2.5%)", "Dose", "Curve", "X32.TYW5", 
"X52.CNOT1", "X55.MAP1LC3B.MAP1LC3B2", "X73.METTL15", "X107.ESPN", 
"X120.WASH3P.WASH2P.WASH1", "X266.KIF2A", "X267.KIF2A", "X269.DDX39A", 
"X273.ACOT7", "X274.MYO1C", "X276.MYO1F", "X279.HAX1", "X283.STXBP3", 
"X288.AP3B1", "X290.APBB1", "X291.LGALS8", "X294.NFKBIE", "X295.PSMD11", 
"X296.PSMD12", "X297.PSMD9", "X298.ATOX1", "X299.MEN1", "X302.SUPT5H", 
"X308.HIP1", "X311.WWP2", "X312.CDC7", "X313.PIK3CD", "X318.DCTN6", 
"X320.WASL", "X321.IPO5", "X322.IPO5", "X327.IGF2BP3"))
# frozen pre-migration golden: the real-L1 status payload key order and
# the analyst_feature_general entry
PANEL_W7_FG_GOLDEN_L1_STATUS <- list(
  names = c("analyst_active_tab", "analyst_feature_general", "analyst_sample_general", 
"analyst_gene_shot", "analyst_fgsea", "analyst_ora", "analyst_ptm", 
"analyst_dose_response", "analyst_additional"),
  fg = list(xax = NULL, showRegLine = FALSE, attr4 = list(selectColor = NULL, 
    selectShape = NULL, selectSize = NULL, selectTooltip = NULL, 
    searchOnCol = NULL, searchValue = NULL, xcut = NULL, ycut = NULL, 
    acorner = NULL), plotType = NULL))

# ------------------------------------------------- boot factories ----
# synthetic data: the agentWidgets feature_general shape (group/score
# pheno columns + a Surv column the cascade must EXCLUDE); seeded for
# determinism (values never enter the store surface, but pin anyway)
set.seed(7)
.w7_pd <- data.frame(
  `General|All|group` = factor(rep(c("A", "B"), each = 4)),
  `General|All|score` = rnorm(8),
  `Surv|All|time` = rexp(8),
  row.names = paste0("S", 1:8), check.names = FALSE)
.w7_ex <- matrix(rnorm(16 * 8), nrow = 16,
                 dimnames = list(paste0("F", 1:16), paste0("S", 1:8)))
.w7_fd <- data.frame(`General|All|Gene.name` = paste0("g", 1:16),
                     row.names = paste0("F", 1:16), check.names = FALSE)
# context-built (the WP-7 L1 call-site shape; holder receives the
# context so the keep closure is inspectable) and legacy direct
# (formals only - the pre-WP-7 shape, still supported)
.w7_mkctx <- function(holder, pd = .w7_pd, ex = .w7_ex, fd = .w7_fd,
                      feat = 3L, status_payload = NULL) {
  force(holder)
  function(input, output, session) {
    .pd <- reactive(pd); .ex <- reactive(ex); .fd <- reactive(fd)
    .ft <- reactive(feat); .st <- reactive(status_payload)
    holder$store <- widget_store_new()
    holder$child <- widget_store_child(holder$store,
                                       "resultspace.feature_general")
    holder$ctx <- panel_context(
      id = "feature_general", tab_value = "Feature", store = holder$child,
      data = list(expr = .ex, pdata = .pd, fdata = .fd),
      selection = list(feature = .ft, sample = reactive(NULL)),
      status = .st, session = session)
    holder$ret <- feature_general_module(
      "feature_general", reactive_expr = .ex, reactive_i = .ft,
      reactive_highlight = reactive(NULL), reactive_phenoData = .pd,
      reactive_featureData = .fd, reactive_status = .st,
      store = holder$child, context = holder$ctx)
  }
}
.w7_mklegacy <- function(pd = .w7_pd, ex = .w7_ex, fd = .w7_fd,
                          feat = 3L, status_payload = NULL) {
  function(input, output, session) {
    feature_general_module(
      "feature_general", reactive_expr = reactive(ex),
      reactive_i = reactive(feat), reactive_highlight = reactive(NULL),
      reactive_phenoData = reactive(pd),
      reactive_featureData = reactive(fd),
      reactive_status = reactive(status_payload))
  }
}
# store-backed legacy boot (formals only - the pre-WP-7 shape; the
# module return lands in env$ret for the payload assertions)
.w7_boot <- function(env, pd = .w7_pd, ex = .w7_ex, fd = .w7_fd,
                     feat = 3L, status_payload = NULL) {
  force(env)
  function(input, output, session) {
    env$store <- widget_store_new()
    env$child <- widget_store_child(env$store,
                                    "resultspace.feature_general")
    env$ret <- feature_general_module(
      "feature_general", reactive_expr = reactive(ex),
      reactive_i = reactive(feat), reactive_highlight = reactive(NULL),
      reactive_phenoData = reactive(pd),
      reactive_featureData = reactive(fd),
      reactive_status = reactive(status_payload),
      store = env$child)
  }
}
.w7_meta <- function(store) {
  lapply(grep("^resultspace\\.feature_general(\\.|$)", names(store$bindings),
              value = TRUE),
         function(k) {
           b <- store$bindings[[k]]
           list(id = k, kind = b$kind, label = b$label, help = b$help,
                depends_on = b$depends_on, agent_writable = b$agent_writable)
         })
}

# ------------------- registration: context boot = the frozen golden ----
w7 <- new.env()
w7_err <- NULL; w7_view <- NULL; w7_meta <- NULL
shiny::testServer(.w7_mkctx(w7), {
  w7_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  # registry view INSIDE the session (WP-3 lore: choices providers read
  # live reactives; outside a reactive context they degrade to empty)
  w7_view <<- store_registry_view(w7$store, "resultspace.feature_general")
  w7_meta <<- .w7_meta(w7$store)
})
ok(is.null(w7_err),
   "WP-7 fg: module boots headless context-built (the new L1 call-site shape)")
ok(ut_cmp_identical(sort(vapply(w7_view, function(r) r$id, "")),
                    PANEL_W7_FG_IDS),
   "WP-7 fg: exactly the 24 frozen bindings register (5 own + 19 attr4 through the child store; D12 untouched)")
ok(ut_cmp_identical(w7_view, PANEL_W7_FG_GOLDEN_VIEW),
   "WP-7 fg: registry view byte-identical to the frozen pre-migration golden (registration order preserved: xax trio, plot_type, regression_line, then attr4)")
ok(ut_cmp_identical(w7_meta, PANEL_W7_FG_GOLDEN_META),
   "WP-7 fg: raw binding meta byte-identical incl. agent_writable TRUE x24")

# ------------------- cold values + module return (context boot) ----
w7c <- new.env()
w7c_err <- NULL; w7c_cold <- NULL; w7c_ret <- NULL
shiny::testServer(.w7_mkctx(w7c), {
  w7c_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w7c_cold <<- store_read(w7c$store)
  w7c_ret <<- tryCatch(shiny::isolate(w7c$ret()), error = function(e) NULL)
})
ok(is.null(w7c_err), "WP-7 fg: cold-boot capture stays alive")
ok(ut_cmp_identical(w7c_cold, PANEL_W7_FG_GOLDEN_COLD),
   "WP-7 fg: cold store values byte-identical (cascade unset; boolean default FALSE; attr4 defaults untouched)")
ok(ut_cmp_identical(w7c_ret, PANEL_W7_FG_GOLDEN_RET),
   "WP-7 fg: module return (the analyst_feature_general payload source) byte-identical in a cold boot")

# ------------------- legacy boot registers identically (D13) ----
w7l <- new.env()
w7l_err <- NULL; w7l_view <- NULL; w7l_meta <- NULL; w7l_cold <- NULL
shiny::testServer(.w7_boot(w7l), {
  w7l_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w7l_view <<- store_registry_view(w7l$store, "resultspace.feature_general")
  w7l_meta <<- .w7_meta(w7l$store)
  w7l_cold <<- store_read(w7l$store)
})
ok(is.null(w7l_err),
   "WP-7 fg: legacy direct instantiation (store formal only, no context) still boots")
ok(ut_cmp_identical(w7l_view, PANEL_W7_FG_GOLDEN_VIEW) &&
   ut_cmp_identical(w7l_meta, PANEL_W7_FG_GOLDEN_META) &&
   ut_cmp_identical(w7l_cold, PANEL_W7_FG_GOLDEN_COLD),
   "WP-7 fg: legacy store boot registers and seeds byte-identically to the golden (both boot paths agree)")
ok(ut_cmp_identical(w7l_view, w7_view),
   "WP-7 fg: context boot and legacy boot register identically (D13 argument coalescing)")

# ------------------- cascade user picks + choices (synthetic) ----
w7p <- new.env()
w7p_err <- NULL; w7p_vals <- NULL; w7p_ret <- NULL; w7p_ch <- NULL
shiny::testServer(.w7_mkctx(w7p), {
  w7p_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  session$setInputs(`feature_general-tris_feature_general-analysis` = "General")
  session$flushReact()
  session$setInputs(`feature_general-tris_feature_general-subset` = "All")
  session$flushReact()
  session$setInputs(`feature_general-tris_feature_general-variable` = "group")
  session$flushReact(); session$flushReact()
  # plot-type radio: warm the ignoreInit swallow, then a user edit
  session$setInputs(`feature_general-internal_radio` = "Bees")
  session$flushReact()
  session$setInputs(`feature_general-internal_radio` = "Curve")
  session$flushReact(); session$flushReact()
  # agent regression-line write (the mirrored-reactiveVal ack path)
  store_apply(w7p$child, list(regression_line = TRUE), origin = "agent")
  session$flushReact(); session$flushReact()
  w7p_vals <<- store_read(w7p$store)
  w7p_ret <<- tryCatch(shiny::isolate(w7p$ret()), error = function(e) NULL)
  v <- store_registry_view(w7p$store, "resultspace.feature_general")
  ids5 <- paste0("resultspace.feature_general.", c(
    "xax_analysis", "xax_subset", "xax_variable", "plot_type",
    "regression_line"))
  w7p_ch <<- stats::setNames(
    lapply(v[match(ids5, vapply(v, function(r) r$id, ""))],
           `[[`, "allowed_values"), ids5)
})
ok(is.null(w7p_err), "WP-7 fg: user-pick boot stays alive")
ok(ut_cmp_identical(w7p_vals, PANEL_W7_FG_GOLDEN_PICK_VALS),
   "WP-7 fg: user cascade pick + radio edit + agent regline write reproduce the frozen store values (24 keys, byte-identical)")
ok(ut_cmp_identical(w7p_ret, PANEL_W7_FG_GOLDEN_PICK_RET),
   "WP-7 fg: module return reports the settled triple, the regline mirror and the radio (analyst_feature_general payload content)")
ok(ut_cmp_identical(w7p_ch, PANEL_W7_FG_GOLDEN_PICK_CH),
   "WP-7 fg: cascade choices stay state-dependent (subset narrows to All, variable filtered to the General branch; Surv excluded)")

# ------------------- demo-data boot (the L1 data shape) ----
w7d <- new.env()
w7d_err <- NULL; w7d_view <- NULL; w7d_cold <- NULL; w7d_ret <- NULL
w7d_attr4 <- NULL
shiny::testServer(.w7_boot(
  w7d, pd = Biobase::pData(dat), ex = Biobase::exprs(dat),
  fd = fd0, feat = 1L), {
  w7d_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  v <- store_registry_view(w7d$store, "resultspace.feature_general")
  vids <- vapply(v, function(r) r$id, "")
  ids5 <- paste0("resultspace.feature_general.", c(
    "xax_analysis", "xax_subset", "xax_variable", "plot_type",
    "regression_line"))
  w7d_view <<- v[match(ids5, vids)]
  w7d_attr4 <<- v[match(c("resultspace.feature_general.attr4.color_analysis",
                           "resultspace.feature_general.attr4.search_analysis"),
                         vids)]
  w7d_cold <<- store_read(w7d$store)[ids5]
  w7d_ret <<- tryCatch(shiny::isolate(w7d$ret()), error = function(e) NULL)
})
ok(is.null(w7d_err), "WP-7 fg: demo-data boot stays alive")
ok(ut_cmp_identical(w7d_ret, PANEL_W7_FG_DEMO_RET),
   "WP-7 fg: demo module return byte-identical in a cold boot")
ok(ut_cmp_identical(w7d_view, PANEL_W7_FG_DEMO_OWN_VIEW),
   "WP-7 fg: demo own-binding registry view byte-identical (General/PCA/Feature analysis levels; Surv excluded; head-50 variable order)")
ok(ut_cmp_identical(w7d_cold, PANEL_W7_FG_DEMO_OWN_COLD),
   "WP-7 fg: demo own-binding cold values byte-identical")
ok(ut_cmp_identical(w7d_attr4, PANEL_W7_FG_DEMO_ATTR4_CH),
   "WP-7 fg: the attr4 registrations share the demo triset source through this panel's child store (color/search analysis levels = the own cascade's)")

# demo user pick: PCA | All | PC1(10.5%)
w7dp <- new.env()
w7dp_err <- NULL; w7dp_vals <- NULL; w7dp_ret <- NULL; w7dp_ch <- NULL
shiny::testServer(.w7_boot(
  w7dp, pd = Biobase::pData(dat), ex = Biobase::exprs(dat),
  fd = fd0, feat = 1L), {
  w7dp_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  session$setInputs(`feature_general-tris_feature_general-analysis` = "PCA")
  session$flushReact()
  session$setInputs(`feature_general-tris_feature_general-subset` = "All")
  session$flushReact()
  session$setInputs(`feature_general-tris_feature_general-variable` = "PC1(10.5%)")
  session$flushReact(); session$flushReact()
  w7dp_vals <<- store_read(w7dp$store)[paste0(
    "resultspace.feature_general.", c("xax_analysis", "xax_subset",
    "xax_variable", "plot_type", "regression_line"))]
  w7dp_ret <<- tryCatch(shiny::isolate(w7dp$ret()), error = function(e) NULL)
  v <- store_registry_view(w7dp$store, "resultspace.feature_general")
  vids <- vapply(v, function(r) r$id, "")
  ids5 <- paste0("resultspace.feature_general.", c(
    "xax_analysis", "xax_subset", "xax_variable", "plot_type",
    "regression_line"))
  w7dp_ch <<- stats::setNames(
    lapply(v[match(ids5, vids)], `[[`, "allowed_values"), ids5)
})
ok(is.null(w7dp_err), "WP-7 fg: demo user-pick boot stays alive")
ok(ut_cmp_identical(w7dp_vals, PANEL_W7_FG_DEMO_PICK_VALS),
   "WP-7 fg: demo PCA|All|PC1(10.5%) pick reproduces the frozen store values")
ok(ut_cmp_identical(w7dp_ret, PANEL_W7_FG_DEMO_PICK_RET),
   "WP-7 fg: demo module return reports the settled PCA triple (analyst_feature_general payload content)")
ok(ut_cmp_identical(w7dp_ch, PANEL_W7_FG_DEMO_PICK_CH),
   "WP-7 fg: demo cascade choices stay state-dependent (variable narrows to the 8 PC columns)")

# ------------------- legacy status restore paths ----
w7r <- new.env()
w7r_err <- NULL; w7r_vals <- NULL; w7r_ret <- NULL
shiny::testServer(.w7_boot(w7r, status_payload = list(
  xax = c("General", "All", "score"), showRegLine = TRUE,
  plotType = "Curve")), {
  w7r_err <<- tryCatch({
    for (i in 1:6) session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w7r_vals <<- store_read(w7r$store)
  w7r_ret <<- tryCatch(shiny::isolate(w7r$ret()), error = function(e) NULL)
})
ok(is.null(w7r_err),
   "WP-7 fg: boot with a legacy analyst_feature_general payload stays alive")
ok(ut_cmp_identical(w7r_vals, PANEL_W7_FG_GOLDEN_RESTORE_VALS),
   "WP-7 fg: the legacy restore path re-applies the payload triple + showRegLine + plotType (24 keys byte-identical)")
ok(ut_cmp_identical(w7r_ret, PANEL_W7_FG_GOLDEN_RESTORE_RET),
   "WP-7 fg: module return stays cold after a restore-only boot except the mirrored showRegLine (v1 commits through the UI, not the store)")

# a saved triple with an UNSET variable pushes the canonical unset
w7u <- new.env()
w7u_err <- NULL; w7u_vals <- NULL
shiny::testServer(.w7_boot(w7u, status_payload = list(
  xax = list("General", "All", NULL), showRegLine = TRUE,
  plotType = "Curve")), {
  w7u_err <<- tryCatch({
    for (i in 1:6) session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w7u_vals <<- store_read(w7u$store)[paste0(
    "resultspace.feature_general.", c("xax_analysis", "xax_subset",
    "xax_variable", "plot_type", "regression_line"))]
})
ok(is.null(w7u_err), "WP-7 fg: unset-variable restore boot stays alive")
ok(ut_cmp_identical(w7u_vals, PANEL_W7_FG_GOLDEN_UNSET_VALS),
   "WP-7 fg: an unset saved variable pushes the canonical unset --select-- (todo 2.3 semantics preserved)")

# demo restore triple
w7dr <- new.env()
w7dr_err <- NULL; w7dr_vals <- NULL
shiny::testServer(.w7_boot(
  w7dr, pd = Biobase::pData(dat), ex = Biobase::exprs(dat),
  fd = fd0, feat = 1L,
  status_payload = list(xax = c("General", "All", "Cell.line"),
                        showRegLine = TRUE, plotType = "Curve")), {
  w7dr_err <<- tryCatch({
    for (i in 1:6) session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w7dr_vals <<- store_read(w7dr$store)[paste0(
    "resultspace.feature_general.", c("xax_analysis", "xax_subset",
    "xax_variable", "plot_type", "regression_line"))]
})
ok(is.null(w7dr_err), "WP-7 fg: demo restore boot stays alive")
ok(ut_cmp_identical(w7dr_vals, PANEL_W7_FG_DEMO_RESTORE_VALS),
   "WP-7 fg: the demo restore triple re-applies byte-identically (General|All|Cell.line + plotType + showRegLine)")

# ------------------- the triset derivation (minus Surv) ----
ts_full <- trisetter(expr = Biobase::exprs(dat), meta = Biobase::pData(dat),
                     combine = "pheno")
ts_w7 <- ts_full[ts_full[, 1] != "Surv", , drop = FALSE]
ok(ut_cmp_identical(list(dim = dim(ts_w7),
                         analysis_levels = unique(ts_w7[, 1]),
                         non_feature_rows = ts_w7[ts_w7[, 1] != "Feature", , drop = FALSE]),
                    list(dim = PANEL_W7_FG_TRISET$dim,
                         analysis_levels = PANEL_W7_FG_TRISET$analysis_levels,
                         non_feature_rows = PANEL_W7_FG_TRISET$non_feature_rows)),
   "WP-7 fg: the cascade triset source over demo pdata = trisetter(combine=pheno) MINUS Surv (2719x3; General/PCA/Feature; the 17 pheno rows)")
ok(ut_cmp_identical(ts_full[ts_full[, 1] == "Surv", , drop = FALSE],
                    PANEL_W7_FG_TRISET$surv_rows_in_full) &&
   ut_cmp_identical(head(unique(ts_w7[, 3]), 50),
                    PANEL_W7_FG_TRISET$variable_head50),
   "WP-7 fg: the 2 Surv rows exist in the unfiltered triset and are excluded; variable head-50 order frozen")

# ------------------- retention (D3 unification) ----
# structural (deterministic, WP-0 MINOR-2 pattern): the module's ONE
# panel_keep closure holds the helper's bind-sync observer, the
# view-type carrier observer (previously the SPLIT .fg_view_keep list)
# and the three plot_type/regression sync/seed/push observers; the
# beeswarm htest observer keeps its own module-local list (out of the
# plan's unification scope). Functional (retention direction, PD-8):
# after gc() x2 a user cascade pick, a radio edit and an agent regline
# write still reach the store.
w7k <- new.env()
w7k_err <- NULL; w7k_nobs <- NULL; w7k_vals <- NULL
shiny::testServer(.w7_mkctx(w7k), {
  w7k_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  obs <- get("obs", envir = environment(w7k$ctx$keep))
  w7k_nobs <<- length(obs)
  gc(); gc()
  session$setInputs(`feature_general-tris_feature_general-analysis` = "General")
  session$flushReact()
  session$setInputs(`feature_general-tris_feature_general-subset` = "All")
  session$flushReact()
  session$setInputs(`feature_general-tris_feature_general-variable` = "score")
  session$flushReact()
  session$setInputs(`feature_general-internal_radio` = "Bees")
  session$flushReact()
  session$setInputs(`feature_general-internal_radio` = "Curve")
  session$flushReact(); session$flushReact()
  store_apply(w7k$child, list(regression_line = TRUE), origin = "agent")
  session$flushReact(); session$flushReact()
  w7k_vals <<- store_read(w7k$child)[paste0(
    "resultspace.feature_general.", c("xax_analysis", "xax_subset",
    "xax_variable", "plot_type", "regression_line"))]
})
ok(is.null(w7k_err), "WP-7 fg: retention boot stays alive")
ok(ut_cmp_identical(w7k_nobs, 6L),
   "WP-7 fg: exactly the bind-sync + view-type + radio-sync + regline-sync + seed + push observers land in ONE panel_keep list (D3 unification of .fg_keep + .fg_view_keep)")
ok(ut_cmp_identical(
  unname(unlist(w7k_vals[c("resultspace.feature_general.xax_analysis",
                    "resultspace.feature_general.xax_subset",
                    "resultspace.feature_general.xax_variable",
                    "resultspace.feature_general.plot_type",
                    "resultspace.feature_general.regression_line")])),
  c("General", "All", "score", "Curve", TRUE)),
   "WP-7 fg: the kept sync observers still fire after gc() x2 (cascade pick, radio edit and agent regline write all reach the store)")

# ------------------- the store-less legacy branch stays alive ----
w7n_err <- NULL
shiny::testServer(.w7_mklegacy(), {
  w7n_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
})
ok(is.null(w7n_err),
   "WP-7 fg: the no-store legacy branch (reactiveVal cascade path) still boots")

# ------------------- app level: the REAL L1 composition ----
w7l1 <- new.env()
w7l1$store <- widget_store_new()
app_w7l1 <- function(input, output, session) {
  w7l1$status <- omicsViewer:::L1_result_space_module(
    "rs",
    reactive_expr = reactive(Biobase::exprs(dat)),
    reactive_phenoData = reactive(Biobase::pData(dat)),
    reactive_featureData = reactive(fd0),
    status = reactive(NULL),
    store = w7l1$store)
}
w7l1_err <- NULL; w7l1_payload <- NULL
shiny::testServer(app_w7l1, {
  w7l1_err <<- tryCatch({
    for (i in 1:3) session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w7l1_payload <<- tryCatch(shiny::isolate(w7l1$status()),
                            error = function(e) conditionMessage(e))
})
ok(is.null(w7l1_err),
   "WP-7 fg: the real L1 composition boots headless post-migration (context built at the call site)")
ok(ut_cmp_identical(names(w7l1_payload), PANEL_W7_FG_GOLDEN_L1_STATUS$names),
   "WP-7 fg: L1 status payload keeps its exact key order (snapshot schema untouched)")
ok(ut_cmp_identical(w7l1_payload$analyst_feature_general,
                    PANEL_W7_FG_GOLDEN_L1_STATUS$fg),
   "WP-7 fg: analyst_feature_general payload byte-identical pre/post migration (cold boot: xax/plotType NULL, showRegLine FALSE, attr4 nine-field NULL shape)")
ok(ut_cmp_identical(
  sort(intersect(vapply(store_registry_view(w7l1$store), function(r) r$id, ""),
                 PANEL_W7_FG_IDS)),
  PANEL_W7_FG_IDS),
   "WP-7 fg: the real L1 boot carries all 24 feature_general keys (the D12 freeze above pins the full 72-id vector)")

## --------------- WP-8: Sample-general migration (D3/D4/D13) ----
# The second-largest pure migration (the WP-5/WP-6/WP-7 pattern):
# sample-general composes attr4selector, beeswarm, contingency-table and
# SURVIVAL submodules, and its store surface is 27 keys (3 own cascade +
# 19 attr4 + 1 survival_censor + 4 batch_* registered through this
# panel's child store). SUBMODULES ARE UNTOUCHED (the plan's WP-8 scope
# is panel glue only; attr4/plotly/contingency/survival/batch internals
# are front-run boundaries of later work, not this WP). The migration
# unifies retention (.sg_keep + .sg_view_keep -> ONE panel_keep
# instance, D3), moves the cascade boilerplate into panel_triset_bind
# (D4: triset_source = trisetter(combine="pheno") KEEPING the Surv rows
# - Kaplan-Meier, the deliberate opposite of feature-general's
# minus-Surv cascade; allow_unset = TRUE) and threads the context
# (D13). The .sg_apply_triple dual-path helper (store / legacy xax)
# stays verbatim in the module (its internal legacy fallback, not L1
# glue).
#
# Golden-fixture provenance (PD-29/PD-33/PD-38 school): all goldens
# below were captured from PRE-MIGRATION HEAD 74924e3 BEFORE any edit
# (before the WP-7 NIT cleanup commit dd69868, for strict provenance),
# by booting the then-current sample_general_module exactly like the L1
# call site / the agentWidgets S4 shape and dput()ing the outputs:
# the FULL 27-binding registry view + raw meta (incl. agent_writable) +
# cold store values + module return, user-pick/restore triples over two
# datasets (synthetic + demo.RDS), the real-L1 analyst_sample_general
# payload + status key order, and the KEEPING-Surv triset derivation
# over demo pdata. The IDENTICAL capture script was re-run
# post-migration: all 21 outputs byte-identical (the PD-33 pre/post
# proof). The fixtures are FROZEN: never regenerate them from the
# migrated code (that would make the equivalence claim circular).
sample_general_module <- omicsViewer:::sample_general_module

# -------------------------------------------------- frozen goldens ----
# the 27 ids (3 own cascade + 19 attr4 + 1 survival + 4 batch, all
# under the panel's child store)

PANEL_W8_SG_IDS <- sort(c(
  paste0("resultspace.sample_general.", c(
    "xax_analysis", "xax_subset", "xax_variable", "survival_censor",
    "batch_show_phenotype", "batch_show_features",
    "batch_phenotype_selected_row", "batch_feature_selected_row")),
  paste0("resultspace.sample_general.attr4.", c(
    outer(c("color", "shape", "size", "tooltip", "search"),
          c("_analysis", "_subset", "_variable"), FUN = paste0),
    "search_value", "xcut", "ycut", "scorner"))))


PANEL_W8_SG_GOLDEN_VIEW <- list(list(id = "resultspace.sample_general.xax_analysis", kind = "select", 
    label = "Link variable category", help = "Annotation category of the variable the selected samples are analyzed against", 
    depends_on = character(0), min = NULL, max = NULL, allowed_values = c("General", 
    "Surv", "Feature")), list(id = "resultspace.sample_general.xax_subset", 
    kind = "select", label = "Link variable subcategory", help = "Subcategory within the link-variable category", 
    depends_on = "resultspace.sample_general.xax_analysis", min = NULL, 
    max = NULL, allowed_values = character(0)), list(id = "resultspace.sample_general.xax_variable", 
    kind = "select_cascaded", label = "Link variable", help = "Variable driving the analysis view: numeric gives a beeswarm, categorical a contingency table, Surv a Kaplan-Meier curve", 
    depends_on = c("resultspace.sample_general.xax_analysis", 
    "resultspace.sample_general.xax_subset"), min = NULL, max = NULL, 
    allowed_values = c("group", "score", "time", "F1", "F2", 
    "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12", 
    "F13", "F14", "F15", "F16")), list(id = "resultspace.sample_general.attr4.color_analysis", 
    kind = "select", label = "Color attribute category", help = "First level of the color attribute cascade (annotation category) to color the figure's points by a group or annotation column", 
    depends_on = character(0), min = NULL, max = NULL, allowed_values = c("General", 
    "Surv", "Feature")), list(id = "resultspace.sample_general.attr4.color_subset", 
    kind = "select", label = "Color attribute subcategory", help = "Second level of the color attribute cascade to color the figure's points by a group or annotation column", 
    depends_on = "resultspace.sample_general.attr4.color_analysis", 
    min = NULL, max = NULL, allowed_values = character(0)), list(
    id = "resultspace.sample_general.attr4.color_variable", kind = "select_cascaded", 
    label = "Color attribute variable", help = "Annotation variable to color the figure's points by a group or annotation column", 
    depends_on = c("resultspace.sample_general.attr4.color_analysis", 
    "resultspace.sample_general.attr4.color_subset"), min = NULL, 
    max = NULL, allowed_values = c("group", "score", "time", 
    "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", 
    "F11", "F12", "F13", "F14", "F15", "F16")), list(id = "resultspace.sample_general.attr4.shape_analysis", 
    kind = "select", label = "Shape attribute category", help = "First level of the shape attribute cascade (annotation category) to give the figure's points a different shape per group", 
    depends_on = character(0), min = NULL, max = NULL, allowed_values = c("General", 
    "Surv", "Feature")), list(id = "resultspace.sample_general.attr4.shape_subset", 
    kind = "select", label = "Shape attribute subcategory", help = "Second level of the shape attribute cascade to give the figure's points a different shape per group", 
    depends_on = "resultspace.sample_general.attr4.shape_analysis", 
    min = NULL, max = NULL, allowed_values = character(0)), list(
    id = "resultspace.sample_general.attr4.shape_variable", kind = "select_cascaded", 
    label = "Shape attribute variable", help = "Annotation variable to give the figure's points a different shape per group", 
    depends_on = c("resultspace.sample_general.attr4.shape_analysis", 
    "resultspace.sample_general.attr4.shape_subset"), min = NULL, 
    max = NULL, allowed_values = c("group", "score", "time", 
    "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", 
    "F11", "F12", "F13", "F14", "F15", "F16")), list(id = "resultspace.sample_general.attr4.size_analysis", 
    kind = "select", label = "Size attribute category", help = "First level of the size attribute cascade (annotation category) to scale the figure's point size by a numeric column", 
    depends_on = character(0), min = NULL, max = NULL, allowed_values = c("General", 
    "Surv", "Feature")), list(id = "resultspace.sample_general.attr4.size_subset", 
    kind = "select", label = "Size attribute subcategory", help = "Second level of the size attribute cascade to scale the figure's point size by a numeric column", 
    depends_on = "resultspace.sample_general.attr4.size_analysis", 
    min = NULL, max = NULL, allowed_values = character(0)), list(
    id = "resultspace.sample_general.attr4.size_variable", kind = "select_cascaded", 
    label = "Size attribute variable", help = "Annotation variable to scale the figure's point size by a numeric column", 
    depends_on = c("resultspace.sample_general.attr4.size_analysis", 
    "resultspace.sample_general.attr4.size_subset"), min = NULL, 
    max = NULL, allowed_values = c("group", "score", "time", 
    "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", 
    "F11", "F12", "F13", "F14", "F15", "F16")), list(id = "resultspace.sample_general.attr4.tooltip_analysis", 
    kind = "select", label = "Tooltips attribute category", help = "First level of the tooltips attribute cascade (annotation category) shown in the figure's hover tooltip", 
    depends_on = character(0), min = NULL, max = NULL, allowed_values = c("General", 
    "Surv", "Feature")), list(id = "resultspace.sample_general.attr4.tooltip_subset", 
    kind = "select", label = "Tooltips attribute subcategory", 
    help = "Second level of the tooltips attribute cascade shown in the figure's hover tooltip", 
    depends_on = "resultspace.sample_general.attr4.tooltip_analysis", 
    min = NULL, max = NULL, allowed_values = character(0)), list(
    id = "resultspace.sample_general.attr4.tooltip_variable", 
    kind = "select_cascaded", label = "Tooltips attribute variable", 
    help = "Annotation variable shown in the figure's hover tooltip", 
    depends_on = c("resultspace.sample_general.attr4.tooltip_analysis", 
    "resultspace.sample_general.attr4.tooltip_subset"), min = NULL, 
    max = NULL, allowed_values = c("group", "score", "time", 
    "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", 
    "F11", "F12", "F13", "F14", "F15", "F16")), list(id = "resultspace.sample_general.attr4.search_analysis", 
    kind = "select", label = "Search attribute category", help = "First level of the search attribute cascade (annotation category) whose values the gene highlight picker circles as open-ring markers on the figure", 
    depends_on = character(0), min = NULL, max = NULL, allowed_values = c("General", 
    "Surv", "Feature")), list(id = "resultspace.sample_general.attr4.search_subset", 
    kind = "select", label = "Search attribute subcategory", 
    help = "Second level of the search attribute cascade whose values the gene highlight picker circles as open-ring markers on the figure", 
    depends_on = "resultspace.sample_general.attr4.search_analysis", 
    min = NULL, max = NULL, allowed_values = character(0)), list(
    id = "resultspace.sample_general.attr4.search_variable", 
    kind = "select_cascaded", label = "Search attribute variable", 
    help = "Annotation variable whose values the gene highlight picker circles as open-ring markers on the figure", 
    depends_on = c("resultspace.sample_general.attr4.search_analysis", 
    "resultspace.sample_general.attr4.search_subset"), min = NULL, 
    max = NULL, allowed_values = c("group", "score", "time", 
    "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", 
    "F11", "F12", "F13", "F14", "F15", "F16")), list(id = "resultspace.sample_general.attr4.search_value", 
    kind = "multi_select", label = "Search highlight values", 
    help = "Gene name or other value of the selected Search column to draw as open-ring highlights on the figure; set the Search column first; empty clears the highlight", 
    depends_on = "resultspace.sample_general.attr4.search_variable", 
    min = NULL, max = NULL, allowed_values = character(0)), list(
    id = "resultspace.sample_general.attr4.xcut", kind = "string", 
    label = "X cutoff", help = "X-axis cutoff for area selection; a number or an R expression such as -log10(0.05)", 
    depends_on = character(0), min = NULL, max = NULL, allowed_values = NULL), 
    list(id = "resultspace.sample_general.attr4.ycut", kind = "string", 
        label = "Y cutoff", help = "Y-axis cutoff for area selection; a number or an R expression such as log10(2)", 
        depends_on = character(0), min = NULL, max = NULL, allowed_values = NULL), 
    list(id = "resultspace.sample_general.attr4.scorner", kind = "select", 
        label = "Selected area", help = "Plot region selected by the cutoffs (corner or volcano mode)", 
        depends_on = character(0), min = NULL, max = NULL, allowed_values = "None"), 
    list(id = "resultspace.sample_general.survival_censor", kind = "numeric", 
        label = "Censor time", help = "Time at which the Kaplan-Meier analysis right-censors observations. PRECONDITION: requires a survival-type (Surv) sample variable selected as the x-axis variable - the censor slider exists only in the survival view; a write that lands before the view renders is applied when it does", 
        depends_on = character(0), min = NULL, max = NULL, allowed_values = NULL), 
    list(id = "resultspace.sample_general.batch_show_phenotype", 
        kind = "boolean", label = "Compare phenotype variables", 
        help = "Test all phenotype variables between selected and unselected samples", 
        depends_on = character(0), min = NULL, max = NULL, allowed_values = NULL), 
    list(id = "resultspace.sample_general.batch_show_features", 
        kind = "boolean", label = "Compare features", help = "Test all molecular features between selected and unselected samples", 
        depends_on = character(0), min = NULL, max = NULL, allowed_values = NULL), 
    list(id = "resultspace.sample_general.batch_phenotype_selected_row", 
        kind = "select", label = "Selected phenotype variable", 
        help = "Phenotype-comparison row selected in the batch results table; drives the sample link. PRECONDITION: requires the Compare-phenotype toggle on, a sample selection defining two groups, and the batch panel visible - the comparison computes only then, and the row choices are its computed rows", 
        depends_on = character(0), min = NULL, max = NULL, allowed_values = character(0)), 
    list(id = "resultspace.sample_general.batch_feature_selected_row", 
        kind = "select", label = "Selected feature", help = "Feature-comparison row selected in the batch results table; drives the sample link. PRECONDITION: requires the Compare-features toggle on, a sample selection defining two groups, and the batch panel visible - the comparison computes only then, and the row choices are its computed rows", 
        depends_on = character(0), min = NULL, max = NULL, allowed_values = character(0)))
PANEL_W8_SG_GOLDEN_META <- list(list(id = "resultspace.sample_general.xax_analysis", kind = "select", 
    label = "Link variable category", help = "Annotation category of the variable the selected samples are analyzed against", 
    depends_on = character(0), agent_writable = TRUE), list(id = "resultspace.sample_general.xax_subset", 
    kind = "select", label = "Link variable subcategory", help = "Subcategory within the link-variable category", 
    depends_on = "resultspace.sample_general.xax_analysis", agent_writable = TRUE), 
    list(id = "resultspace.sample_general.xax_variable", kind = "select_cascaded", 
        label = "Link variable", help = "Variable driving the analysis view: numeric gives a beeswarm, categorical a contingency table, Surv a Kaplan-Meier curve", 
        depends_on = c("resultspace.sample_general.xax_analysis", 
        "resultspace.sample_general.xax_subset"), agent_writable = TRUE), 
    list(id = "resultspace.sample_general.attr4.color_analysis", 
        kind = "select", label = "Color attribute category", 
        help = "First level of the color attribute cascade (annotation category) to color the figure's points by a group or annotation column", 
        depends_on = character(0), agent_writable = TRUE), list(
        id = "resultspace.sample_general.attr4.color_subset", 
        kind = "select", label = "Color attribute subcategory", 
        help = "Second level of the color attribute cascade to color the figure's points by a group or annotation column", 
        depends_on = "resultspace.sample_general.attr4.color_analysis", 
        agent_writable = TRUE), list(id = "resultspace.sample_general.attr4.color_variable", 
        kind = "select_cascaded", label = "Color attribute variable", 
        help = "Annotation variable to color the figure's points by a group or annotation column", 
        depends_on = c("resultspace.sample_general.attr4.color_analysis", 
        "resultspace.sample_general.attr4.color_subset"), agent_writable = TRUE), 
    list(id = "resultspace.sample_general.attr4.shape_analysis", 
        kind = "select", label = "Shape attribute category", 
        help = "First level of the shape attribute cascade (annotation category) to give the figure's points a different shape per group", 
        depends_on = character(0), agent_writable = TRUE), list(
        id = "resultspace.sample_general.attr4.shape_subset", 
        kind = "select", label = "Shape attribute subcategory", 
        help = "Second level of the shape attribute cascade to give the figure's points a different shape per group", 
        depends_on = "resultspace.sample_general.attr4.shape_analysis", 
        agent_writable = TRUE), list(id = "resultspace.sample_general.attr4.shape_variable", 
        kind = "select_cascaded", label = "Shape attribute variable", 
        help = "Annotation variable to give the figure's points a different shape per group", 
        depends_on = c("resultspace.sample_general.attr4.shape_analysis", 
        "resultspace.sample_general.attr4.shape_subset"), agent_writable = TRUE), 
    list(id = "resultspace.sample_general.attr4.size_analysis", 
        kind = "select", label = "Size attribute category", help = "First level of the size attribute cascade (annotation category) to scale the figure's point size by a numeric column", 
        depends_on = character(0), agent_writable = TRUE), list(
        id = "resultspace.sample_general.attr4.size_subset", 
        kind = "select", label = "Size attribute subcategory", 
        help = "Second level of the size attribute cascade to scale the figure's point size by a numeric column", 
        depends_on = "resultspace.sample_general.attr4.size_analysis", 
        agent_writable = TRUE), list(id = "resultspace.sample_general.attr4.size_variable", 
        kind = "select_cascaded", label = "Size attribute variable", 
        help = "Annotation variable to scale the figure's point size by a numeric column", 
        depends_on = c("resultspace.sample_general.attr4.size_analysis", 
        "resultspace.sample_general.attr4.size_subset"), agent_writable = TRUE), 
    list(id = "resultspace.sample_general.attr4.tooltip_analysis", 
        kind = "select", label = "Tooltips attribute category", 
        help = "First level of the tooltips attribute cascade (annotation category) shown in the figure's hover tooltip", 
        depends_on = character(0), agent_writable = TRUE), list(
        id = "resultspace.sample_general.attr4.tooltip_subset", 
        kind = "select", label = "Tooltips attribute subcategory", 
        help = "Second level of the tooltips attribute cascade shown in the figure's hover tooltip", 
        depends_on = "resultspace.sample_general.attr4.tooltip_analysis", 
        agent_writable = TRUE), list(id = "resultspace.sample_general.attr4.tooltip_variable", 
        kind = "select_cascaded", label = "Tooltips attribute variable", 
        help = "Annotation variable shown in the figure's hover tooltip", 
        depends_on = c("resultspace.sample_general.attr4.tooltip_analysis", 
        "resultspace.sample_general.attr4.tooltip_subset"), agent_writable = TRUE), 
    list(id = "resultspace.sample_general.attr4.search_analysis", 
        kind = "select", label = "Search attribute category", 
        help = "First level of the search attribute cascade (annotation category) whose values the gene highlight picker circles as open-ring markers on the figure", 
        depends_on = character(0), agent_writable = TRUE), list(
        id = "resultspace.sample_general.attr4.search_subset", 
        kind = "select", label = "Search attribute subcategory", 
        help = "Second level of the search attribute cascade whose values the gene highlight picker circles as open-ring markers on the figure", 
        depends_on = "resultspace.sample_general.attr4.search_analysis", 
        agent_writable = TRUE), list(id = "resultspace.sample_general.attr4.search_variable", 
        kind = "select_cascaded", label = "Search attribute variable", 
        help = "Annotation variable whose values the gene highlight picker circles as open-ring markers on the figure", 
        depends_on = c("resultspace.sample_general.attr4.search_analysis", 
        "resultspace.sample_general.attr4.search_subset"), agent_writable = TRUE), 
    list(id = "resultspace.sample_general.attr4.search_value", 
        kind = "multi_select", label = "Search highlight values", 
        help = "Gene name or other value of the selected Search column to draw as open-ring highlights on the figure; set the Search column first; empty clears the highlight", 
        depends_on = "resultspace.sample_general.attr4.search_variable", 
        agent_writable = TRUE), list(id = "resultspace.sample_general.attr4.xcut", 
        kind = "string", label = "X cutoff", help = "X-axis cutoff for area selection; a number or an R expression such as -log10(0.05)", 
        depends_on = character(0), agent_writable = TRUE), list(
        id = "resultspace.sample_general.attr4.ycut", kind = "string", 
        label = "Y cutoff", help = "Y-axis cutoff for area selection; a number or an R expression such as log10(2)", 
        depends_on = character(0), agent_writable = TRUE), list(
        id = "resultspace.sample_general.attr4.scorner", kind = "select", 
        label = "Selected area", help = "Plot region selected by the cutoffs (corner or volcano mode)", 
        depends_on = character(0), agent_writable = TRUE), list(
        id = "resultspace.sample_general.survival_censor", kind = "numeric", 
        label = "Censor time", help = "Time at which the Kaplan-Meier analysis right-censors observations. PRECONDITION: requires a survival-type (Surv) sample variable selected as the x-axis variable - the censor slider exists only in the survival view; a write that lands before the view renders is applied when it does", 
        depends_on = character(0), agent_writable = TRUE), list(
        id = "resultspace.sample_general.batch_show_phenotype", 
        kind = "boolean", label = "Compare phenotype variables", 
        help = "Test all phenotype variables between selected and unselected samples", 
        depends_on = character(0), agent_writable = TRUE), list(
        id = "resultspace.sample_general.batch_show_features", 
        kind = "boolean", label = "Compare features", help = "Test all molecular features between selected and unselected samples", 
        depends_on = character(0), agent_writable = TRUE), list(
        id = "resultspace.sample_general.batch_phenotype_selected_row", 
        kind = "select", label = "Selected phenotype variable", 
        help = "Phenotype-comparison row selected in the batch results table; drives the sample link. PRECONDITION: requires the Compare-phenotype toggle on, a sample selection defining two groups, and the batch panel visible - the comparison computes only then, and the row choices are its computed rows", 
        depends_on = character(0), agent_writable = TRUE), list(
        id = "resultspace.sample_general.batch_feature_selected_row", 
        kind = "select", label = "Selected feature", help = "Feature-comparison row selected in the batch results table; drives the sample link. PRECONDITION: requires the Compare-features toggle on, a sample selection defining two groups, and the batch panel visible - the comparison computes only then, and the row choices are its computed rows", 
        depends_on = character(0), agent_writable = TRUE))
PANEL_W8_SG_GOLDEN_COLD <- list(resultspace.sample_general.xax_analysis = NULL, resultspace.sample_general.xax_subset = NULL, 
    resultspace.sample_general.xax_variable = NULL, resultspace.sample_general.attr4.color_analysis = NULL, 
    resultspace.sample_general.attr4.color_subset = NULL, resultspace.sample_general.attr4.color_variable = NULL, 
    resultspace.sample_general.attr4.shape_analysis = NULL, resultspace.sample_general.attr4.shape_subset = NULL, 
    resultspace.sample_general.attr4.shape_variable = NULL, resultspace.sample_general.attr4.size_analysis = NULL, 
    resultspace.sample_general.attr4.size_subset = NULL, resultspace.sample_general.attr4.size_variable = NULL, 
    resultspace.sample_general.attr4.tooltip_analysis = NULL, 
    resultspace.sample_general.attr4.tooltip_subset = NULL, resultspace.sample_general.attr4.tooltip_variable = NULL, 
    resultspace.sample_general.attr4.search_analysis = NULL, 
    resultspace.sample_general.attr4.search_subset = NULL, resultspace.sample_general.attr4.search_variable = NULL, 
    resultspace.sample_general.attr4.search_value = character(0), 
    resultspace.sample_general.attr4.xcut = NULL, resultspace.sample_general.attr4.ycut = NULL, 
    resultspace.sample_general.attr4.scorner = "None", resultspace.sample_general.survival_censor = NULL, 
    resultspace.sample_general.batch_show_phenotype = NULL, resultspace.sample_general.batch_show_features = NULL, 
    resultspace.sample_general.batch_phenotype_selected_row = NULL, 
    resultspace.sample_general.batch_feature_selected_row = NULL)
PANEL_W8_SG_GOLDEN_RET <- list(xax = NULL, attr4 = list(selectColor = NULL, selectShape = NULL, 
    selectSize = NULL, selectTooltip = NULL, searchOnCol = NULL, 
    searchValue = NULL, xcut = NULL, ycut = NULL, acorner = NULL))
PANEL_W8_SG_GOLDEN_PICK_VALS <- list(resultspace.sample_general.xax_analysis = "Surv", resultspace.sample_general.xax_subset = "All", 
    resultspace.sample_general.xax_variable = "time", resultspace.sample_general.attr4.color_analysis = NULL, 
    resultspace.sample_general.attr4.color_subset = NULL, resultspace.sample_general.attr4.color_variable = NULL, 
    resultspace.sample_general.attr4.shape_analysis = NULL, resultspace.sample_general.attr4.shape_subset = NULL, 
    resultspace.sample_general.attr4.shape_variable = NULL, resultspace.sample_general.attr4.size_analysis = NULL, 
    resultspace.sample_general.attr4.size_subset = NULL, resultspace.sample_general.attr4.size_variable = NULL, 
    resultspace.sample_general.attr4.tooltip_analysis = NULL, 
    resultspace.sample_general.attr4.tooltip_subset = NULL, resultspace.sample_general.attr4.tooltip_variable = NULL, 
    resultspace.sample_general.attr4.search_analysis = NULL, 
    resultspace.sample_general.attr4.search_subset = NULL, resultspace.sample_general.attr4.search_variable = NULL, 
    resultspace.sample_general.attr4.search_value = character(0), 
    resultspace.sample_general.attr4.xcut = NULL, resultspace.sample_general.attr4.ycut = NULL, 
    resultspace.sample_general.attr4.scorner = "None", resultspace.sample_general.survival_censor = NULL, 
    resultspace.sample_general.batch_show_phenotype = NULL, resultspace.sample_general.batch_show_features = NULL, 
    resultspace.sample_general.batch_phenotype_selected_row = NULL, 
    resultspace.sample_general.batch_feature_selected_row = NULL)
PANEL_W8_SG_GOLDEN_PICK_RET <- list(xax = list(analysis = "Surv", subset = "All", variable = "time"), 
    attr4 = list(selectColor = NULL, selectShape = NULL, selectSize = NULL, 
        selectTooltip = NULL, searchOnCol = NULL, searchValue = NULL, 
        xcut = NULL, ycut = NULL, acorner = NULL))
PANEL_W8_SG_GOLDEN_PICK_CH <- list(resultspace.sample_general.xax_analysis = c("General", "Surv", 
"Feature"), resultspace.sample_general.xax_subset = "All", resultspace.sample_general.xax_variable = "time")
PANEL_W8_SG_GOLDEN_RESTORE_VALS <- list(resultspace.sample_general.xax_analysis = "Surv", resultspace.sample_general.xax_subset = "All", 
    resultspace.sample_general.xax_variable = "time", resultspace.sample_general.attr4.color_analysis = NULL, 
    resultspace.sample_general.attr4.color_subset = NULL, resultspace.sample_general.attr4.color_variable = NULL, 
    resultspace.sample_general.attr4.shape_analysis = NULL, resultspace.sample_general.attr4.shape_subset = NULL, 
    resultspace.sample_general.attr4.shape_variable = NULL, resultspace.sample_general.attr4.size_analysis = NULL, 
    resultspace.sample_general.attr4.size_subset = NULL, resultspace.sample_general.attr4.size_variable = NULL, 
    resultspace.sample_general.attr4.tooltip_analysis = NULL, 
    resultspace.sample_general.attr4.tooltip_subset = NULL, resultspace.sample_general.attr4.tooltip_variable = NULL, 
    resultspace.sample_general.attr4.search_analysis = NULL, 
    resultspace.sample_general.attr4.search_subset = NULL, resultspace.sample_general.attr4.search_variable = NULL, 
    resultspace.sample_general.attr4.search_value = character(0), 
    resultspace.sample_general.attr4.xcut = NULL, resultspace.sample_general.attr4.ycut = NULL, 
    resultspace.sample_general.attr4.scorner = "None", resultspace.sample_general.survival_censor = NULL, 
    resultspace.sample_general.batch_show_phenotype = NULL, resultspace.sample_general.batch_show_features = NULL, 
    resultspace.sample_general.batch_phenotype_selected_row = NULL, 
    resultspace.sample_general.batch_feature_selected_row = NULL)
PANEL_W8_SG_GOLDEN_RESTORE_RET <- list(xax = NULL, attr4 = list(selectColor = NULL, selectShape = NULL, 
    selectSize = NULL, selectTooltip = NULL, searchOnCol = NULL, 
    searchValue = NULL, xcut = NULL, ycut = NULL, acorner = NULL))
PANEL_W8_SG_GOLDEN_UNSET_VALS <- list(resultspace.sample_general.xax_analysis = "General", resultspace.sample_general.xax_subset = "All", 
    resultspace.sample_general.xax_variable = "--select--")
PANEL_W8_SG_DEMO_OWN_VIEW <- list(list(id = "resultspace.sample_general.xax_analysis", kind = "select", 
    label = "Link variable category", help = "Annotation category of the variable the selected samples are analyzed against", 
    depends_on = character(0), min = NULL, max = NULL, allowed_values = c("General", 
    "PCA", "Surv", "Feature")), list(id = "resultspace.sample_general.xax_subset", 
    kind = "select", label = "Link variable subcategory", help = "Subcategory within the link-variable category", 
    depends_on = "resultspace.sample_general.xax_analysis", min = NULL, 
    max = NULL, allowed_values = character(0)), list(id = "resultspace.sample_general.xax_variable", 
    kind = "select_cascaded", label = "Link variable", help = "Variable driving the analysis view: numeric gives a beeswarm, categorical a contingency table, Surv a Kaplan-Meier curve", 
    depends_on = c("resultspace.sample_general.xax_analysis", 
    "resultspace.sample_general.xax_subset"), min = NULL, max = NULL, 
    allowed_values = c("Cell.line", "MDR", "Doubleing.Time", 
    "Origin", "TP53.Status", "Gender", "numberOfFeatures", "PC1(10.5%)", 
    "PC2(7.2%)", "PC3(5.4%)", "PC4(4.8%)", "PC5(3.6%)", "PC6(2.9%)", 
    "PC7(2.8%)", "PC8(2.5%)", "OS", "DFS", "Dose", "Curve", "X32.TYW5", 
    "X52.CNOT1", "X55.MAP1LC3B.MAP1LC3B2", "X73.METTL15", "X107.ESPN", 
    "X120.WASH3P.WASH2P.WASH1", "X266.KIF2A", "X267.KIF2A", "X269.DDX39A", 
    "X273.ACOT7", "X274.MYO1C", "X276.MYO1F", "X279.HAX1", "X283.STXBP3", 
    "X288.AP3B1", "X290.APBB1", "X291.LGALS8", "X294.NFKBIE", 
    "X295.PSMD11", "X296.PSMD12", "X297.PSMD9", "X298.ATOX1", 
    "X299.MEN1", "X302.SUPT5H", "X308.HIP1", "X311.WWP2", "X312.CDC7", 
    "X313.PIK3CD", "X318.DCTN6", "X320.WASL", "X321.IPO5")))
PANEL_W8_SG_DEMO_OWN_COLD <- list(resultspace.sample_general.xax_analysis = NULL, resultspace.sample_general.xax_subset = NULL, 
    resultspace.sample_general.xax_variable = NULL)
PANEL_W8_SG_DEMO_RET <- list(xax = NULL, attr4 = list(selectColor = NULL, selectShape = NULL, 
    selectSize = NULL, selectTooltip = NULL, searchOnCol = NULL, 
    searchValue = NULL, xcut = NULL, ycut = NULL, acorner = NULL))
PANEL_W8_SG_DEMO_PICK_VALS <- list(resultspace.sample_general.xax_analysis = "Surv", resultspace.sample_general.xax_subset = "all", 
    resultspace.sample_general.xax_variable = "OS")
PANEL_W8_SG_DEMO_PICK_RET <- list(xax = list(analysis = "Surv", subset = "all", variable = "OS"), 
    attr4 = list(selectColor = NULL, selectShape = NULL, selectSize = NULL, 
        selectTooltip = NULL, searchOnCol = NULL, searchValue = NULL, 
        xcut = NULL, ycut = NULL, acorner = NULL))
PANEL_W8_SG_DEMO_PICK_CH <- list(resultspace.sample_general.xax_analysis = c("General", "PCA", 
"Surv", "Feature"), resultspace.sample_general.xax_subset = "all", 
    resultspace.sample_general.xax_variable = c("OS", "DFS"))
PANEL_W8_SG_DEMO_RESTORE_VALS <- list(resultspace.sample_general.xax_analysis = "General", resultspace.sample_general.xax_subset = "All", 
    resultspace.sample_general.xax_variable = "Cell.line")
PANEL_W8_SG_TRISET <- list(dim = c(2721L, 3L), analysis_levels = c("General", "PCA", 
"Surv", "Feature"), surv_rows = structure(c("Surv", "Surv", "all", 
"all", "OS", "DFS"), dim = 2:3), non_feature_pca_rows = structure(c("General", 
"General", "General", "General", "General", "General", "General", 
"Surv", "Surv", "General", "General", "All", "All", "All", "All", 
"All", "All", "All", "all", "all", "All", "All", "Cell.line", 
"MDR", "Doubleing.Time", "Origin", "TP53.Status", "Gender", "numberOfFeatures", 
"OS", "DFS", "Dose", "Curve"), dim = c(11L, 3L)), variable_head50 = c("Cell.line", 
"MDR", "Doubleing.Time", "Origin", "TP53.Status", "Gender", "numberOfFeatures", 
"PC1(10.5%)", "PC2(7.2%)", "PC3(5.4%)", "PC4(4.8%)", "PC5(3.6%)", 
"PC6(2.9%)", "PC7(2.8%)", "PC8(2.5%)", "OS", "DFS", "Dose", "Curve", 
"X32.TYW5", "X52.CNOT1", "X55.MAP1LC3B.MAP1LC3B2", "X73.METTL15", 
"X107.ESPN", "X120.WASH3P.WASH2P.WASH1", "X266.KIF2A", "X267.KIF2A", 
"X269.DDX39A", "X273.ACOT7", "X274.MYO1C", "X276.MYO1F", "X279.HAX1", 
"X283.STXBP3", "X288.AP3B1", "X290.APBB1", "X291.LGALS8", "X294.NFKBIE", 
"X295.PSMD11", "X296.PSMD12", "X297.PSMD9", "X298.ATOX1", "X299.MEN1", 
"X302.SUPT5H", "X308.HIP1", "X311.WWP2", "X312.CDC7", "X313.PIK3CD", 
"X318.DCTN6", "X320.WASL", "X321.IPO5"))
PANEL_W8_SG_GOLDEN_L1_NAMES <- c("analyst_active_tab", "analyst_feature_general", "analyst_sample_general", 
"analyst_gene_shot", "analyst_fgsea", "analyst_ora", "analyst_ptm", 
"analyst_dose_response", "analyst_additional")
PANEL_W8_SG_GOLDEN_L1_SG <- list(xax = NULL, attr4 = list(selectColor = NULL, selectShape = NULL, 
    selectSize = NULL, selectTooltip = NULL, searchOnCol = NULL, 
    searchValue = NULL, xcut = NULL, ycut = NULL, acorner = NULL))


# -------------------------------------------------- boot factories ----
# synthetic data: the agentWidgets S4 sample shape (group/score pheno
# columns + a Surv column the cascade must KEEP); seeded for determinism
# (values never enter the store surface, but pin anyway)
set.seed(7)
.w8_pd <- data.frame(
  `General|All|group` = factor(rep(c("A", "B"), each = 4)),
  `General|All|score` = rnorm(8),
  `Surv|All|time` = rexp(8),
  row.names = paste0("S", 1:8), check.names = FALSE)
.w8_ex <- matrix(rnorm(16 * 8), nrow = 16,
                 dimnames = list(paste0("F", 1:16), paste0("S", 1:8)))
# context-built (the WP-8 L1 call-site shape; holder receives the
# context so the keep closure is inspectable) and legacy direct
# (formals only - the pre-WP-8 shape, still supported)
.w8_mkctx <- function(holder, pd = .w8_pd, ex = .w8_ex,
                      j = c("S1", "S3"), status_payload = NULL) {
  force(holder)
  function(input, output, session) {
    .pd <- reactive(pd); .ex <- reactive(ex); .jj <- reactive(j)
    .st <- reactive(status_payload)
    holder$store <- widget_store_new()
    holder$child <- widget_store_child(holder$store,
                                       "resultspace.sample_general")
    holder$ctx <- panel_context(
      id = "sample_general", tab_value = "Sample", store = holder$child,
      data = list(expr = .ex, pdata = .pd),
      selection = list(sample = .jj),
      status = .st, session = session)
    holder$ret <- sample_general_module(
      "sg", reactive_phenoData = .pd, reactive_expr = .ex,
      reactive_j = .jj, reactive_status = .st,
      store = holder$child, context = holder$ctx)
  }
}
# no-context legacy boot (store formal only - the pre-WP-8 shape)
.w8_boot <- function(env, pd = .w8_pd, ex = .w8_ex,
                     j = c("S1", "S3"), status_payload = NULL) {
  force(env)
  function(input, output, session) {
    env$store <- widget_store_new()
    env$child <- widget_store_child(env$store,
                                    "resultspace.sample_general")
    env$ret <- sample_general_module(
      "sg", reactive_phenoData = reactive(pd), reactive_expr = reactive(ex),
      reactive_j = reactive(j),
      reactive_status = reactive(status_payload),
      store = env$child)
  }
}
# store-less legacy boot (no store, no context - the reactiveVal path)
.w8_mklegacy <- function(pd = .w8_pd, ex = .w8_ex, j = c("S1", "S3"),
                         status_payload = NULL) {
  function(input, output, session) {
    sample_general_module(
      "sg", reactive_phenoData = reactive(pd), reactive_expr = reactive(ex),
      reactive_j = reactive(j), reactive_status = reactive(status_payload))
  }
}
.w8_meta <- function(store) {
  lapply(grep("^resultspace\\.sample_general(\\.|$)", names(store$bindings),
              value = TRUE),
         function(k) {
           b <- store$bindings[[k]]
           list(id = k, kind = b$kind, label = b$label, help = b$help,
                depends_on = b$depends_on, agent_writable = b$agent_writable)
         })
}
.w8_TRIO <- paste0("resultspace.sample_general.",
                   c("xax_analysis", "xax_subset", "xax_variable"))

# ------------------- registration: context boot = the frozen golden ----
w8 <- new.env()
w8_err <- NULL; w8_view <- NULL; w8_meta <- NULL
shiny::testServer(.w8_mkctx(w8), {
  w8_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  # registry view INSIDE the session (WP-3 lore: choices providers read
  # live reactives; outside a reactive context they degrade to empty)
  w8_view <<- store_registry_view(w8$store, "resultspace.sample_general")
  w8_meta <<- .w8_meta(w8$store)
})
ok(is.null(w8_err),
   "WP-8 sg: module boots headless context-built (the new L1 call-site shape)")
ok(ut_cmp_identical(sort(vapply(w8_view, function(r) r$id, "")),
                    PANEL_W8_SG_IDS),
   "WP-8 sg: exactly the 27 frozen bindings register (3 own + 19 attr4 + 1 survival + 4 batch through the child store; D12 untouched)")
ok(ut_cmp_identical(w8_view, PANEL_W8_SG_GOLDEN_VIEW),
   "WP-8 sg: registry view byte-identical to the frozen pre-migration golden (registration order preserved: xax trio, attr4, survival_censor, batch)")
ok(ut_cmp_identical(w8_meta, PANEL_W8_SG_GOLDEN_META),
   "WP-8 sg: raw binding meta byte-identical incl. agent_writable TRUE x27 (survival_censor + batch_* ride this child store)")

# ------------------- cold values + module return (context boot) ----
w8c <- new.env()
w8c_err <- NULL; w8c_cold <- NULL; w8c_ret <- NULL
shiny::testServer(.w8_mkctx(w8c), {
  w8c_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w8c_cold <<- store_read(w8c$store)
  w8c_ret <<- tryCatch(shiny::isolate(w8c$ret()), error = function(e) NULL)
})
ok(is.null(w8c_err), "WP-8 sg: cold-boot capture stays alive")
ok(ut_cmp_identical(w8c_cold, PANEL_W8_SG_GOLDEN_COLD),
   "WP-8 sg: cold store values byte-identical (cascade unset; attr4 search_value/scorner defaults; survival + batch NULL)")
ok(ut_cmp_identical(w8c_ret, PANEL_W8_SG_GOLDEN_RET),
   "WP-8 sg: module return (the analyst_sample_general payload source) byte-identical in a cold boot")

# ------------------- legacy boot registers identically (D13) ----
w8l <- new.env()
w8l_err <- NULL; w8l_view <- NULL; w8l_meta <- NULL; w8l_cold <- NULL
shiny::testServer(.w8_boot(w8l), {
  w8l_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w8l_view <<- store_registry_view(w8l$store, "resultspace.sample_general")
  w8l_meta <<- .w8_meta(w8l$store)
  w8l_cold <<- store_read(w8l$store)
})
ok(is.null(w8l_err),
   "WP-8 sg: legacy direct instantiation (store formal only, no context) still boots")
ok(ut_cmp_identical(w8l_view, PANEL_W8_SG_GOLDEN_VIEW) &&
   ut_cmp_identical(w8l_meta, PANEL_W8_SG_GOLDEN_META) &&
   ut_cmp_identical(w8l_cold, PANEL_W8_SG_GOLDEN_COLD),
   "WP-8 sg: legacy store boot registers and seeds byte-identically to the golden (both boot paths agree)")
ok(ut_cmp_identical(w8l_view, w8_view),
   "WP-8 sg: context boot and legacy boot register identically (D13 argument coalescing)")

# ------------------- cascade user picks + choices (synthetic) ----
w8p <- new.env()
w8p_err <- NULL; w8p_vals <- NULL; w8p_ret <- NULL; w8p_ch <- NULL
shiny::testServer(.w8_mkctx(w8p), {
  w8p_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  session$setInputs(`sg-tris_sample_general-analysis` = "Surv")
  session$flushReact()
  session$setInputs(`sg-tris_sample_general-subset` = "All")
  session$flushReact()
  session$setInputs(`sg-tris_sample_general-variable` = "time")
  session$flushReact(); session$flushReact()
  w8p_vals <<- store_read(w8p$store)
  w8p_ret <<- tryCatch(shiny::isolate(w8p$ret()), error = function(e) NULL)
  v <- store_registry_view(w8p$store, "resultspace.sample_general")
  w8p_ch <<- stats::setNames(
    lapply(v[match(.w8_TRIO, vapply(v, function(r) r$id, ""))],
           `[[`, "allowed_values"), .w8_TRIO)
})
ok(is.null(w8p_err), "WP-8 sg: user-pick boot stays alive")
ok(ut_cmp_identical(w8p_vals, PANEL_W8_SG_GOLDEN_PICK_VALS),
   "WP-8 sg: user Surv|All|time pick reproduces the frozen store values (27 keys, byte-identical)")
ok(ut_cmp_identical(w8p_ret, PANEL_W8_SG_GOLDEN_PICK_RET),
   "WP-8 sg: module return reports the settled Surv triple (analyst_sample_general payload content)")
ok(ut_cmp_identical(w8p_ch, PANEL_W8_SG_GOLDEN_PICK_CH),
   "WP-8 sg: cascade choices stay state-dependent (subset narrows to All, variable to time)")
ok("Surv" %in% PANEL_W8_SG_GOLDEN_VIEW[[1]]$allowed_values &&
   !("Surv" %in% PANEL_W7_FG_GOLDEN_VIEW[[1]]$allowed_values),
   "WP-8 sg: KEEPING-Surv is the panel's signature - the analysis level list carries Surv where feature-general's frozen golden excludes it (deliberate opposites)")

# ------------------- demo-data boot (the L1 data shape) ----
w8d <- new.env()
w8d_err <- NULL; w8d_view <- NULL; w8d_cold <- NULL; w8d_ret <- NULL
shiny::testServer(.w8_boot(
  w8d, pd = Biobase::pData(dat), ex = Biobase::exprs(dat),
  j = colnames(Biobase::exprs(dat))[1]), {
  w8d_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  v <- store_registry_view(w8d$store, "resultspace.sample_general")
  vids <- vapply(v, function(r) r$id, "")
  w8d_view <<- v[match(.w8_TRIO, vids)]
  w8d_cold <<- store_read(w8d$store)[.w8_TRIO]
  w8d_ret <<- tryCatch(shiny::isolate(w8d$ret()), error = function(e) NULL)
})
ok(is.null(w8d_err), "WP-8 sg: demo-data boot stays alive")
ok(ut_cmp_identical(w8d_ret, PANEL_W8_SG_DEMO_RET),
   "WP-8 sg: demo module return byte-identical in a cold boot")
ok(ut_cmp_identical(w8d_view, PANEL_W8_SG_DEMO_OWN_VIEW),
   "WP-8 sg: demo own-binding registry view byte-identical (General/PCA/Surv/Feature analysis levels; head-50 variable order)")
ok(ut_cmp_identical(w8d_cold, PANEL_W8_SG_DEMO_OWN_COLD),
   "WP-8 sg: demo own-binding cold values byte-identical")

# demo user pick: Surv | all | OS (the Kaplan-Meier signature pick)
w8dp <- new.env()
w8dp_err <- NULL; w8dp_vals <- NULL; w8dp_ret <- NULL; w8dp_ch <- NULL
shiny::testServer(.w8_boot(
  w8dp, pd = Biobase::pData(dat), ex = Biobase::exprs(dat),
  j = colnames(Biobase::exprs(dat))[1]), {
  w8dp_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  session$setInputs(`sg-tris_sample_general-analysis` = "Surv")
  session$flushReact()
  session$setInputs(`sg-tris_sample_general-subset` = "all")
  session$flushReact()
  session$setInputs(`sg-tris_sample_general-variable` = "OS")
  session$flushReact(); session$flushReact()
  w8dp_vals <<- store_read(w8dp$store)[.w8_TRIO]
  w8dp_ret <<- tryCatch(shiny::isolate(w8dp$ret()), error = function(e) NULL)
  v <- store_registry_view(w8dp$store, "resultspace.sample_general")
  vids <- vapply(v, function(r) r$id, "")
  w8dp_ch <<- stats::setNames(
    lapply(v[match(.w8_TRIO, vids)], `[[`, "allowed_values"), .w8_TRIO)
})
ok(is.null(w8dp_err), "WP-8 sg: demo user-pick boot stays alive")
ok(ut_cmp_identical(w8dp_vals, PANEL_W8_SG_DEMO_PICK_VALS),
   "WP-8 sg: demo Surv|all|OS pick reproduces the frozen store values")
ok(ut_cmp_identical(w8dp_ret, PANEL_W8_SG_DEMO_PICK_RET),
   "WP-8 sg: demo module return reports the settled Surv triple")
ok(ut_cmp_identical(w8dp_ch, PANEL_W8_SG_DEMO_PICK_CH),
   "WP-8 sg: demo cascade choices stay state-dependent (variable narrows to OS/DFS - the two demo survival variables)")

# ------------------- legacy status restore paths ----
# the Surv triple re-applies through .sg_apply_triple(origin="restore")
w8r <- new.env()
w8r_err <- NULL; w8r_vals <- NULL; w8r_ret <- NULL
shiny::testServer(.w8_boot(w8r, status_payload = list(
  xax = c("Surv", "All", "time"))), {
  w8r_err <<- tryCatch({
    for (i in 1:6) session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w8r_vals <<- store_read(w8r$store)
  w8r_ret <<- tryCatch(shiny::isolate(w8r$ret()), error = function(e) NULL)
})
ok(is.null(w8r_err),
   "WP-8 sg: boot with a legacy analyst_sample_general payload stays alive")
ok(ut_cmp_identical(w8r_vals, PANEL_W8_SG_GOLDEN_RESTORE_VALS),
   "WP-8 sg: the legacy restore path re-applies the Surv triple through the store (27 keys byte-identical)")
ok(ut_cmp_identical(w8r_ret, PANEL_W8_SG_GOLDEN_RESTORE_RET),
   "WP-8 sg: module return stays cold after a restore-only boot (the triple lands in the store; v1 has not been user-driven)")

# a saved triple with an UNSET variable pushes the canonical unset
w8u <- new.env()
w8u_err <- NULL; w8u_vals <- NULL
shiny::testServer(.w8_boot(w8u, status_payload = list(
  xax = list("General", "All", NULL))), {
  w8u_err <<- tryCatch({
    for (i in 1:6) session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w8u_vals <<- store_read(w8u$store)[.w8_TRIO]
})
ok(is.null(w8u_err), "WP-8 sg: unset-variable restore boot stays alive")
ok(ut_cmp_identical(w8u_vals, PANEL_W8_SG_GOLDEN_UNSET_VALS),
   "WP-8 sg: an unset saved variable pushes the canonical unset --select-- (todo 2.3 semantics preserved)")

# demo restore triple
w8dr <- new.env()
w8dr_err <- NULL; w8dr_vals <- NULL
shiny::testServer(.w8_boot(
  w8dr, pd = Biobase::pData(dat), ex = Biobase::exprs(dat),
  j = colnames(Biobase::exprs(dat))[1],
  status_payload = list(xax = c("General", "All", "Cell.line"))), {
  w8dr_err <<- tryCatch({
    for (i in 1:6) session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w8dr_vals <<- store_read(w8dr$store)[.w8_TRIO]
})
ok(is.null(w8dr_err), "WP-8 sg: demo restore boot stays alive")
ok(ut_cmp_identical(w8dr_vals, PANEL_W8_SG_DEMO_RESTORE_VALS),
   "WP-8 sg: the demo restore triple re-applies byte-identically (General|All|Cell.line)")

# ------------------- the triset derivation (KEEPING Surv) ----
ts_full8 <- trisetter(expr = Biobase::exprs(dat), meta = Biobase::pData(dat),
                      combine = "pheno")
ok(ut_cmp_identical(list(dim = dim(ts_full8),
                         analysis_levels = unique(ts_full8[, 1]),
                         non_feature_pca_rows =
                           ts_full8[!ts_full8[, 1] %in% c("Feature", "PCA"), , drop = FALSE]),
                    list(dim = PANEL_W8_SG_TRISET$dim,
                         analysis_levels = PANEL_W8_SG_TRISET$analysis_levels,
                         non_feature_pca_rows = PANEL_W8_SG_TRISET$non_feature_pca_rows)),
   "WP-8 sg: the cascade triset source over demo pdata = trisetter(combine=pheno) KEEPING Surv (2721x3; General/PCA/Surv/Feature; the 11 pheno rows)")
ok(ut_cmp_identical(PANEL_W8_SG_TRISET$surv_rows,
                    PANEL_W7_FG_TRISET$surv_rows_in_full) &&
   ut_cmp_identical(PANEL_W8_SG_TRISET$dim - PANEL_W7_FG_TRISET$dim, c(2L, 0L)) &&
   ut_cmp_identical(PANEL_W8_SG_TRISET$variable_head50[1:15],
                    PANEL_W7_FG_TRISET$variable_head50[1:15]) &&
   ut_cmp_identical(setdiff(PANEL_W8_SG_TRISET$variable_head50,
                            PANEL_W7_FG_TRISET$variable_head50),
                    c("OS", "DFS")),
   "WP-8 sg: cross-panel contrast - the same 2 Surv rows feature-general excluded (2719) are exactly the delta here (2721); the shared variable order agrees through PC8(2.5%) and diverges only by OS/DFS entering at positions 16-17")

# ------------------- retention (D3 unification) ----
# structural (deterministic, WP-0 MINOR-2 pattern): the module's ONE
# panel_keep closure holds the helper's bind-sync observer and the
# view-type carrier observer (previously the SPLIT .sg_view_keep list);
# the beeswarm htest observer keeps its own module-local list (out of
# the plan's unification scope), and the survival/batch/attr4
# submodules own their retention internally (untouched). Functional
# (retention direction, PD-8): after gc() x2 a user Surv pick still
# reaches the store.
w8k <- new.env()
w8k_err <- NULL; w8k_nobs <- NULL; w8k_vals <- NULL
w8k_vt_before <- NULL; w8k_vt_after <- NULL
shiny::testServer(.w8_mkctx(w8k), {
  w8k_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  obs <- get("obs", envir = environment(w8k$ctx$keep))
  w8k_nobs <<- length(obs)
  w8k_vt_before <<- is.null(session$output$`sg-sample_general_plot`)
  gc(); gc()
  session$setInputs(`sg-tris_sample_general-analysis` = "Surv")
  session$flushReact()
  session$setInputs(`sg-tris_sample_general-subset` = "All")
  session$flushReact()
  session$setInputs(`sg-tris_sample_general-variable` = "time")
  session$flushReact(); session$flushReact()
  w8k_vals <<- store_read(w8k$child)[.w8_TRIO]
  # WP-8 review NIT-2: the view-type WRITE direction after gc() x2 - the
  # kept carrier observer still flips .sg_view_type to "surv", visible
  # as the renderUI swapping NULL for the survival container (mirrors
  # the cascade-pick direction asserted below)
  vt <- session$output$`sg-sample_general_plot`
  w8k_vt_after <<- if (is.null(vt)) "NULL" else
    if (grepl("sample_general_surv",
              paste(as.character(vt), collapse = ""), fixed = TRUE))
      "surv" else "other"
})
ok(is.null(w8k_err), "WP-8 sg: retention boot stays alive")
ok(ut_cmp_identical(w8k_nobs, 2L),
   "WP-8 sg: exactly the bind-sync + view-type observers land in ONE panel_keep list (D3 unification of .sg_keep + .sg_view_keep; htest stays module-local)")
ok(ut_cmp_identical(
  unname(unlist(w8k_vals)),
  c("Surv", "All", "time")),
   "WP-8 sg: the kept bind observer still fires after gc() x2 (the Surv cascade pick reaches the store)")
ok(ut_cmp_identical(w8k_vt_before, TRUE) && ut_cmp_identical(w8k_vt_after, "surv"),
   "WP-8 sg: the kept view-type carrier still fires after gc() x2 (the Surv pick flips the rendered view from NULL to the Kaplan-Meier container; WP-8 review NIT-2)")

# ------------------- submodule bindings riding the store ----
# the survival submodule registers survival_censor through THIS panel's
# child store (store formal wiring verbatim); the batch comparison
# module the four batch_* keys - all pinned by the goldens above; these
# spot-checks make the group explicit
w8sv <- Filter(function(r) identical(r$id, "resultspace.sample_general.survival_censor"),
               PANEL_W8_SG_GOLDEN_VIEW)[[1]]
ok(ut_cmp_identical(w8sv$kind, "numeric") &&
   ut_cmp_identical(w8sv$label, "Censor time") &&
   ut_cmp_identical(w8sv$allowed_values, NULL),
   "WP-8 sg: the survival submodule's censor slider rides this child store (kind numeric, open bounds, agent-writable)")
ok(ut_cmp_identical(
  unname(vapply(Filter(function(r) grepl("\\.batch_", r$id),
                        PANEL_W8_SG_GOLDEN_VIEW),
                function(r) r$kind, "")),
  c("boolean", "boolean", "select", "select")),
   "WP-8 sg: the batch comparison module's two switches + two selected-row bindings ride this child store")

# ------------------- the store-less legacy branch stays alive ----
w8n_err <- NULL
shiny::testServer(.w8_mklegacy(), {
  w8n_err <<- tryCatch({
    session$flushReact(); session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
})
ok(is.null(w8n_err),
   "WP-8 sg: the no-store legacy branch (reactiveVal cascade path) still boots")

# ------------------- app level: the REAL L1 composition ----
w8l1 <- new.env()
w8l1$store <- widget_store_new()
app_w8l1 <- function(input, output, session) {
  w8l1$status <- omicsViewer:::L1_result_space_module(
    "rs",
    reactive_expr = reactive(Biobase::exprs(dat)),
    reactive_phenoData = reactive(Biobase::pData(dat)),
    reactive_featureData = reactive(fd0),
    status = reactive(NULL),
    store = w8l1$store)
}
w8l1_err <- NULL; w8l1_payload <- NULL
shiny::testServer(app_w8l1, {
  w8l1_err <<- tryCatch({
    for (i in 1:3) session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w8l1_payload <<- tryCatch(shiny::isolate(w8l1$status()),
                            error = function(e) conditionMessage(e))
})
ok(is.null(w8l1_err),
   "WP-8 sg: the real L1 composition boots headless post-migration (context built at the call site)")
ok(ut_cmp_identical(names(w8l1_payload), PANEL_W8_SG_GOLDEN_L1_NAMES),
   "WP-8 sg: L1 status payload keeps its exact key order (snapshot schema untouched)")
ok(ut_cmp_identical(w8l1_payload$analyst_sample_general,
                    PANEL_W8_SG_GOLDEN_L1_SG),
   "WP-8 sg: analyst_sample_general payload byte-identical pre/post migration (cold boot: xax NULL, attr4 nine-field NULL shape)")
ok(ut_cmp_identical(
  sort(intersect(vapply(store_registry_view(w8l1$store), function(r) r$id, ""),
                 PANEL_W8_SG_IDS)),
  PANEL_W8_SG_IDS),
   "WP-8 sg: the real L1 boot carries all 27 sample_general keys (the D12 freeze above pins the full 72-id vector)")

## --------------- WP-9: ORA migration (D3/D4/D13) ----
# The FIRST agent-coupled migration: ORA is the plan's reference citizen,
# and three agent surfaces address this panel's store keys - the
# set_enrichment_parameters tool (apply_agent_enrichment: validate ->
# store_apply over the resultspace.ora.* keys), the panel_results
# excerpt reader (the module return's table attr, WP-7/F4) and the
# enrichment normalizer (agent_normalize_enrichment_update). NONE of
# them are rewritten here; the migration is the WP-5..WP-8
# pure-migration pattern (retention -> panel_keep, D3; cascade
# boilerplate -> panel_triset_bind, D4; context threading, D13) and the
# agent coupling is proven UNCHANGED by goldens: the L0 apply-callback
# receipt over a real GS-attached dataset is byte-identical pre/post.
# The restore observer and the return-with-table-attr stay in the
# module verbatim (D4 panel semantics; the run-preconditioned
# selected_row semantics of the results table ride the untouched
# dataTableDownload_module wiring).
#
# The WP-0 ORA-cascade FROZEN GOLDEN FIXTURE earlier in this suite (the
# verbatim pre-WP-0 copy vs panel_triset_bind) now doubles as THE
# pre-migration record of the hand-rolled block: the real module adopted
# the helper in this WP, so that equivalence assertion (helper == frozen
# verbatim copy) is exactly the claim that the migration is
# behavior-preserving; it stays green and untouched.
#
# Golden-fixture provenance (PD-29/PD-33/PD-38/PD-43 school): all
# goldens below were captured from PRE-MIGRATION HEAD 7213ffa BEFORE
# any WP-9 edit (before the WP-8 NIT-2 test commit fe248b7, for strict
# provenance - it touches no R/ code), by booting the then-current
# enrichment_analysis_module exactly like the L1 call site / the
# agentWidgets S4 shape and dput()ing the outputs: the 4-binding
# registry view + raw meta (incl. agent_writable) + cold values +
# module return + pre-run selected_row semantics over demo fd; the
# user-pick and legacy-restore triples; the GS-attached run surface
# (boot table, collapsed-run table attr, live selected_row choices);
# the set_enrichment_parameters round-trip through the L0 apply
# callback path over the REAL L1 composition with a GS-attached
# dataset; the real-L1 analyst_ora payload + status key order; and the
# triset derivation over demo fd (trisetter combine="none"). The
# IDENTICAL capture script was re-run POST-migration: all 20 outputs
# byte-identical (the PD-33/PD-38 pre/post proof, stash-free). The
# fixtures are FROZEN: never regenerate them from the migrated code
# (that would make the equivalence claim circular).
enrichment_analysis_module <- omicsViewer:::enrichment_analysis_module
agent_normalize_enrichment_update <- omicsViewer:::agent_normalize_enrichment_update

# -------------------------------------------------- frozen goldens ----
# the 4 ids (3 cascade + the run-preconditioned results-row selection)
PANEL_W9_ORA_IDS <- sort(paste0("resultspace.ora.",
  c("xax_analysis", "xax_subset", "xax_variable", "selected_row")))

PANEL_W9_ORA_GOLDEN_VIEW <- 
list(list(id = "resultspace.ora.xax_analysis", kind = "select", 
    label = "Collapse category", help = "Annotation category whose values the input features are collapsed on before testing over-representation", 
    depends_on = character(0), min = NULL, max = NULL, allowed_values = c("General", 
    "PCA", "mean", "n value", "quantile", "ttest", "Rank", "GS", 
    "StringDB", "Cor", "SeqLogo", "ResponseCurve")), list(id = "resultspace.ora.xax_subset", 
    kind = "select", label = "Collapse subcategory", help = "Subcategory within the collapse category", 
    depends_on = "resultspace.ora.xax_analysis", min = NULL, 
    max = NULL, allowed_values = character(0)), list(id = "resultspace.ora.xax_variable", 
    kind = "select_cascaded", label = "Collapse variable", help = "Variable the input features are collapsed on (feature annotation column)", 
    depends_on = c("resultspace.ora.xax_analysis", "resultspace.ora.xax_subset"
    ), min = NULL, max = NULL, allowed_values = c("Gene.name", 
    "Protein.ID", "Intensity", "PC1(10.5%)", "PC2(7.2%)", "PC3(5.4%)", 
    "PC4(4.8%)", "PC5(3.6%)", "PC6(2.9%)", "PC7(2.8%)", "PC8(2.5%)", 
    "RE", "ME", "LE", "MT", "WT", "pvalue", "log.pvalue", "fdr", 
    "log.fdr", "mean.diff", "X786O_NCI60", "A498_NCI60", "A549_NCI60", 
    "ACHN_NCI60", "BT549_NCI60", "CAKI1_NCI60", "CCRFCEM_NCI60", 
    "COLO205_NCI60", "DU145_NCI60", "EKVX_NCI60", "HCC2998_NCI60", 
    "HCT116_NCI60", "HCT15_NCI60", "HL60_NCI60", "HOP62_NCI60", 
    "HOP92_NCI60", "HS578T_NCI60", "HT29_NCI60", "IGROV1_NCI60", 
    "K562_NCI60", "KM12_NCI60", "LOXIMVI_NCI60", "M14_NCI60", 
    "MALME3M_NCI60", "MCF7_NCI60", "MDAMB231_NCI60", "MDAMB435_NCI60", 
    "MDAMB468_NCI60", "MOLT4_NCI60")), list(id = "resultspace.ora.selected_row", 
    kind = "select", label = "Selected pathway", help = "Gene set whose row is selected in the results table; drives the overlap-genes table below", 
    depends_on = character(0), min = NULL, max = NULL, allowed_values = character(0)))

PANEL_W9_ORA_GOLDEN_META <- 
list(list(id = "resultspace.ora.xax_analysis", kind = "select", 
    label = "Collapse category", help = "Annotation category whose values the input features are collapsed on before testing over-representation", 
    depends_on = character(0), agent_writable = TRUE), list(id = "resultspace.ora.xax_subset", 
    kind = "select", label = "Collapse subcategory", help = "Subcategory within the collapse category", 
    depends_on = "resultspace.ora.xax_analysis", agent_writable = TRUE), 
    list(id = "resultspace.ora.xax_variable", kind = "select_cascaded", 
        label = "Collapse variable", help = "Variable the input features are collapsed on (feature annotation column)", 
        depends_on = c("resultspace.ora.xax_analysis", "resultspace.ora.xax_subset"
        ), agent_writable = TRUE), list(id = "resultspace.ora.selected_row", 
        kind = "select", label = "Selected pathway", help = "Gene set whose row is selected in the results table; drives the overlap-genes table below", 
        depends_on = character(0), agent_writable = TRUE))

PANEL_W9_ORA_GOLDEN_COLD <- 
list(resultspace.ora.xax_analysis = NULL, resultspace.ora.xax_subset = NULL, 
    resultspace.ora.xax_variable = NULL, resultspace.ora.selected_row = NULL)

PANEL_W9_ORA_GOLDEN_ALLOWED_COLD <- 
list(resultspace.ora.xax_analysis = c("General", "PCA", "mean", 
"n value", "quantile", "ttest", "Rank", "GS", "StringDB", "Cor", 
"SeqLogo", "ResponseCurve"), resultspace.ora.xax_subset = character(0), 
    resultspace.ora.xax_variable = c("Gene.name", "Protein.ID", 
    "Intensity", "PC1(10.5%)", "PC2(7.2%)", "PC3(5.4%)", "PC4(4.8%)", 
    "PC5(3.6%)", "PC6(2.9%)", "PC7(2.8%)", "PC8(2.5%)", "RE", 
    "ME", "LE", "MT", "WT", "pvalue", "log.pvalue", "fdr", "log.fdr", 
    "mean.diff", "X786O_NCI60", "A498_NCI60", "A549_NCI60", "ACHN_NCI60", 
    "BT549_NCI60", "CAKI1_NCI60", "CCRFCEM_NCI60", "COLO205_NCI60", 
    "DU145_NCI60", "EKVX_NCI60", "HCC2998_NCI60", "HCT116_NCI60", 
    "HCT15_NCI60", "HL60_NCI60", "HOP62_NCI60", "HOP92_NCI60", 
    "HS578T_NCI60", "HT29_NCI60", "IGROV1_NCI60", "K562_NCI60", 
    "KM12_NCI60", "LOXIMVI_NCI60", "M14_NCI60", "MALME3M_NCI60", 
    "MCF7_NCI60", "MDAMB231_NCI60", "MDAMB435_NCI60", "MDAMB468_NCI60", 
    "MOLT4_NCI60"), resultspace.ora.selected_row = character(0))

PANEL_W9_ORA_GOLDEN_RET <- 
list(xax = NULL, table = NULL)

PANEL_W9_ORA_GOLDEN_PICK_VALS <- 
list(resultspace.ora.xax_analysis = "General", resultspace.ora.xax_subset = "All", 
    resultspace.ora.xax_variable = "Gene.name", resultspace.ora.selected_row = NULL)

PANEL_W9_ORA_GOLDEN_PICK_RET <- 
list(xax = list(analysis = "General", subset = "All", variable = "Gene.name"), 
    table = NULL)

PANEL_W9_ORA_GOLDEN_PICK_ALLOWED <- 
list(resultspace.ora.xax_analysis = c("General", "PCA", "mean", 
"n value", "quantile", "ttest", "Rank", "GS", "StringDB", "Cor", 
"SeqLogo", "ResponseCurve"), resultspace.ora.xax_subset = "All", 
    resultspace.ora.xax_variable = c("Gene.name", "Protein.ID", 
    "Intensity"), resultspace.ora.selected_row = character(0))

PANEL_W9_ORA_GOLDEN_RESTORE_VALS <- 
list(resultspace.ora.xax_analysis = "General", resultspace.ora.xax_subset = "All", 
    resultspace.ora.xax_variable = "Gene.name", resultspace.ora.selected_row = NULL)

PANEL_W9_ORA_GOLDEN_RESTORE_RET <- 
list(xax = NULL, table = NULL)

PANEL_W9_ORA_TRISET <- 
list(dim = c(713L, 3L), analysis_levels = c("General", "PCA", 
"mean", "n value", "quantile", "ttest", "Rank", "GS", "StringDB", 
"Cor", "SeqLogo", "ResponseCurve"), first_rows = structure(c("General", 
"General", "General", "PCA", "PCA", "PCA", "All", "All", "All", 
"All", "All", "All", "Gene.name", "Protein.ID", "Intensity", 
"PC1(10.5%)", "PC2(7.2%)", "PC3(5.4%)"), dim = c(6L, 3L)), variable_head50 = c("Gene.name", 
"Protein.ID", "Intensity", "PC1(10.5%)", "PC2(7.2%)", "PC3(5.4%)", 
"PC4(4.8%)", "PC5(3.6%)", "PC6(2.9%)", "PC7(2.8%)", "PC8(2.5%)", 
"RE", "ME", "LE", "MT", "WT", "pvalue", "log.pvalue", "fdr", 
"log.fdr", "mean.diff", "X786O_NCI60", "A498_NCI60", "A549_NCI60", 
"ACHN_NCI60", "BT549_NCI60", "CAKI1_NCI60", "CCRFCEM_NCI60", 
"COLO205_NCI60", "DU145_NCI60", "EKVX_NCI60", "HCC2998_NCI60", 
"HCT116_NCI60", "HCT15_NCI60", "HL60_NCI60", "HOP62_NCI60", "HOP92_NCI60", 
"HS578T_NCI60", "HT29_NCI60", "IGROV1_NCI60", "K562_NCI60", "KM12_NCI60", 
"LOXIMVI_NCI60", "M14_NCI60", "MALME3M_NCI60", "MCF7_NCI60", 
"MDAMB231_NCI60", "MDAMB435_NCI60", "MDAMB468_NCI60", "MOLT4_NCI60"
))

PANEL_W9_ORA_AGENT_RECEIPT <- 
list(method = "ora", panel_tab = "ORA", applied = c("resultspace.ora.xax_analysis", 
"resultspace.ora.selected_row", "resultspace.analyst_tab", "resultspace.ora.xax_subset", 
"resultspace.ora.xax_variable"), applied_values = list(resultspace.ora.xax_analysis = "General", 
    resultspace.ora.selected_row = "GS_A", resultspace.analyst_tab = "ORA", 
    resultspace.ora.xax_subset = "All", resultspace.ora.xax_variable = "Gene.name"), 
    unchanged = character(0), rejected = list(), note = "The panel opens and the enrichment recomputes from the currently selected features (ORA) or the full ranking (fGSEA); read the top rows back through get_omics_viewer_state(sections='panel_results').")

PANEL_W9_ORA_AGENT_VALS <- 
list(resultspace.analyst_tab = "ORA", resultspace.feature_general.xax_analysis = NULL, 
    resultspace.feature_general.xax_subset = NULL, resultspace.feature_general.xax_variable = NULL, 
    resultspace.feature_general.plot_type = NULL, resultspace.feature_general.regression_line = FALSE, 
    resultspace.feature_general.attr4.color_analysis = NULL, 
    resultspace.feature_general.attr4.color_subset = NULL, resultspace.feature_general.attr4.color_variable = NULL, 
    resultspace.feature_general.attr4.shape_analysis = NULL, 
    resultspace.feature_general.attr4.shape_subset = NULL, resultspace.feature_general.attr4.shape_variable = NULL, 
    resultspace.feature_general.attr4.size_analysis = NULL, resultspace.feature_general.attr4.size_subset = NULL, 
    resultspace.feature_general.attr4.size_variable = NULL, resultspace.feature_general.attr4.tooltip_analysis = NULL, 
    resultspace.feature_general.attr4.tooltip_subset = NULL, 
    resultspace.feature_general.attr4.tooltip_variable = NULL, 
    resultspace.feature_general.attr4.search_analysis = NULL, 
    resultspace.feature_general.attr4.search_subset = NULL, resultspace.feature_general.attr4.search_variable = NULL, 
    resultspace.feature_general.attr4.search_value = character(0), 
    resultspace.feature_general.attr4.xcut = NULL, resultspace.feature_general.attr4.ycut = NULL, 
    resultspace.feature_general.attr4.scorner = "None", resultspace.fgsea.xax_analysis = NULL, 
    resultspace.fgsea.xax_subset = NULL, resultspace.fgsea.xax_variable = NULL, 
    resultspace.fgsea.selected_row = NULL, resultspace.ora.xax_analysis = "General", 
    resultspace.ora.xax_subset = "All", resultspace.ora.xax_variable = "Gene.name", 
    resultspace.ora.selected_row = "GS_A", resultspace.stringdb.taxonomy = NULL, 
    resultspace.stringdb.show_labels = NULL, resultspace.stringdb.selected_row = NULL, 
    resultspace.sample_general.xax_analysis = NULL, resultspace.sample_general.xax_subset = NULL, 
    resultspace.sample_general.xax_variable = NULL, resultspace.sample_general.attr4.color_analysis = NULL, 
    resultspace.sample_general.attr4.color_subset = NULL, resultspace.sample_general.attr4.color_variable = NULL, 
    resultspace.sample_general.attr4.shape_analysis = NULL, resultspace.sample_general.attr4.shape_subset = NULL, 
    resultspace.sample_general.attr4.shape_variable = NULL, resultspace.sample_general.attr4.size_analysis = NULL, 
    resultspace.sample_general.attr4.size_subset = NULL, resultspace.sample_general.attr4.size_variable = NULL, 
    resultspace.sample_general.attr4.tooltip_analysis = NULL, 
    resultspace.sample_general.attr4.tooltip_subset = NULL, resultspace.sample_general.attr4.tooltip_variable = NULL, 
    resultspace.sample_general.attr4.search_analysis = NULL, 
    resultspace.sample_general.attr4.search_subset = NULL, resultspace.sample_general.attr4.search_variable = NULL, 
    resultspace.sample_general.attr4.search_value = character(0), 
    resultspace.sample_general.attr4.xcut = NULL, resultspace.sample_general.attr4.ycut = NULL, 
    resultspace.sample_general.attr4.scorner = "None", resultspace.sample_general.survival_censor = NULL, 
    resultspace.sample_general.batch_show_phenotype = NULL, resultspace.sample_general.batch_show_features = NULL, 
    resultspace.sample_general.batch_phenotype_selected_row = NULL, 
    resultspace.sample_general.batch_feature_selected_row = NULL, 
    resultspace.geneshot.term = NULL, resultspace.geneshot.xax_analysis = NULL, 
    resultspace.geneshot.xax_subset = NULL, resultspace.geneshot.xax_variable = NULL, 
    resultspace.ptm.xax_analysis = NULL, resultspace.ptm.xax_subset = NULL, 
    resultspace.ptm.xax_variable = NULL, resultspace.dose_response.feature_table = NULL, 
    resultspace.dose_response.param_table = NULL)

PANEL_W9_ORA_GS_BOOT_RET <- 
list(xax = NULL, table = list(class = "data.frame", names = c("pathway", 
"desc", "p.value", "p.adjusted", "OR", "size_overlap", "size_geneset", 
"size_input", "size_backgroung", "overlap_ids"), nrow = 2L, pathways = c("GS_A", 
"GS_D")))

PANEL_W9_ORA_GS_VIEW <- 
list(list(id = "resultspace.ora.xax_analysis", kind = "select", 
    label = "Collapse category", help = "Annotation category whose values the input features are collapsed on before testing over-representation", 
    depends_on = character(0), min = NULL, max = NULL, allowed_values = "General"), 
    list(id = "resultspace.ora.xax_subset", kind = "select", 
        label = "Collapse subcategory", help = "Subcategory within the collapse category", 
        depends_on = "resultspace.ora.xax_analysis", min = NULL, 
        max = NULL, allowed_values = character(0)), list(id = "resultspace.ora.xax_variable", 
        kind = "select_cascaded", label = "Collapse variable", 
        help = "Variable the input features are collapsed on (feature annotation column)", 
        depends_on = c("resultspace.ora.xax_analysis", "resultspace.ora.xax_subset"
        ), min = NULL, max = NULL, allowed_values = "Gene.name"), 
    list(id = "resultspace.ora.selected_row", kind = "select", 
        label = "Selected pathway", help = "Gene set whose row is selected in the results table; drives the overlap-genes table below", 
        depends_on = character(0), min = NULL, max = NULL, allowed_values = c("GS_A", 
        "GS_D")))

PANEL_W9_ORA_GS_PICK_VALS <- 
list(resultspace.ora.xax_analysis = "General", resultspace.ora.xax_subset = "All", 
    resultspace.ora.xax_variable = "Gene.name", resultspace.ora.selected_row = NULL)

PANEL_W9_ORA_GS_PICK_RET <- 
list(xax = list(analysis = "General", subset = "All", variable = "Gene.name"), 
    table = list(class = "data.frame", names = c("pathway", "desc", 
    "p.value", "p.adjusted", "OR", "size_overlap", "size_geneset", 
    "size_input", "size_backgroung", "overlap_ids"), nrow = 2L, 
        pathways = c("GS_A", "GS_D")))

PANEL_W9_ORA_GS_PICK_ALLOWED <- 
list(resultspace.ora.xax_analysis = "General", resultspace.ora.xax_subset = "All", 
    resultspace.ora.xax_variable = "Gene.name", resultspace.ora.selected_row = c("GS_A", 
    "GS_D"))

PANEL_W9_ORA_GOLDEN_L1_NAMES <- 
c("analyst_active_tab", "analyst_feature_general", "analyst_sample_general", 
"analyst_gene_shot", "analyst_fgsea", "analyst_ora", "analyst_ptm", 
"analyst_dose_response", "analyst_additional")

PANEL_W9_ORA_GOLDEN_L1_ORA <- 
list(xax = NULL)


# -------------------------------------------------- boot factories ----
# GS-attached synthetic fd (the agentWidgets S4 ora shape,
# deterministic) drives the run/table/agent surfaces; demo fd0 (no GS)
# drives the registration/pick/restore surfaces
.w9_gs_fd <- data.frame(
  `General|All|Gene.name` = paste0("g", 1:30),
  row.names = paste0("F", 1:30), check.names = FALSE)
attr(.w9_gs_fd, "GS") <- data.frame(
  featureId = factor(paste0("F", c(1:10, 6:20, 21:30, 1:5, 25:28))),
  gsId = factor(rep(c("GS_A", "GS_B", "GS_C", "GS_D"),
                    c(10, 15, 10, 9))))
set.seed(9)
.w9_pd <- data.frame(group = c("A", "A", "B", "B"),
                     row.names = paste0("S", 1:4))
.w9_ex <- matrix(rnorm(120), nrow = 30,
                 dimnames = list(paste0("F", 1:30), paste0("S", 1:4)))
.w9_sel <- paste0("F", 1:5)
# context-built (the WP-9 L1 call-site shape; holder receives the
# context so the keep closure is inspectable)
.w9_mkctx <- function(holder, fd = fd0, i = NULL, status_payload = NULL) {
  force(holder)
  function(input, output, session) {
    .fd <- reactive(fd); .ii <- reactive(i); .st <- reactive(status_payload)
    holder$store <- widget_store_new()
    holder$child <- widget_store_child(holder$store, "resultspace.ora")
    holder$ctx <- panel_context(
      id = "ora", tab_value = "ORA", store = holder$child,
      data = list(fdata = .fd), selection = list(feature = .ii),
      status = .st, session = session)
    holder$ret <- enrichment_analysis_module(
      "ora", reactive_featureData = .fd, reactive_i = .ii,
      reactive_status = .st, store = holder$child, context = holder$ctx)
  }
}
# legacy direct boot (formals only - the pre-WP-9 shape, still supported)
.w9_boot <- function(env, fd = fd0, i = NULL, status_payload = NULL) {
  force(env)
  function(input, output, session) {
    env$store <- widget_store_new()
    env$child <- widget_store_child(env$store, "resultspace.ora")
    env$ret <- enrichment_analysis_module(
      "ora", reactive_featureData = reactive(fd), reactive_i = reactive(i),
      reactive_status = reactive(status_payload), store = env$child)
  }
}
# store-less legacy boot (no store, no context - the reactiveVal path)
.w9_mklegacy <- function(fd = fd0, i = NULL, status_payload = NULL) {
  function(input, output, session) {
    enrichment_analysis_module(
      "ora", reactive_featureData = reactive(fd), reactive_i = reactive(i),
      reactive_status = reactive(status_payload))
  }
}
.w9_meta <- function(root) {
  ks <- grep("^resultspace\\.ora(\\.|$)", names(root$bindings), value = TRUE)
  lapply(ks, function(k) {
    b <- root$bindings[[k]]
    list(id = k, kind = b$kind, label = b$label, help = b$help,
         depends_on = b$depends_on, agent_writable = b$agent_writable)
  })
}
.w9_allowed <- function(root) {
  v <- store_registry_view(root, "resultspace.ora")
  stats::setNames(lapply(v, `[[`, "allowed_values"),
                  vapply(v, function(r) r$id, ""))
}
# the module return as the status builder / excerpt reader see it (the
# table attr summarized - content pins, not the whole data frame)
.w9_tab_sum <- function(tab) {
  if (!is.data.frame(tab)) return(tab)
  list(class = class(tab)[[1]], names = colnames(tab), nrow = nrow(tab),
       pathways = as.character(tab$pathway))
}
.w9_ret_shape <- function(ret) {
  r <- tryCatch(shiny::isolate(ret()), error = function(e) NULL)
  list(xax = r$xax, table = .w9_tab_sum(attr(r, "table")))
}
.w9_flush <- function(session, n = 6)
  for (k in seq_len(n)) tryCatch(session$flushReact(), error = function(e) NULL)
.w9_TRIO <- paste0("resultspace.ora.",
                   c("xax_analysis", "xax_subset", "xax_variable"))

# ------------------- registration: context boot = the frozen golden ----
w9 <- new.env()
w9_err <- NULL; w9_view <- NULL; w9_meta <- NULL
shiny::testServer(.w9_mkctx(w9), {
  w9_err <<- tryCatch({
    .w9_flush(session)
    NULL
  }, error = function(e) conditionMessage(e))
  # registry view INSIDE the session (WP-3 lore: choices providers read
  # live reactives; outside a reactive context they degrade to empty)
  w9_view <<- store_registry_view(w9$store, "resultspace.ora")
  w9_meta <<- .w9_meta(w9$store)
})
ok(is.null(w9_err),
   "WP-9 ora: module boots headless context-built (the new L1 call-site shape)")
ok(ut_cmp_identical(sort(vapply(w9_view, function(r) r$id, "")),
                    PANEL_W9_ORA_IDS),
   "WP-9 ora: exactly the four frozen bindings register (3 cascade + run-preconditioned selected_row; D12 untouched)")
ok(ut_cmp_identical(w9_view, PANEL_W9_ORA_GOLDEN_VIEW),
   "WP-9 ora: registry view byte-identical to the frozen pre-migration golden (over demo fd: 12 analysis levels, full variable order)")
ok(ut_cmp_identical(w9_meta, PANEL_W9_ORA_GOLDEN_META),
   "WP-9 ora: raw binding meta byte-identical incl. agent_writable TRUE x4 (the set_enrichment_parameters surface)")

# ------------------- cold values + module return + pre-run table attr ----
w9c <- new.env()
w9c_err <- NULL; w9c_cold <- NULL; w9c_ret <- NULL; w9c_allow <- NULL
shiny::testServer(.w9_mkctx(w9c), {
  w9c_err <<- tryCatch({
    .w9_flush(session)
    NULL
  }, error = function(e) conditionMessage(e))
  w9c_cold <<- store_read(w9c$store)
  w9c_ret <<- .w9_ret_shape(w9c$ret)
  w9c_allow <<- .w9_allowed(w9c$store)
})
ok(is.null(w9c_err), "WP-9 ora: cold-boot capture stays alive")
ok(ut_cmp_identical(w9c_cold, PANEL_W9_ORA_GOLDEN_COLD),
   "WP-9 ora: cold store values byte-identical (cascade + selected_row all NULL without a run)")
ok(ut_cmp_identical(w9c_ret, PANEL_W9_ORA_GOLDEN_RET),
   "WP-9 ora: module return (the analyst_ora payload source + panel_results excerpt source) byte-identical: list(xax = NULL), table attr NULL pre-run")
ok(ut_cmp_identical(w9c_allow, PANEL_W9_ORA_GOLDEN_ALLOWED_COLD),
   "WP-9 ora: pre-run selected_row choices are empty (run-preconditioned semantics; only tested pathways become rows)")

# ------------------- user pick + state-dependent choices (demo fd) ----
w9p <- new.env()
w9p_err <- NULL; w9p_vals <- NULL; w9p_ret <- NULL; w9p_allow <- NULL
shiny::testServer(.w9_mkctx(w9p), {
  w9p_err <<- tryCatch({
    .w9_flush(session)
    NULL
  }, error = function(e) conditionMessage(e))
  session$setInputs(`ora-tris_ora-analysis` = "General")
  session$flushReact()
  session$setInputs(`ora-tris_ora-subset` = "All")
  session$flushReact()
  session$setInputs(`ora-tris_ora-variable` = "Gene.name")
  .w9_flush(session, 3)
  w9p_vals <<- store_read(w9p$store)
  w9p_ret <<- .w9_ret_shape(w9p$ret)
  w9p_allow <<- .w9_allowed(w9p$store)
})
ok(is.null(w9p_err), "WP-9 ora: user-pick boot stays alive")
ok(ut_cmp_identical(w9p_vals, PANEL_W9_ORA_GOLDEN_PICK_VALS),
   "WP-9 ora: user General|All|Gene.name pick reproduces the frozen store values (4 keys, byte-identical)")
ok(ut_cmp_identical(w9p_ret, PANEL_W9_ORA_GOLDEN_PICK_RET),
   "WP-9 ora: module return reports the settled triple (analyst_ora payload content)")
ok(ut_cmp_identical(w9p_allow, PANEL_W9_ORA_GOLDEN_PICK_ALLOWED),
   "WP-9 ora: cascade choices stay state-dependent (subset narrows to All, variable to the 3 General columns)")

# ------------------- legacy status restore path ----
w9r <- new.env()
w9r_err <- NULL; w9r_vals <- NULL; w9r_ret <- NULL
shiny::testServer(.w9_boot(w9r, status_payload = list(
  xax = c("General", "All", "Gene.name"))), {
  w9r_err <<- tryCatch({
    .w9_flush(session)
    NULL
  }, error = function(e) conditionMessage(e))
  w9r_vals <<- store_read(w9r$store)
  w9r_ret <<- .w9_ret_shape(w9r$ret)
})
ok(is.null(w9r_err),
   "WP-9 ora: boot with a legacy analyst_ora payload stays alive")
ok(ut_cmp_identical(w9r_vals, PANEL_W9_ORA_GOLDEN_RESTORE_VALS),
   "WP-9 ora: the legacy restore path re-applies the triple through the store (observer untouched, 4 keys byte-identical)")
ok(ut_cmp_identical(w9r_ret, PANEL_W9_ORA_GOLDEN_RESTORE_RET),
   "WP-9 ora: module return stays cold after a restore-only boot (the triple lands in the store; v1 has not been user-driven)")

# ------------------- legacy boot registers identically (D13) ----
w9l <- new.env()
w9l_err <- NULL; w9l_view <- NULL; w9l_meta <- NULL; w9l_cold <- NULL
shiny::testServer(.w9_boot(w9l), {
  w9l_err <<- tryCatch({
    .w9_flush(session)
    NULL
  }, error = function(e) conditionMessage(e))
  w9l_view <<- store_registry_view(w9l$store, "resultspace.ora")
  w9l_meta <<- .w9_meta(w9l$store)
  w9l_cold <<- store_read(w9l$store)
})
ok(is.null(w9l_err),
   "WP-9 ora: legacy direct instantiation (store formal only, no context) still boots")
ok(ut_cmp_identical(w9l_view, PANEL_W9_ORA_GOLDEN_VIEW) &&
   ut_cmp_identical(w9l_meta, PANEL_W9_ORA_GOLDEN_META) &&
   ut_cmp_identical(w9l_cold, PANEL_W9_ORA_GOLDEN_COLD),
   "WP-9 ora: legacy store boot registers and seeds byte-identically to the golden (both boot paths agree)")
ok(ut_cmp_identical(w9l_view, w9_view),
   "WP-9 ora: context boot and legacy boot register identically (D13 argument coalescing)")

# ------------------- the triset derivation (combine = "none") ----
tsd9 <- trisetter(meta = fd0, combine = "none")
ok(ut_cmp_identical(list(dim = dim(tsd9),
                         analysis_levels = unique(tsd9[, 1]),
                         first_rows = utils::head(tsd9, 6L),
                         variable_head50 = utils::head(unique(tsd9[, 3]), 50)),
                    list(dim = PANEL_W9_ORA_TRISET$dim,
                         analysis_levels = PANEL_W9_ORA_TRISET$analysis_levels,
                         first_rows = PANEL_W9_ORA_TRISET$first_rows,
                         variable_head50 = PANEL_W9_ORA_TRISET$variable_head50)),
   "WP-9 ora: the cascade triset source over demo fd = trisetter(combine=none) (713x3; the raw fd column inventory incl. GS/StringDB/Cor/SeqLogo/ResponseCurve levels)")

# ------------------- GS dataset: the run / table-attr surface ----
# module-direct boot (legacy formals) with a real GS attr: the unset
# triple runs ORA on the raw selection by historical contract, so the
# table attr is live at boot; the user collapse pick re-runs on the
# collapsed pathway; selected_row choices follow the run
w9g <- new.env()
w9g_err <- NULL; w9g_boot <- NULL; w9g_view <- NULL
w9g_vals <- NULL; w9g_ret <- NULL; w9g_allow <- NULL
shiny::testServer(.w9_boot(w9g, fd = .w9_gs_fd, i = .w9_sel), {
  w9g_err <<- tryCatch({
    .w9_flush(session, 8)
    NULL
  }, error = function(e) conditionMessage(e))
  w9g_boot <<- .w9_ret_shape(w9g$ret)
  w9g_view <<- store_registry_view(w9g$store, "resultspace.ora")
  session$setInputs(`ora-tris_ora-analysis` = "General")
  session$flushReact()
  session$setInputs(`ora-tris_ora-subset` = "All")
  session$flushReact()
  session$setInputs(`ora-tris_ora-variable` = "Gene.name")
  .w9_flush(session, 8)
  w9g_vals <<- store_read(w9g$store)
  w9g_ret <<- .w9_ret_shape(w9g$ret)
  w9g_allow <<- .w9_allowed(w9g$store)
})
ok(is.null(w9g_err), "WP-9 ora: GS-attached run boot stays alive")
ok(ut_cmp_identical(w9g_boot, PANEL_W9_ORA_GS_BOOT_RET),
   "WP-9 ora: the unset-triple boot run puts the real results table on the return's table attr (panel_results excerpt source; xax still NULL)")
ok(ut_cmp_identical(w9g_view, PANEL_W9_ORA_GS_VIEW),
   "WP-9 ora: GS registry view byte-identical (single General branch; selected_row choices live from the boot run)")
ok(ut_cmp_identical(w9g_vals, PANEL_W9_ORA_GS_PICK_VALS),
   "WP-9 ora: the collapsed run keeps the same settled triple + NULL selected_row (byte-identical)")
ok(ut_cmp_identical(w9g_ret, PANEL_W9_ORA_GS_PICK_RET),
   "WP-9 ora: the collapsed run's table attr is the gene-name-collapsed table (same shape/pathway set pinned)")
ok(ut_cmp_identical(w9g_allow, PANEL_W9_ORA_GS_PICK_ALLOWED),
   "WP-9 ora: post-run selected_row choices are the tested gene sets (run-preconditioned semantics)")

# ------------------- agent round-trip: the L0 apply-callback path ----
# set_enrichment_parameters lands through apply_agent_enrichment's
# exact body (validate -> store_apply on the ROOT store -> receipt,
# note verbatim) over the REAL L1 composition booted with the
# GS-attached dataset (so resultspace.analyst_tab is registered exactly
# as in the app). The frozen receipt is the pre-migration capture: the
# first agent-coupled migration proves the tool surface byte-identical.
w9a <- new.env()
w9a$store <- widget_store_new()
app_w9a <- function(input, output, session) {
  w9a$status <- omicsViewer:::L1_result_space_module(
    "rs",
    reactive_expr = reactive(.w9_ex),
    reactive_phenoData = reactive(.w9_pd),
    reactive_featureData = reactive(.w9_gs_fd),
    reactive_i = reactive(.w9_sel),
    reactive_highlight = reactive(NULL),
    status = reactive(NULL),
    store = w9a$store)
}
w9a_err <- NULL; w9a_receipt <- NULL; w9a_vals <- NULL
shiny::testServer(app_w9a, {
  w9a_err <<- tryCatch({
    .w9_flush(session, 8)
    NULL
  }, error = function(e) conditionMessage(e))
  update <- list(method = "ora", collapse = "General|All|Gene.name",
                 selected_pathway = "GS_A")
  validated <- agent_normalize_enrichment_update(update, .w9_gs_fd)
  receipt <- shiny::isolate(
    store_apply(w9a$store, validated$patch, origin = "agent",
                strict = FALSE))
  w9a_receipt <<- list(
    method = validated$method,
    panel_tab = validated$tab,
    applied = receipt$applied,
    applied_values = receipt$diff,
    unchanged = receipt$skipped,
    rejected = receipt$rejected %||% list(),
    note = paste(
      "The panel opens and the enrichment recomputes from the currently",
      "selected features (ORA) or the full ranking (fGSEA); read the top",
      "rows back through get_omics_viewer_state(sections='panel_results')."))
  .w9_flush(session, 8)
  w9a_vals <<- store_read(w9a$store)
})
ok(is.null(w9a_err),
   "WP-9 ora: the real L1 composition with a GS-attached dataset boots for the agent round-trip")
ok(ut_cmp_identical(w9a_receipt, PANEL_W9_ORA_AGENT_RECEIPT),
   "WP-9 ora: set_enrichment_parameters receipt byte-identical pre/post migration (all five patch keys applied incl. the run-preconditioned selected_row + the tab switch)")
ok(ut_cmp_identical(w9a_vals, PANEL_W9_ORA_AGENT_VALS),
   "WP-9 ora: the full result-space store state after the agent write byte-identical (zero diff outside the migration surface)")

# ------------------- retention (D3 unification) ----
# structural (deterministic, WP-0 MINOR-2 pattern): the module's ONE
# panel_keep closure holds exactly the helper's bind-sync observer (the
# hand-rolled .ora_keep list held the same single observer; the restore
# observer was never kept and stays untouched; the results tables own
# their retention internally via dataTableDownload_module's .dtd_keep).
# Functional (retention direction, PD-8): after gc() x2 BOTH write
# directions still land - the user cascade pick (UI -> store sync) and
# the agent cascade write (store_apply -> store with a GENUINELY
# DIFFERENT value - store_apply skips same-value keys, so a vacuous
# write would arm nothing and prove nothing - acknowledged through
# the same kept observer once the user input follows; agentWidgets ORA
# pattern).
w9k <- new.env()
w9k_err <- NULL; w9k_nobs <- NULL; w9k_pick <- NULL; w9k_agent <- NULL
w9k_pending <- NULL
shiny::testServer(.w9_mkctx(w9k), {
  w9k_err <<- tryCatch({
    .w9_flush(session)
    NULL
  }, error = function(e) conditionMessage(e))
  obs <- get("obs", envir = environment(w9k$ctx$keep))
  w9k_nobs <<- length(obs)
  gc(); gc()
  # write direction 1: user picks (UI -> store)
  session$setInputs(`ora-tris_ora-analysis` = "General")
  session$flushReact()
  session$setInputs(`ora-tris_ora-subset` = "All")
  session$flushReact()
  session$setInputs(`ora-tris_ora-variable` = "Gene.name")
  session$flushReact(); session$flushReact()
  w9k_pick <<- store_read(w9k$child)[.w9_TRIO]
  # write direction 2: agent write (store_apply -> store + ack). The
  # written value must GENUINELY DIFFER from the user pick ("Gene.name"
  # -> "Protein.ID"): store_apply skips same-value keys (diff against
  # the store), so a vacuous write arms no pending entry and the
  # identical setInputs fires no input event - both assertions would
  # hold trivially and a broken kept observer would go undetected.
  store_apply(w9k$child, list(xax_variable = "Protein.ID"),
              origin = "agent")
  session$flushReact()
  session$setInputs(`ora-tris_ora-variable` = "Protein.ID")
  session$flushReact(); session$flushReact()
  w9k_agent <<- store_read(w9k$child)[.w9_TRIO]
  w9k_pending <<- length(w9k$store$pending)
})
ok(is.null(w9k_err), "WP-9 ora: retention boot stays alive")
ok(ut_cmp_identical(w9k_nobs, 1L),
   "WP-9 ora: exactly the bind-sync observer lands in the panel_keep list (D3; the pre-WP-9 .ora_keep list held the same single observer)")
ok(ut_cmp_identical(unname(unlist(w9k_pick)), c("General", "All", "Gene.name")),
   "WP-9 ora: kept bind observer still fires after gc() x2 (user cascade pick reaches the store)")
ok(ut_cmp_identical(unname(unlist(w9k_agent)), c("General", "All", "Protein.ID")) &&
   ut_cmp_identical(w9k_pending, 0L),
   "WP-9 ora: agent cascade write with a genuinely changed value still lands and acknowledges after gc() x2 (store write direction; pending drains through the kept observer)")

# ------------------- the store-less legacy branch stays alive ----
w9n_err <- NULL
shiny::testServer(.w9_mklegacy(), {
  w9n_err <<- tryCatch({
    .w9_flush(session)
    NULL
  }, error = function(e) conditionMessage(e))
})
ok(is.null(w9n_err),
   "WP-9 ora: the no-store legacy branch (reactiveVal cascade path) still boots")

# ------------------- app level: the REAL L1 composition ----
w9l1 <- new.env()
w9l1$store <- widget_store_new()
app_w9l1 <- function(input, output, session) {
  w9l1$status <- omicsViewer:::L1_result_space_module(
    "rs",
    reactive_expr = reactive(Biobase::exprs(dat)),
    reactive_phenoData = reactive(Biobase::pData(dat)),
    reactive_featureData = reactive(fd0),
    status = reactive(NULL),
    store = w9l1$store)
}
w9l1_err <- NULL; w9l1_payload <- NULL
shiny::testServer(app_w9l1, {
  w9l1_err <<- tryCatch({
    .w9_flush(session, 3)
    NULL
  }, error = function(e) conditionMessage(e))
  w9l1_payload <<- tryCatch(shiny::isolate(w9l1$status()),
                            error = function(e) conditionMessage(e))
})
ok(is.null(w9l1_err),
   "WP-9 ora: the real L1 composition boots headless post-migration (context built at the call site)")
ok(ut_cmp_identical(names(w9l1_payload), PANEL_W9_ORA_GOLDEN_L1_NAMES),
   "WP-9 ora: L1 status payload keeps its exact key order (snapshot schema untouched)")
ok(ut_cmp_identical(w9l1_payload$analyst_ora, PANEL_W9_ORA_GOLDEN_L1_ORA),
   "WP-9 ora: analyst_ora payload byte-identical pre/post migration (cold boot: sanitized list(xax = NULL); the table attr never rides the legacy transport)")
ok(ut_cmp_identical(
  sort(intersect(vapply(store_registry_view(w9l1$store), function(r) r$id, ""),
                 PANEL_W9_ORA_IDS)),
  PANEL_W9_ORA_IDS),
   "WP-9 ora: the real L1 boot carries all four ora keys (the D12 freeze above pins the full 72-id vector)")
## --------------- WP-10: fGSEA migration (D3/D4/D13) ----
# The LAST panel migration, identical shape to WP-9 ("as ORA" per the
# plan appendix): fGSEA is agent-coupled exactly like ORA - the
# set_enrichment_parameters tool (method = "fgsea" through
# agent_normalize_enrichment_update -> apply_agent_enrichment), the
# panel_results excerpt reader (the module return's table attr,
# WP-7/F4) and the L0/L1 wiring all address this panel's store keys,
# NONE of them rewritten here. The migration is the WP-5..WP-9
# pure-migration pattern: retention -> panel_keep (D3), cascade
# boilerplate -> panel_triset_bind (D4), context threading (D13); the
# restore observer and the return-with-table-attr stay in the module
# verbatim (panel semantics).
#
# Golden-fixture provenance (PD-29/PD-33/PD-38/PD-43/PD-48 school):
# all goldens below were captured from PRE-MIGRATION HEAD 40de15e
# BEFORE any WP-10 edit, by booting the then-current
# enrichment_fgsea_module exactly like the L1 call site (legacy
# formals) and dput()ing the outputs: the 4-binding registry view +
# raw meta (incl. agent_writable) + cold values + module return +
# pre-run selected_row semantics over demo fd; the legacy-restore
# triple; the GS-attached run surface (no-run boot - fGSEA req()s a
# settled ranking variable, the deliberate contrast to ORA's
# unset-triple no-collapse boot run; the picked ranking's table attr;
# the live selected_row choices); the set_enrichment_parameters
# round-trip through the L0 apply-callback path over the REAL L1
# composition with a GS-attached dataset (two phases: the cold apply
# rejects the run-preconditioned selected_row with per-key
# resilience, the acknowledged push settles the ranking, runs fGSEA
# and the pathway row applies); the real-L1 analyst_fgsea payload +
# status key order; and the panel's OWN triset derivation over demo
# fd. The IDENTICAL capture script was re-run POST-migration: all
# outputs byte-identical (the pre/post proof). The fixtures are
# FROZEN: never regenerate them from the migrated code.
#
# Two panel-specific notes (recorded as decisions in
# PANEL_DEBT_STATUS.md):
# - The appendix's triset source "trisetter(fd, combine='none') /
#   as ORA" refers to the WP SHAPE, not the literal source: the real
#   module derives its cascade from the NUMERIC non-GS feature
#   annotations (a ranking statistic must be numeric). Swapping in
#   trisetter(combine="none") would list non-numeric and GS columns -
#   a user-visible behavior change the plan forbids (goal 1: zero
#   user-visible change). The helper carries the panel's own
#   derivation verbatim (141x3 over demo fd vs ORA's 713x3).
# - No demo-fd (GS-less) user pick: the fGSEA tab only exists with a
#   GS attr (availability predicate), and a headless pick over GS-less
#   fd drives tab() into its pre-existing no-GS error ("argument is of
#   length zero" at the nrow check - the error predates WP-10 and
#   survives it byte-identically; headless it also destroys the module
#   session, so no golden may be built on it). The GS-attached fd is
#   the realistic pick surface; demo fd covers registration/cold/
#   restore, where no run is attempted.
enrichment_fgsea_module <- omicsViewer:::enrichment_fgsea_module

# -------------------------------------------------- frozen goldens ----
# the 4 ids (3 cascade + the run-preconditioned results-row selection)
PANEL_W10_FGSEA_IDS <- sort(paste0("resultspace.fgsea.",
  c("xax_analysis", "xax_subset", "xax_variable", "selected_row")))

PANEL_W10_FGSEA_GOLDEN_VIEW <- 
list(list(id = "resultspace.fgsea.xax_analysis", kind = "select", 
    label = "Ranking category", help = "Annotation category of the ranking statistic the enrichment scores are computed from", 
    depends_on = character(0), min = NULL, max = NULL, allowed_values = c("General", 
    "PCA", "mean", "n value", "quantile", "ttest", "Rank", "Cor", 
    "ResponseCurve")), list(id = "resultspace.fgsea.xax_subset", 
    kind = "select", label = "Ranking subcategory", help = "Subcategory within the ranking category", 
    depends_on = "resultspace.fgsea.xax_analysis", min = NULL, 
    max = NULL, allowed_values = character(0)), list(id = "resultspace.fgsea.xax_variable", 
    kind = "select_cascaded", label = "Ranking variable", help = "Numeric feature annotation the features are ranked by for the enrichment analysis", 
    depends_on = c("resultspace.fgsea.xax_analysis", "resultspace.fgsea.xax_subset"
    ), min = NULL, max = NULL, allowed_values = c("Intensity", 
    "PC1(10.5%)", "PC2(7.2%)", "PC3(5.4%)", "PC4(4.8%)", "PC5(3.6%)", 
    "PC6(2.9%)", "PC7(2.8%)", "PC8(2.5%)", "RE", "ME", "LE", 
    "MT", "WT", "pvalue", "log.pvalue", "fdr", "log.fdr", "mean.diff", 
    "X786O_NCI60", "A498_NCI60", "A549_NCI60", "ACHN_NCI60", 
    "BT549_NCI60", "CAKI1_NCI60", "CCRFCEM_NCI60", "COLO205_NCI60", 
    "DU145_NCI60", "EKVX_NCI60", "HCC2998_NCI60", "HCT116_NCI60", 
    "HCT15_NCI60", "HL60_NCI60", "HOP62_NCI60", "HOP92_NCI60", 
    "HS578T_NCI60", "HT29_NCI60", "IGROV1_NCI60", "K562_NCI60", 
    "KM12_NCI60", "LOXIMVI_NCI60", "M14_NCI60", "MALME3M_NCI60", 
    "MCF7_NCI60", "MDAMB231_NCI60", "MDAMB435_NCI60", "MDAMB468_NCI60", 
    "MOLT4_NCI60", "NCIADRES_NCI60", "NCIH226_NCI60")), list(
    id = "resultspace.fgsea.selected_row", kind = "select", label = "Selected pathway", 
    help = "Gene set whose row is selected in the results table; highlights its leading edge in the bar plot above", 
    depends_on = character(0), min = NULL, max = NULL, allowed_values = character(0)))

PANEL_W10_FGSEA_GOLDEN_META <- 
list(list(id = "resultspace.fgsea.xax_analysis", kind = "select", 
    label = "Ranking category", help = "Annotation category of the ranking statistic the enrichment scores are computed from", 
    depends_on = character(0), agent_writable = TRUE), list(id = "resultspace.fgsea.xax_subset", 
    kind = "select", label = "Ranking subcategory", help = "Subcategory within the ranking category", 
    depends_on = "resultspace.fgsea.xax_analysis", agent_writable = TRUE), 
    list(id = "resultspace.fgsea.xax_variable", kind = "select_cascaded", 
        label = "Ranking variable", help = "Numeric feature annotation the features are ranked by for the enrichment analysis", 
        depends_on = c("resultspace.fgsea.xax_analysis", "resultspace.fgsea.xax_subset"
        ), agent_writable = TRUE), list(id = "resultspace.fgsea.selected_row", 
        kind = "select", label = "Selected pathway", help = "Gene set whose row is selected in the results table; highlights its leading edge in the bar plot above", 
        depends_on = character(0), agent_writable = TRUE))

PANEL_W10_FGSEA_GOLDEN_COLD <- 
list(resultspace.fgsea.xax_analysis = NULL, resultspace.fgsea.xax_subset = NULL, 
    resultspace.fgsea.xax_variable = NULL, resultspace.fgsea.selected_row = NULL)

PANEL_W10_FGSEA_GOLDEN_ALLOWED_COLD <- 
list(resultspace.fgsea.xax_analysis = c("General", "PCA", "mean", 
"n value", "quantile", "ttest", "Rank", "Cor", "ResponseCurve"
), resultspace.fgsea.xax_subset = character(0), resultspace.fgsea.xax_variable = c("Intensity", 
"PC1(10.5%)", "PC2(7.2%)", "PC3(5.4%)", "PC4(4.8%)", "PC5(3.6%)", 
"PC6(2.9%)", "PC7(2.8%)", "PC8(2.5%)", "RE", "ME", "LE", "MT", 
"WT", "pvalue", "log.pvalue", "fdr", "log.fdr", "mean.diff", 
"X786O_NCI60", "A498_NCI60", "A549_NCI60", "ACHN_NCI60", "BT549_NCI60", 
"CAKI1_NCI60", "CCRFCEM_NCI60", "COLO205_NCI60", "DU145_NCI60", 
"EKVX_NCI60", "HCC2998_NCI60", "HCT116_NCI60", "HCT15_NCI60", 
"HL60_NCI60", "HOP62_NCI60", "HOP92_NCI60", "HS578T_NCI60", "HT29_NCI60", 
"IGROV1_NCI60", "K562_NCI60", "KM12_NCI60", "LOXIMVI_NCI60", 
"M14_NCI60", "MALME3M_NCI60", "MCF7_NCI60", "MDAMB231_NCI60", 
"MDAMB435_NCI60", "MDAMB468_NCI60", "MOLT4_NCI60", "NCIADRES_NCI60", 
"NCIH226_NCI60"), resultspace.fgsea.selected_row = character(0))

PANEL_W10_FGSEA_GOLDEN_RET <- 
list(xax = NULL, table = NULL)

PANEL_W10_FGSEA_GOLDEN_RESTORE_VALS <- 
list(resultspace.fgsea.xax_analysis = "ttest", resultspace.fgsea.xax_subset = "RE_vs_ME", 
    resultspace.fgsea.xax_variable = "log.fdr", resultspace.fgsea.selected_row = NULL)

PANEL_W10_FGSEA_GOLDEN_RESTORE_RET <- 
list(xax = NULL, table = NULL)

PANEL_W10_FGSEA_TRISET <- 
list(dim = c(141L, 3L), analysis_levels = c("General", "PCA", 
"mean", "n value", "quantile", "ttest", "Rank", "Cor", "ResponseCurve"
), first_rows = structure(c("General", "PCA", "PCA", "PCA", "PCA", 
"PCA", "All", "All", "All", "All", "All", "All", "Intensity", 
"PC1(10.5%)", "PC2(7.2%)", "PC3(5.4%)", "PC4(4.8%)", "PC5(3.6%)"
), dim = c(6L, 3L)), variable_head50 = c("Intensity", "PC1(10.5%)", 
"PC2(7.2%)", "PC3(5.4%)", "PC4(4.8%)", "PC5(3.6%)", "PC6(2.9%)", 
"PC7(2.8%)", "PC8(2.5%)", "RE", "ME", "LE", "MT", "WT", "pvalue", 
"log.pvalue", "fdr", "log.fdr", "mean.diff", "X786O_NCI60", "A498_NCI60", 
"A549_NCI60", "ACHN_NCI60", "BT549_NCI60", "CAKI1_NCI60", "CCRFCEM_NCI60", 
"COLO205_NCI60", "DU145_NCI60", "EKVX_NCI60", "HCC2998_NCI60", 
"HCT116_NCI60", "HCT15_NCI60", "HL60_NCI60", "HOP62_NCI60", "HOP92_NCI60", 
"HS578T_NCI60", "HT29_NCI60", "IGROV1_NCI60", "K562_NCI60", "KM12_NCI60", 
"LOXIMVI_NCI60", "M14_NCI60", "MALME3M_NCI60", "MCF7_NCI60", 
"MDAMB231_NCI60", "MDAMB435_NCI60", "MDAMB468_NCI60", "MOLT4_NCI60", 
"NCIADRES_NCI60", "NCIH226_NCI60"))

PANEL_W10_FGSEA_AGENT_RECEIPT <- 
list(method = "fgsea", panel_tab = "fGSEA", applied = c("resultspace.fgsea.xax_analysis", 
"resultspace.analyst_tab", "resultspace.fgsea.xax_subset", "resultspace.fgsea.xax_variable"
), applied_values = list(resultspace.fgsea.xax_analysis = "ttest", 
    resultspace.analyst_tab = "fGSEA", resultspace.fgsea.xax_subset = "A_vs_B", 
    resultspace.fgsea.xax_variable = "stat"), unchanged = character(0), 
    rejected = list(list(id = "resultspace.fgsea.selected_row", 
        reason = "Choices are currently unavailable for resultspace.fgsea.selected_row: its choices are not computed yet (the owning panel may be inactive, its data not ready, or a cascade parent unset). Open or activate the owning panel, then query the control again for its current allowed values.")), 
    note = "The panel opens and the enrichment recomputes from the currently selected features (ORA) or the full ranking (fGSEA); read the top rows back through get_omics_viewer_state(sections='panel_results').")

PANEL_W10_FGSEA_AGENT_RECEIPT2 <- 
list(method = "fgsea", panel_tab = "fGSEA", applied = "resultspace.fgsea.selected_row", 
    applied_values = list(resultspace.fgsea.selected_row = "GS_A"), 
    unchanged = "resultspace.analyst_tab", rejected = list(), 
    note = "The panel opens and the enrichment recomputes from the currently selected features (ORA) or the full ranking (fGSEA); read the top rows back through get_omics_viewer_state(sections='panel_results').")

PANEL_W10_FGSEA_AGENT_VALS <- 
list(resultspace.analyst_tab = "fGSEA", resultspace.feature_general.xax_analysis = NULL, 
    resultspace.feature_general.xax_subset = NULL, resultspace.feature_general.xax_variable = NULL, 
    resultspace.feature_general.plot_type = NULL, resultspace.feature_general.regression_line = FALSE, 
    resultspace.feature_general.attr4.color_analysis = NULL, 
    resultspace.feature_general.attr4.color_subset = NULL, resultspace.feature_general.attr4.color_variable = NULL, 
    resultspace.feature_general.attr4.shape_analysis = NULL, 
    resultspace.feature_general.attr4.shape_subset = NULL, resultspace.feature_general.attr4.shape_variable = NULL, 
    resultspace.feature_general.attr4.size_analysis = NULL, resultspace.feature_general.attr4.size_subset = NULL, 
    resultspace.feature_general.attr4.size_variable = NULL, resultspace.feature_general.attr4.tooltip_analysis = NULL, 
    resultspace.feature_general.attr4.tooltip_subset = NULL, 
    resultspace.feature_general.attr4.tooltip_variable = NULL, 
    resultspace.feature_general.attr4.search_analysis = NULL, 
    resultspace.feature_general.attr4.search_subset = NULL, resultspace.feature_general.attr4.search_variable = NULL, 
    resultspace.feature_general.attr4.search_value = character(0), 
    resultspace.feature_general.attr4.xcut = NULL, resultspace.feature_general.attr4.ycut = NULL, 
    resultspace.feature_general.attr4.scorner = "None", resultspace.fgsea.xax_analysis = "ttest", 
    resultspace.fgsea.xax_subset = "A_vs_B", resultspace.fgsea.xax_variable = "stat", 
    resultspace.fgsea.selected_row = "GS_A", resultspace.ora.xax_analysis = NULL, 
    resultspace.ora.xax_subset = NULL, resultspace.ora.xax_variable = NULL, 
    resultspace.ora.selected_row = NULL, resultspace.stringdb.taxonomy = NULL, 
    resultspace.stringdb.show_labels = NULL, resultspace.stringdb.selected_row = NULL, 
    resultspace.sample_general.xax_analysis = NULL, resultspace.sample_general.xax_subset = NULL, 
    resultspace.sample_general.xax_variable = NULL, resultspace.sample_general.attr4.color_analysis = NULL, 
    resultspace.sample_general.attr4.color_subset = NULL, resultspace.sample_general.attr4.color_variable = NULL, 
    resultspace.sample_general.attr4.shape_analysis = NULL, resultspace.sample_general.attr4.shape_subset = NULL, 
    resultspace.sample_general.attr4.shape_variable = NULL, resultspace.sample_general.attr4.size_analysis = NULL, 
    resultspace.sample_general.attr4.size_subset = NULL, resultspace.sample_general.attr4.size_variable = NULL, 
    resultspace.sample_general.attr4.tooltip_analysis = NULL, 
    resultspace.sample_general.attr4.tooltip_subset = NULL, resultspace.sample_general.attr4.tooltip_variable = NULL, 
    resultspace.sample_general.attr4.search_analysis = NULL, 
    resultspace.sample_general.attr4.search_subset = NULL, resultspace.sample_general.attr4.search_variable = NULL, 
    resultspace.sample_general.attr4.search_value = character(0), 
    resultspace.sample_general.attr4.xcut = NULL, resultspace.sample_general.attr4.ycut = NULL, 
    resultspace.sample_general.attr4.scorner = "None", resultspace.sample_general.survival_censor = NULL, 
    resultspace.sample_general.batch_show_phenotype = NULL, resultspace.sample_general.batch_show_features = NULL, 
    resultspace.sample_general.batch_phenotype_selected_row = NULL, 
    resultspace.sample_general.batch_feature_selected_row = NULL, 
    resultspace.geneshot.term = NULL, resultspace.geneshot.xax_analysis = NULL, 
    resultspace.geneshot.xax_subset = NULL, resultspace.geneshot.xax_variable = NULL, 
    resultspace.ptm.xax_analysis = NULL, resultspace.ptm.xax_subset = NULL, 
    resultspace.ptm.xax_variable = NULL, resultspace.dose_response.feature_table = NULL, 
    resultspace.dose_response.param_table = NULL)

PANEL_W10_FGSEA_GS_BOOT_RET <- 
list(xax = NULL, table = NULL)

PANEL_W10_FGSEA_GS_VIEW <- 
list(list(id = "resultspace.fgsea.xax_analysis", kind = "select", 
    label = "Ranking category", help = "Annotation category of the ranking statistic the enrichment scores are computed from", 
    depends_on = character(0), min = NULL, max = NULL, allowed_values = "ttest"), 
    list(id = "resultspace.fgsea.xax_subset", kind = "select", 
        label = "Ranking subcategory", help = "Subcategory within the ranking category", 
        depends_on = "resultspace.fgsea.xax_analysis", min = NULL, 
        max = NULL, allowed_values = character(0)), list(id = "resultspace.fgsea.xax_variable", 
        kind = "select_cascaded", label = "Ranking variable", 
        help = "Numeric feature annotation the features are ranked by for the enrichment analysis", 
        depends_on = c("resultspace.fgsea.xax_analysis", "resultspace.fgsea.xax_subset"
        ), min = NULL, max = NULL, allowed_values = c("pval", 
        "stat")), list(id = "resultspace.fgsea.selected_row", 
        kind = "select", label = "Selected pathway", help = "Gene set whose row is selected in the results table; highlights its leading edge in the bar plot above", 
        depends_on = character(0), min = NULL, max = NULL, allowed_values = character(0)))

PANEL_W10_FGSEA_GS_PICK_VALS <- 
list(resultspace.fgsea.xax_analysis = "ttest", resultspace.fgsea.xax_subset = "A_vs_B", 
    resultspace.fgsea.xax_variable = "stat", resultspace.fgsea.selected_row = NULL)

PANEL_W10_FGSEA_GS_PICK_RET <- 
list(xax = list(analysis = "ttest", subset = "A_vs_B", variable = "stat"), 
    table = list(class = "data.frame", names = c("pathway", "pval", 
    "padj", "log2err", "enrichment score (ES)", "normalized ES", 
    "size", "leadingEdge", "description"), nrow = 4L, pathways = c("GS_A", 
    "GS_D", "GS_B", "GS_C")))

PANEL_W10_FGSEA_GS_PICK_ALLOWED <- 
list(resultspace.fgsea.xax_analysis = "ttest", resultspace.fgsea.xax_subset = "A_vs_B", 
    resultspace.fgsea.xax_variable = c("pval", "stat"), resultspace.fgsea.selected_row = c("GS_A", 
    "GS_D", "GS_B", "GS_C"))

PANEL_W10_FGSEA_GOLDEN_L1_NAMES <- 
c("analyst_active_tab", "analyst_feature_general", "analyst_sample_general", 
"analyst_gene_shot", "analyst_fgsea", "analyst_ora", "analyst_ptm", 
"analyst_dose_response", "analyst_additional")

PANEL_W10_FGSEA_GOLDEN_L1_FGSEA <- 
list(xax = NULL)


# -------------------------------------------------- boot factories ----
# GS-attached synthetic fd (the agentWidgets S4 fgsea shape, two
# NUMERIC ranking columns - the retention case needs a GENUINELY
# DIFFERENT valid agent write target in the same cascade branch,
# stat -> pval; deterministic under set.seed)
set.seed(10)
.w10_gs_fd <- data.frame(
  `ttest|A_vs_B|pval` = c(rnorm(10), rnorm(20)),
  `ttest|A_vs_B|stat` = c(rnorm(10, mean = 2), rnorm(20)),
  row.names = paste0("F", 1:30), check.names = FALSE)
attr(.w10_gs_fd, "GS") <- data.frame(
  featureId = factor(paste0("F", c(1:10, 6:20, 21:30, 1:5, 25:28))),
  gsId = factor(rep(c("GS_A", "GS_B", "GS_C", "GS_D"),
                    c(10, 15, 10, 9))))
.w10_pd <- data.frame(group = c("A", "A", "B", "B"),
                     row.names = paste0("S", 1:4))
.w10_ex <- matrix(rnorm(120), nrow = 30,
                 dimnames = list(paste0("F", 1:30), paste0("S", 1:4)))
.w10_sel <- paste0("F", 1:5)
# context-built (the WP-10 L1 call-site shape; holder receives the
# context so the keep closure is inspectable)
.w10_mkctx <- function(holder, fd = fd0, status_payload = NULL) {
  force(holder)
  function(input, output, session) {
    .fd <- reactive(fd); .st <- reactive(status_payload)
    holder$store <- widget_store_new()
    holder$child <- widget_store_child(holder$store, "resultspace.fgsea")
    holder$ctx <- panel_context(
      id = "fgsea", tab_value = "fGSEA", store = holder$child,
      data = list(fdata = .fd), status = .st, session = session)
    holder$ret <- enrichment_fgsea_module(
      "fgsea", reactive_featureData = .fd,
      reactive_status = .st, store = holder$child, context = holder$ctx)
  }
}
# legacy direct boot (formals only - the pre-WP-10 shape, still supported)
.w10_boot <- function(env, fd = fd0, status_payload = NULL) {
  force(env)
  function(input, output, session) {
    env$store <- widget_store_new()
    env$child <- widget_store_child(env$store, "resultspace.fgsea")
    env$ret <- enrichment_fgsea_module(
      "fgsea", reactive_featureData = reactive(fd),
      reactive_status = reactive(status_payload), store = env$child)
  }
}
# store-less legacy boot (no store, no context - the reactiveVal path)
.w10_mklegacy <- function(fd = fd0, status_payload = NULL) {
  function(input, output, session) {
    enrichment_fgsea_module(
      "fgsea", reactive_featureData = reactive(fd),
      reactive_status = reactive(status_payload))
  }
}
.w10_meta <- function(root) {
  ks <- grep("^resultspace\\.fgsea(\\.|$)", names(root$bindings), value = TRUE)
  lapply(ks, function(k) {
    b <- root$bindings[[k]]
    list(id = k, kind = b$kind, label = b$label, help = b$help,
         depends_on = b$depends_on, agent_writable = b$agent_writable)
  })
}
.w10_allowed <- function(root) {
  v <- store_registry_view(root, "resultspace.fgsea")
  stats::setNames(lapply(v, `[[`, "allowed_values"),
                  vapply(v, function(r) r$id, ""))
}
# the module return as the status builder / excerpt reader see it (the
# table attr summarized - content pins, not the whole data frame)
.w10_tab_sum <- function(tab) {
  if (!is.data.frame(tab)) return(tab)
  list(class = class(tab)[[1]], names = colnames(tab), nrow = nrow(tab),
       pathways = as.character(tab$pathway))
}
.w10_ret_shape <- function(ret) {
  r <- tryCatch(shiny::isolate(ret()), error = function(e) NULL)
  list(xax = r$xax, table = .w10_tab_sum(attr(r, "table")))
}
.w10_flush <- function(session, n = 6)
  for (k in seq_len(n)) tryCatch(session$flushReact(), error = function(e) NULL)
.w10_TRIO <- paste0("resultspace.fgsea.",
                    c("xax_analysis", "xax_subset", "xax_variable"))

# ------------------- registration: context boot = the frozen golden ----
w10 <- new.env()
w10_err <- NULL; w10_view <- NULL; w10_meta <- NULL
shiny::testServer(.w10_mkctx(w10), {
  w10_err <<- tryCatch({
    .w10_flush(session)
    NULL
  }, error = function(e) conditionMessage(e))
  # registry view INSIDE the session (WP-3 lore: choices providers read
  # live reactives; outside a reactive context they degrade to empty)
  w10_view <<- store_registry_view(w10$store, "resultspace.fgsea")
  w10_meta <<- .w10_meta(w10$store)
})
ok(is.null(w10_err),
   "WP-10 fgsea: module boots headless context-built (the new L1 call-site shape)")
ok(ut_cmp_identical(sort(vapply(w10_view, function(r) r$id, "")),
                    PANEL_W10_FGSEA_IDS),
   "WP-10 fgsea: exactly the four frozen bindings register (3 cascade + run-preconditioned selected_row; D12 untouched)")
ok(ut_cmp_identical(w10_view, PANEL_W10_FGSEA_GOLDEN_VIEW),
   "WP-10 fgsea: registry view byte-identical to the frozen pre-migration golden (over demo fd: the 9 NUMERIC analysis levels - the ranking filter; full variable order)")
ok(ut_cmp_identical(w10_meta, PANEL_W10_FGSEA_GOLDEN_META),
   "WP-10 fgsea: raw binding meta byte-identical incl. agent_writable TRUE x4 (the set_enrichment_parameters surface)")

# ------------------- cold values + module return + pre-run table attr ----
w10c <- new.env()
w10c_err <- NULL; w10c_cold <- NULL; w10c_ret <- NULL; w10c_allow <- NULL
shiny::testServer(.w10_mkctx(w10c), {
  w10c_err <<- tryCatch({
    .w10_flush(session)
    NULL
  }, error = function(e) conditionMessage(e))
  w10c_cold <<- store_read(w10c$store)
  w10c_ret <<- .w10_ret_shape(w10c$ret)
  w10c_allow <<- .w10_allowed(w10c$store)
})
ok(is.null(w10c_err), "WP-10 fgsea: cold-boot capture stays alive")
ok(ut_cmp_identical(w10c_cold, PANEL_W10_FGSEA_GOLDEN_COLD),
   "WP-10 fgsea: cold store values byte-identical (cascade + selected_row all NULL without a run)")
ok(ut_cmp_identical(w10c_ret, PANEL_W10_FGSEA_GOLDEN_RET),
   "WP-10 fgsea: module return (the analyst_fgsea payload source + panel_results excerpt source) byte-identical: list(xax = NULL), table attr NULL pre-run")
ok(ut_cmp_identical(w10c_allow, PANEL_W10_FGSEA_GOLDEN_ALLOWED_COLD),
   "WP-10 fgsea: pre-run selected_row choices are empty (run-preconditioned semantics; only tested pathways become rows)")

# ------------------- legacy status restore path ----
w10r <- new.env()
w10r_err <- NULL; w10r_vals <- NULL; w10r_ret <- NULL
shiny::testServer(.w10_boot(w10r, status_payload = list(
  xax = c("ttest", "RE_vs_ME", "log.fdr"))), {
  w10r_err <<- tryCatch({
    .w10_flush(session)
    NULL
  }, error = function(e) conditionMessage(e))
  w10r_vals <<- store_read(w10r$store)
  w10r_ret <<- .w10_ret_shape(w10r$ret)
})
ok(is.null(w10r_err),
   "WP-10 fgsea: boot with a legacy analyst_fgsea payload stays alive")
ok(ut_cmp_identical(w10r_vals, PANEL_W10_FGSEA_GOLDEN_RESTORE_VALS),
   "WP-10 fgsea: the legacy restore path re-applies the triple through the store (observer untouched, 4 keys byte-identical)")
ok(ut_cmp_identical(w10r_ret, PANEL_W10_FGSEA_GOLDEN_RESTORE_RET),
   "WP-10 fgsea: module return stays cold after a restore-only boot (the triple lands in the store; v1 has not been user-driven)")

# ------------------- legacy boot registers identically (D13) ----
w10l <- new.env()
w10l_err <- NULL; w10l_view <- NULL; w10l_meta <- NULL; w10l_cold <- NULL
shiny::testServer(.w10_boot(w10l), {
  w10l_err <<- tryCatch({
    .w10_flush(session)
    NULL
  }, error = function(e) conditionMessage(e))
  w10l_view <<- store_registry_view(w10l$store, "resultspace.fgsea")
  w10l_meta <<- .w10_meta(w10l$store)
  w10l_cold <<- store_read(w10l$store)
})
ok(is.null(w10l_err),
   "WP-10 fgsea: legacy direct instantiation (store formal only, no context) still boots")
ok(ut_cmp_identical(w10l_view, PANEL_W10_FGSEA_GOLDEN_VIEW) &&
   ut_cmp_identical(w10l_meta, PANEL_W10_FGSEA_GOLDEN_META) &&
   ut_cmp_identical(w10l_cold, PANEL_W10_FGSEA_GOLDEN_COLD),
   "WP-10 fgsea: legacy store boot registers and seeds byte-identically to the golden (both boot paths agree)")
ok(ut_cmp_identical(w10l_view, w10_view),
   "WP-10 fgsea: context boot and legacy boot register identically (D13 argument coalescing)")

# ------------------- the triset derivation (the panel's OWN filter) ----
cn10 <- colnames(fd0)[vapply(fd0, is.numeric, logical(1)) &
                        !grepl("^GS\\|", colnames(fd0))]
tsd10 <- stringr::str_split_fixed(cn10, "\\|", n = 3)
ok(ut_cmp_identical(list(dim = dim(tsd10),
                         analysis_levels = unique(tsd10[, 1]),
                         first_rows = utils::head(tsd10, 6L),
                         variable_head50 = utils::head(unique(tsd10[, 3]), 50)),
                    list(dim = PANEL_W10_FGSEA_TRISET$dim,
                         analysis_levels = PANEL_W10_FGSEA_TRISET$analysis_levels,
                         first_rows = PANEL_W10_FGSEA_TRISET$first_rows,
                         variable_head50 = PANEL_W10_FGSEA_TRISET$variable_head50)),
   "WP-10 fgsea: the cascade triset source over demo fd = the numeric non-GS column inventory (141x3; NOT trisetter(combine=none)'s 713x3 - a ranking statistic must be numeric; the appendix's 'as ORA' is the WP shape, not the source)")

# ------------------- GS dataset: the run / table-attr surface ----
# module-direct boot (legacy formals) with a real GS attr: fGSEA does
# NOT run at boot (tab() req()s a settled ranking variable - the
# deliberate contrast to ORA's unset-triple no-collapse boot run); the
# ranking pick runs fGSEA and the table attr + selected_row choices go
# live
w10g <- new.env()
w10g_err <- NULL; w10g_boot <- NULL; w10g_view <- NULL
w10g_vals <- NULL; w10g_ret <- NULL; w10g_allow <- NULL
shiny::testServer(.w10_boot(w10g, fd = .w10_gs_fd), {
  w10g_err <<- tryCatch({
    .w10_flush(session, 8)
    NULL
  }, error = function(e) conditionMessage(e))
  w10g_boot <<- .w10_ret_shape(w10g$ret)
  w10g_view <<- store_registry_view(w10g$store, "resultspace.fgsea")
  session$setInputs(`fgsea-tris_fgsea-analysis` = "ttest")
  session$flushReact()
  session$setInputs(`fgsea-tris_fgsea-subset` = "A_vs_B")
  session$flushReact()
  session$setInputs(`fgsea-tris_fgsea-variable` = "stat")
  .w10_flush(session, 8)
  w10g_vals <<- store_read(w10g$store)
  w10g_ret <<- .w10_ret_shape(w10g$ret)
  w10g_allow <<- .w10_allowed(w10g$store)
})
ok(is.null(w10g_err), "WP-10 fgsea: GS-attached run boot stays alive")
ok(ut_cmp_identical(w10g_boot, PANEL_W10_FGSEA_GS_BOOT_RET),
   "WP-10 fgsea: the unset-triple boot does NOT run (table attr NULL until a ranking settles - fGSEA's run-preconditioned semantics, the contrast to ORA's boot run)")
ok(ut_cmp_identical(w10g_view, PANEL_W10_FGSEA_GS_VIEW),
   "WP-10 fgsea: GS registry view byte-identical (single ttest branch; both ranking variables; selected_row choices cold at boot)")
ok(ut_cmp_identical(w10g_vals, PANEL_W10_FGSEA_GS_PICK_VALS),
   "WP-10 fgsea: the ranking pick keeps the settled triple + NULL selected_row (byte-identical)")
ok(ut_cmp_identical(w10g_ret, PANEL_W10_FGSEA_GS_PICK_RET),
   "WP-10 fgsea: the ranking pick's table attr is the real fGSEA table (panel_results excerpt source; NES-ordered pathways pinned)")
ok(ut_cmp_identical(w10g_allow, PANEL_W10_FGSEA_GS_PICK_ALLOWED),
   "WP-10 fgsea: post-run selected_row choices are the tested gene sets (run-preconditioned semantics)")

# ------------------- agent round-trip: the L0 apply-callback path ----
# set_enrichment_parameters(method = "fgsea") lands through
# apply_agent_enrichment's exact body (validate -> store_apply on the
# ROOT store -> receipt, note verbatim) over the REAL L1 composition
# booted with the GS-attached dataset. TWO phases, both frozen from
# the pre-migration capture: (1) the cold full-five-key apply - the
# ranking has not settled, so the run-preconditioned selected_row
# REJECTS with per-key resilience (the L0 comment's own example),
# exactly like ORA's cold-path rejection semantics; (2) the store->UI
# push acknowledged headless (setInputs with the pushed triple - the
# retention-ack pattern; a real browser applies the input messages),
# which settles the ranking, runs fGSEA and the pathway row APPLIES.
w10a <- new.env()
w10a$store <- widget_store_new()
app_w10a <- function(input, output, session) {
  w10a$status <- omicsViewer:::L1_result_space_module(
    "rs",
    reactive_expr = reactive(.w10_ex),
    reactive_phenoData = reactive(.w10_pd),
    reactive_featureData = reactive(.w10_gs_fd),
    reactive_i = reactive(.w10_sel),
    reactive_highlight = reactive(NULL),
    status = reactive(NULL),
    store = w10a$store)
}
w10a_err <- NULL; w10a_receipt <- NULL; w10a_receipt2 <- NULL; w10a_vals <- NULL
shiny::testServer(app_w10a, {
  w10a_err <<- tryCatch({
    .w10_flush(session, 8)
    update <- list(method = "fgsea", collapse = "ttest|A_vs_B|stat",
                   selected_pathway = "GS_A")
    validated <- agent_normalize_enrichment_update(update, .w10_gs_fd)
    receipt <- shiny::isolate(
      store_apply(w10a$store, validated$patch, origin = "agent",
                  strict = FALSE))
    w10a_receipt <<- list(
      method = validated$method,
      panel_tab = validated$tab,
      applied = receipt$applied,
      applied_values = receipt$diff,
      unchanged = receipt$skipped,
      rejected = receipt$rejected %||% list(),
      note = paste(
        "The panel opens and the enrichment recomputes from the currently",
        "selected features (ORA) or the full ranking (fGSEA); read the top",
        "rows back through get_omics_viewer_state(sections='panel_results')."))
    .w10_flush(session, 8)
    # acknowledge the pushed triple headless (the browser's input events)
    session$setInputs(`rs-fgsea-tris_fgsea-analysis` = "ttest")
    session$flushReact()
    session$setInputs(`rs-fgsea-tris_fgsea-subset` = "A_vs_B")
    session$flushReact()
    session$setInputs(`rs-fgsea-tris_fgsea-variable` = "stat")
    .w10_flush(session, 8)
    update2 <- list(method = "fgsea", selected_pathway = "GS_A")
    validated2 <- agent_normalize_enrichment_update(update2, .w10_gs_fd)
    receipt2 <- shiny::isolate(
      store_apply(w10a$store, validated2$patch, origin = "agent",
                  strict = FALSE))
    w10a_receipt2 <<- list(
      method = validated2$method,
      panel_tab = validated2$tab,
      applied = receipt2$applied,
      applied_values = receipt2$diff,
      unchanged = receipt2$skipped,
      rejected = receipt2$rejected %||% list(),
      note = paste(
        "The panel opens and the enrichment recomputes from the currently",
        "selected features (ORA) or the full ranking (fGSEA); read the top",
        "rows back through get_omics_viewer_state(sections='panel_results')."))
    .w10_flush(session, 8)
    NULL
  }, error = function(e) conditionMessage(e))
  w10a_vals <<- store_read(w10a$store)
})
ok(is.null(w10a_err),
   "WP-10 fgsea: the real L1 composition with a GS-attached dataset boots for the agent round-trip")
ok(ut_cmp_identical(w10a_receipt, PANEL_W10_FGSEA_AGENT_RECEIPT),
   "WP-10 fgsea: cold set_enrichment_parameters receipt byte-identical pre/post migration (the four patch keys apply incl. the tab switch; the run-preconditioned selected_row rejects with the suggestion-bearing reason)")
ok(ut_cmp_identical(w10a_receipt2, PANEL_W10_FGSEA_AGENT_RECEIPT2),
   "WP-10 fgsea: after the acknowledged push settles the ranking and fGSEA runs, the pathway row applies (receipt 2 byte-identical; the agent-coupling proof, PD-49 school)")
ok(ut_cmp_identical(w10a_vals, PANEL_W10_FGSEA_AGENT_VALS),
   "WP-10 fgsea: the full result-space store state after both writes byte-identical (zero diff outside the migration surface)")

# ------------------- retention (D3 unification) ----
# structural (deterministic, WP-0 MINOR-2 pattern): the module's ONE
# panel_keep closure holds exactly TWO observers - the helper's
# bind-sync observer and the fgs_tab visibility-gated mirror (the
# pre-WP-10 hand-rolled .fgs_keep list held the same two via its two
# keep() call sites; the restore observer was never kept and stays
# untouched; the results table owns its retention internally via
# dataTableDownload_module's .dtd_keep).
# Functional (retention direction, PD-8): after gc() x2 BOTH write
# directions still land - the user ranking pick (UI -> store sync)
# and the agent cascade write (store_apply -> store with a GENUINELY
# DIFFERENT value - store_apply skips same-value keys, so a vacuous
# write would arm nothing and prove nothing; the WP-9 fix-round shape,
# 40de15e): pending ARMED (1L) after the apply, then acknowledged
# through the same kept observer once the user input follows, and
# DRAINED (0L).
w10k <- new.env()
w10k_err <- NULL; w10k_nobs <- NULL; w10k_pick <- NULL; w10k_agent <- NULL
w10k_armed <- NULL; w10k_preack <- NULL; w10k_pending <- NULL
shiny::testServer(.w10_mkctx(w10k, fd = .w10_gs_fd), {
  w10k_err <<- tryCatch({
    .w10_flush(session)
    NULL
  }, error = function(e) conditionMessage(e))
  obs <- get("obs", envir = environment(w10k$ctx$keep))
  w10k_nobs <<- length(obs)
  gc(); gc()
  # write direction 1: user picks (UI -> store)
  session$setInputs(`fgsea-tris_fgsea-analysis` = "ttest")
  session$flushReact()
  session$setInputs(`fgsea-tris_fgsea-subset` = "A_vs_B")
  session$flushReact()
  session$setInputs(`fgsea-tris_fgsea-variable` = "stat")
  session$flushReact(); session$flushReact()
  w10k_pick <<- store_read(w10k$child)[.w10_TRIO]
  # write direction 2: agent write (store_apply -> store + ack). The
  # written value must GENUINELY DIFFER from the user pick ("stat" ->
  # "pval", the second ranking column of .w10_gs_fd): store_apply
  # skips same-value keys (diff against the store), so a vacuous write
  # arms no pending entry and the identical setInputs fires no input
  # event - both assertions would hold trivially and a broken kept
  # observer would go undetected.
  store_apply(w10k$child, list(xax_variable = "pval"),
              origin = "agent")
  session$flushReact()
  w10k_armed <<- length(w10k$store$pending)
  w10k_preack <<- store_read(w10k$child)$resultspace.fgsea.xax_variable
  session$setInputs(`fgsea-tris_fgsea-variable` = "pval")
  session$flushReact(); session$flushReact()
  w10k_agent <<- store_read(w10k$child)[.w10_TRIO]
  w10k_pending <<- length(w10k$store$pending)
})
ok(is.null(w10k_err), "WP-10 fgsea: retention boot stays alive")
ok(ut_cmp_identical(w10k_nobs, 2L),
   "WP-10 fgsea: exactly the bind-sync + fgs_tab-mirror observers land in the panel_keep list (D3; the pre-WP-10 .fgs_keep list held the same two)")
ok(ut_cmp_identical(unname(unlist(w10k_pick)), c("ttest", "A_vs_B", "stat")),
   "WP-10 fgsea: kept bind observer still fires after gc() x2 (user ranking pick reaches the store)")
ok(ut_cmp_identical(w10k_armed, 1L) &&
   ut_cmp_identical(w10k_preack, "pval"),
   "WP-10 fgsea: the agent write with a genuinely changed value arms pending (1) and lands in the store pre-ack (non-vacuous write)")
ok(ut_cmp_identical(unname(unlist(w10k_agent)), c("ttest", "A_vs_B", "pval")) &&
   ut_cmp_identical(w10k_pending, 0L),
   "WP-10 fgsea: the agent ranking write still acknowledges after gc() x2 (store write direction; pending armed then drained through the kept observer)")

# ------------------- the store-less legacy branch stays alive ----
w10n_err <- NULL
shiny::testServer(.w10_mklegacy(), {
  w10n_err <<- tryCatch({
    .w10_flush(session)
    NULL
  }, error = function(e) conditionMessage(e))
})
ok(is.null(w10n_err),
   "WP-10 fgsea: the no-store legacy branch (reactiveVal cascade path) still boots")

# ------------------- app level: the REAL L1 composition ----
w10l1 <- new.env()
w10l1$store <- widget_store_new()
app_w10l1 <- function(input, output, session) {
  w10l1$status <- omicsViewer:::L1_result_space_module(
    "rs",
    reactive_expr = reactive(Biobase::exprs(dat)),
    reactive_phenoData = reactive(Biobase::pData(dat)),
    reactive_featureData = reactive(fd0),
    status = reactive(NULL),
    store = w10l1$store)
}
w10l1_err <- NULL; w10l1_payload <- NULL
shiny::testServer(app_w10l1, {
  w10l1_err <<- tryCatch({
    .w10_flush(session, 3)
    NULL
  }, error = function(e) conditionMessage(e))
  w10l1_payload <<- tryCatch(shiny::isolate(w10l1$status()),
                            error = function(e) conditionMessage(e))
})
ok(is.null(w10l1_err),
   "WP-10 fgsea: the real L1 composition boots headless post-migration (context built at the call site)")
ok(ut_cmp_identical(names(w10l1_payload), PANEL_W10_FGSEA_GOLDEN_L1_NAMES),
   "WP-10 fgsea: L1 status payload keeps its exact key order (snapshot schema untouched)")
ok(ut_cmp_identical(w10l1_payload$analyst_fgsea, PANEL_W10_FGSEA_GOLDEN_L1_FGSEA),
   "WP-10 fgsea: analyst_fgsea payload byte-identical pre/post migration (cold boot: sanitized list(xax = NULL); the table attr never rides the legacy transport)")
ok(ut_cmp_identical(
  sort(intersect(vapply(store_registry_view(w10l1$store), function(r) r$id, ""),
                 PANEL_W10_FGSEA_IDS)),
  PANEL_W10_FGSEA_IDS),
   "WP-10 fgsea: the real L1 boot carries all four fgsea keys (the D12 freeze above pins the full 72-id vector)")

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


## ------------------- WP-11: L1 server composition from the registry ----
# The structural centerpiece: L1's SERVER-side composition collapses from
# the eight hand-wired per-panel call sites (WP-3..WP-10 context builds +
# bespoke formals + the seven hardcoded widget_store_child lines) into
# ONE loop over panel_records (per record: child store, panel_context,
# the record's server with its module id and legacy formals; results
# collected into a map keyed by status_key). Zero behavior change by
# construction and by proof.
#
# GOLDEN PROVENANCE (strict, the PD-29/PD-33/PD-38/PD-43/PD-48/PD-55
# school): every literal below was captured from PRE-WP-11 HEAD b0f920a
# BEFORE any edit (2026-11-06), by booting the real L1 composition
# exactly like the freeze harness above (demo.RDS fd0, fresh store, 3
# flushes) and dput()ing the outputs: the FULL store snapshot
# (values+epochs+actions, names AND order), the store_registry_view ids
# in registration order, the status payload (reader attr stripped for
# the literal; its behavior is its own golden), the agent state overview
# shape (agent_compact_state over the real payload + store), the
# capability-record ids in order with their per-tab counts, and the
# enrichment reader's cold behavior. The fixtures are FROZEN: never
# regenerate them from the registry-composed code.
#
# Record-order DATA fix: the pre-WP-11 L1 literal call order was
# feature_general, fgsea, ora, stringdb, sample_general, geneshotTab,
# ptm, rescurve - and call order IS binding registration order (snapshot
# bytes, capability-record order). The registry LIST was reordered to
# exactly that order; the tab-composition order lives in the order
# field, and rs_tab_set sorts by it (output unchanged - the truth table
# at the top of this suite pins it). The D12 freeze vector and these
# goldens are untouched by definition of the invariant.
#
# The per-WP real-L1 payload tests earlier in this suite (w5l1, w6l1,
# the WP-7/8/9/10 blocks) all boot the registry-composed L1 now and stay
# green - the per-panel golden harnesses are the regression net.

PANEL_W11_SNAP <- list(values = list(resultspace.analyst_tab = NULL, resultspace.feature_general.xax_analysis = NULL, 
    resultspace.feature_general.xax_subset = NULL, resultspace.feature_general.xax_variable = NULL, 
    resultspace.feature_general.plot_type = NULL, resultspace.feature_general.regression_line = FALSE, 
    resultspace.feature_general.attr4.color_analysis = NULL, 
    resultspace.feature_general.attr4.color_subset = NULL, resultspace.feature_general.attr4.color_variable = NULL, 
    resultspace.feature_general.attr4.shape_analysis = NULL, 
    resultspace.feature_general.attr4.shape_subset = NULL, resultspace.feature_general.attr4.shape_variable = NULL, 
    resultspace.feature_general.attr4.size_analysis = NULL, resultspace.feature_general.attr4.size_subset = NULL, 
    resultspace.feature_general.attr4.size_variable = NULL, resultspace.feature_general.attr4.tooltip_analysis = NULL, 
    resultspace.feature_general.attr4.tooltip_subset = NULL, 
    resultspace.feature_general.attr4.tooltip_variable = NULL, 
    resultspace.feature_general.attr4.search_analysis = NULL, 
    resultspace.feature_general.attr4.search_subset = NULL, resultspace.feature_general.attr4.search_variable = NULL, 
    resultspace.feature_general.attr4.search_value = character(0), 
    resultspace.feature_general.attr4.xcut = NULL, resultspace.feature_general.attr4.ycut = NULL, 
    resultspace.feature_general.attr4.scorner = "None", resultspace.fgsea.xax_analysis = NULL, 
    resultspace.fgsea.xax_subset = NULL, resultspace.fgsea.xax_variable = NULL, 
    resultspace.fgsea.selected_row = NULL, resultspace.ora.xax_analysis = NULL, 
    resultspace.ora.xax_subset = NULL, resultspace.ora.xax_variable = NULL, 
    resultspace.ora.selected_row = NULL, resultspace.stringdb.taxonomy = NULL, 
    resultspace.stringdb.show_labels = NULL, resultspace.stringdb.selected_row = NULL, 
    resultspace.sample_general.xax_analysis = NULL, resultspace.sample_general.xax_subset = NULL, 
    resultspace.sample_general.xax_variable = NULL, resultspace.sample_general.attr4.color_analysis = NULL, 
    resultspace.sample_general.attr4.color_subset = NULL, resultspace.sample_general.attr4.color_variable = NULL, 
    resultspace.sample_general.attr4.shape_analysis = NULL, resultspace.sample_general.attr4.shape_subset = NULL, 
    resultspace.sample_general.attr4.shape_variable = NULL, resultspace.sample_general.attr4.size_analysis = NULL, 
    resultspace.sample_general.attr4.size_subset = NULL, resultspace.sample_general.attr4.size_variable = NULL, 
    resultspace.sample_general.attr4.tooltip_analysis = NULL, 
    resultspace.sample_general.attr4.tooltip_subset = NULL, resultspace.sample_general.attr4.tooltip_variable = NULL, 
    resultspace.sample_general.attr4.search_analysis = NULL, 
    resultspace.sample_general.attr4.search_subset = NULL, resultspace.sample_general.attr4.search_variable = NULL, 
    resultspace.sample_general.attr4.search_value = character(0), 
    resultspace.sample_general.attr4.xcut = NULL, resultspace.sample_general.attr4.ycut = NULL, 
    resultspace.sample_general.attr4.scorner = "None", resultspace.sample_general.survival_censor = NULL, 
    resultspace.sample_general.batch_show_phenotype = NULL, resultspace.sample_general.batch_show_features = NULL, 
    resultspace.sample_general.batch_phenotype_selected_row = NULL, 
    resultspace.sample_general.batch_feature_selected_row = NULL, 
    resultspace.geneshot.term = NULL, resultspace.geneshot.xax_analysis = NULL, 
    resultspace.geneshot.xax_subset = NULL, resultspace.geneshot.xax_variable = NULL, 
    resultspace.ptm.xax_analysis = "SeqLogo", resultspace.ptm.xax_subset = "All", 
    resultspace.ptm.xax_variable = "seq.window", resultspace.dose_response.feature_table = NULL, 
    resultspace.dose_response.param_table = NULL), epochs = list(
    resultspace.analyst_tab = 0L, resultspace.feature_general.xax_analysis = 0L, 
    resultspace.feature_general.xax_subset = 0L, resultspace.feature_general.xax_variable = 0L, 
    resultspace.feature_general.plot_type = 0L, resultspace.feature_general.regression_line = 0L, 
    resultspace.feature_general.attr4.color_analysis = 0L, resultspace.feature_general.attr4.color_subset = 0L, 
    resultspace.feature_general.attr4.color_variable = 0L, resultspace.feature_general.attr4.shape_analysis = 0L, 
    resultspace.feature_general.attr4.shape_subset = 0L, resultspace.feature_general.attr4.shape_variable = 0L, 
    resultspace.feature_general.attr4.size_analysis = 0L, resultspace.feature_general.attr4.size_subset = 0L, 
    resultspace.feature_general.attr4.size_variable = 0L, resultspace.feature_general.attr4.tooltip_analysis = 0L, 
    resultspace.feature_general.attr4.tooltip_subset = 0L, resultspace.feature_general.attr4.tooltip_variable = 0L, 
    resultspace.feature_general.attr4.search_analysis = 0L, resultspace.feature_general.attr4.search_subset = 0L, 
    resultspace.feature_general.attr4.search_variable = 0L, resultspace.feature_general.attr4.search_value = 0L, 
    resultspace.feature_general.attr4.xcut = 0L, resultspace.feature_general.attr4.ycut = 0L, 
    resultspace.feature_general.attr4.scorner = 1L, resultspace.fgsea.xax_analysis = 0L, 
    resultspace.fgsea.xax_subset = 0L, resultspace.fgsea.xax_variable = 0L, 
    resultspace.fgsea.selected_row = 0L, resultspace.ora.xax_analysis = 0L, 
    resultspace.ora.xax_subset = 0L, resultspace.ora.xax_variable = 0L, 
    resultspace.ora.selected_row = 0L, resultspace.stringdb.taxonomy = 0L, 
    resultspace.stringdb.show_labels = 0L, resultspace.stringdb.selected_row = 0L, 
    resultspace.sample_general.xax_analysis = 0L, resultspace.sample_general.xax_subset = 0L, 
    resultspace.sample_general.xax_variable = 0L, resultspace.sample_general.attr4.color_analysis = 0L, 
    resultspace.sample_general.attr4.color_subset = 0L, resultspace.sample_general.attr4.color_variable = 0L, 
    resultspace.sample_general.attr4.shape_analysis = 0L, resultspace.sample_general.attr4.shape_subset = 0L, 
    resultspace.sample_general.attr4.shape_variable = 0L, resultspace.sample_general.attr4.size_analysis = 0L, 
    resultspace.sample_general.attr4.size_subset = 0L, resultspace.sample_general.attr4.size_variable = 0L, 
    resultspace.sample_general.attr4.tooltip_analysis = 0L, resultspace.sample_general.attr4.tooltip_subset = 0L, 
    resultspace.sample_general.attr4.tooltip_variable = 0L, resultspace.sample_general.attr4.search_analysis = 0L, 
    resultspace.sample_general.attr4.search_subset = 0L, resultspace.sample_general.attr4.search_variable = 0L, 
    resultspace.sample_general.attr4.search_value = 0L, resultspace.sample_general.attr4.xcut = 0L, 
    resultspace.sample_general.attr4.ycut = 0L, resultspace.sample_general.attr4.scorner = 1L, 
    resultspace.sample_general.survival_censor = 0L, resultspace.sample_general.batch_show_phenotype = 0L, 
    resultspace.sample_general.batch_show_features = 0L, resultspace.sample_general.batch_phenotype_selected_row = 0L, 
    resultspace.sample_general.batch_feature_selected_row = 0L, 
    resultspace.geneshot.term = 0L, resultspace.geneshot.xax_analysis = 0L, 
    resultspace.geneshot.xax_subset = 0L, resultspace.geneshot.xax_variable = 0L, 
    resultspace.ptm.xax_analysis = 1L, resultspace.ptm.xax_subset = 1L, 
    resultspace.ptm.xax_variable = 1L, resultspace.dose_response.feature_table = 0L, 
    resultspace.dose_response.param_table = 0L), actions = list(
    resultspace.stringdb.run = 0L, resultspace.geneshot.submit = 0L))
PANEL_W11_REG_ORDER <- c("resultspace.analyst_tab", "resultspace.feature_general.xax_analysis", 
"resultspace.feature_general.xax_subset", "resultspace.feature_general.xax_variable", 
"resultspace.feature_general.plot_type", "resultspace.feature_general.regression_line", 
"resultspace.feature_general.attr4.color_analysis", "resultspace.feature_general.attr4.color_subset", 
"resultspace.feature_general.attr4.color_variable", "resultspace.feature_general.attr4.shape_analysis", 
"resultspace.feature_general.attr4.shape_subset", "resultspace.feature_general.attr4.shape_variable", 
"resultspace.feature_general.attr4.size_analysis", "resultspace.feature_general.attr4.size_subset", 
"resultspace.feature_general.attr4.size_variable", "resultspace.feature_general.attr4.tooltip_analysis", 
"resultspace.feature_general.attr4.tooltip_subset", "resultspace.feature_general.attr4.tooltip_variable", 
"resultspace.feature_general.attr4.search_analysis", "resultspace.feature_general.attr4.search_subset", 
"resultspace.feature_general.attr4.search_variable", "resultspace.feature_general.attr4.search_value", 
"resultspace.feature_general.attr4.xcut", "resultspace.feature_general.attr4.ycut", 
"resultspace.feature_general.attr4.scorner", "resultspace.fgsea.xax_analysis", 
"resultspace.fgsea.xax_subset", "resultspace.fgsea.xax_variable", 
"resultspace.fgsea.selected_row", "resultspace.ora.xax_analysis", 
"resultspace.ora.xax_subset", "resultspace.ora.xax_variable", 
"resultspace.ora.selected_row", "resultspace.stringdb.taxonomy", 
"resultspace.stringdb.show_labels", "resultspace.stringdb.selected_row", 
"resultspace.sample_general.xax_analysis", "resultspace.sample_general.xax_subset", 
"resultspace.sample_general.xax_variable", "resultspace.sample_general.attr4.color_analysis", 
"resultspace.sample_general.attr4.color_subset", "resultspace.sample_general.attr4.color_variable", 
"resultspace.sample_general.attr4.shape_analysis", "resultspace.sample_general.attr4.shape_subset", 
"resultspace.sample_general.attr4.shape_variable", "resultspace.sample_general.attr4.size_analysis", 
"resultspace.sample_general.attr4.size_subset", "resultspace.sample_general.attr4.size_variable", 
"resultspace.sample_general.attr4.tooltip_analysis", "resultspace.sample_general.attr4.tooltip_subset", 
"resultspace.sample_general.attr4.tooltip_variable", "resultspace.sample_general.attr4.search_analysis", 
"resultspace.sample_general.attr4.search_subset", "resultspace.sample_general.attr4.search_variable", 
"resultspace.sample_general.attr4.search_value", "resultspace.sample_general.attr4.xcut", 
"resultspace.sample_general.attr4.ycut", "resultspace.sample_general.attr4.scorner", 
"resultspace.sample_general.survival_censor", "resultspace.sample_general.batch_show_phenotype", 
"resultspace.sample_general.batch_show_features", "resultspace.sample_general.batch_phenotype_selected_row", 
"resultspace.sample_general.batch_feature_selected_row", "resultspace.geneshot.term", 
"resultspace.geneshot.xax_analysis", "resultspace.geneshot.xax_subset", 
"resultspace.geneshot.xax_variable", "resultspace.ptm.xax_analysis", 
"resultspace.ptm.xax_subset", "resultspace.ptm.xax_variable", 
"resultspace.dose_response.feature_table", "resultspace.dose_response.param_table"
)
PANEL_W11_STATUS_NAMES <- c("analyst_active_tab", "analyst_feature_general",
  "analyst_sample_general", "analyst_gene_shot", "analyst_fgsea",
  "analyst_ora", "analyst_ptm", "analyst_dose_response",
  "analyst_additional")
PANEL_W11_STATUS_PAYLOAD <- list(analyst_active_tab = NULL, analyst_feature_general = list(
    xax = NULL, showRegLine = FALSE, attr4 = list(selectColor = NULL, 
        selectShape = NULL, selectSize = NULL, selectTooltip = NULL, 
        searchOnCol = NULL, searchValue = NULL, xcut = NULL, 
        ycut = NULL, acorner = NULL), plotType = NULL), analyst_sample_general = list(
    xax = NULL, attr4 = list(selectColor = NULL, selectShape = NULL, 
        selectSize = NULL, selectTooltip = NULL, searchOnCol = NULL, 
        searchValue = NULL, xcut = NULL, ycut = NULL, acorner = NULL)), 
    analyst_gene_shot = list(xax = NULL, term = NULL), analyst_fgsea = list(
        xax = NULL), analyst_ora = list(xax = NULL), analyst_ptm = list(
        xax = NULL), analyst_dose_response = list(state_support = FALSE), 
    analyst_additional = list())
PANEL_W11_OVERVIEW_NAMES <- c("dataset", "active_tabs", "selection", "available_tabs", "quick_views", 
"available_sections", "state_policy", "capabilities", "panels"
)
PANEL_W11_OVERVIEW_TABS <- list(data_space = c("Feature", "Sample"), analysis_space = c("Feature", 
"Response", "StringDB", "SeqLogo", "Geneshot", "Sample"))
PANEL_W11_OVERVIEW_CAPS <- list(capability_count = 72L, panels = "Analysis space", semantic_tools = c("set_omics_viewer_state", 
"set_scatter_view", "set_enrichment_parameters", "set_table_view", 
"create_figure", "update_figure"))
PANEL_W11_OVERVIEW_PANELS <- list(data_space = list(eset_active_tab = "Feature"), result_space = list(
    analyst_active_tab = NULL, analyst_feature_general = list(
        xax = NULL, showRegLine = FALSE, attr4 = list(selectColor = NULL, 
            selectShape = NULL, selectSize = NULL, selectTooltip = NULL, 
            searchOnCol = NULL, searchValue = NULL, xcut = NULL, 
            ycut = NULL, acorner = NULL), plotType = NULL), analyst_sample_general = list(
        xax = NULL, attr4 = list(selectColor = NULL, selectShape = NULL, 
            selectSize = NULL, selectTooltip = NULL, searchOnCol = NULL, 
            searchValue = NULL, xcut = NULL, ycut = NULL, acorner = NULL)), 
    analyst_gene_shot = list(xax = NULL, term = NULL), analyst_fgsea = list(
        xax = NULL), analyst_ora = list(xax = NULL), analyst_ptm = list(
        xax = NULL), analyst_dose_response = list(state_support = FALSE), 
    analyst_additional = list()))
PANEL_W11_CAP_IDS <- c("resultspace.analyst_tab", "resultspace.feature_general.xax_analysis", 
"resultspace.feature_general.xax_subset", "resultspace.feature_general.xax_variable", 
"resultspace.feature_general.plot_type", "resultspace.feature_general.regression_line", 
"resultspace.feature_general.attr4.color_analysis", "resultspace.feature_general.attr4.color_subset", 
"resultspace.feature_general.attr4.color_variable", "resultspace.feature_general.attr4.shape_analysis", 
"resultspace.feature_general.attr4.shape_subset", "resultspace.feature_general.attr4.shape_variable", 
"resultspace.feature_general.attr4.size_analysis", "resultspace.feature_general.attr4.size_subset", 
"resultspace.feature_general.attr4.size_variable", "resultspace.feature_general.attr4.tooltip_analysis", 
"resultspace.feature_general.attr4.tooltip_subset", "resultspace.feature_general.attr4.tooltip_variable", 
"resultspace.feature_general.attr4.search_analysis", "resultspace.feature_general.attr4.search_subset", 
"resultspace.feature_general.attr4.search_variable", "resultspace.feature_general.attr4.search_value", 
"resultspace.feature_general.attr4.xcut", "resultspace.feature_general.attr4.ycut", 
"resultspace.feature_general.attr4.scorner", "resultspace.fgsea.xax_analysis", 
"resultspace.fgsea.xax_subset", "resultspace.fgsea.xax_variable", 
"resultspace.fgsea.selected_row", "resultspace.ora.xax_analysis", 
"resultspace.ora.xax_subset", "resultspace.ora.xax_variable", 
"resultspace.ora.selected_row", "resultspace.stringdb.taxonomy", 
"resultspace.stringdb.show_labels", "resultspace.stringdb.selected_row", 
"resultspace.sample_general.xax_analysis", "resultspace.sample_general.xax_subset", 
"resultspace.sample_general.xax_variable", "resultspace.sample_general.attr4.color_analysis", 
"resultspace.sample_general.attr4.color_subset", "resultspace.sample_general.attr4.color_variable", 
"resultspace.sample_general.attr4.shape_analysis", "resultspace.sample_general.attr4.shape_subset", 
"resultspace.sample_general.attr4.shape_variable", "resultspace.sample_general.attr4.size_analysis", 
"resultspace.sample_general.attr4.size_subset", "resultspace.sample_general.attr4.size_variable", 
"resultspace.sample_general.attr4.tooltip_analysis", "resultspace.sample_general.attr4.tooltip_subset", 
"resultspace.sample_general.attr4.tooltip_variable", "resultspace.sample_general.attr4.search_analysis", 
"resultspace.sample_general.attr4.search_subset", "resultspace.sample_general.attr4.search_variable", 
"resultspace.sample_general.attr4.search_value", "resultspace.sample_general.attr4.xcut", 
"resultspace.sample_general.attr4.ycut", "resultspace.sample_general.attr4.scorner", 
"resultspace.sample_general.survival_censor", "resultspace.sample_general.batch_show_phenotype", 
"resultspace.sample_general.batch_show_features", "resultspace.sample_general.batch_phenotype_selected_row", 
"resultspace.sample_general.batch_feature_selected_row", "resultspace.geneshot.term", 
"resultspace.geneshot.xax_analysis", "resultspace.geneshot.xax_subset", 
"resultspace.geneshot.xax_variable", "resultspace.ptm.xax_analysis", 
"resultspace.ptm.xax_subset", "resultspace.ptm.xax_variable", 
"resultspace.dose_response.feature_table", "resultspace.dose_response.param_table", 
"action:resultspace.stringdb.run", "action:resultspace.geneshot.submit", 
"panel:Feature Scatter Plot", "panel:Feature Metadata Table", 
"panel:Sample Scatter Plot", "panel:Sample Metadata Table", "panel:Correlation Heatmap", 
"panel:Expression Heatmap", "panel:Dynamic Heatmap with Selection", 
"panel:Expression Matrix Table", "panel:Gene Set Membership Table", 
"panel:Feature Analysis", "panel:Sample Analysis", "panel:Geneshot Literature Search", 
"panel:Over-Representation Analysis", "panel:Fast Gene Set Enrichment Analysis", 
"panel:Dose-Response Curves", "panel:STRING Protein Interaction Network", 
"panel:Post-Translational Modification Motifs", "tool:get_omics_viewer_state", 
"tool:search_annotations", "tool:summarize_annotation", "tool:set_omics_viewer_state", 
"tool:set_scatter_view", "tool:set_enrichment_parameters", "tool:set_table_view", 
"tool:create_figure", "tool:update_figure", "tool:find_controls", 
"tool:describe_control", "tool:set_widgets", "tool:run_control"
)
PANEL_W11_CAP_TABS <- structure(c(Cor = 1L, `Dynamic heatmap` = 1L, Expression = 1L, 
`Feature table` = 1L, GSList = 1L, Heatmap = 1L, `Sample table` = 1L, 
`Analysis space` = 2L, Annotations = 2L, `Data space` = 2L, Figures = 2L, 
`Global state` = 2L, Response = 3L, Controls = 4L, SeqLogo = 4L, 
fGSEA = 5L, ORA = 5L, StringDB = 5L, Geneshot = 6L, Feature = 26L, 
Sample = 29L), dim = 21L, dimnames = structure(list(c("Cor", 
"Dynamic heatmap", "Expression", "Feature table", "GSList", "Heatmap", 
"Sample table", "Analysis space", "Annotations", "Data space", 
"Figures", "Global state", "Response", "Controls", "SeqLogo", 
"fGSEA", "ORA", "StringDB", "Geneshot", "Feature", "Sample")), names = ""), class = "table")
PANEL_W11_READER_COLD <- list(ora = list(method = "ORA", status = "not_run", message = "No results are computed yet. Run the analysis first (set_enrichment_parameters or the panel itself), then request this section again."), 
    fgsea = list(method = "fGSEA", status = "not_run", message = "No results are computed yet. Run the analysis first (set_enrichment_parameters or the panel itself), then request this section again."))

store_snapshot <- omicsViewer:::store_snapshot
build_app_state <- omicsViewer:::build_app_state
agent_compact_state <- omicsViewer:::agent_compact_state
store_action_view <- omicsViewer:::store_action_view

# --- the registry-composed L1 reproduces EVERY golden byte-for-byte
# (this IS the loop-vs-old-call-site equivalence proof: the literals are
# pre-change captures of the hand-wired composition, booted the same way)
w11 <- new.env()
w11$store <- widget_store_new()
w11$ret <- NULL
app_w11 <- function(input, output, session) {
  w11$ret <- omicsViewer:::L1_result_space_module(
    "rs",
    reactive_expr = reactive(Biobase::exprs(dat)),
    reactive_phenoData = reactive(Biobase::pData(dat)),
    reactive_featureData = reactive(fd0),
    status = reactive(NULL),
    store = w11$store)
}
w11_err <- NULL; w11_snap <- NULL; w11_reg_order <- NULL
w11_payload <- NULL; w11_cap <- NULL; w11_reader <- NULL
shiny::testServer(app_w11, {
  w11_err <<- tryCatch({
    for (i in 1:3) session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w11_snap <<- store_snapshot(w11$store)
  w11_reg_order <<- unname(vapply(store_registry_view(w11$store),
                                  function(r) r$id, ""))
  w11_payload <<- tryCatch(shiny::isolate(w11$ret()),
                           error = function(e) NULL)
  w11_cap <<- agent_capability_records(w11$store)
  w11_reader <<- tryCatch(
    attr(shiny::isolate(w11$ret()), "read_enrichment_results")(),
    error = function(e) conditionMessage(e))
})
ok(is.null(w11_err),
   "WP-11: the registry-composed L1 boots headless (one loop, eight records)")
ok(ut_cmp_identical(w11_snap, PANEL_W11_SNAP),
   "WP-11: FULL store snapshot byte-identical (values+epochs+actions, names AND order = the pre-change call order)")
ok(ut_cmp_identical(w11_reg_order, PANEL_W11_REG_ORDER),
   "WP-11: store_registry_view ids in the frozen registration order (the capability-record order)")
ok(ut_cmp_identical(sort(unique(w11_reg_order)), PANEL_REGISTRY_FREEZE_IDS),
   "WP-11: the D12 enumeration freeze stays green through the loop (72 ids)")
p11 <- w11_payload
attributes(p11) <- list(names = names(p11))
ok(ut_cmp_identical(names(w11_payload), PANEL_W11_STATUS_NAMES),
   "WP-11: status payload keeps its exact key order (analyst_active_tab first, analyst_additional last; the plan's 'registry order == literal order' parenthetical was wrong - the goldens pin the truth)")
ok(ut_cmp_identical(p11, PANEL_W11_STATUS_PAYLOAD),
   "WP-11: status payload content byte-identical (registry-keyed map, emission-order-pinned builder; dose marker, string unemitted)")
ok(is.function(attr(w11_payload, "read_enrichment_results")),
   "WP-11: the enrichment reader still rides the payload attribute (L1-owned, map-fed)")
ok(ut_cmp_identical(w11_reader, PANEL_W11_READER_COLD),
   "WP-11: the reader's cold behavior over the demo boot is byte-identical (ora/fgsea not_run notes)")
ok(ut_cmp_identical(vapply(w11_cap, function(r) r$id, ""), PANEL_W11_CAP_IDS),
   "WP-11: agent capability records keep their exact id order (104 records)")
ok(ut_cmp_identical(sort(table(vapply(w11_cap, function(r) r$tab, ""))),
                    PANEL_W11_CAP_TABS),
   "WP-11: agent capability per-tab counts unchanged")

# --- agent state overview shape over the WP-11 composition (the
# get_omics_viewer_state output shape: top-level section names, the
# available-tabs block with rs_tab_set's analysis_space, the capability
# summary over the real store, and the panels section over the payload)
w11_state <- build_app_state(
  dataset = dat, dataset_id = "demo.RDS",
  data_status = list(eset_active_tab = "Feature"),
  result_status = w11_payload,
  selected_features = "EGFR", selected_samples = character(),
  label = "wp11")
w11_ov <- agent_compact_state(
  state = w11_state,
  available_tabs = list(
    data_space = c("Feature", "Sample"),
    analysis_space = rs_tab_set(fd0)),
  store = w11$store,
  sections = "panels")
ok(ut_cmp_identical(names(w11_ov), PANEL_W11_OVERVIEW_NAMES),
   "WP-11: agent state overview section names/order unchanged")
ok(ut_cmp_identical(w11_ov$available_tabs, PANEL_W11_OVERVIEW_TABS),
   "WP-11: overview available_tabs unchanged (analysis_space = rs_tab_set output; sort-by-order keeps it byte-identical)")
ok(ut_cmp_identical(w11_ov$capabilities, PANEL_W11_OVERVIEW_CAPS),
   "WP-11: overview capability summary unchanged (count 72, panel labels, semantic tools)")
ok(ut_cmp_identical(w11_ov$panels, PANEL_W11_OVERVIEW_PANELS),
   "WP-11: overview panels section byte-identical (the bounded view of the status payload)")

# --- loop data threading, dose side: ctx$data$object rides the loop
# exactly as at the old call site (the S6.6_drc derivation goes live;
# the WP-3 module-level golden re-proven at L1 level)
w11o <- new.env()
w11o$store <- widget_store_new()
w11o_err <- NULL; w11o_par <- NULL; w11o_feat <- NULL
shiny::testServer(function(input, output, session) {
  omicsViewer:::L1_result_space_module(
    "rs",
    reactive_expr = reactive(Biobase::exprs(dat)),
    reactive_phenoData = reactive(Biobase::pData(dat)),
    reactive_featureData = reactive(fd0),
    reactive_i = reactive(1L),
    object = reactive(dat),
    status = reactive(NULL),
    store = w11o$store)
}, {
  w11o_err <<- tryCatch({
    for (i in 1:3) session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  v <- store_registry_view(w11o$store, "resultspace.dose_response")
  w11o_par <<- Filter(function(r) grepl("param_table$", r$id), v)[[1]]$allowed_values
  w11o_feat <<- Filter(function(r) grepl("feature_table$", r$id), v)[[1]]$allowed_values
})
ok(is.null(w11o_err),
   "WP-11: the loop boots L1 with object threaded (dose derivation active)")
ok(ut_cmp_identical(w11o_par, .w3_params),
   "WP-11: dose param choices live through the loop (identical to the WP-3 module golden rows)")
ok(ut_cmp_identical(w11o_feat, rownames(fd0)[1]),
   "WP-11: dose feature-table choice is the selected feature id through the loop")

# --- loop data threading, string side: ctx$selection$feature rides the
# loop (the moved-in ids derivation feeds the run action's availability;
# the WP-4 module-level golden re-proven at L1 level)
w11s <- new.env()
w11s$store <- widget_store_new()
w11s$sel <- shiny::reactiveVal(NULL)
w11s_err <- NULL; w11s_act0 <- NULL; w11s_act1 <- NULL
shiny::testServer(function(input, output, session) {
  omicsViewer:::L1_result_space_module(
    "rs",
    reactive_expr = reactive(Biobase::exprs(dat)),
    reactive_phenoData = reactive(Biobase::pData(dat)),
    reactive_featureData = reactive(fd0),
    reactive_i = w11s$sel,
    status = reactive(NULL),
    store = w11s$store)
}, {
  w11s_err <<- tryCatch({
    for (i in 1:3) session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  w11s_act0 <<- store_action_view(
    w11s$store, "resultspace.stringdb")[[1]]$available
  w11s$sel(1L)
  for (i in 1:3) session$flushReact()
  w11s_act1 <<- store_action_view(
    w11s$store, "resultspace.stringdb")[[1]]$available
})
ok(is.null(w11s_err),
   "WP-11: the loop boots L1 with the string selection threaded")
ok(ut_cmp_identical(
  w11s_act0,
  "select features first (set_omics_viewer_state features)"),
   "WP-11: no selection -> the loop-built string context derives no ids (run action gated)")
ok(ut_cmp_identical(w11s_act1, TRUE),
   "WP-11: a feature selection -> the derived ids go live through the loop (identical to the WP-4 module golden)")

# --- store-less L1 (exported/embedded use): the loop passes NULL child
# stores, every panel boots on its legacy path, the status payload still
# emits (the 4.6 regression shape, now with the full key-order pin)
w11n_err <- NULL; w11n_names <- NULL; w11n_dose <- NULL
shiny::testServer(function(input, output, session) {
  v <- omicsViewer:::L1_result_space_module(
    "rs",
    reactive_expr = reactive(Biobase::exprs(dat)),
    reactive_phenoData = reactive(Biobase::pData(dat)),
    reactive_featureData = reactive(fd0),
    status = reactive(list(analyst_active_tab = "Sample")),
    store = NULL)
  w11n_names <<- tryCatch({
    for (i in 1:10) session$flushReact()
    names(shiny::isolate(v()))
  }, error = function(e) conditionMessage(e))
  w11n_dose <<- tryCatch(
    shiny::isolate(v())$analyst_dose_response, error = function(e) NULL)
}, { })
ok(ut_cmp_identical(w11n_names, PANEL_W11_STATUS_NAMES),
   "WP-11: store-less L1 boots through the loop and emits the frozen key set in order")
ok(ut_cmp_identical(w11n_dose, list(state_support = FALSE)),
   "WP-11: store-less L1 keeps the analyst_dose_response state_support marker")
## ------------------- WP-12: L1 client composition (DOM parity) ----
# Plan section 4 WP-12: the inline tabPanel list and the renderUI
# conditional visibility became registry-driven (records' order field
# drives tab order, rs_tab_set drives visibility, records'
# presentation fields drive the tabPanels; the additionalTabs UI
# loop and the pinned Geneshot/Sample tail stay byte-verbatim,
# D15/Q1). Parity method (recorded as the WP-12 decision):
#
# GOLDENS (frozen, captured 2026-11-06 at PRE-CHANGE HEAD f7fd24e
# BEFORE any WP-12 edit - the PD-29/PD-33/PD-48/PD-55/PD-66 school;
# never regenerate from the migrated code): the renderUI output of
# the REAL L1 composition, rendered through the identical renderUI
# pipeline (testServer getOutput), under every headless availability
# variant. Determinism was proven by a double pre-change capture run
# (byte-identical, separate processes) and the IDENTICAL capture
# script re-run post-change reproduced every golden byte-identically.
#
# NORMALIZATION RULES (the only nondeterminism; bslib/shinyWidgets
# random ids per render - shiny keeps a PRIVATE seed, so set.seed
# cannot pin them):
#   data-tabsetid="<N>" -> data-tabsetid="TSID"   (bslib tabset ids,
#                                                     outer + nested)
#   tab-<N>-            -> tab-TSID-               (pane ids / hrefs)
#   btn-<N>             -> btn-RID                 (shinyWidgets ids)
# Everything else must be byte-identical. The full normalized demo
# HTML is frozen below (WP12_DEMO_HTML); the other five variants are
# frozen as md5 digests over their full normalized HTML plus complete
# parsed truth tables (per-tab value/tooltip/icon/label/sr-only h2/
# aria-label).
#
# LEGACY FIXTURE: the pre-change renderUI expression kept verbatim
# below (the plan's "old expression kept one commit as fixture"),
# wrapped in function(ns, reactive_featureData, additionalTabs) and
# run through its own renderUI output in the SAME session - the
# same-commit old-vs-new HTML diff, asserted per variant.

# the verbatim pre-WP-12 composition (extracted mechanically from the
# f7fd24e blob; its output is pinned against the frozen goldens, so a
# corrupted paste cannot pass silently)
.wp12_legacy_tabset <- function(ns, reactive_featureData, additionalTabs) {
      titleTabs <- list(
        title = tagList(
          tags$span(class = "omicsviewer-nav-title", "Analysis"),
          tags$span(class = "omicsviewer-nav-sub", " \u2014 analyse the current selection on the fly")
        ),
        windowTitle = "Analysis",
        id = ns("analyst"),
        theme = shinytheme("spacelab"),
        tabPanel(
          tags$span(icon("chart-column"), "Feature", title = "Selected-feature analysis.\nCompares the abundance of the currently selected feature(s) across sample groups:\nbeeswarm/box plots with group-comparison statistics (t-test or Wilcoxon for two groups,\nANOVA/Kruskal-Wallis for more), or ROC / precision-recall curves when a binary\noutcome is chosen.\nInterpretation: the box spans the interquartile range with the median as center line.\nA small p-value (e.g. < 0.05) indicates a statistically significant group difference -\nalways judge effect size and overlap on the plot too. For ROC/PR curves: AUC near 1.0\nmeans the feature separates the two classes well; 0.5 is no better than chance."),
          value = "Feature",
          tags$h2("Feature Analysis", class = "sr-only", `aria-label` = "Statistical analysis and visualization of selected features including boxplots with group comparisons and ROC curves for binary outcomes"),
          feature_general_ui(ns("feature_general"))
        )
      )
      sampleAnalyst <- list(
        tabPanel(
          tags$span(icon("arrows-down-to-people"), "Sample", title = "Selected-sample analysis.\nAnalyses the currently selected samples against the rest of the cohort:\nbeeswarm plots for continuous variables, contingency tables of categorical\nannotations (chi-square / Fisher's exact test), and Kaplan-Meier survival\ncurves with the log-rank test when time-to-event data are available.\nInterpretation: a small contingency p-value means the selected samples differ in\nthat annotation more than expected by chance. Kaplan-Meier curves show the\nevent-free proportion over time per group; the log-rank p-value tests whether the\ncurves differ. Check group sizes (n) - small groups make tests unreliable."),
          value = "Sample",
          tags$h2("Sample Analysis", class = "sr-only", `aria-label` = "Statistical analysis of selected samples including group comparisons, contingency tables, and Kaplan-Meier survival curves"),
          sample_general_ui(ns("sample_general"))
        )
      )

      ### geneshot
      geneshot <- list(
        tabPanel(
          tags$span(icon("book-open"), "Geneshot", title = "Literature-based gene discovery.\nSearches PubMed-associated literature (AutoRIF) for genes co-mentioned with your\nsearch term(s); separate multiple terms with ';'. Useful for hypothesis generation.\nInterpretation: the scatter plot shows publication count (x) versus the share of a\ngene's publications mentioning your term (y). Genes in the upper right are both\nwell studied and specifically associated with your topic. A '+' marks genes that\nare present in your current feature selection. Ranking reflects literature bias,\nnot statistical evidence from your data."),
          value = "Geneshot",
          tags$h2("Geneshot Literature Search", class = "sr-only", `aria-label` = "Literature-based gene discovery using PubMed co-mention analysis to find genes associated with search terms"),
          geneshot_ui(ns("geneshotTab"))
        )
      )
      ### end

      optionalTabs <- list()

      if (!is.null(attr(reactive_featureData(), "GS"))) {
        optionalTabs <- c(optionalTabs, list(tabPanel(
          tags$span(icon("object-group"), "ORA", title = "Over-representation analysis (hypergeometric test).\nTests whether your currently selected features belong to a gene set (pathway, GO\ncategory, custom set) more often than expected by chance. Needs a well-defined\nselection, e.g. significant hits from the Feature tab.\nInterpretation: the p-value is the chance probability of the observed overlap; FDR\ncorrects for testing many gene sets (FDR < 0.05 is a common significance cutoff);\nthe odds ratio / gene ratio measures enrichment strength. The second table lists\nthe overlapping features per gene set. Results depend on selection size and the\nbackground gene set collection."),
          value = "ORA",
          tags$h2("Over-Representation Analysis", class = "sr-only", `aria-label` = "Gene set over-representation analysis using hypergeometric test to identify enriched pathways and functional categories"),
          enrichment_analysis_ui(ns("ora"))
        )))
        optionalTabs <- c(optionalTabs, list(tabPanel(
          tags$span(icon("arrow-down-short-wide"), "fGSEA", title = "Fast gene set enrichment analysis (GSEA on a ranked list).\nUnlike ORA, fGSEA uses ALL features, ranked by a statistic you choose (e.g. fold\nchange or signed p-value), and asks whether a gene set's members cluster at the\ntop or bottom of that ranking - no hard significance cutoff needed.\nInterpretation: NES (normalized enrichment score) > 0 means the set is enriched\namong up-ranked features, NES < 0 depleted/enriched at the other end; padj < 0.05\nis a common significance threshold; the leading edge lists the members driving\nthe enrichment signal."),
          value = "fGSEA",
          tags$h2("Fast Gene Set Enrichment Analysis", class = "sr-only", `aria-label` = "Ranked gene set enrichment analysis computing normalized enrichment scores and identifying leading edge genes"),
          enrichment_fgsea_ui(ns("fgsea"))
        )))
      }

      if (any(grepl("^ResponseCurve\\|", colnames(reactive_featureData())))) {
        optionalTabs <- c(optionalTabs, list(tabPanel(
          tags$span(icon("capsules"), "Response", title = "Dose-response curves (4-parameter logistic fit).\nFits abundance vs concentration for the selected features and estimates EC50/IC50\n(concentration giving half-maximal effect), Hill slope (curve steepness) and the\ntop/bottom response plateaus.\nInterpretation: a lower EC50/IC50 means higher potency; the Hill slope shows how\nabruptly the response switches; the parameter table gives fitted values with\nconfidence intervals. Points are observed data, the line the model fit."),
          value = "Response",
          tags$h2("Dose-Response Curves", class = "sr-only", `aria-label` = "Dose-response curve fitting with EC50 and IC50 estimation using 4-parameter logistic regression model"),
          dose_response_ui(ns("rescurve"))
        )))
      }

      if (any(grepl("^StringDB\\|", colnames(reactive_featureData())))) {
        optionalTabs <- c(optionalTabs, list(tabPanel(
          tags$span(icon("share-nodes"), "StringDB", title = "Protein-protein interaction network (STRING database).\nQueries STRING for known and predicted interactions among the currently selected\nfeatures (max ~300; set the taxonomy code, e.g. 9606 = human) and draws the\nnetwork - nodes are proteins, edges are interactions. Requires internet access.\nInterpretation: edge confidence reflects combined evidence (experiments,\ndatabases, text mining); dense clusters often correspond to functional modules or\ncomplexes; highly connected hub proteins may be key regulators."),
          value = "StringDB",
          tags$h2("STRING Protein Interaction Network", class = "sr-only", `aria-label` = "Protein-protein interaction network from STRING database showing physical and functional associations"),
          string_ui(ns("stringdb"))
        )))
      }

      if (any(grepl("^SeqLogo\\|", colnames(reactive_featureData())))) {
        optionalTabs <- c(optionalTabs, list(tabPanel(
          tags$span(icon("font"), "SeqLogo", title = "PTM motif enrichment (sequence logo).\nTakes the amino-acid sequence windows around the selected modification sites and\ntests which residues are over-represented at each position, shown as a sequence\nlogo around the modification site (position 0).\nInterpretation: taller letters mean the amino acid (or property class) is more\nconserved/enriched at that position; significant enrichment versus all sites in\nthe dataset hints at the kinase or recognition machinery involved."),
          value = "SeqLogo",
          tags$h2("Post-Translational Modification Motifs", class = "sr-only", `aria-label` = "PTM motif enrichment analysis identifying overrepresented amino acid sequence patterns around modification sites"),
          ptmotif_ui(ns("ptm"))
        )))
      }

      ######
      if (length(additionalTabs) > 0) {
        for (lo in additionalTabs) {
          optionalTabs <- c(optionalTabs, list(tabPanel(lo$tabName, lo$moduleUi(ns(lo$moduleName)))))
        }
      }
      do.call(navbarPage, c(titleTabs, optionalTabs, geneshot, sampleAnalyst))
}
environment(.wp12_legacy_tabset) <- asNamespace("omicsViewer")

WP12_MD5 <- c(demo = "e7fb101f6d820a2d9545f8c9e80585dd", gs = "a53fbc828a4d8d953e1271e67613b777", 
noopt = "0ea56978968da2493bd8dbe8edc6481b", nullfd = "0ea56978968da2493bd8dbe8edc6481b", 
demo_addl = "6b14a9d69155b379746947d80cbb94e6", gs_addl = "95f1d59896aba91634e388e44ed8c8cd"
)

WP12_TABS_DEMO <- list(list(value = "Feature", title = "Selected-feature analysis.&#10;Compares the abundance of the currently selected feature(s) across sample groups:&#10;beeswarm/box plots with group-comparison statistics (t-test or Wilcoxon for two groups,&#10;ANOVA/Kruskal-Wallis for more), or ROC / precision-recall curves when a binary&#10;outcome is chosen.&#10;Interpretation: the box spans the interquartile range with the median as center line.&#10;A small p-value (e.g. &lt; 0.05) indicates a statistically significant group difference -&#10;always judge effect size and overlap on the plot too. For ROC/PR curves: AUC near 1.0&#10;means the feature separates the two classes well; 0.5 is no better than chance.", 
    icon = "chart-column", label = "Feature", sr_title = "Feature Analysis", 
    aria = "Statistical analysis and visualization of selected features including boxplots with group comparisons and ROC curves for binary outcomes"), 
    list(value = "Response", title = "Dose-response curves (4-parameter logistic fit).&#10;Fits abundance vs concentration for the selected features and estimates EC50/IC50&#10;(concentration giving half-maximal effect), Hill slope (curve steepness) and the&#10;top/bottom response plateaus.&#10;Interpretation: a lower EC50/IC50 means higher potency; the Hill slope shows how&#10;abruptly the response switches; the parameter table gives fitted values with&#10;confidence intervals. Points are observed data, the line the model fit.", 
        icon = "capsules", label = "Response", sr_title = "Dose-Response Curves", 
        aria = "Dose-response curve fitting with EC50 and IC50 estimation using 4-parameter logistic regression model"), 
    list(value = "StringDB", title = "Protein-protein interaction network (STRING database).&#10;Queries STRING for known and predicted interactions among the currently selected&#10;features (max ~300; set the taxonomy code, e.g. 9606 = human) and draws the&#10;network - nodes are proteins, edges are interactions. Requires internet access.&#10;Interpretation: edge confidence reflects combined evidence (experiments,&#10;databases, text mining); dense clusters often correspond to functional modules or&#10;complexes; highly connected hub proteins may be key regulators.", 
        icon = "share-nodes", label = "StringDB", sr_title = "STRING Protein Interaction Network", 
        aria = "Protein-protein interaction network from STRING database showing physical and functional associations"), 
    list(value = "SeqLogo", title = "PTM motif enrichment (sequence logo).&#10;Takes the amino-acid sequence windows around the selected modification sites and&#10;tests which residues are over-represented at each position, shown as a sequence&#10;logo around the modification site (position 0).&#10;Interpretation: taller letters mean the amino acid (or property class) is more&#10;conserved/enriched at that position; significant enrichment versus all sites in&#10;the dataset hints at the kinase or recognition machinery involved.", 
        icon = "font", label = "SeqLogo", sr_title = "Post-Translational Modification Motifs", 
        aria = "PTM motif enrichment analysis identifying overrepresented amino acid sequence patterns around modification sites"), 
    list(value = "Geneshot", title = "Literature-based gene discovery.&#10;Searches PubMed-associated literature (AutoRIF) for genes co-mentioned with your&#10;search term(s); separate multiple terms with &#39;;&#39;. Useful for hypothesis generation.&#10;Interpretation: the scatter plot shows publication count (x) versus the share of a&#10;gene&#39;s publications mentioning your term (y). Genes in the upper right are both&#10;well studied and specifically associated with your topic. A &#39;+&#39; marks genes that&#10;are present in your current feature selection. Ranking reflects literature bias,&#10;not statistical evidence from your data.", 
        icon = "book-open", label = "Geneshot", sr_title = "Geneshot Literature Search", 
        aria = "Literature-based gene discovery using PubMed co-mention analysis to find genes associated with search terms"), 
    list(value = "Sample", title = "Selected-sample analysis.&#10;Analyses the currently selected samples against the rest of the cohort:&#10;beeswarm plots for continuous variables, contingency tables of categorical&#10;annotations (chi-square / Fisher&#39;s exact test), and Kaplan-Meier survival&#10;curves with the log-rank test when time-to-event data are available.&#10;Interpretation: a small contingency p-value means the selected samples differ in&#10;that annotation more than expected by chance. Kaplan-Meier curves show the&#10;event-free proportion over time per group; the log-rank p-value tests whether the&#10;curves differ. Check group sizes (n) - small groups make tests unreliable.", 
        icon = "arrows-down-to-people", label = "Sample", sr_title = "Sample Analysis", 
        aria = "Statistical analysis of selected samples including group comparisons, contingency tables, and Kaplan-Meier survival curves"))
WP12_TABS_GS <- list(list(value = "Feature", title = "Selected-feature analysis.&#10;Compares the abundance of the currently selected feature(s) across sample groups:&#10;beeswarm/box plots with group-comparison statistics (t-test or Wilcoxon for two groups,&#10;ANOVA/Kruskal-Wallis for more), or ROC / precision-recall curves when a binary&#10;outcome is chosen.&#10;Interpretation: the box spans the interquartile range with the median as center line.&#10;A small p-value (e.g. &lt; 0.05) indicates a statistically significant group difference -&#10;always judge effect size and overlap on the plot too. For ROC/PR curves: AUC near 1.0&#10;means the feature separates the two classes well; 0.5 is no better than chance.", 
    icon = "chart-column", label = "Feature", sr_title = "Feature Analysis", 
    aria = "Statistical analysis and visualization of selected features including boxplots with group comparisons and ROC curves for binary outcomes"), 
    list(value = "ORA", title = "Over-representation analysis (hypergeometric test).&#10;Tests whether your currently selected features belong to a gene set (pathway, GO&#10;category, custom set) more often than expected by chance. Needs a well-defined&#10;selection, e.g. significant hits from the Feature tab.&#10;Interpretation: the p-value is the chance probability of the observed overlap; FDR&#10;corrects for testing many gene sets (FDR &lt; 0.05 is a common significance cutoff);&#10;the odds ratio / gene ratio measures enrichment strength. The second table lists&#10;the overlapping features per gene set. Results depend on selection size and the&#10;background gene set collection.", 
        icon = "object-group", label = "ORA", sr_title = "Over-Representation Analysis", 
        aria = "Gene set over-representation analysis using hypergeometric test to identify enriched pathways and functional categories"), 
    list(value = "fGSEA", title = "Fast gene set enrichment analysis (GSEA on a ranked list).&#10;Unlike ORA, fGSEA uses ALL features, ranked by a statistic you choose (e.g. fold&#10;change or signed p-value), and asks whether a gene set&#39;s members cluster at the&#10;top or bottom of that ranking - no hard significance cutoff needed.&#10;Interpretation: NES (normalized enrichment score) &gt; 0 means the set is enriched&#10;among up-ranked features, NES &lt; 0 depleted/enriched at the other end; padj &lt; 0.05&#10;is a common significance threshold; the leading edge lists the members driving&#10;the enrichment signal.", 
        icon = "arrow-down-short-wide", label = "fGSEA", sr_title = "Fast Gene Set Enrichment Analysis", 
        aria = "Ranked gene set enrichment analysis computing normalized enrichment scores and identifying leading edge genes"), 
    list(value = "Response", title = "Dose-response curves (4-parameter logistic fit).&#10;Fits abundance vs concentration for the selected features and estimates EC50/IC50&#10;(concentration giving half-maximal effect), Hill slope (curve steepness) and the&#10;top/bottom response plateaus.&#10;Interpretation: a lower EC50/IC50 means higher potency; the Hill slope shows how&#10;abruptly the response switches; the parameter table gives fitted values with&#10;confidence intervals. Points are observed data, the line the model fit.", 
        icon = "capsules", label = "Response", sr_title = "Dose-Response Curves", 
        aria = "Dose-response curve fitting with EC50 and IC50 estimation using 4-parameter logistic regression model"), 
    list(value = "StringDB", title = "Protein-protein interaction network (STRING database).&#10;Queries STRING for known and predicted interactions among the currently selected&#10;features (max ~300; set the taxonomy code, e.g. 9606 = human) and draws the&#10;network - nodes are proteins, edges are interactions. Requires internet access.&#10;Interpretation: edge confidence reflects combined evidence (experiments,&#10;databases, text mining); dense clusters often correspond to functional modules or&#10;complexes; highly connected hub proteins may be key regulators.", 
        icon = "share-nodes", label = "StringDB", sr_title = "STRING Protein Interaction Network", 
        aria = "Protein-protein interaction network from STRING database showing physical and functional associations"), 
    list(value = "SeqLogo", title = "PTM motif enrichment (sequence logo).&#10;Takes the amino-acid sequence windows around the selected modification sites and&#10;tests which residues are over-represented at each position, shown as a sequence&#10;logo around the modification site (position 0).&#10;Interpretation: taller letters mean the amino acid (or property class) is more&#10;conserved/enriched at that position; significant enrichment versus all sites in&#10;the dataset hints at the kinase or recognition machinery involved.", 
        icon = "font", label = "SeqLogo", sr_title = "Post-Translational Modification Motifs", 
        aria = "PTM motif enrichment analysis identifying overrepresented amino acid sequence patterns around modification sites"), 
    list(value = "Geneshot", title = "Literature-based gene discovery.&#10;Searches PubMed-associated literature (AutoRIF) for genes co-mentioned with your&#10;search term(s); separate multiple terms with &#39;;&#39;. Useful for hypothesis generation.&#10;Interpretation: the scatter plot shows publication count (x) versus the share of a&#10;gene&#39;s publications mentioning your term (y). Genes in the upper right are both&#10;well studied and specifically associated with your topic. A &#39;+&#39; marks genes that&#10;are present in your current feature selection. Ranking reflects literature bias,&#10;not statistical evidence from your data.", 
        icon = "book-open", label = "Geneshot", sr_title = "Geneshot Literature Search", 
        aria = "Literature-based gene discovery using PubMed co-mention analysis to find genes associated with search terms"), 
    list(value = "Sample", title = "Selected-sample analysis.&#10;Analyses the currently selected samples against the rest of the cohort:&#10;beeswarm plots for continuous variables, contingency tables of categorical&#10;annotations (chi-square / Fisher&#39;s exact test), and Kaplan-Meier survival&#10;curves with the log-rank test when time-to-event data are available.&#10;Interpretation: a small contingency p-value means the selected samples differ in&#10;that annotation more than expected by chance. Kaplan-Meier curves show the&#10;event-free proportion over time per group; the log-rank p-value tests whether the&#10;curves differ. Check group sizes (n) - small groups make tests unreliable.", 
        icon = "arrows-down-to-people", label = "Sample", sr_title = "Sample Analysis", 
        aria = "Statistical analysis of selected samples including group comparisons, contingency tables, and Kaplan-Meier survival curves"))
WP12_TABS_NOOPT <- list(list(value = "Feature", title = "Selected-feature analysis.&#10;Compares the abundance of the currently selected feature(s) across sample groups:&#10;beeswarm/box plots with group-comparison statistics (t-test or Wilcoxon for two groups,&#10;ANOVA/Kruskal-Wallis for more), or ROC / precision-recall curves when a binary&#10;outcome is chosen.&#10;Interpretation: the box spans the interquartile range with the median as center line.&#10;A small p-value (e.g. &lt; 0.05) indicates a statistically significant group difference -&#10;always judge effect size and overlap on the plot too. For ROC/PR curves: AUC near 1.0&#10;means the feature separates the two classes well; 0.5 is no better than chance.", 
    icon = "chart-column", label = "Feature", sr_title = "Feature Analysis", 
    aria = "Statistical analysis and visualization of selected features including boxplots with group comparisons and ROC curves for binary outcomes"), 
    list(value = "Geneshot", title = "Literature-based gene discovery.&#10;Searches PubMed-associated literature (AutoRIF) for genes co-mentioned with your&#10;search term(s); separate multiple terms with &#39;;&#39;. Useful for hypothesis generation.&#10;Interpretation: the scatter plot shows publication count (x) versus the share of a&#10;gene&#39;s publications mentioning your term (y). Genes in the upper right are both&#10;well studied and specifically associated with your topic. A &#39;+&#39; marks genes that&#10;are present in your current feature selection. Ranking reflects literature bias,&#10;not statistical evidence from your data.", 
        icon = "book-open", label = "Geneshot", sr_title = "Geneshot Literature Search", 
        aria = "Literature-based gene discovery using PubMed co-mention analysis to find genes associated with search terms"), 
    list(value = "Sample", title = "Selected-sample analysis.&#10;Analyses the currently selected samples against the rest of the cohort:&#10;beeswarm plots for continuous variables, contingency tables of categorical&#10;annotations (chi-square / Fisher&#39;s exact test), and Kaplan-Meier survival&#10;curves with the log-rank test when time-to-event data are available.&#10;Interpretation: a small contingency p-value means the selected samples differ in&#10;that annotation more than expected by chance. Kaplan-Meier curves show the&#10;event-free proportion over time per group; the log-rank p-value tests whether the&#10;curves differ. Check group sizes (n) - small groups make tests unreliable.", 
        icon = "arrows-down-to-people", label = "Sample", sr_title = "Sample Analysis", 
        aria = "Statistical analysis of selected samples including group comparisons, contingency tables, and Kaplan-Meier survival curves"))
WP12_TABS_NULLFD <- list(list(value = "Feature", title = "Selected-feature analysis.&#10;Compares the abundance of the currently selected feature(s) across sample groups:&#10;beeswarm/box plots with group-comparison statistics (t-test or Wilcoxon for two groups,&#10;ANOVA/Kruskal-Wallis for more), or ROC / precision-recall curves when a binary&#10;outcome is chosen.&#10;Interpretation: the box spans the interquartile range with the median as center line.&#10;A small p-value (e.g. &lt; 0.05) indicates a statistically significant group difference -&#10;always judge effect size and overlap on the plot too. For ROC/PR curves: AUC near 1.0&#10;means the feature separates the two classes well; 0.5 is no better than chance.", 
    icon = "chart-column", label = "Feature", sr_title = "Feature Analysis", 
    aria = "Statistical analysis and visualization of selected features including boxplots with group comparisons and ROC curves for binary outcomes"), 
    list(value = "Geneshot", title = "Literature-based gene discovery.&#10;Searches PubMed-associated literature (AutoRIF) for genes co-mentioned with your&#10;search term(s); separate multiple terms with &#39;;&#39;. Useful for hypothesis generation.&#10;Interpretation: the scatter plot shows publication count (x) versus the share of a&#10;gene&#39;s publications mentioning your term (y). Genes in the upper right are both&#10;well studied and specifically associated with your topic. A &#39;+&#39; marks genes that&#10;are present in your current feature selection. Ranking reflects literature bias,&#10;not statistical evidence from your data.", 
        icon = "book-open", label = "Geneshot", sr_title = "Geneshot Literature Search", 
        aria = "Literature-based gene discovery using PubMed co-mention analysis to find genes associated with search terms"), 
    list(value = "Sample", title = "Selected-sample analysis.&#10;Analyses the currently selected samples against the rest of the cohort:&#10;beeswarm plots for continuous variables, contingency tables of categorical&#10;annotations (chi-square / Fisher&#39;s exact test), and Kaplan-Meier survival&#10;curves with the log-rank test when time-to-event data are available.&#10;Interpretation: a small contingency p-value means the selected samples differ in&#10;that annotation more than expected by chance. Kaplan-Meier curves show the&#10;event-free proportion over time per group; the log-rank p-value tests whether the&#10;curves differ. Check group sizes (n) - small groups make tests unreliable.", 
        icon = "arrows-down-to-people", label = "Sample", sr_title = "Sample Analysis", 
        aria = "Statistical analysis of selected samples including group comparisons, contingency tables, and Kaplan-Meier survival curves"))
WP12_TABS_DEMO_ADDL <- list(list(value = "Feature", title = "Selected-feature analysis.&#10;Compares the abundance of the currently selected feature(s) across sample groups:&#10;beeswarm/box plots with group-comparison statistics (t-test or Wilcoxon for two groups,&#10;ANOVA/Kruskal-Wallis for more), or ROC / precision-recall curves when a binary&#10;outcome is chosen.&#10;Interpretation: the box spans the interquartile range with the median as center line.&#10;A small p-value (e.g. &lt; 0.05) indicates a statistically significant group difference -&#10;always judge effect size and overlap on the plot too. For ROC/PR curves: AUC near 1.0&#10;means the feature separates the two classes well; 0.5 is no better than chance.", 
    icon = "chart-column", label = "Feature", sr_title = "Feature Analysis", 
    aria = "Statistical analysis and visualization of selected features including boxplots with group comparisons and ROC curves for binary outcomes"), 
    list(value = "Response", title = "Dose-response curves (4-parameter logistic fit).&#10;Fits abundance vs concentration for the selected features and estimates EC50/IC50&#10;(concentration giving half-maximal effect), Hill slope (curve steepness) and the&#10;top/bottom response plateaus.&#10;Interpretation: a lower EC50/IC50 means higher potency; the Hill slope shows how&#10;abruptly the response switches; the parameter table gives fitted values with&#10;confidence intervals. Points are observed data, the line the model fit.", 
        icon = "capsules", label = "Response", sr_title = "Dose-Response Curves", 
        aria = "Dose-response curve fitting with EC50 and IC50 estimation using 4-parameter logistic regression model"), 
    list(value = "StringDB", title = "Protein-protein interaction network (STRING database).&#10;Queries STRING for known and predicted interactions among the currently selected&#10;features (max ~300; set the taxonomy code, e.g. 9606 = human) and draws the&#10;network - nodes are proteins, edges are interactions. Requires internet access.&#10;Interpretation: edge confidence reflects combined evidence (experiments,&#10;databases, text mining); dense clusters often correspond to functional modules or&#10;complexes; highly connected hub proteins may be key regulators.", 
        icon = "share-nodes", label = "StringDB", sr_title = "STRING Protein Interaction Network", 
        aria = "Protein-protein interaction network from STRING database showing physical and functional associations"), 
    list(value = "SeqLogo", title = "PTM motif enrichment (sequence logo).&#10;Takes the amino-acid sequence windows around the selected modification sites and&#10;tests which residues are over-represented at each position, shown as a sequence&#10;logo around the modification site (position 0).&#10;Interpretation: taller letters mean the amino acid (or property class) is more&#10;conserved/enriched at that position; significant enrichment versus all sites in&#10;the dataset hints at the kinase or recognition machinery involved.", 
        icon = "font", label = "SeqLogo", sr_title = "Post-Translational Modification Motifs", 
        aria = "PTM motif enrichment analysis identifying overrepresented amino acid sequence patterns around modification sites"), 
    list(value = "Custom1", title = "<li>\n          <a href=\"#tab-TSID-5\" data-toggle=\"tab\" data-bs-toggle=\"tab\" data-value=\"Custom1\">Custom1</a>\n        </li>", 
        icon = "<li>\n          <a href=\"#tab-TSID-5\" data-toggle=\"tab\" data-bs-toggle=\"tab\" data-value=\"Custom1\">Custom1</a>\n        </li>", 
        label = "Custom1", sr_title = NA_character_, aria = NA_character_), list(
        value = "Custom2", title = "<li>\n          <a href=\"#tab-TSID-6\" data-toggle=\"tab\" data-bs-toggle=\"tab\" data-value=\"Custom2\">Custom2</a>\n        </li>", 
        icon = "<li>\n          <a href=\"#tab-TSID-6\" data-toggle=\"tab\" data-bs-toggle=\"tab\" data-value=\"Custom2\">Custom2</a>\n        </li>", 
        label = "Custom2", sr_title = NA_character_, aria = NA_character_), list(
        value = "Geneshot", title = "Literature-based gene discovery.&#10;Searches PubMed-associated literature (AutoRIF) for genes co-mentioned with your&#10;search term(s); separate multiple terms with &#39;;&#39;. Useful for hypothesis generation.&#10;Interpretation: the scatter plot shows publication count (x) versus the share of a&#10;gene&#39;s publications mentioning your term (y). Genes in the upper right are both&#10;well studied and specifically associated with your topic. A &#39;+&#39; marks genes that&#10;are present in your current feature selection. Ranking reflects literature bias,&#10;not statistical evidence from your data.", 
        icon = "book-open", label = "Geneshot", sr_title = "Geneshot Literature Search", 
        aria = "Literature-based gene discovery using PubMed co-mention analysis to find genes associated with search terms"), 
    list(value = "Sample", title = "Selected-sample analysis.&#10;Analyses the currently selected samples against the rest of the cohort:&#10;beeswarm plots for continuous variables, contingency tables of categorical&#10;annotations (chi-square / Fisher&#39;s exact test), and Kaplan-Meier survival&#10;curves with the log-rank test when time-to-event data are available.&#10;Interpretation: a small contingency p-value means the selected samples differ in&#10;that annotation more than expected by chance. Kaplan-Meier curves show the&#10;event-free proportion over time per group; the log-rank p-value tests whether the&#10;curves differ. Check group sizes (n) - small groups make tests unreliable.", 
        icon = "arrows-down-to-people", label = "Sample", sr_title = "Sample Analysis", 
        aria = "Statistical analysis of selected samples including group comparisons, contingency tables, and Kaplan-Meier survival curves"))
WP12_TABS_GS_ADDL <- list(list(value = "Feature", title = "Selected-feature analysis.&#10;Compares the abundance of the currently selected feature(s) across sample groups:&#10;beeswarm/box plots with group-comparison statistics (t-test or Wilcoxon for two groups,&#10;ANOVA/Kruskal-Wallis for more), or ROC / precision-recall curves when a binary&#10;outcome is chosen.&#10;Interpretation: the box spans the interquartile range with the median as center line.&#10;A small p-value (e.g. &lt; 0.05) indicates a statistically significant group difference -&#10;always judge effect size and overlap on the plot too. For ROC/PR curves: AUC near 1.0&#10;means the feature separates the two classes well; 0.5 is no better than chance.", 
    icon = "chart-column", label = "Feature", sr_title = "Feature Analysis", 
    aria = "Statistical analysis and visualization of selected features including boxplots with group comparisons and ROC curves for binary outcomes"), 
    list(value = "ORA", title = "Over-representation analysis (hypergeometric test).&#10;Tests whether your currently selected features belong to a gene set (pathway, GO&#10;category, custom set) more often than expected by chance. Needs a well-defined&#10;selection, e.g. significant hits from the Feature tab.&#10;Interpretation: the p-value is the chance probability of the observed overlap; FDR&#10;corrects for testing many gene sets (FDR &lt; 0.05 is a common significance cutoff);&#10;the odds ratio / gene ratio measures enrichment strength. The second table lists&#10;the overlapping features per gene set. Results depend on selection size and the&#10;background gene set collection.", 
        icon = "object-group", label = "ORA", sr_title = "Over-Representation Analysis", 
        aria = "Gene set over-representation analysis using hypergeometric test to identify enriched pathways and functional categories"), 
    list(value = "fGSEA", title = "Fast gene set enrichment analysis (GSEA on a ranked list).&#10;Unlike ORA, fGSEA uses ALL features, ranked by a statistic you choose (e.g. fold&#10;change or signed p-value), and asks whether a gene set&#39;s members cluster at the&#10;top or bottom of that ranking - no hard significance cutoff needed.&#10;Interpretation: NES (normalized enrichment score) &gt; 0 means the set is enriched&#10;among up-ranked features, NES &lt; 0 depleted/enriched at the other end; padj &lt; 0.05&#10;is a common significance threshold; the leading edge lists the members driving&#10;the enrichment signal.", 
        icon = "arrow-down-short-wide", label = "fGSEA", sr_title = "Fast Gene Set Enrichment Analysis", 
        aria = "Ranked gene set enrichment analysis computing normalized enrichment scores and identifying leading edge genes"), 
    list(value = "Response", title = "Dose-response curves (4-parameter logistic fit).&#10;Fits abundance vs concentration for the selected features and estimates EC50/IC50&#10;(concentration giving half-maximal effect), Hill slope (curve steepness) and the&#10;top/bottom response plateaus.&#10;Interpretation: a lower EC50/IC50 means higher potency; the Hill slope shows how&#10;abruptly the response switches; the parameter table gives fitted values with&#10;confidence intervals. Points are observed data, the line the model fit.", 
        icon = "capsules", label = "Response", sr_title = "Dose-Response Curves", 
        aria = "Dose-response curve fitting with EC50 and IC50 estimation using 4-parameter logistic regression model"), 
    list(value = "StringDB", title = "Protein-protein interaction network (STRING database).&#10;Queries STRING for known and predicted interactions among the currently selected&#10;features (max ~300; set the taxonomy code, e.g. 9606 = human) and draws the&#10;network - nodes are proteins, edges are interactions. Requires internet access.&#10;Interpretation: edge confidence reflects combined evidence (experiments,&#10;databases, text mining); dense clusters often correspond to functional modules or&#10;complexes; highly connected hub proteins may be key regulators.", 
        icon = "share-nodes", label = "StringDB", sr_title = "STRING Protein Interaction Network", 
        aria = "Protein-protein interaction network from STRING database showing physical and functional associations"), 
    list(value = "SeqLogo", title = "PTM motif enrichment (sequence logo).&#10;Takes the amino-acid sequence windows around the selected modification sites and&#10;tests which residues are over-represented at each position, shown as a sequence&#10;logo around the modification site (position 0).&#10;Interpretation: taller letters mean the amino acid (or property class) is more&#10;conserved/enriched at that position; significant enrichment versus all sites in&#10;the dataset hints at the kinase or recognition machinery involved.", 
        icon = "font", label = "SeqLogo", sr_title = "Post-Translational Modification Motifs", 
        aria = "PTM motif enrichment analysis identifying overrepresented amino acid sequence patterns around modification sites"), 
    list(value = "Custom1", title = "<li>\n          <a href=\"#tab-TSID-7\" data-toggle=\"tab\" data-bs-toggle=\"tab\" data-value=\"Custom1\">Custom1</a>\n        </li>", 
        icon = "<li>\n          <a href=\"#tab-TSID-7\" data-toggle=\"tab\" data-bs-toggle=\"tab\" data-value=\"Custom1\">Custom1</a>\n        </li>", 
        label = "Custom1", sr_title = NA_character_, aria = NA_character_), list(
        value = "Custom2", title = "<li>\n          <a href=\"#tab-TSID-8\" data-toggle=\"tab\" data-bs-toggle=\"tab\" data-value=\"Custom2\">Custom2</a>\n        </li>", 
        icon = "<li>\n          <a href=\"#tab-TSID-8\" data-toggle=\"tab\" data-bs-toggle=\"tab\" data-value=\"Custom2\">Custom2</a>\n        </li>", 
        label = "Custom2", sr_title = NA_character_, aria = NA_character_), list(
        value = "Geneshot", title = "Literature-based gene discovery.&#10;Searches PubMed-associated literature (AutoRIF) for genes co-mentioned with your&#10;search term(s); separate multiple terms with &#39;;&#39;. Useful for hypothesis generation.&#10;Interpretation: the scatter plot shows publication count (x) versus the share of a&#10;gene&#39;s publications mentioning your term (y). Genes in the upper right are both&#10;well studied and specifically associated with your topic. A &#39;+&#39; marks genes that&#10;are present in your current feature selection. Ranking reflects literature bias,&#10;not statistical evidence from your data.", 
        icon = "book-open", label = "Geneshot", sr_title = "Geneshot Literature Search", 
        aria = "Literature-based gene discovery using PubMed co-mention analysis to find genes associated with search terms"), 
    list(value = "Sample", title = "Selected-sample analysis.&#10;Analyses the currently selected samples against the rest of the cohort:&#10;beeswarm plots for continuous variables, contingency tables of categorical&#10;annotations (chi-square / Fisher&#39;s exact test), and Kaplan-Meier survival&#10;curves with the log-rank test when time-to-event data are available.&#10;Interpretation: a small contingency p-value means the selected samples differ in&#10;that annotation more than expected by chance. Kaplan-Meier curves show the&#10;event-free proportion over time per group; the log-rank p-value tests whether the&#10;curves differ. Check group sizes (n) - small groups make tests unreliable.", 
        icon = "arrows-down-to-people", label = "Sample", sr_title = "Sample Analysis", 
        aria = "Statistical analysis of selected samples including group comparisons, contingency tables, and Kaplan-Meier survival curves"))

WP12_DEMO_HTML <- structure("<head>\n  <title>Analysis</title>\n</head>\n<head>\n  <link rel=\"stylesheet\" type=\"text/css\" href=\"shinythemes/css/spacelab.min.css\"/>\n</head>\n<body class=\"bslib-page-navbar\">\n  <nav class=\"navbar navbar-default navbar-static-top\" role=\"navigation\" data-bs-theme=\"light\">\n    <div class=\"container-fluid\">\n      <div class=\"navbar-header\">\n        <span class=\"navbar-brand\">\n          <span class=\"omicsviewer-nav-title\">Analysis</span>\n          <span class=\"omicsviewer-nav-sub\"> — analyse the current selection on the fly</span>\n        </span>\n      </div>\n      <ul class=\"nav navbar-nav nav-underline shiny-tab-input\" id=\"rs-analyst\" data-tabsetid=\"TSID\">\n        <li class=\"active\">\n          <a href=\"#tab-TSID-1\" data-toggle=\"tab\" data-bs-toggle=\"tab\" data-value=\"Feature\">\n            <span title=\"Selected-feature analysis.&#10;Compares the abundance of the currently selected feature(s) across sample groups:&#10;beeswarm/box plots with group-comparison statistics (t-test or Wilcoxon for two groups,&#10;ANOVA/Kruskal-Wallis for more), or ROC / precision-recall curves when a binary&#10;outcome is chosen.&#10;Interpretation: the box spans the interquartile range with the median as center line.&#10;A small p-value (e.g. &lt; 0.05) indicates a statistically significant group difference -&#10;always judge effect size and overlap on the plot too. For ROC/PR curves: AUC near 1.0&#10;means the feature separates the two classes well; 0.5 is no better than chance.\">\n              <i class=\"fas fa-chart-column\" role=\"presentation\" aria-label=\"chart-column icon\"></i>\n              Feature\n            </span>\n          </a>\n        </li>\n        <li>\n          <a href=\"#tab-TSID-2\" data-toggle=\"tab\" data-bs-toggle=\"tab\" data-value=\"Response\">\n            <span title=\"Dose-response curves (4-parameter logistic fit).&#10;Fits abundance vs concentration for the selected features and estimates EC50/IC50&#10;(concentration giving half-maximal effect), Hill slope (curve steepness) and the&#10;top/bottom response plateaus.&#10;Interpretation: a lower EC50/IC50 means higher potency; the Hill slope shows how&#10;abruptly the response switches; the parameter table gives fitted values with&#10;confidence intervals. Points are observed data, the line the model fit.\">\n              <i class=\"fas fa-capsules\" role=\"presentation\" aria-label=\"capsules icon\"></i>\n              Response\n            </span>\n          </a>\n        </li>\n        <li>\n          <a href=\"#tab-TSID-3\" data-toggle=\"tab\" data-bs-toggle=\"tab\" data-value=\"StringDB\">\n            <span title=\"Protein-protein interaction network (STRING database).&#10;Queries STRING for known and predicted interactions among the currently selected&#10;features (max ~300; set the taxonomy code, e.g. 9606 = human) and draws the&#10;network - nodes are proteins, edges are interactions. Requires internet access.&#10;Interpretation: edge confidence reflects combined evidence (experiments,&#10;databases, text mining); dense clusters often correspond to functional modules or&#10;complexes; highly connected hub proteins may be key regulators.\">\n              <i class=\"fas fa-share-nodes\" role=\"presentation\" aria-label=\"share-nodes icon\"></i>\n              StringDB\n            </span>\n          </a>\n        </li>\n        <li>\n          <a href=\"#tab-TSID-4\" data-toggle=\"tab\" data-bs-toggle=\"tab\" data-value=\"SeqLogo\">\n            <span title=\"PTM motif enrichment (sequence logo).&#10;Takes the amino-acid sequence windows around the selected modification sites and&#10;tests which residues are over-represented at each position, shown as a sequence&#10;logo around the modification site (position 0).&#10;Interpretation: taller letters mean the amino acid (or property class) is more&#10;conserved/enriched at that position; significant enrichment versus all sites in&#10;the dataset hints at the kinase or recognition machinery involved.\">\n              <i class=\"fas fa-font\" role=\"presentation\" aria-label=\"font icon\"></i>\n              SeqLogo\n            </span>\n          </a>\n        </li>\n        <li>\n          <a href=\"#tab-TSID-5\" data-toggle=\"tab\" data-bs-toggle=\"tab\" data-value=\"Geneshot\">\n            <span title=\"Literature-based gene discovery.&#10;Searches PubMed-associated literature (AutoRIF) for genes co-mentioned with your&#10;search term(s); separate multiple terms with &#39;;&#39;. Useful for hypothesis generation.&#10;Interpretation: the scatter plot shows publication count (x) versus the share of a&#10;gene&#39;s publications mentioning your term (y). Genes in the upper right are both&#10;well studied and specifically associated with your topic. A &#39;+&#39; marks genes that&#10;are present in your current feature selection. Ranking reflects literature bias,&#10;not statistical evidence from your data.\">\n              <i class=\"fas fa-book-open\" role=\"presentation\" aria-label=\"book-open icon\"></i>\n              Geneshot\n            </span>\n          </a>\n        </li>\n        <li>\n          <a href=\"#tab-TSID-6\" data-toggle=\"tab\" data-bs-toggle=\"tab\" data-value=\"Sample\">\n            <span title=\"Selected-sample analysis.&#10;Analyses the currently selected samples against the rest of the cohort:&#10;beeswarm plots for continuous variables, contingency tables of categorical&#10;annotations (chi-square / Fisher&#39;s exact test), and Kaplan-Meier survival&#10;curves with the log-rank test when time-to-event data are available.&#10;Interpretation: a small contingency p-value means the selected samples differ in&#10;that annotation more than expected by chance. Kaplan-Meier curves show the&#10;event-free proportion over time per group; the log-rank p-value tests whether the&#10;curves differ. Check group sizes (n) - small groups make tests unreliable.\">\n              <i class=\"fas fa-arrows-down-to-people\" role=\"presentation\" aria-label=\"arrows-down-to-people icon\"></i>\n              Sample\n            </span>\n          </a>\n        </li>\n      </ul>\n    </div>\n  </nav>\n  <div class=\"container-fluid\">\n    <div class=\"tab-content\" data-tabsetid=\"TSID\">\n      <div class=\"tab-pane active\" data-value=\"Feature\" id=\"tab-TSID-1\">\n        <h2 class=\"sr-only\" aria-label=\"Statistical analysis and visualization of selected features including boxplots with group comparisons and ROC curves for binary outcomes\">Feature Analysis</h2>\n        <div class=\"sr-only\" id=\"rs-feature_general-module-help\">\n          <h4>About Feature Statistical Analysis</h4>\n          <p>This module performs comprehensive statistical analysis to compare selected features (genes/proteins) across different sample groups. It provides multiple visualization options including beeswarm boxplots for distribution comparison and ROC/PR curves for classification performance evaluation.</p>\n          <h4>When to use this analysis</h4>\n          <p>Use this module when you want to test if a feature shows significant differences between groups (e.g., disease vs. control, treatment vs. untreated). The beeswarm plot is ideal for visualizing distributions and identifying outliers, while ROC/PR curves help assess how well a feature can distinguish between two classes.</p>\n          <h4>How to interpret results</h4>\n          <p>For boxplots: The box shows the interquartile range (25th to 75th percentile), the line inside is the median, and individual points show data distributions. P-values indicate statistical significance of group differences. For ROC curves: Area Under Curve (AUC) near 1.0 indicates excellent classification, 0.5 is random chance. The curve shows trade-off between sensitivity and specificity at different thresholds.</p>\n        </div>\n        <h3 class=\"sr-only\" aria-label=\"Controls for selecting statistical test results and choosing between beeswarm boxplot or ROC curve visualization\">Variable Selection and Plot Type</h3>\n        <div class=\"row\">\n          <div class=\"col-sm-12\" style=\"margin-top: 0px;\">\n            <div class=\"row\">\n              <div class=\"col-sm-2\" align=\"right\" style=\"padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px\">\n                <div id=\"rs-feature_general-tris_feature_general-groupLabel\" class=\"shiny-html-output\"></div>\n              </div>\n              <div class=\"col-sm-3\" style=\"padding:2px;\">\n                <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-feature_general-tris_feature_general-analysis-selector\">\n                  <label class=\"control-label shiny-label-null\" for=\"rs-feature_general-tris_feature_general-analysis\" id=\"rs-feature_general-tris_feature_general-analysis-label\"></label>\n                  <div>\n                    <select id=\"rs-feature_general-tris_feature_general-analysis\" class=\"shiny-input-select\"></select>\n                    <script type=\"application/json\" data-for=\"rs-feature_general-tris_feature_general-analysis\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                  </div>\n                </div>\n              </div>\n              <div class=\"col-sm-4\" style=\"padding:2px;\">\n                <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-feature_general-tris_feature_general-subset-selector\">\n                  <label class=\"control-label shiny-label-null\" for=\"rs-feature_general-tris_feature_general-subset\" id=\"rs-feature_general-tris_feature_general-subset-label\"></label>\n                  <div>\n                    <select id=\"rs-feature_general-tris_feature_general-subset\" class=\"shiny-input-select\"></select>\n                    <script type=\"application/json\" data-for=\"rs-feature_general-tris_feature_general-subset\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                  </div>\n                </div>\n              </div>\n              <div class=\"col-sm-3\" style=\"padding-left:2px; padding-right:5px; padding-top:2px; padding-bottom:2px\">\n                <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-feature_general-tris_feature_general-variable-selector\">\n                  <label class=\"control-label shiny-label-null\" for=\"rs-feature_general-tris_feature_general-variable\" id=\"rs-feature_general-tris_feature_general-variable-label\"></label>\n                  <div>\n                    <select id=\"rs-feature_general-tris_feature_general-variable\" class=\"shiny-input-select\"></select>\n                    <script type=\"application/json\" data-for=\"rs-feature_general-tris_feature_general-variable\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                  </div>\n                </div>\n              </div>\n            </div>\n          </div>\n          <div class=\"col-sm-11\">\n            <div id=\"rs-feature_general-feature_general_plot\" class=\"shiny-html-output\"></div>\n          </div>\n          <div class=\"col-sm-1\">\n            <div class=\"sw-dropdown\" id=\"sw-drop-btn-RID\">\n              <button class=\"btn btn-default  action-button\" data-bs-html=\"false\" data-bs-placement=\"right\" data-bs-title=\"Click to modify figure!\" data-bs-toggle=\"tooltip\" id=\"btn-RID\" type=\"button\">\n                <i class=\"fas fa-gear\" role=\"presentation\" aria-label=\"gear icon\"></i>\n                <span class=\"glyphicon glyphicon-triangle-bottom\"></span>\n              </button>\n              <div class=\"sw-dropdown-content animated sw-dropright-content\" id=\"sw-content-btn-RID\" style=\"width:788px;\">\n                <div class=\"sw-dropdown-in\" margin=\"25px\" circle=\"FALSE\">\n                  <br/>\n                  <div class=\"row\">\n                    <div class=\"col-sm-2\" align=\"right\" style=\"padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px\">\n                      <div id=\"rs-feature_general-a4_gf-selectColorUI-groupLabel\" class=\"shiny-html-output\"></div>\n                    </div>\n                    <div class=\"col-sm-3\" style=\"padding:2px;\">\n                      <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-feature_general-a4_gf-selectColorUI-analysis-selector\">\n                        <label class=\"control-label shiny-label-null\" for=\"rs-feature_general-a4_gf-selectColorUI-analysis\" id=\"rs-feature_general-a4_gf-selectColorUI-analysis-label\"></label>\n                        <div>\n                          <select id=\"rs-feature_general-a4_gf-selectColorUI-analysis\" class=\"shiny-input-select\"></select>\n                          <script type=\"application/json\" data-for=\"rs-feature_general-a4_gf-selectColorUI-analysis\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                        </div>\n                      </div>\n                    </div>\n                    <div class=\"col-sm-4\" style=\"padding:2px;\">\n                      <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-feature_general-a4_gf-selectColorUI-subset-selector\">\n                        <label class=\"control-label shiny-label-null\" for=\"rs-feature_general-a4_gf-selectColorUI-subset\" id=\"rs-feature_general-a4_gf-selectColorUI-subset-label\"></label>\n                        <div>\n                          <select id=\"rs-feature_general-a4_gf-selectColorUI-subset\" class=\"shiny-input-select\"></select>\n                          <script type=\"application/json\" data-for=\"rs-feature_general-a4_gf-selectColorUI-subset\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                        </div>\n                      </div>\n                    </div>\n                    <div class=\"col-sm-3\" style=\"padding-left:2px; padding-right:20px; padding-top:2px; padding-bottom:2px\">\n                      <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-feature_general-a4_gf-selectColorUI-variable-selector\">\n                        <label class=\"control-label shiny-label-null\" for=\"rs-feature_general-a4_gf-selectColorUI-variable\" id=\"rs-feature_general-a4_gf-selectColorUI-variable-label\"></label>\n                        <div>\n                          <select id=\"rs-feature_general-a4_gf-selectColorUI-variable\" class=\"shiny-input-select\"></select>\n                          <script type=\"application/json\" data-for=\"rs-feature_general-a4_gf-selectColorUI-variable\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                        </div>\n                      </div>\n                    </div>\n                  </div>\n                  <div class=\"row\">\n                    <div class=\"col-sm-2\" align=\"right\" style=\"padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px\">\n                      <div id=\"rs-feature_general-a4_gf-selectShapeUI-groupLabel\" class=\"shiny-html-output\"></div>\n                    </div>\n                    <div class=\"col-sm-3\" style=\"padding:2px;\">\n                      <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-feature_general-a4_gf-selectShapeUI-analysis-selector\">\n                        <label class=\"control-label shiny-label-null\" for=\"rs-feature_general-a4_gf-selectShapeUI-analysis\" id=\"rs-feature_general-a4_gf-selectShapeUI-analysis-label\"></label>\n                        <div>\n                          <select id=\"rs-feature_general-a4_gf-selectShapeUI-analysis\" class=\"shiny-input-select\"></select>\n                          <script type=\"application/json\" data-for=\"rs-feature_general-a4_gf-selectShapeUI-analysis\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                        </div>\n                      </div>\n                    </div>\n                    <div class=\"col-sm-4\" style=\"padding:2px;\">\n                      <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-feature_general-a4_gf-selectShapeUI-subset-selector\">\n                        <label class=\"control-label shiny-label-null\" for=\"rs-feature_general-a4_gf-selectShapeUI-subset\" id=\"rs-feature_general-a4_gf-selectShapeUI-subset-label\"></label>\n                        <div>\n                          <select id=\"rs-feature_general-a4_gf-selectShapeUI-subset\" class=\"shiny-input-select\"></select>\n                          <script type=\"application/json\" data-for=\"rs-feature_general-a4_gf-selectShapeUI-subset\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                        </div>\n                      </div>\n                    </div>\n                    <div class=\"col-sm-3\" style=\"padding-left:2px; padding-right:20px; padding-top:2px; padding-bottom:2px\">\n                      <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-feature_general-a4_gf-selectShapeUI-variable-selector\">\n                        <label class=\"control-label shiny-label-null\" for=\"rs-feature_general-a4_gf-selectShapeUI-variable\" id=\"rs-feature_general-a4_gf-selectShapeUI-variable-label\"></label>\n                        <div>\n                          <select id=\"rs-feature_general-a4_gf-selectShapeUI-variable\" class=\"shiny-input-select\"></select>\n                          <script type=\"application/json\" data-for=\"rs-feature_general-a4_gf-selectShapeUI-variable\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                        </div>\n                      </div>\n                    </div>\n                  </div>\n                  <div class=\"row\">\n                    <div class=\"col-sm-2\" align=\"right\" style=\"padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px\">\n                      <div id=\"rs-feature_general-a4_gf-selectSizeUI-groupLabel\" class=\"shiny-html-output\"></div>\n                    </div>\n                    <div class=\"col-sm-3\" style=\"padding:2px;\">\n                      <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-feature_general-a4_gf-selectSizeUI-analysis-selector\">\n                        <label class=\"control-label shiny-label-null\" for=\"rs-feature_general-a4_gf-selectSizeUI-analysis\" id=\"rs-feature_general-a4_gf-selectSizeUI-analysis-label\"></label>\n                        <div>\n                          <select id=\"rs-feature_general-a4_gf-selectSizeUI-analysis\" class=\"shiny-input-select\"></select>\n                          <script type=\"application/json\" data-for=\"rs-feature_general-a4_gf-selectSizeUI-analysis\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                        </div>\n                      </div>\n                    </div>\n                    <div class=\"col-sm-4\" style=\"padding:2px;\">\n                      <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-feature_general-a4_gf-selectSizeUI-subset-selector\">\n                        <label class=\"control-label shiny-label-null\" for=\"rs-feature_general-a4_gf-selectSizeUI-subset\" id=\"rs-feature_general-a4_gf-selectSizeUI-subset-label\"></label>\n                        <div>\n                          <select id=\"rs-feature_general-a4_gf-selectSizeUI-subset\" class=\"shiny-input-select\"></select>\n                          <script type=\"application/json\" data-for=\"rs-feature_general-a4_gf-selectSizeUI-subset\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                        </div>\n                      </div>\n                    </div>\n                    <div class=\"col-sm-3\" style=\"padding-left:2px; padding-right:20px; padding-top:2px; padding-bottom:2px\">\n                      <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-feature_general-a4_gf-selectSizeUI-variable-selector\">\n                        <label class=\"control-label shiny-label-null\" for=\"rs-feature_general-a4_gf-selectSizeUI-variable\" id=\"rs-feature_general-a4_gf-selectSizeUI-variable-label\"></label>\n                        <div>\n                          <select id=\"rs-feature_general-a4_gf-selectSizeUI-variable\" class=\"shiny-input-select\"></select>\n                          <script type=\"application/json\" data-for=\"rs-feature_general-a4_gf-selectSizeUI-variable\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                        </div>\n                      </div>\n                    </div>\n                  </div>\n                  <div class=\"row\">\n                    <div class=\"col-sm-2\" align=\"right\" style=\"padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px\">\n                      <div id=\"rs-feature_general-a4_gf-selectTooltipUI-groupLabel\" class=\"shiny-html-output\"></div>\n                    </div>\n                    <div class=\"col-sm-3\" style=\"padding:2px;\">\n                      <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-feature_general-a4_gf-selectTooltipUI-analysis-selector\">\n                        <label class=\"control-label shiny-label-null\" for=\"rs-feature_general-a4_gf-selectTooltipUI-analysis\" id=\"rs-feature_general-a4_gf-selectTooltipUI-analysis-label\"></label>\n                        <div>\n                          <select id=\"rs-feature_general-a4_gf-selectTooltipUI-analysis\" class=\"shiny-input-select\"></select>\n                          <script type=\"application/json\" data-for=\"rs-feature_general-a4_gf-selectTooltipUI-analysis\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                        </div>\n                      </div>\n                    </div>\n                    <div class=\"col-sm-4\" style=\"padding:2px;\">\n                      <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-feature_general-a4_gf-selectTooltipUI-subset-selector\">\n                        <label class=\"control-label shiny-label-null\" for=\"rs-feature_general-a4_gf-selectTooltipUI-subset\" id=\"rs-feature_general-a4_gf-selectTooltipUI-subset-label\"></label>\n                        <div>\n                          <select id=\"rs-feature_general-a4_gf-selectTooltipUI-subset\" class=\"shiny-input-select\"></select>\n                          <script type=\"application/json\" data-for=\"rs-feature_general-a4_gf-selectTooltipUI-subset\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                        </div>\n                      </div>\n                    </div>\n                    <div class=\"col-sm-3\" style=\"padding-left:2px; padding-right:20px; padding-top:2px; padding-bottom:2px\">\n                      <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-feature_general-a4_gf-selectTooltipUI-variable-selector\">\n                        <label class=\"control-label shiny-label-null\" for=\"rs-feature_general-a4_gf-selectTooltipUI-variable\" id=\"rs-feature_general-a4_gf-selectTooltipUI-variable-label\"></label>\n                        <div>\n                          <select id=\"rs-feature_general-a4_gf-selectTooltipUI-variable\" class=\"shiny-input-select\"></select>\n                          <script type=\"application/json\" data-for=\"rs-feature_general-a4_gf-selectTooltipUI-variable\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                        </div>\n                      </div>\n                    </div>\n                  </div>\n                  <div class=\"row\">\n                    <div class=\"col-sm-2\" align=\"right\" style=\"padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px\">\n                      <div id=\"rs-feature_general-a4_gf-selectSearchCol-groupLabel\" class=\"shiny-html-output\"></div>\n                    </div>\n                    <div class=\"col-sm-3\" style=\"padding:2px;\">\n                      <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-feature_general-a4_gf-selectSearchCol-analysis-selector\">\n                        <label class=\"control-label shiny-label-null\" for=\"rs-feature_general-a4_gf-selectSearchCol-analysis\" id=\"rs-feature_general-a4_gf-selectSearchCol-analysis-label\"></label>\n                        <div>\n                          <select id=\"rs-feature_general-a4_gf-selectSearchCol-analysis\" class=\"shiny-input-select\"></select>\n                          <script type=\"application/json\" data-for=\"rs-feature_general-a4_gf-selectSearchCol-analysis\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                        </div>\n                      </div>\n                    </div>\n                    <div class=\"col-sm-4\" style=\"padding:2px;\">\n                      <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-feature_general-a4_gf-selectSearchCol-subset-selector\">\n                        <label class=\"control-label shiny-label-null\" for=\"rs-feature_general-a4_gf-selectSearchCol-subset\" id=\"rs-feature_general-a4_gf-selectSearchCol-subset-label\"></label>\n                        <div>\n                          <select id=\"rs-feature_general-a4_gf-selectSearchCol-subset\" class=\"shiny-input-select\"></select>\n                          <script type=\"application/json\" data-for=\"rs-feature_general-a4_gf-selectSearchCol-subset\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                        </div>\n                      </div>\n                    </div>\n                    <div class=\"col-sm-3\" style=\"padding-left:2px; padding-right:20px; padding-top:2px; padding-bottom:2px\">\n                      <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-feature_general-a4_gf-selectSearchCol-variable-selector\">\n                        <label class=\"control-label shiny-label-null\" for=\"rs-feature_general-a4_gf-selectSearchCol-variable\" id=\"rs-feature_general-a4_gf-selectSearchCol-variable-label\"></label>\n                        <div>\n                          <select id=\"rs-feature_general-a4_gf-selectSearchCol-variable\" class=\"shiny-input-select\"></select>\n                          <script type=\"application/json\" data-for=\"rs-feature_general-a4_gf-selectSearchCol-variable\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                        </div>\n                      </div>\n                    </div>\n                  </div>\n                  <div class=\"shiny-panel-conditional\" data-display-if=\"1 == 2\" data-ns-prefix=\"\">\n                    <div class=\"form-group shiny-input-container\">\n                      <div class=\"checkbox\">\n                        <label>\n                          <input id=\"rs-feature_general-a4_gf-showSearchBox\" type=\"checkbox\" class=\"shiny-input-checkbox\"/>\n                          <span>show</span>\n                        </label>\n                      </div>\n                    </div>\n                  </div>\n                  <div class=\"shiny-panel-conditional\" data-display-if=\"input.showSearchBox == true\" data-ns-prefix=\"rs-feature_general-a4_gf-\">\n                    <div style=\"padding-left:100px; padding-right:0px; padding-top:0px; padding-bottom:0px\">\n                      <div class=\"form-group shiny-input-container\" style=\"width:100%;\">\n                        <label class=\"control-label shiny-label-null\" for=\"rs-feature_general-a4_gf-searchon\" id=\"rs-feature_general-a4_gf-searchon-label\"></label>\n                        <div>\n                          <select id=\"rs-feature_general-a4_gf-searchon\" class=\"shiny-input-select\" multiple=\"multiple\"></select>\n                          <script type=\"application/json\" data-for=\"rs-feature_general-a4_gf-searchon\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                        </div>\n                      </div>\n                    </div>\n                  </div>\n                  <div class=\"row\">\n                    <div class=\"col-sm-2\"></div>\n                    <div class=\"col-sm-4\" style=\"padding:2px;\">\n                      <div class=\"form-group shiny-input-container\">\n                        <label class=\"control-label\" id=\"rs-feature_general-a4_gf-xcut-label\" for=\"rs-feature_general-a4_gf-xcut\">Select points by x/y cutoffs</label>\n                        <div class=\"input-group\">\n                          <div class=\"input-group-addon sw-input-icon input-group-prepend\">\n                            <span class=\"input-group-text\">x-cut</span>\n                          </div>\n                          <input id=\"rs-feature_general-a4_gf-xcut\" type=\"text\" class=\"form-control text-input-icon\" value=\"log10(2)\" placeholder=\"e.g. -1 or -log10(2)\"/>\n                        </div>\n                      </div>\n                    </div>\n                    <div class=\"col-sm-4\" style=\"padding-left:5px; padding-right:2px; padding-top:27px; padding-bottom:2px;\">\n                      <div class=\"form-group shiny-input-container\">\n                        <label class=\"control-label shiny-label-null\" for=\"rs-feature_general-a4_gf-ycut\" id=\"rs-feature_general-a4_gf-ycut-label\"></label>\n                        <div class=\"input-group\">\n                          <div class=\"input-group-addon sw-input-icon input-group-prepend\">\n                            <span class=\"input-group-text\">y-cut</span>\n                          </div>\n                          <input id=\"rs-feature_general-a4_gf-ycut\" type=\"text\" class=\"form-control text-input-icon\" value=\"-log10(0.05)\" placeholder=\"e.g 2 or -log10(0.05)\"/>\n                        </div>\n                      </div>\n                    </div>\n                    <div class=\"col-sm-2\" style=\"padding-left:5px; padding-right:20px; padding-top:4px; padding-bottom:2px;\">\n                      <div class=\"form-group shiny-input-container\">\n                        <label class=\"control-label\" id=\"rs-feature_general-a4_gf-scorner-label\" for=\"rs-feature_general-a4_gf-scorner\">Area</label>\n                        <div>\n                          <select id=\"rs-feature_general-a4_gf-scorner\" class=\"shiny-input-select\"><option value=\"None\" selected>None</option></select>\n                          <script type=\"application/json\" data-for=\"rs-feature_general-a4_gf-scorner\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                        </div>\n                      </div>\n                    </div>\n                  </div>\n                </div>\n              </div>\n              <script>$('#btn-RID').tooltip({ placement: 'right', title: 'Click to modify figure!', html: false });</script>\n              <script>$(function() {swDrop('btn-RID', 'sw-content-btn-RID', 'sw-drop-btn-RID', 'sw-none', 'sw-none', '1');});</script>\n            </div>\n            <div class=\"form-group shiny-input-container shiny-input-radiogroup shiny-input-container-inline\" data-testid=\"rs-feature_general-plot-type-selector\">\n              <label class=\"control-label\" id=\"rs-feature_general-internal_radio-label\" for=\"rs-feature_general-internal_radio\"> </label>\n              <br/>\n              <div id=\"rs-feature_general-internal_radio\" class=\"radio-group-buttons\">\n                <div aria-labelledby=\"rs-feature_general-internal_radio-label\" class=\"btn-group-vertical btn-group-xs btn-group-container-sw\" data-toggle=\"buttons\" role=\"group\">\n                  <div class=\"btn-group btn-group-toggle btn-group-xs\" role=\"group\">\n                    <button class=\"btn radiobtn btn-success active\">\n                      <input type=\"radio\" autocomplete=\"off\" name=\"rs-feature_general-internal_radio\" value=\"Bees\" checked=\"checked\"/>\n                      Bees\n                    </button>\n                  </div>\n                  <div class=\"btn-group btn-group-toggle btn-group-xs\" role=\"group\">\n                    <button class=\"btn radiobtn btn-success\">\n                      <input type=\"radio\" autocomplete=\"off\" name=\"rs-feature_general-internal_radio\" value=\"Curve\"/>\n                      Curve\n                    </button>\n                  </div>\n                </div>\n              </div>\n            </div>\n          </div>\n        </div>\n        <h3 class=\"sr-only\" aria-label=\"Table showing statistical test results including p-values, effect sizes, and group comparisons for selected features\">Statistical Results Table</h3>\n        <div id=\"rs-feature_general-mtab-showButton\" class=\"shiny-html-output\"></div>\n        <div class=\"datatables html-widget html-widget-output shiny-report-size html-fill-item\" id=\"rs-feature_general-mtab-table\" style=\"width:100%;height:auto;\"></div>\n        <div class=\"sr-only\" aria-live=\"polite\" aria-atomic=\"true\">\n          <div id=\"rs-feature_general-plotSummary\" class=\"shiny-html-output\"></div>\n        </div>\n      </div>\n      <div class=\"tab-pane\" data-value=\"Response\" id=\"tab-TSID-2\">\n        <h2 class=\"sr-only\" aria-label=\"Dose-response curve fitting with EC50 and IC50 estimation using 4-parameter logistic regression model\">Dose-Response Curves</h2>\n        <div class=\"sr-only\" id=\"rs-rescurve-module-help\">\n          <h4>About Dose-Response Curve Analysis</h4>\n          <p>Dose-response analysis models the relationship between drug/compound concentration and biological response using a 4-parameter logistic curve. This analysis fits sigmoid curves to your data and estimates key pharmacological parameters including EC50 (half-maximal effective concentration), IC50 (half-maximal inhibitory concentration), Hill slope, and maximum/minimum response values.</p>\n          <h4>When to use dose-response analysis</h4>\n          <p>Use this analysis when you have measured responses across multiple concentrations of a compound or treatment. This is essential for drug discovery, toxicology studies, and understanding concentration-dependent biological effects. It helps determine the potency of compounds and compare their effectiveness.</p>\n          <h4>How to interpret results</h4>\n          <p>The curve plot shows the fitted sigmoid relationship between concentration (x-axis, usually log-scale) and response (y-axis). EC50/IC50 represents the concentration producing 50% of maximum effect - lower values indicate higher potency. The Hill slope indicates the steepness of the curve. The parameter table shows fitted values with their confidence intervals. Top and bottom parameters represent the maximum and minimum response plateaus.</p>\n        </div>\n        <div class=\"row\">\n          <div class=\"col-sm-7\">\n            <div class=\"shiny-plot-output html-fill-item\" id=\"rs-rescurve-plot\" style=\"width:100%;height:550px;\"></div>\n          </div>\n          <div class=\"col-sm-5\" style=\"margin-top: 30px;\">\n            <div id=\"rs-rescurve-param-showButton\" class=\"shiny-html-output\"></div>\n            <div class=\"datatables html-widget html-widget-output shiny-report-size html-fill-item\" id=\"rs-rescurve-param-table\" style=\"width:100%;height:auto;\"></div>\n          </div>\n          <div class=\"col-sm-12\">\n            <div id=\"rs-rescurve-feature-showButton\" class=\"shiny-html-output\"></div>\n            <div class=\"datatables html-widget html-widget-output shiny-report-size html-fill-item\" id=\"rs-rescurve-feature-table\" style=\"width:100%;height:auto;\"></div>\n          </div>\n        </div>\n      </div>\n      <div class=\"tab-pane\" data-value=\"StringDB\" id=\"tab-TSID-3\">\n        <h2 class=\"sr-only\" aria-label=\"Protein-protein interaction network from STRING database showing physical and functional associations\">STRING Protein Interaction Network</h2>\n        <div class=\"sr-only\" id=\"rs-stringdb-module-help\">\n          <h4>About STRING Protein-Protein Interaction Network</h4>\n          <p>STRING (Search Tool for the Retrieval of Interacting Genes/Proteins) is a database of known and predicted protein-protein interactions. This module queries the STRING database to build an interaction network for your selected proteins, showing how they might work together in biological processes.</p>\n          <h4>When to use STRING</h4>\n          <p>Use this analysis when you want to understand how your proteins of interest interact with each other, identify protein complexes, or discover hub proteins that connect many others. This is especially useful for interpreting proteomics data or understanding the molecular machinery involved in a biological process.</p>\n          <h4>How to interpret results</h4>\n          <p>The network visualization shows proteins as nodes and interactions as edges connecting them. The interaction table lists all protein pairs with confidence scores (0-1, higher is better) and evidence types (experimental, database, text mining, etc.). Connected clusters of proteins often represent functional modules or complexes. You can also see pathway enrichment results to understand what biological processes the network is involved in.</p>\n        </div>\n        <h3 class=\"sr-only\" aria-label=\"Input for organism taxonomy code and run button to query STRING database. Maximum 300 proteins allowed\">Query Configuration</h3>\n        <div class=\"row\">\n          <div class=\"col-sm-4\" style=\"padding-left:15px; padding-right:2px; padding-top:0px; padding-bottom:0px\">\n            <div class=\"form-group shiny-input-container\" data-testid=\"rs-stringdb-taxonomy-input\">\n              <label class=\"control-label shiny-label-null\" for=\"rs-stringdb-tax\" id=\"rs-stringdb-tax-label\"></label>\n              <div class=\"input-group\">\n                <div class=\"input-group-addon sw-input-icon input-group-prepend\">\n                  <span class=\"input-group-text\">Taxomony Code</span>\n                </div>\n                <input id=\"rs-stringdb-tax\" type=\"text\" class=\"form-control text-input-icon\" value=\"9606\"/>\n              </div>\n            </div>\n          </div>\n          <div class=\"col-sm-6\" style=\"padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px\">\n            <pre aria-live=\"assertive\" class=\"shiny-text-output noplaceholder\" id=\"rs-stringdb-error.msg\"></pre>\n          </div>\n          <div class=\"col-sm-2\" style=\"padding-left:15px; padding-right:2px; padding-top:0px; padding-bottom:0px\">\n            <button id=\"rs-stringdb-run\" type=\"button\" class=\"btn btn-default action-button\" title=\"Query the STRING database for interactions among the currently selected features (max ~300; requires internet)\" data-testid=\"rs-stringdb-run-analysis-button\"><span class=\"action-label\">Run!</span></button>\n          </div>\n        </div>\n        <div id=\"rs-stringdb-noresRet\" class=\"shiny-html-output\"></div>\n        <h3 class=\"sr-only\" aria-label=\"Table of protein-protein interactions from STRING database with confidence scores and evidence types\">Interaction Results Table</h3>\n        <div id=\"rs-stringdb-strtab-showButton\" class=\"shiny-html-output\"></div>\n        <div class=\"datatables html-widget html-widget-output shiny-report-size html-fill-item\" id=\"rs-stringdb-strtab-table\" style=\"width:100%;height:auto;\"></div>\n        <h3 class=\"sr-only\" aria-label=\"Interactive force-directed network graph showing protein interactions with node labels toggle\">Network Visualization</h3>\n        <div class=\"form-group shiny-input-container\" data-testid=\"rs-stringdb-show-labels-toggle\">\n          <div class=\"checkbox\">\n            <label>\n              <input id=\"rs-stringdb-showLabel\" type=\"checkbox\" class=\"shiny-input-checkbox\"/>\n              <span>Show labels</span>\n            </label>\n          </div>\n        </div>\n        <div class=\"forceNetwork html-widget html-widget-output shiny-report-size html-fill-item\" id=\"rs-stringdb-network\" style=\"width:100%;height:500px;\"></div>\n        <div class=\"sr-only\" aria-live=\"polite\" aria-atomic=\"true\">\n          <div id=\"rs-stringdb-networkSummary\" class=\"shiny-html-output\"></div>\n        </div>\n      </div>\n      <div class=\"tab-pane\" data-value=\"SeqLogo\" id=\"tab-TSID-4\">\n        <h2 class=\"sr-only\" aria-label=\"PTM motif enrichment analysis identifying overrepresented amino acid sequence patterns around modification sites\">Post-Translational Modification Motifs</h2>\n        <div class=\"sr-only\" id=\"rs-ptm-module-help\">\n          <h4>About Post-Translational Modification (PTM) Motif Analysis</h4>\n          <p>PTM motif analysis identifies enriched sequence patterns surrounding post-translational modification sites (e.g., phosphorylation, acetylation, ubiquitination) in your selected proteins. This analysis uses sequence logo visualization to show the consensus amino acid sequence motif and performs statistical testing to determine if certain motifs are over-represented in your dataset compared to background.</p>\n          <h4>When to use PTM motif analysis</h4>\n          <p>Use this analysis when you have identified PTM sites (typically from mass spectrometry) and want to understand which kinases, enzymes, or regulatory mechanisms might be responsible for those modifications. This helps identify active signaling pathways and regulatory networks. It's particularly valuable for phosphoproteomics data to predict kinase-substrate relationships.</p>\n          <h4>How to interpret results</h4>\n          <p>The sequence logo shows the consensus motif with letter height indicating amino acid frequency at each position relative to the PTM site (position 0). Larger letters indicate more conserved positions. The enrichment statistics show whether your selected PTM sites have a different motif pattern compared to all PTM sites in the dataset. Significant enrichment suggests specific kinases or regulatory proteins are active in your experimental condition.</p>\n        </div>\n        <div class=\"row\">\n          <div class=\"col-sm-2\" align=\"right\" style=\"padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px\">\n            <div id=\"rs-ptm-tris_seqlogo-groupLabel\" class=\"shiny-html-output\"></div>\n          </div>\n          <div class=\"col-sm-3\" style=\"padding:2px;\">\n            <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-ptm-tris_seqlogo-analysis-selector\">\n              <label class=\"control-label shiny-label-null\" for=\"rs-ptm-tris_seqlogo-analysis\" id=\"rs-ptm-tris_seqlogo-analysis-label\"></label>\n              <div>\n                <select id=\"rs-ptm-tris_seqlogo-analysis\" class=\"shiny-input-select\"></select>\n                <script type=\"application/json\" data-for=\"rs-ptm-tris_seqlogo-analysis\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n              </div>\n            </div>\n          </div>\n          <div class=\"col-sm-4\" style=\"padding:2px;\">\n            <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-ptm-tris_seqlogo-subset-selector\">\n              <label class=\"control-label shiny-label-null\" for=\"rs-ptm-tris_seqlogo-subset\" id=\"rs-ptm-tris_seqlogo-subset-label\"></label>\n              <div>\n                <select id=\"rs-ptm-tris_seqlogo-subset\" class=\"shiny-input-select\"></select>\n                <script type=\"application/json\" data-for=\"rs-ptm-tris_seqlogo-subset\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n              </div>\n            </div>\n          </div>\n          <div class=\"col-sm-3\" style=\"padding-left:2px; padding-right:20px; padding-top:2px; padding-bottom:2px\">\n            <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-ptm-tris_seqlogo-variable-selector\">\n              <label class=\"control-label shiny-label-null\" for=\"rs-ptm-tris_seqlogo-variable\" id=\"rs-ptm-tris_seqlogo-variable-label\"></label>\n              <div>\n                <select id=\"rs-ptm-tris_seqlogo-variable\" class=\"shiny-input-select\"></select>\n                <script type=\"application/json\" data-for=\"rs-ptm-tris_seqlogo-variable\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n              </div>\n            </div>\n          </div>\n        </div>\n        <div id=\"rs-ptm-msg_ui\" class=\"shiny-html-output\"></div>\n      </div>\n      <div class=\"tab-pane\" data-value=\"Geneshot\" id=\"tab-TSID-5\">\n        <h2 class=\"sr-only\" aria-label=\"Literature-based gene discovery using PubMed co-mention analysis to find genes associated with search terms\">Geneshot Literature Search</h2>\n        <div class=\"sr-only\" id=\"rs-geneshotTab-module-help\">\n          <h4>About Geneshot Literature Search</h4>\n          <p>Geneshot is a search engine that finds genes and proteins associated with arbitrary text queries by mining the scientific literature. It ranks genes by how frequently they co-occur with your search terms in PubMed abstracts and other biomedical texts.</p>\n          <h4>When to use Geneshot</h4>\n          <p>Use Geneshot when you want to discover genes related to a disease, phenotype, biological process, or any concept of interest. It's particularly useful for literature-based hypothesis generation, finding candidate genes for follow-up studies, or understanding which genes are most studied in relation to your topic of interest.</p>\n          <h4>How to interpret results</h4>\n          <p>Genes are ranked by their association score with your search terms. The scatter plot shows publication count (x-axis) versus percentage of publications mentioning your search terms (y-axis). Genes in the upper right are both well-studied and highly associated with your query. The plus symbol (+) in the 'selected' column indicates genes that overlap with your currently selected features in the app.</p>\n        </div>\n        <div class=\"omicsviewer-toolbar\" style=\"display:flex; flex-wrap:wrap; align-items:center;justify-content:flex-end; gap:8px; padding:2px 0 6px 0;\">\n          <div class=\"form-group shiny-input-container\" data-testid=\"rs-geneshotTab-search-input\" style=\"width:100%; flex:1 1 360px; margin-bottom:0;\">\n            <label class=\"control-label shiny-label-null\" for=\"rs-geneshotTab-term\" id=\"rs-geneshotTab-term-label\"></label>\n            <input id=\"rs-geneshotTab-term\" type=\"text\" class=\"shiny-input-text form-control\" value=\"\" placeholder=\"Search any term here, multiple separated by (;) ...\" data-update-on=\"change\"/>\n          </div>\n          <button id=\"rs-geneshotTab-submit\" type=\"button\" class=\"btn btn-default action-button\" title=\"Mine the literature for genes associated with the search term(s); separate multiple terms with &#39;;&#39;\" data-testid=\"rs-geneshotTab-search-button\"><span class=\"action-label\">Search related genes!</span></button>\n          <div class=\"sw-dropdown\" id=\"sw-drop-rs-geneshotTab-mapIdPanel\" data-testid=\"rs-geneshotTab-map-id-button\">\n            <button class=\"btn btn-default  action-button\" data-bs-html=\"false\" data-bs-placement=\"right\" data-bs-title=\"Choose the feature annotation column that Geneshot genes are matched against!\" data-bs-toggle=\"tooltip\" id=\"rs-geneshotTab-mapIdPanel\" type=\"button\">\n              <i class=\"far fa-id-card\" role=\"presentation\" aria-label=\"id-card icon\"></i>\n              Map ID\n              <span class=\"glyphicon glyphicon-triangle-bottom\"></span>\n            </button>\n            <div class=\"sw-dropdown-content animated sw-dropright-content\" id=\"sw-content-rs-geneshotTab-mapIdPanel\" style=\"width:620px;\">\n              <div class=\"sw-dropdown-in\" circle=\"FALSE\" margin=\"10px\">\n                <div class=\"row\">\n                  <div class=\"col-sm-2\" align=\"right\" style=\"padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px\">\n                    <div id=\"rs-geneshotTab-geneNameCol-groupLabel\" class=\"shiny-html-output\"></div>\n                  </div>\n                  <div class=\"col-sm-3\" style=\"padding:2px;\">\n                    <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-geneshotTab-geneNameCol-analysis-selector\">\n                      <label class=\"control-label shiny-label-null\" for=\"rs-geneshotTab-geneNameCol-analysis\" id=\"rs-geneshotTab-geneNameCol-analysis-label\"></label>\n                      <div>\n                        <select id=\"rs-geneshotTab-geneNameCol-analysis\" class=\"shiny-input-select\"></select>\n                        <script type=\"application/json\" data-for=\"rs-geneshotTab-geneNameCol-analysis\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                      </div>\n                    </div>\n                  </div>\n                  <div class=\"col-sm-4\" style=\"padding:2px;\">\n                    <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-geneshotTab-geneNameCol-subset-selector\">\n                      <label class=\"control-label shiny-label-null\" for=\"rs-geneshotTab-geneNameCol-subset\" id=\"rs-geneshotTab-geneNameCol-subset-label\"></label>\n                      <div>\n                        <select id=\"rs-geneshotTab-geneNameCol-subset\" class=\"shiny-input-select\"></select>\n                        <script type=\"application/json\" data-for=\"rs-geneshotTab-geneNameCol-subset\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                      </div>\n                    </div>\n                  </div>\n                  <div class=\"col-sm-3\" style=\"padding-left:2px; padding-right:20px; padding-top:2px; padding-bottom:2px\">\n                    <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-geneshotTab-geneNameCol-variable-selector\">\n                      <label class=\"control-label shiny-label-null\" for=\"rs-geneshotTab-geneNameCol-variable\" id=\"rs-geneshotTab-geneNameCol-variable-label\"></label>\n                      <div>\n                        <select id=\"rs-geneshotTab-geneNameCol-variable\" class=\"shiny-input-select\"></select>\n                        <script type=\"application/json\" data-for=\"rs-geneshotTab-geneNameCol-variable\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                      </div>\n                    </div>\n                  </div>\n                </div>\n              </div>\n            </div>\n            <script>$('#rs-geneshotTab-mapIdPanel').tooltip({ placement: 'right', title: 'Choose the feature annotation column that Geneshot genes are matched against!', html: false });</script>\n            <script>$(function() {swDrop('rs-geneshotTab-mapIdPanel', 'sw-content-rs-geneshotTab-mapIdPanel', 'sw-drop-rs-geneshotTab-mapIdPanel', 'sw-none', 'sw-none', '1');});</script>\n          </div>\n        </div>\n        <div class=\"plotly html-widget html-widget-output shiny-report-size shiny-report-theme html-fill-item\" id=\"rs-geneshotTab-plt\" style=\"width:100%;height:400px;\"></div>\n        <br/>\n        <div id=\"rs-geneshotTab-autorif-showButton\" class=\"shiny-html-output\"></div>\n        <div class=\"datatables html-widget html-widget-output shiny-report-size html-fill-item\" id=\"rs-geneshotTab-autorif-table\" style=\"width:100%;height:auto;\"></div>\n        <div class=\"sr-only\" aria-live=\"polite\" aria-atomic=\"true\">\n          <div id=\"rs-geneshotTab-searchSummary\" class=\"shiny-html-output\"></div>\n        </div>\n      </div>\n      <div class=\"tab-pane\" data-value=\"Sample\" id=\"tab-TSID-6\">\n        <h2 class=\"sr-only\" aria-label=\"Statistical analysis of selected samples including group comparisons, contingency tables, and Kaplan-Meier survival curves\">Sample Analysis</h2>\n        <div class=\"sr-only\" id=\"rs-sample_general-module-help\">\n          <h4>About Sample Group Analysis</h4>\n          <p>This module provides comprehensive statistical analysis to compare and visualize sample groups or patient cohorts. It offers multiple analysis types including beeswarm plots for continuous variables, contingency table analysis for categorical variables, and survival analysis for time-to-event data. This is the sample-centric complement to the feature analysis module.</p>\n          <h4>When to use sample analysis</h4>\n          <p>Use this module when you want to analyze relationships between sample metadata variables, compare patient groups, test associations between clinical variables, or perform survival analysis stratified by sample characteristics. This is essential for understanding sample clustering, identifying clinical correlations, and validating patient groupings.</p>\n          <h4>How to interpret results</h4>\n          <p>The visualization type depends on your selected variables: beeswarm plots show distributions of continuous variables across groups with statistical tests, contingency tables show relationships between categorical variables with chi-square or Fisher's exact test, and survival curves show time-to-event outcomes stratified by groups with log-rank test. The results table provides detailed statistical metrics appropriate for each analysis type.</p>\n        </div>\n        <div class=\"row\">\n          <div class=\"col-sm-12\">\n            <div class=\"box\" style=\"height: 490px\">\n              <div class=\"box-body\">\n                <div class=\"sr-only\" id=\"rs-sample_general-batch_comp-module-help\">\n                  <h4>About Batch Comparison Analysis</h4>\n                  <p>Batch comparison analysis automatically tests all phenotype variables or all molecular features for differences between your selected and unselected sample groups. This comprehensive approach applies appropriate statistical tests based on variable types (numeric, categorical, or survival data) and corrects for multiple testing using FDR. It provides a systematic way to identify which clinical variables or molecular features significantly differ between your sample groups.</p>\n                  <h4>When to use batch comparison</h4>\n                  <p>Use this analysis when you have selected a group of samples (e.g., based on clustering, clinical criteria, or molecular signature) and want to systematically identify: (1) which clinical/phenotype variables characterize your selected group, or (2) which genes/proteins are differentially expressed in your selected group. This is essential for biomarker discovery, patient stratification validation, and comprehensive group characterization.</p>\n                  <h4>How to interpret results</h4>\n                  <p>The results table shows one row per variable/feature tested. P-values indicate statistical significance of group differences. FDR-corrected p-values (q-values) control for false discoveries when testing many variables - use FDR &lt; 0.05 as a stringent threshold. For phenotype comparisons, summary statistics show distributions in each group. For feature comparisons, mean differences indicate effect sizes. Sort by FDR to identify the most significant differences.</p>\n                </div>\n                <div class=\"row\">\n                  <div class=\"col-sm-6\">\n                    <div class=\"form-group shiny-input-container\">\n                      <div class=\"checkbox\">\n                        <label>\n                          <input id=\"rs-sample_general-batch_comp-show_phenotype\" type=\"checkbox\" class=\"shiny-input-checkbox\" checked=\"checked\"/>\n                          <span>Compare All Phenotype Variables</span>\n                        </label>\n                      </div>\n                    </div>\n                  </div>\n                  <div class=\"col-sm-6\">\n                    <div class=\"form-group shiny-input-container\">\n                      <div class=\"checkbox\">\n                        <label>\n                          <input id=\"rs-sample_general-batch_comp-show_features\" type=\"checkbox\" class=\"shiny-input-checkbox\"/>\n                          <span>Compare All Features/Expression</span>\n                        </label>\n                      </div>\n                    </div>\n                  </div>\n                  <div class=\"tabbable\">\n                    <ul class=\"nav nav-tabs\" data-tabsetid=\"TSID\">\n                      <li class=\"active\">\n                        <a href=\"#tab-TSID-1\" data-toggle=\"tab\" data-bs-toggle=\"tab\" data-value=\"Phenotype variable\">Phenotype variable</a>\n                      </li>\n                      <li>\n                        <a href=\"#tab-TSID-2\" data-toggle=\"tab\" data-bs-toggle=\"tab\" data-value=\"Feature differential expression\">Feature differential expression</a>\n                      </li>\n                    </ul>\n                    <div class=\"tab-content\" data-tabsetid=\"TSID\">\n                      <div class=\"tab-pane active\" data-value=\"Phenotype variable\" id=\"tab-TSID-1\">\n                        <div id=\"rs-sample_general-batch_comp-phenotype_results\" class=\"shiny-html-output\"></div>\n                      </div>\n                      <div class=\"tab-pane\" data-value=\"Feature differential expression\" id=\"tab-TSID-2\">\n                        <div id=\"rs-sample_general-batch_comp-features_results\" class=\"shiny-html-output\"></div>\n                      </div>\n                    </div>\n                  </div>\n                </div>\n              </div>\n            </div>\n          </div>\n          <div class=\"col-sm-12\" style=\"margin-top: 0px;\">\n            <div class=\"row\">\n              <div class=\"col-sm-2\" align=\"right\" style=\"padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px\">\n                <div id=\"rs-sample_general-tris_sample_general-groupLabel\" class=\"shiny-html-output\"></div>\n              </div>\n              <div class=\"col-sm-3\" style=\"padding:2px;\">\n                <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-sample_general-tris_sample_general-analysis-selector\">\n                  <label class=\"control-label shiny-label-null\" for=\"rs-sample_general-tris_sample_general-analysis\" id=\"rs-sample_general-tris_sample_general-analysis-label\"></label>\n                  <div>\n                    <select id=\"rs-sample_general-tris_sample_general-analysis\" class=\"shiny-input-select\"></select>\n                    <script type=\"application/json\" data-for=\"rs-sample_general-tris_sample_general-analysis\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                  </div>\n                </div>\n              </div>\n              <div class=\"col-sm-4\" style=\"padding:2px;\">\n                <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-sample_general-tris_sample_general-subset-selector\">\n                  <label class=\"control-label shiny-label-null\" for=\"rs-sample_general-tris_sample_general-subset\" id=\"rs-sample_general-tris_sample_general-subset-label\"></label>\n                  <div>\n                    <select id=\"rs-sample_general-tris_sample_general-subset\" class=\"shiny-input-select\"></select>\n                    <script type=\"application/json\" data-for=\"rs-sample_general-tris_sample_general-subset\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                  </div>\n                </div>\n              </div>\n              <div class=\"col-sm-3\" style=\"padding-left:2px; padding-right:5px; padding-top:2px; padding-bottom:2px\">\n                <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-sample_general-tris_sample_general-variable-selector\">\n                  <label class=\"control-label shiny-label-null\" for=\"rs-sample_general-tris_sample_general-variable\" id=\"rs-sample_general-tris_sample_general-variable-label\"></label>\n                  <div>\n                    <select id=\"rs-sample_general-tris_sample_general-variable\" class=\"shiny-input-select\"></select>\n                    <script type=\"application/json\" data-for=\"rs-sample_general-tris_sample_general-variable\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                  </div>\n                </div>\n              </div>\n            </div>\n          </div>\n          <div class=\"col-sm-5\">\n            <div class=\"box\" style=\"height: 500px\">\n              <div class=\"box-body\">\n                <div id=\"rs-sample_general-msatab-showButton\" class=\"shiny-html-output\"></div>\n                <div class=\"datatables html-widget html-widget-output shiny-report-size html-fill-item\" id=\"rs-sample_general-msatab-table\" style=\"width:100%;height:auto;\"></div>\n              </div>\n            </div>\n          </div>\n          <div class=\"col-sm-7\">\n            <div class=\"box\" style=\"height: 500px\">\n              <div class=\"box-body\">\n                <div class=\"col-sm-11\">\n                  <div id=\"rs-sample_general-sample_general_plot\" class=\"shiny-html-output\"></div>\n                </div>\n                <div class=\"col-sm-1\">\n                  <div class=\"sw-dropdown\" id=\"sw-drop-btn-RID\">\n                    <button class=\"btn btn-default  action-button\" data-bs-html=\"false\" data-bs-placement=\"right\" data-bs-title=\"Click to modify figure!\" data-bs-toggle=\"tooltip\" id=\"btn-RID\" type=\"button\">\n                      <i class=\"fas fa-gear\" role=\"presentation\" aria-label=\"gear icon\"></i>\n                      <span class=\"glyphicon glyphicon-triangle-bottom\"></span>\n                    </button>\n                    <div class=\"sw-dropdown-content animated sw-dropright-content\" id=\"sw-content-btn-RID\" style=\"width:788px;\">\n                      <div class=\"sw-dropdown-in\" margin=\"25px\" circle=\"FALSE\">\n                        <br/>\n                        <div class=\"row\">\n                          <div class=\"col-sm-2\" align=\"right\" style=\"padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px\">\n                            <div id=\"rs-sample_general-a4_gp-selectColorUI-groupLabel\" class=\"shiny-html-output\"></div>\n                          </div>\n                          <div class=\"col-sm-3\" style=\"padding:2px;\">\n                            <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-sample_general-a4_gp-selectColorUI-analysis-selector\">\n                              <label class=\"control-label shiny-label-null\" for=\"rs-sample_general-a4_gp-selectColorUI-analysis\" id=\"rs-sample_general-a4_gp-selectColorUI-analysis-label\"></label>\n                              <div>\n                                <select id=\"rs-sample_general-a4_gp-selectColorUI-analysis\" class=\"shiny-input-select\"></select>\n                                <script type=\"application/json\" data-for=\"rs-sample_general-a4_gp-selectColorUI-analysis\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                              </div>\n                            </div>\n                          </div>\n                          <div class=\"col-sm-4\" style=\"padding:2px;\">\n                            <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-sample_general-a4_gp-selectColorUI-subset-selector\">\n                              <label class=\"control-label shiny-label-null\" for=\"rs-sample_general-a4_gp-selectColorUI-subset\" id=\"rs-sample_general-a4_gp-selectColorUI-subset-label\"></label>\n                              <div>\n                                <select id=\"rs-sample_general-a4_gp-selectColorUI-subset\" class=\"shiny-input-select\"></select>\n                                <script type=\"application/json\" data-for=\"rs-sample_general-a4_gp-selectColorUI-subset\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                              </div>\n                            </div>\n                          </div>\n                          <div class=\"col-sm-3\" style=\"padding-left:2px; padding-right:20px; padding-top:2px; padding-bottom:2px\">\n                            <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-sample_general-a4_gp-selectColorUI-variable-selector\">\n                              <label class=\"control-label shiny-label-null\" for=\"rs-sample_general-a4_gp-selectColorUI-variable\" id=\"rs-sample_general-a4_gp-selectColorUI-variable-label\"></label>\n                              <div>\n                                <select id=\"rs-sample_general-a4_gp-selectColorUI-variable\" class=\"shiny-input-select\"></select>\n                                <script type=\"application/json\" data-for=\"rs-sample_general-a4_gp-selectColorUI-variable\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                              </div>\n                            </div>\n                          </div>\n                        </div>\n                        <div class=\"row\">\n                          <div class=\"col-sm-2\" align=\"right\" style=\"padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px\">\n                            <div id=\"rs-sample_general-a4_gp-selectShapeUI-groupLabel\" class=\"shiny-html-output\"></div>\n                          </div>\n                          <div class=\"col-sm-3\" style=\"padding:2px;\">\n                            <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-sample_general-a4_gp-selectShapeUI-analysis-selector\">\n                              <label class=\"control-label shiny-label-null\" for=\"rs-sample_general-a4_gp-selectShapeUI-analysis\" id=\"rs-sample_general-a4_gp-selectShapeUI-analysis-label\"></label>\n                              <div>\n                                <select id=\"rs-sample_general-a4_gp-selectShapeUI-analysis\" class=\"shiny-input-select\"></select>\n                                <script type=\"application/json\" data-for=\"rs-sample_general-a4_gp-selectShapeUI-analysis\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                              </div>\n                            </div>\n                          </div>\n                          <div class=\"col-sm-4\" style=\"padding:2px;\">\n                            <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-sample_general-a4_gp-selectShapeUI-subset-selector\">\n                              <label class=\"control-label shiny-label-null\" for=\"rs-sample_general-a4_gp-selectShapeUI-subset\" id=\"rs-sample_general-a4_gp-selectShapeUI-subset-label\"></label>\n                              <div>\n                                <select id=\"rs-sample_general-a4_gp-selectShapeUI-subset\" class=\"shiny-input-select\"></select>\n                                <script type=\"application/json\" data-for=\"rs-sample_general-a4_gp-selectShapeUI-subset\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                              </div>\n                            </div>\n                          </div>\n                          <div class=\"col-sm-3\" style=\"padding-left:2px; padding-right:20px; padding-top:2px; padding-bottom:2px\">\n                            <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-sample_general-a4_gp-selectShapeUI-variable-selector\">\n                              <label class=\"control-label shiny-label-null\" for=\"rs-sample_general-a4_gp-selectShapeUI-variable\" id=\"rs-sample_general-a4_gp-selectShapeUI-variable-label\"></label>\n                              <div>\n                                <select id=\"rs-sample_general-a4_gp-selectShapeUI-variable\" class=\"shiny-input-select\"></select>\n                                <script type=\"application/json\" data-for=\"rs-sample_general-a4_gp-selectShapeUI-variable\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                              </div>\n                            </div>\n                          </div>\n                        </div>\n                        <div class=\"row\">\n                          <div class=\"col-sm-2\" align=\"right\" style=\"padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px\">\n                            <div id=\"rs-sample_general-a4_gp-selectSizeUI-groupLabel\" class=\"shiny-html-output\"></div>\n                          </div>\n                          <div class=\"col-sm-3\" style=\"padding:2px;\">\n                            <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-sample_general-a4_gp-selectSizeUI-analysis-selector\">\n                              <label class=\"control-label shiny-label-null\" for=\"rs-sample_general-a4_gp-selectSizeUI-analysis\" id=\"rs-sample_general-a4_gp-selectSizeUI-analysis-label\"></label>\n                              <div>\n                                <select id=\"rs-sample_general-a4_gp-selectSizeUI-analysis\" class=\"shiny-input-select\"></select>\n                                <script type=\"application/json\" data-for=\"rs-sample_general-a4_gp-selectSizeUI-analysis\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                              </div>\n                            </div>\n                          </div>\n                          <div class=\"col-sm-4\" style=\"padding:2px;\">\n                            <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-sample_general-a4_gp-selectSizeUI-subset-selector\">\n                              <label class=\"control-label shiny-label-null\" for=\"rs-sample_general-a4_gp-selectSizeUI-subset\" id=\"rs-sample_general-a4_gp-selectSizeUI-subset-label\"></label>\n                              <div>\n                                <select id=\"rs-sample_general-a4_gp-selectSizeUI-subset\" class=\"shiny-input-select\"></select>\n                                <script type=\"application/json\" data-for=\"rs-sample_general-a4_gp-selectSizeUI-subset\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                              </div>\n                            </div>\n                          </div>\n                          <div class=\"col-sm-3\" style=\"padding-left:2px; padding-right:20px; padding-top:2px; padding-bottom:2px\">\n                            <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-sample_general-a4_gp-selectSizeUI-variable-selector\">\n                              <label class=\"control-label shiny-label-null\" for=\"rs-sample_general-a4_gp-selectSizeUI-variable\" id=\"rs-sample_general-a4_gp-selectSizeUI-variable-label\"></label>\n                              <div>\n                                <select id=\"rs-sample_general-a4_gp-selectSizeUI-variable\" class=\"shiny-input-select\"></select>\n                                <script type=\"application/json\" data-for=\"rs-sample_general-a4_gp-selectSizeUI-variable\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                              </div>\n                            </div>\n                          </div>\n                        </div>\n                        <div class=\"row\">\n                          <div class=\"col-sm-2\" align=\"right\" style=\"padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px\">\n                            <div id=\"rs-sample_general-a4_gp-selectTooltipUI-groupLabel\" class=\"shiny-html-output\"></div>\n                          </div>\n                          <div class=\"col-sm-3\" style=\"padding:2px;\">\n                            <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-sample_general-a4_gp-selectTooltipUI-analysis-selector\">\n                              <label class=\"control-label shiny-label-null\" for=\"rs-sample_general-a4_gp-selectTooltipUI-analysis\" id=\"rs-sample_general-a4_gp-selectTooltipUI-analysis-label\"></label>\n                              <div>\n                                <select id=\"rs-sample_general-a4_gp-selectTooltipUI-analysis\" class=\"shiny-input-select\"></select>\n                                <script type=\"application/json\" data-for=\"rs-sample_general-a4_gp-selectTooltipUI-analysis\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                              </div>\n                            </div>\n                          </div>\n                          <div class=\"col-sm-4\" style=\"padding:2px;\">\n                            <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-sample_general-a4_gp-selectTooltipUI-subset-selector\">\n                              <label class=\"control-label shiny-label-null\" for=\"rs-sample_general-a4_gp-selectTooltipUI-subset\" id=\"rs-sample_general-a4_gp-selectTooltipUI-subset-label\"></label>\n                              <div>\n                                <select id=\"rs-sample_general-a4_gp-selectTooltipUI-subset\" class=\"shiny-input-select\"></select>\n                                <script type=\"application/json\" data-for=\"rs-sample_general-a4_gp-selectTooltipUI-subset\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                              </div>\n                            </div>\n                          </div>\n                          <div class=\"col-sm-3\" style=\"padding-left:2px; padding-right:20px; padding-top:2px; padding-bottom:2px\">\n                            <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-sample_general-a4_gp-selectTooltipUI-variable-selector\">\n                              <label class=\"control-label shiny-label-null\" for=\"rs-sample_general-a4_gp-selectTooltipUI-variable\" id=\"rs-sample_general-a4_gp-selectTooltipUI-variable-label\"></label>\n                              <div>\n                                <select id=\"rs-sample_general-a4_gp-selectTooltipUI-variable\" class=\"shiny-input-select\"></select>\n                                <script type=\"application/json\" data-for=\"rs-sample_general-a4_gp-selectTooltipUI-variable\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                              </div>\n                            </div>\n                          </div>\n                        </div>\n                        <div class=\"row\">\n                          <div class=\"col-sm-2\" align=\"right\" style=\"padding-left:2px; padding-right:2px; padding-top:0px; padding-bottom:0px\">\n                            <div id=\"rs-sample_general-a4_gp-selectSearchCol-groupLabel\" class=\"shiny-html-output\"></div>\n                          </div>\n                          <div class=\"col-sm-3\" style=\"padding:2px;\">\n                            <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-sample_general-a4_gp-selectSearchCol-analysis-selector\">\n                              <label class=\"control-label shiny-label-null\" for=\"rs-sample_general-a4_gp-selectSearchCol-analysis\" id=\"rs-sample_general-a4_gp-selectSearchCol-analysis-label\"></label>\n                              <div>\n                                <select id=\"rs-sample_general-a4_gp-selectSearchCol-analysis\" class=\"shiny-input-select\"></select>\n                                <script type=\"application/json\" data-for=\"rs-sample_general-a4_gp-selectSearchCol-analysis\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                              </div>\n                            </div>\n                          </div>\n                          <div class=\"col-sm-4\" style=\"padding:2px;\">\n                            <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-sample_general-a4_gp-selectSearchCol-subset-selector\">\n                              <label class=\"control-label shiny-label-null\" for=\"rs-sample_general-a4_gp-selectSearchCol-subset\" id=\"rs-sample_general-a4_gp-selectSearchCol-subset-label\"></label>\n                              <div>\n                                <select id=\"rs-sample_general-a4_gp-selectSearchCol-subset\" class=\"shiny-input-select\"></select>\n                                <script type=\"application/json\" data-for=\"rs-sample_general-a4_gp-selectSearchCol-subset\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                              </div>\n                            </div>\n                          </div>\n                          <div class=\"col-sm-3\" style=\"padding-left:2px; padding-right:20px; padding-top:2px; padding-bottom:2px\">\n                            <div class=\"form-group shiny-input-container\" style=\"width:100%;\" data-testid=\"rs-sample_general-a4_gp-selectSearchCol-variable-selector\">\n                              <label class=\"control-label shiny-label-null\" for=\"rs-sample_general-a4_gp-selectSearchCol-variable\" id=\"rs-sample_general-a4_gp-selectSearchCol-variable-label\"></label>\n                              <div>\n                                <select id=\"rs-sample_general-a4_gp-selectSearchCol-variable\" class=\"shiny-input-select\"></select>\n                                <script type=\"application/json\" data-for=\"rs-sample_general-a4_gp-selectSearchCol-variable\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                              </div>\n                            </div>\n                          </div>\n                        </div>\n                        <div class=\"shiny-panel-conditional\" data-display-if=\"1 == 2\" data-ns-prefix=\"\">\n                          <div class=\"form-group shiny-input-container\">\n                            <div class=\"checkbox\">\n                              <label>\n                                <input id=\"rs-sample_general-a4_gp-showSearchBox\" type=\"checkbox\" class=\"shiny-input-checkbox\"/>\n                                <span>show</span>\n                              </label>\n                            </div>\n                          </div>\n                        </div>\n                        <div class=\"shiny-panel-conditional\" data-display-if=\"input.showSearchBox == true\" data-ns-prefix=\"rs-sample_general-a4_gp-\">\n                          <div style=\"padding-left:100px; padding-right:0px; padding-top:0px; padding-bottom:0px\">\n                            <div class=\"form-group shiny-input-container\" style=\"width:100%;\">\n                              <label class=\"control-label shiny-label-null\" for=\"rs-sample_general-a4_gp-searchon\" id=\"rs-sample_general-a4_gp-searchon-label\"></label>\n                              <div>\n                                <select id=\"rs-sample_general-a4_gp-searchon\" class=\"shiny-input-select\" multiple=\"multiple\"></select>\n                                <script type=\"application/json\" data-for=\"rs-sample_general-a4_gp-searchon\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                              </div>\n                            </div>\n                          </div>\n                        </div>\n                        <div class=\"row\">\n                          <div class=\"col-sm-2\"></div>\n                          <div class=\"col-sm-4\" style=\"padding:2px;\">\n                            <div class=\"form-group shiny-input-container\">\n                              <label class=\"control-label\" id=\"rs-sample_general-a4_gp-xcut-label\" for=\"rs-sample_general-a4_gp-xcut\">Select points by x/y cutoffs</label>\n                              <div class=\"input-group\">\n                                <div class=\"input-group-addon sw-input-icon input-group-prepend\">\n                                  <span class=\"input-group-text\">x-cut</span>\n                                </div>\n                                <input id=\"rs-sample_general-a4_gp-xcut\" type=\"text\" class=\"form-control text-input-icon\" value=\"log10(2)\" placeholder=\"e.g. -1 or -log10(2)\"/>\n                              </div>\n                            </div>\n                          </div>\n                          <div class=\"col-sm-4\" style=\"padding-left:5px; padding-right:2px; padding-top:27px; padding-bottom:2px;\">\n                            <div class=\"form-group shiny-input-container\">\n                              <label class=\"control-label shiny-label-null\" for=\"rs-sample_general-a4_gp-ycut\" id=\"rs-sample_general-a4_gp-ycut-label\"></label>\n                              <div class=\"input-group\">\n                                <div class=\"input-group-addon sw-input-icon input-group-prepend\">\n                                  <span class=\"input-group-text\">y-cut</span>\n                                </div>\n                                <input id=\"rs-sample_general-a4_gp-ycut\" type=\"text\" class=\"form-control text-input-icon\" value=\"-log10(0.05)\" placeholder=\"e.g 2 or -log10(0.05)\"/>\n                              </div>\n                            </div>\n                          </div>\n                          <div class=\"col-sm-2\" style=\"padding-left:5px; padding-right:20px; padding-top:4px; padding-bottom:2px;\">\n                            <div class=\"form-group shiny-input-container\">\n                              <label class=\"control-label\" id=\"rs-sample_general-a4_gp-scorner-label\" for=\"rs-sample_general-a4_gp-scorner\">Area</label>\n                              <div>\n                                <select id=\"rs-sample_general-a4_gp-scorner\" class=\"shiny-input-select\"><option value=\"None\" selected>None</option></select>\n                                <script type=\"application/json\" data-for=\"rs-sample_general-a4_gp-scorner\" data-nonempty=\"\">{\"plugins\":[\"selectize-plugin-a11y\"]}</script>\n                              </div>\n                            </div>\n                          </div>\n                        </div>\n                      </div>\n                    </div>\n                    <script>$('#btn-RID').tooltip({ placement: 'right', title: 'Click to modify figure!', html: false });</script>\n                    <script>$(function() {swDrop('btn-RID', 'sw-content-btn-RID', 'sw-drop-btn-RID', 'sw-none', 'sw-none', '1');});</script>\n                  </div>\n                </div>\n              </div>\n            </div>\n          </div>\n        </div>\n      </div>\n    </div>\n  </div>\n</body>", html = TRUE, class = c("html", 
"character"))


## runtime helpers (byte-mirrors of the capture script's logic)
.wp12_norm <- function(h) {
  # WP-12 review NIT-2 (anchored patterns, tightened in the WP-13 round):
  # every rule is anchored to the exact syntactic position the
  # nondeterministic id occupies, so a content difference confined to
  # digits inside a "tab-<n>-" / "btn-<n>" substring of PANEL CONTENT
  # (tooltip, label, prose) can no longer be masked. tabset ids appear as
  # pane ids and navbar hrefs (tab-<seed>-<n>; the trailing pane index is
  # deterministic and stays); the shinyWidgets btn ids ride attribute
  # values AND the same ids echoed inside the widget <script> (jQuery
  # selector + swDrop args), so the rule anchors on the surrounding
  # quote and covers the sw-drop-/sw-content- prefixed forms - the frozen
  # goldens normalize identically under the anchored set (byte-checked
  # by the md5/full-literal assertions below).
  h <- gsub('data-tabsetid="[0-9]+"', 'data-tabsetid="TSID"', h)
  h <- gsub('id="tab-[0-9]+-', 'id="tab-TSID-', h)
  h <- gsub('href="#tab-[0-9]+-', 'href="#tab-TSID-', h)
  gsub("(['\"])(#?)(sw-drop-|sw-content-|)btn-[0-9]+(['\"])",
       "\\1\\2\\3btn-RID\\4", h, perl = TRUE)
}
.wp12_tabs <- function(h) {
  ul <- regexpr('shiny-tab-input', h)
  end <- ul + regexpr("</ul>", substring(h, ul)) - 2L
  nav <- substring(h, ul, end)
  lis <- regmatches(nav, gregexpr("(?s)<li[^>]*>.*?</li>", nav, perl = TRUE))[[1]]
  pane_of <- function(v) {
    p <- regexpr(paste0('<div class="tab-pane[^"]*" data-value="', v,
                        '" id="tab-TSID-[0-9]+">\\s*<h2 class="sr-only" aria-label="([^"]*)">([^<]*)</h2>'),
                h, perl = TRUE)
    if (p == -1) return(list(sr_title = NA_character_, aria = NA_character_))
    capt <- attr(p, "capture.start"); lens <- attr(p, "capture.length")
    list(sr_title = paste0(substring(h, capt[2], capt[2] + lens[2] - 1L)),
         aria = paste0(substring(h, capt[1], capt[1] + lens[1] - 1L)))
  }
  lapply(lis, function(li) {
    v <- sub('.*data-value="([^"]*)".*', "\\1", li)
    t <- sub('.*<span title="([^"]*)".*', "\\1", li)
    ic <- sub('.*fa-([a-z0-9-]+).*', "\\1", li)
    span <- sub('(?s).*<span title="[^"]*">(.*?)</span>.*', "\\1", li, perl = TRUE)
    lab <- trimws(gsub("</?[^>]+>", "", span))
    h2 <- pane_of(v)
    list(value = v, title = t, icon = ic, label = lab,
         sr_title = h2$sr_title, aria = h2$aria)
  })
}
.wp12_md5 <- function(h) {
  f <- tempfile(); writeLines(h, f, useBytes = TRUE)
  on.exit(unlink(f), add = TRUE)
  unname(tools::md5sum(f))
}
# the two additional-tab stubs (same shape as the capture's)
.wp12_addl <- list(
  list(tabName = "Custom1", moduleName = "customOne",
       moduleUi = function(id) tags$div(class = "wp12-addl", paste0("content-", id)),
       moduleServer = function(...) invisible(NULL), stateful = TRUE),
  list(tabName = "Custom2", moduleName = "customTwo",
       moduleUi = function(id) tags$div(class = "wp12-addl", paste0("content-", id)),
       moduleServer = function(...) invisible(NULL), stateful = FALSE))
.wp12_mkapp <- function(fd, addl = NULL, store = NULL)
  function(input, output, session) {
    omicsViewer:::L1_result_space_module(
      "rs",
      reactive_expr = reactive(Biobase::exprs(dat)),
      reactive_phenoData = reactive(Biobase::pData(dat)),
      reactive_featureData = reactive(fd),
      additionalTabs = addl,
      status = reactive(NULL), store = store)
    # the frozen legacy fixture through its own renderUI output in the
    # SAME session - the identical pipeline on both sides
    output$wp12_legacy <- renderUI(
      .wp12_legacy_tabset(ns = NS("rs"),
                          reactive_featureData = reactive(fd),
                          additionalTabs = addl))
  }
.wp12_boot <- function(fd, addl = NULL, store = NULL) {
  e <- new.env()
  shiny::testServer(.wp12_mkapp(fd, addl, store), {
    for (i in 1:3) session$flushReact()
    e$boot <- session$getOutput("rs-optTabs")$html
    e$fix <- session$getOutput("wp12_legacy")$html
  })
  list(boot = .wp12_norm(e$boot), fix = .wp12_norm(e$fix))
}

## ------- the parity proof: goldens + same-commit old-vs-new diff ------
.wp12_variants <- list(
  demo = list(fd = fd0, addl = NULL), gs = list(fd = fd0_gs, addl = NULL),
  noopt = list(fd = fd0_no_opt, addl = NULL),
  nullfd = list(fd = NULL, addl = NULL),
  demo_addl = list(fd = fd0, addl = .wp12_addl),
  gs_addl = list(fd = fd0_gs, addl = .wp12_addl))
.wp12_run <- list()
for (.nm in names(.wp12_variants)) {
  .v <- .wp12_variants[[.nm]]
  .wp12_run[[.nm]] <- .wp12_boot(.v$fd, .v$addl)
}
.wp12_tabs_run <- lapply(.wp12_run, function(r) .wp12_tabs(r$boot))

for (.nm in names(.wp12_variants)) {
  ok(.wp12_md5(.wp12_run[[.nm]]$boot) == WP12_MD5[[.nm]],
     paste0("WP-12 DOM parity (", .nm, "): registry composition md5 equals the frozen pre-change golden"))
  ok(identical(.wp12_run[[.nm]]$boot, .wp12_run[[.nm]]$fix),
     paste0("WP-12 DOM parity (", .nm, "): same-commit old-vs-new normalized-HTML diff is empty"))
  ok(ut_cmp_identical(.wp12_tabs_run[[.nm]], get(paste0("WP12_TABS_", toupper(.nm)))),
     paste0("WP-12 DOM parity (", .nm, "): parsed truth table equals the frozen golden"))
  ok(ut_cmp_identical(
       vapply(.wp12_tabs_run[[.nm]], function(t) t$value, ""),
       rs_tab_set(.wp12_variants[[.nm]]$fd,
                  if (is.null(.wp12_variants[[.nm]]$addl)) character(0) else
                    vapply(.wp12_variants[[.nm]]$addl, function(lo) lo$tabName, ""))),
     paste0("WP-12 DOM parity (", .nm, "): the ACTUAL tabset equals rs_tab_set (availability single source of truth)"))
}
ok(identical(.wp12_run$demo$boot, WP12_DEMO_HTML),
   "WP-12 DOM parity (demo): full normalized HTML byte-identical to the frozen literal")
ok(identical(.wp12_md5(.wp12_run$nullfd$boot), .wp12_md5(.wp12_run$noopt$boot)),
   "WP-12 DOM parity: NULL fd and no-optional-column fd render the same tabset (degradation identical)")
.wp12_store_boot <- .wp12_boot(fd0, NULL, store = widget_store_new())
ok(.wp12_md5(.wp12_store_boot$boot) == WP12_MD5[["demo"]],
   "WP-12 DOM parity: the tabset DOM is store-independent (store-backed boot equals the store-less golden)")
ok(identical(.wp12_store_boot$boot, .wp12_store_boot$fix),
   "WP-12 DOM parity: store-backed old-vs-new diff also empty")

## pathological additional tab duplicating a core tab name (PD-11 case):
## the DOM keeps BOTH (the legacy loop appends unconditionally); no golden
## was captured for this shape pre-change, so this is the same-run
## old-vs-new diff alone
.wp12_dup <- list(list(tabName = "Feature", moduleName = "dupFeat",
                       moduleUi = function(id) tags$div("dup"), stateful = FALSE,
                       moduleServer = function(...) invisible(NULL)))
.wp12_dup_run <- .wp12_boot(fd0, .wp12_dup)
ok(identical(.wp12_dup_run$boot, .wp12_dup_run$fix),
   "WP-12 DOM parity: additional tab duplicating a core name composes identically (legacy loop verbatim)")
ok(ut_cmp_identical(vapply(.wp12_tabs(.wp12_dup_run$boot), function(t) t$value, ""),
                    c("Feature", "Response", "StringDB", "SeqLogo", "Feature", "Geneshot", "Sample")),
   "WP-12 DOM parity: duplicated core name renders BOTH tab panels (PD-11's uniquing is choices-only, not DOM)")

## ------------------------------ unit checks ---------------------------
ok(all(vapply(panel_records, function(r)
       is.character(r$tab_title) && is.character(r$sr_title) &&
       is.character(r$aria_label), TRUE)),
   "WP-12 registry: every seed record carries the three presentation fields")
ok(identical(panel_records$sample_general$sr_title, "Sample Analysis") &&
   identical(panel_records$ptm$icon, "font") &&
   startsWith(panel_records$ora$tab_title, "Over-representation analysis (hypergeometric test)."),
   "WP-12 registry: presentation fields spot-check (tooltips/h2 verbatim from the dissolved tabPanels)")
ok(grepl("tab_title", .err(panel_record(
     "a", "x.y", "T", "L", "i", 1, function(fd, pd) TRUE, identity, identity,
     tab_title = ""))),
   "WP-12 panel_record: empty tab_title rejected with the field named")
ok(grepl("sr_title", .err(panel_record(
     "a", "x.y", "T", "L", "i", 1, function(fd, pd) TRUE, identity, identity,
     sr_title = NA))),
   "WP-12 panel_record: NA sr_title rejected")
ok(grepl("aria_label", .err(panel_record(
     "a", "x.y", "T", "L", "i", 1, function(fd, pd) TRUE, identity, identity,
     aria_label = 1))),
   "WP-12 panel_record: non-character aria_label rejected")
r_new <- panel_record(
  "a", "x.y", "T", "L", "i", 2, function(fd, pd) TRUE, identity, identity,
  stateful = FALSE, tab_title = "tip", sr_title = "sr", aria_label = "ar")
ok(identical(r_new$tab_title, "tip") && identical(r_new$sr_title, "sr") &&
   identical(r_new$aria_label, "ar"),
   "WP-12 panel_record: valid presentation fields returned verbatim (all optional for synthetic records)")
# direct structural check: the generic builder over the ora record is the
# dissolved ORA tabPanel (span+title / value / sr-only h2+aria / ui(ns(id)))
.tp <- omicsViewer:::`.rs_tab_panel`(panel_records$ora, NS("rs"))
ok(identical(.tp$name, "div") && identical(.tp$attribs$`data-value`, "ORA") &&
   inherits(.tp$attribs$title, "shiny.tag") &&
   startsWith(.tp$attribs$title$attribs$title, "Over-representation analysis") &&
   identical(.tp$attribs$title$children[[2]], "ORA") &&
   any(vapply(.tp$children, function(ch)
     inherits(ch, "shiny.tag") && identical(ch$name, "h2") &&
     identical(ch$attribs$`aria-label`, panel_records$ora$aria_label), TRUE)),
   "WP-12 tab builder: record fields drive the span/value/h2 of the tabPanel")
# direct composition call (no session): tab order follows the order field
# and availability follows rs_tab_set - the PD-64 server/tab order split
ok(ut_cmp_identical(
     vapply(.wp12_tabs(.wp12_norm(htmltools::renderTags(
       omicsViewer:::`.rs_compose_opt_tabs`(NS("rs"), fd0_gs, NULL))$html)),
       function(t) t$value, ""),
     c("Feature", "ORA", "fGSEA", "Response", "StringDB", "SeqLogo", "Geneshot", "Sample")),
   "WP-12 composition: direct call (session-less) yields the order-field tab order with ORA/fGSEA at 2-3")


## ------------- O5 fix round: per-record status slices in the L1 loops ----
# O5 (PD-79) and its pre-existing additionalTabs sibling. The WP-11
# registry loop built every panel's status slice as
# reactive(status()[[r$status_key]]) - lazily capturing the LOOP
# VARIABLE r, so all eight slices resolved the LAST record's key
# (analyst_dose_response) and the legacy result-space widget re-apply
# through the REAL L1 (the modules' observeEvent(reactive_status())
# restore observers) was a silent no-op since f7fd24e. The
# additionalTabs server loop carried the same trap around
# lo$moduleName (with >= 2 stateful plugins, BOTH received the last
# plugin's analyst_additional payload; sibling present since v1.17.2).
# Both construction sites now capture the loop variable eagerly inside
# local(); the loop bodies are otherwise untouched (the context builds
# evaluate eagerly and were never affected).
#
# MUTATION HONESTY: these cases FAIL on the unfixed code. They drive
# the restore OBSERVERS through the real L1 composition - a status
# reactiveVal delivered exactly like a snapshot restore - and assert
# the MODULE-SIDE consumption (store_apply origin="restore" landing
# under each panel's own store keys), never the store path directly.
# On unfixed L1 the per-panel delivery lands NOTHING (every slice
# reads analyst_dose_response, absent from the payload), the dose-key
# misdelivery CROSS-APPLIES into ora/geneshot, and both plugins
# receive plugB's payload - the WP-13 reviewer's probe shapes.

# --- per-panel slice resolution through the real L1 (three panels:
# ora + geneshot + feature_general; demo.RDS, real store) ---
o5 <- new.env()
o5$store <- widget_store_new()
o5$status <- shiny::reactiveVal(NULL)
o5$app <- function(input, output, session) {
  omicsViewer:::L1_result_space_module(
    "rs",
    reactive_expr = reactive(Biobase::exprs(dat)),
    reactive_phenoData = reactive(Biobase::pData(dat)),
    reactive_featureData = reactive(fd0),
    status = o5$status,
    store = o5$store)
}
o5_err <- NULL; o5_err2 <- NULL; o5_closed <- NULL
v5 <- NULL; v5b <- NULL
shiny::testServer(o5$app, {
  o5_err <<- tryCatch({
    for (i in 1:3) session$flushReact()
    # a legacy-v1-shaped restore payload: each panel's slice under its
    # OWN analyst_* key, delivered at once (the snapshot restore shape)
    o5$status(list(
      analyst_ora = list(xax = list(analysis = "General", subset = "All",
                                    variable = "Gene.name")),
      analyst_gene_shot = list(term = "cell cycle"),
      analyst_feature_general = list(
        xax = list(analysis = "General", subset = "All",
                   variable = "Cell.line"),
        plotType = "Curve")))
    for (i in 1:5) session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  v5 <<- store_snapshot(o5$store)$values
  # misdelivery probe: the SAME payload fields under the LAST record's
  # key (analyst_dose_response - the one key every unfixed slice
  # resolved), with ora/geneshot-valid values DIFFERENT from the landed
  # ones, so on unfixed L1 they flip (the reviewer's probe (b2) shape)
  o5_err2 <<- tryCatch({
    o5$status(list(analyst_dose_response = list(
      xax = list(analysis = "General", subset = "All",
                 variable = "Protein.ID"),
      term = "apoptosis")))
    for (i in 1:5) session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  v5b <<- store_snapshot(o5$store)$values
  o5_closed <<- session$isClosed()
})
ok(is.null(o5_err) && ut_cmp_identical(
     paste(v5[c("resultspace.ora.xax_analysis", "resultspace.ora.xax_subset",
                "resultspace.ora.xax_variable")], collapse = "|"),
     "General|All|Gene.name"),
   "O5: an analyst_ora payload lands its xax triple under the ora store keys through the real L1")
ok(is.null(o5_err) && ut_cmp_identical(v5[["resultspace.geneshot.term"]], "cell cycle"),
   "O5: an analyst_gene_shot payload lands its term under the geneshot store key through the real L1")
ok(is.null(o5_err) && ut_cmp_identical(
     paste(v5[c("resultspace.feature_general.xax_analysis",
                "resultspace.feature_general.xax_subset",
                "resultspace.feature_general.xax_variable")], collapse = "|"),
     "General|All|Cell.line") &&
     ut_cmp_identical(v5[["resultspace.feature_general.plot_type"]], "Curve"),
   "O5: an analyst_feature_general payload lands its xax triple + plot type under the feature-general store keys through the real L1")
ok(is.null(o5_err2) &&
     ut_cmp_identical(v5b[["resultspace.ora.xax_variable"]], "Gene.name") &&
     ut_cmp_identical(v5b[["resultspace.geneshot.term"]], "cell cycle") &&
     ut_cmp_identical(v5b, v5),
   "O5: an analyst_dose_response payload does NOT cross-apply into ora/geneshot/feature (no module consumes the dose slice; the whole store is value-identical after the misdelivery)")
ok(is.null(o5_err) && is.null(o5_err2) && identical(o5_closed, FALSE),
   "O5: the L1 composition session stays alive through both legacy deliveries")

# --- the additionalTabs sibling: two stateful plugins with DISTINCT
# payloads -> each module receives ITS OWN analyst_additional slice
# (the exact pre-fix failure: both resolved the LAST plugin's payload)
.o5p_env <- new.env()
.o5p_env$store <- widget_store_new()
.o5p_env$status <- shiny::reactiveVal(NULL)
.o5p_env$rec <- new.env(); .o5p_env$rec$A <- NULL; .o5p_env$rec$B <- NULL
.o5_plugin <- function(slot) function(id, reactive_status, ...) {
  observe({
    v <- reactive_status()
    if (!is.null(v)) .o5p_env$rec[[slot]] <<- v
  })
  NULL
}
o5p_err <- NULL; o5p_closed <- NULL
shiny::testServer(function(input, output, session) {
  omicsViewer:::L1_result_space_module(
    "rs",
    reactive_expr = reactive(Biobase::exprs(dat)),
    reactive_phenoData = reactive(Biobase::pData(dat)),
    reactive_featureData = reactive(fd0),
    status = .o5p_env$status,
    store = .o5p_env$store,
    additionalTabs = list(
      list(moduleName = "plugA", moduleServer = .o5_plugin("A"),
           stateful = TRUE, tabName = "A", icon = NULL, label = "A"),
      list(moduleName = "plugB", moduleServer = .o5_plugin("B"),
           stateful = TRUE, tabName = "B", icon = NULL, label = "B")))
}, {
  o5p_err <<- tryCatch({
    for (i in 1:3) session$flushReact()
    .o5p_env$status(list(analyst_additional = list(plugA = list(x = "AAA"),
                                                   plugB = list(x = "BBB"))))
    for (i in 1:5) session$flushReact()
    NULL
  }, error = function(e) conditionMessage(e))
  o5p_closed <<- session$isClosed()
})
ok(is.null(o5p_err) && ut_cmp_identical(.o5p_env$rec$A, list(x = "AAA")) &&
     ut_cmp_identical(.o5p_env$rec$B, list(x = "BBB")),
   "O5 sibling: two stateful plugins each receive their OWN analyst_additional slice (not both-get-last)")
ok(is.null(o5p_err) && identical(o5p_closed, FALSE),
   "O5 sibling: the additionalTabs composition session stays alive through the delivery")
