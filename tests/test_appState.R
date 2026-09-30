if (!requireNamespace("unittest", quietly = TRUE)) {
  message("test_appState needs unittest (Suggests-only)")
  quit(save = "no", status = 0)
}
library(omicsViewer)
library(Biobase)
library(unittest, quietly = TRUE)

new_app_state <- omicsViewer:::new_app_state
build_app_state <- omicsViewer:::build_app_state
migrate_app_state <- omicsViewer:::migrate_app_state
validate_app_state <- omicsViewer:::validate_app_state
normalize_selection <- omicsViewer:::normalize_selection
sanitize_snapshot_name <- omicsViewer:::sanitize_snapshot_name
snapshot_file_name <- omicsViewer:::snapshot_file_name
dataset_fingerprint <- omicsViewer:::dataset_fingerprint
write_app_state <- omicsViewer:::write_app_state
data_table_widget_state <- omicsViewer:::data_table_widget_state
restore_table_page_length <- omicsViewer:::restore_table_page_length

ok(
  ut_cmp_identical(
    sanitize_snapshot_name("../evil/name"),
    "evil_name"
  ),
  "unsafe snapshot name characters are replaced"
)
ok(
  ut_cmp_identical(sanitize_snapshot_name("..", fallback = "safe"), "safe"),
  "parent-directory-only names are rejected"
)
ok(
  ut_cmp_identical(sanitize_snapshot_name(c(" a/b ", "second")), "a_b"),
  "only the first name is used"
)
ok(
  ut_cmp_identical(
    snapshot_file_name("../x", dataset_id = "../demo.RDS"),
    "ESVSnapshot_demo.RDS_x.ESS"
  ),
  "snapshot file name cannot contain a path"
)

expr <- matrix(1:12, nrow = 3, dimnames = list(c("f1", "f2", "f3"), c("s1", "s2", "s3", "s4")))
pd <- data.frame(group = c("a", "a", "b", "b"), row.names = colnames(expr))
fd <- data.frame(score = 1:3, row.names = rownames(expr))
es <- ExpressionSet(
  assayData = expr,
  phenoData = AnnotatedDataFrame(pd),
  featureData = AnnotatedDataFrame(fd)
)

ok(
  ut_cmp_identical(
    dataset_fingerprint(es, id = "demo.RDS"),
    dataset_fingerprint(es, id = "demo.RDS")
  ),
  "dataset fingerprint is deterministic"
)
es_changed <- es
featureNames(es_changed) <- c("f1", "f2", "changed")
ok(
  ut_cmp_identical(
    dataset_fingerprint(es, id = "demo.RDS") == dataset_fingerprint(es_changed, id = "demo.RDS"),
    FALSE
  ),
  "dataset fingerprint detects changed feature identity"
)

legacy <- list(
  eset_active_tab = "Sample",
  eset_selected_samples = "s1",
  analyst_active_tab = "ORA",
  active_feature = c("f1", NA, "f2"),
  active_sample = "s1"
)
state <- migrate_app_state(legacy)
ok(
  ut_cmp_identical(state$format, "omicsViewerState"),
  "legacy snapshots are migrated to the versioned format"
)
ok(
  ut_cmp_identical(state$schema_version, omicsViewer:::APP_STATE_SCHEMA_VERSION),
  "migrated snapshot has the current schema version"
)
ok(
  ut_cmp_identical(state$panels$data_space$eset_active_tab, "Sample"),
  "legacy data-space state is namespaced"
)
ok(
  ut_cmp_identical(state$panels$result_space$analyst_active_tab, "ORA"),
  "legacy result-space state is namespaced"
)
ok(
  ut_cmp_identical(state$selection$features, c("f1", "f2")),
  "legacy semantic feature selection is normalized"
)

table_state <- list(
  start = 20L,
  length = 10L,
  order = list(c(1L, "asc")),
  columns = list(
    list(search = list(search = "", regex = FALSE, smart = TRUE)),
    list(visible = FALSE, search = list(search = "abc"))
  ),
  time = 123456,
  childRows = list()
)
compact_table_state <- data_table_widget_state(table_state)
ok(
  ut_cmp_identical(
    compact_table_state,
    list(
      start = 20L,
      length = 10L,
      order = list(c(1L, "asc")),
      columns = list(
        list(search = list(search = "", regex = FALSE, smart = TRUE)),
        list(search = list(search = "abc"))
      )
    )
  ),
  "DataTable state keeps only portable page, order, and filter widgets"
)
ok(
  ut_cmp_identical(restore_table_page_length(NULL, 25L), 25L),
  "missing DataTable page length falls back to the module default"
)
ok(
  ut_cmp_identical(restore_table_page_length(50L, 25L), 50L),
  "saved DataTable page length is restored"
)

state2 <- build_app_state(
  dataset = es,
  dataset_id = "demo.RDS",
  data_status = list(
    eset_active_tab = "Feature table",
    eset_fdata_tab = compact_table_state,
    eset_fdata_fig = list(
      xax = list("General", "All", "score"),
      htestV1 = list(computed = TRUE),
      rowDendrogram = matrix(1)
    )
  ),
  result_status = list(
    analyst_active_tab = "Feature",
    analyst_feature_general = list(
      plotType = "Bees",
      htestV1 = list(computed = TRUE),
      htestV2 = list(computed = TRUE)
    ),
    analyst_gene_shot = list(term = "cell cycle", rif = list(large = TRUE)),
    analyst_fgsea = list(xax = list("ttest", "KO_vs_WT", "mean.diff")),
    analyst_stringdb = list(tax = "9606", showLabel = TRUE)
  ),
  selected_features = c("f2", "f2"),
  selected_samples = c("s1", "s3"),
  label = "widget-only state"
)
ok(
  ut_cmp_identical(state2$dataset$id, "demo.RDS"),
  "new state records the dataset identifier"
)
ok(
  ut_cmp_identical(state2$app$data_active_tab, "Feature table"),
  "new state records the active data tab"
)
ok(
  ut_cmp_identical(state2$app$analysis_active_tab, "Feature"),
  "new state records the active analysis tab"
)
ok(
  ut_cmp_identical(state2$selection$features, "f2"),
  "new state stores unique semantic feature IDs"
)
ok(
  ut_cmp_identical(state2$selection$samples, c("s1", "s3")),
  "new state stores semantic sample IDs"
)

ok(
  ut_cmp_identical(state2$policy, "widget-only"),
  "snapshot records the widget-only state policy"
)
ok(
  ut_cmp_identical(
    state2$gaps$stringdb,
    "STRING network widget/result state is not captured in this phase."
  ),
  "STRING state is explicitly recorded as a gap"
)
ok(
  all(
    vapply(state2$panels, function(panel) {
      txt <- paste(utils::capture.output(str(panel)), collapse = " ")
      !grepl("htestV1|htestV2|rowDendrogram|analyst_stringdb|large = TRUE", txt)
    }, logical(1))
  ),
  "computed analysis payloads and STRING state are omitted"
)
ok(
  ut_cmp_identical(state2$panels$data_space$eset_fdata_tab$start, 20L),
  "table page position is retained in widget state"
)
ok(
  ut_cmp_identical(state2$panels$result_space$analyst_gene_shot$term, "cell cycle"),
  "text widget values are retained in widget state"
)

validated <- validate_app_state(state2, dataset = es, dataset_id = "demo.RDS")
ok(
  ut_cmp_identical(validated$selection$samples, c("s1", "s3")),
  "compatible state validates without changing selections"
)
ok(
  ut_cmp_identical(length(attr(validated, "warnings")), 0L),
  "compatible state collects no warnings"
)
ok(
  ut_cmp_identical(
    attr(validate_app_state(state2, dataset = es_changed, dataset_id = "changed.RDS"), "warnings"),
    c(
      "This snapshot was saved from a different version of the dataset (features, samples or annotation columns changed). Content that no longer exists will be ignored.",
      "Snapshot was saved for dataset 'demo.RDS' and is being restored into 'changed.RDS'."
    )
  ),
  "dataset incompatibility collects ALL warnings instead of blocking"
)
# unknown ids are intersected with the current dataset, not dropped silently
es_sm_changed <- es
sampleNames(es_sm_changed) <- c("s1", "s2", "b1", "b2")
v2 <- validate_app_state(state2, dataset = es_sm_changed, dataset_id = "demo.RDS")
ok(
  ut_cmp_identical(v2$selection$samples, "s1"),
  "unknown selection ids are intersected away after being reported"
)
ok(
  any(grepl("1 of 2 selected samples", attr(v2, "warnings"))),
  "unknown selection ids are reported with counts"
)
res <- tryCatch(migrate_app_state(list(foo = 1, bar = 2)), error = function(e) e)
ok(
  inherits(res, "error") && grepl("not an omicsViewer snapshot", res$message),
  "a stray list .ESS is rejected as not a snapshot"
)

ok(
  ut_cmp_identical(
    normalize_selection(list(features = c(NA, "", "x"), samples = NULL)),
    list(features = "x", samples = character())
  ),
  "selection normalization drops missing and empty IDs"
)

tmpdir <- tempfile()
dir.create(tmpdir)
path <- file.path(tmpdir, "snapshot.ESS")
write_app_state(state2, path)
ok(file.exists(path), "atomic state writer creates the final file")
ok(
  length(list.files(tmpdir, pattern = "\\.tmp$")) == 0,
  "atomic state writer leaves no temporary files"
)
loaded <- readRDS(path)
ok(
  ut_cmp_identical(loaded$format, "omicsViewerState"),
  "written state can be loaded"
)

## ------------- S4 completion: widget-store snapshot round-trip ---------
# The .ESS snapshot embeds the canonical widget store; restoring it must
# reproduce the saved store state exactly (per-key resilient only for
# values the live data invalidates). Driven through the real app_module:
# widgets patched via the agent path, saved through the snapshot modal
# observer, drifted, then restored through the savedSS table selection.
.rt_dir <- file.path(tempdir(), paste0("esv-rt-", Sys.getpid()))
if (dir.exists(.rt_dir)) unlink(.rt_dir, recursive = TRUE)
dir.create(.rt_dir, recursive = TRUE)
Sys.setenv(OMICSVIEWER_TEST_HOOKS = "true")
.rt_demo_path <- system.file("extdata", "demo.RDS", package = "omicsViewer")
if (!nzchar(.rt_demo_path))
  .rt_demo_path <- file.path("inst", "extdata", "demo.RDS")
.rt_dat <- readRDS(.rt_demo_path)
.rt_store_vals <- function() {
  omicsViewer:::store_snapshot(.rt_app_store)$values
}
.rt_run <- new.env(); .rt_run$n <- 0L
.rt_hook <- function(session, output, op, payload = list()) {
  session$setInputs(`app-agentTestHooks-op` = op)
  session$setInputs(`app-agentTestHooks-payload` = jsonlite::toJSON(
    payload, auto_unbox = TRUE))
  .rt_run$n <- .rt_run$n + 1L
  session$setInputs(`app-agentTestHooks-run` = .rt_run$n)
  session$flushReact()
  txt <- as.character(output[["app-agentTestHooks-result"]])
  jsonlite::fromJSON(txt, simplifyVector = FALSE)
}
.rt_app_store <- omicsViewer:::widget_store_new()
app_rt <- function(input, output, session) {
  omicsViewer:::app_module(
    "app", .dir = shiny::reactive(.rt_dir),
    ESVObj = shiny::reactive(.rt_dat), store = .rt_app_store)
}
shiny::testServer(app_rt, {
  # the navbar input must exist before the first flush, otherwise the
  # data-space status assembly errors under the mock session
  session$setInputs(`app-dataspace-eset` = "Feature")
  session$flushReact()
  # agent-path widget writes across modules (scatter, heatmap, table,
  # result space); the data-space navbar is driven through its input like
  # a real browser (a push relays updateNavbarPage, the mock input map
  # would otherwise keep the stale tab and the saved status disagree)
  r1 <- .rt_hook(session, output, "widgets", list(patch = list(
    `dataspace.expr_heatmap.heatmap_colors` = "RdGy",
    `dataspace.expr_heatmap.margin_bottom` = 7,
    `dataspace.tab_feature.multi_selection` = TRUE,
    `resultspace.feature_general.plot_type` = "Curve")))
  ok(length(r1$applied) == 4L && !length(r1$rejected),
     "round-trip setup: four widget writes apply")
  session$setInputs(`app-dataspace-eset` = "Heatmap")
  session$flushReact()
  s1 <- .rt_store_vals()
  ok(
    identical(s1$`dataspace.expr_heatmap.heatmap_colors`, "RdGy") &&
      identical(s1$`dataspace.expr_heatmap.margin_bottom`, 7L) &&
      isTRUE(s1$`dataspace.tab_feature.multi_selection`) &&
      identical(s1$`resultspace.feature_general.plot_type`, "Curve") &&
      identical(s1$`dataspace.active_tab`, "Heatmap"),
    "pre-save store state holds the requested widget values"
  )
  # save through the real snapshot observer (.ESS to disk)
  session$setInputs(`app-snapshot_name` = "rt1")
  session$setInputs(`app-snapshot_save` = 1L)
  session$flushReact()
  ess <- list.files(.rt_dir, pattern = "\\.ESS$", ignore.case = TRUE)
  ok(length(ess) == 1L, "snapshot save writes exactly one .ESS file")
  saved <- if (length(ess)) readRDS(file.path(.rt_dir, ess))
  ok(!is.null(saved) && is.list(saved$widget_store$values),
     "saved .ESS embeds the widget-store snapshot")
  # drift the widgets after saving
  r2 <- .rt_hook(session, output, "widgets", list(patch = list(
    `dataspace.expr_heatmap.heatmap_colors` = "PiYG",
    `dataspace.expr_heatmap.margin_bottom` = 3,
    `dataspace.tab_feature.multi_selection` = FALSE,
    `resultspace.feature_general.plot_type` = "Bees")))
  ok(length(r2$applied) == 4L,
     "post-save drift applies")
  smid <- .rt_store_vals()
  ok(!identical(smid$`dataspace.expr_heatmap.heatmap_colors`,
                s1$`dataspace.expr_heatmap.heatmap_colors`),
     "drift actually changed the store state")
  # restore by selecting the saved row in the snapshot table; the confirm
  # dialog (todo 2.7) must be accepted before the restore runs
  session$setInputs(`app-tab_saveSS_cells_selected` = c(1, 1))
  session$flushReact()
  smid <- .rt_store_vals()
  ok(
    !identical(smid$`dataspace.expr_heatmap.heatmap_colors`,
               s1$`dataspace.expr_heatmap.heatmap_colors`),
    "row click alone does not restore (confirm dialog gates it)"
  )
  session$setInputs(`app-snapshot_restore_confirm` = 1L)
  session$flushReact()
  session$flushReact()
  s2 <- .rt_store_vals()
  diffs <- Filter(function(k) !identical(s1[[k]], s2[[k]]), names(s1))
  ok(
    length(diffs) == 0L,
    sprintf("restore reproduces the saved store state exactly%s",
            if (length(diffs)) paste0(" (differs: ",
                                      paste(head(diffs, 5), collapse = ", "),
                                      ")") else "")
  )
  # WP11 wiring: opting in without a conversation saves no assistant field;
  # the history hooks expose the exact save/restore path headlessly.
  session$setInputs(`app-snapshot_include_chat` = TRUE)
  session$setInputs(`app-snapshot_name` = "rt-chat")
  session$setInputs(`app-snapshot_save` = 2L)
  session$flushReact()
  ess2 <- list.files(.rt_dir, pattern = "rt-chat\\.ESS$", ignore.case = TRUE)
  ok(length(ess2) == 1L, "opt-in snapshot save writes a second .ESS file")
  saved2 <- if (length(ess2)) readRDS(file.path(.rt_dir, ess2))
  ok(is.null(saved2$assistant),
     "opt-in save without a conversation stores no assistant payload")
  hs <- .rt_hook(session, output, "history_save", list())
  ok(is.null(hs) || identical(hs$hook_run > 0L, TRUE),
     "history save hook runs and reports no conversation without a provider")
})
Sys.unsetenv("OMICSVIEWER_TEST_HOOKS")
unlink(.rt_dir, recursive = TRUE)

## ------------- 1.5: save into an unwritable directory must not crash ---
# A read-only data dir is normal on shared servers. The save observer must
# surface the failure as a notification and keep the session alive instead
# of letting the error escape (ShinySession$unhandledError closes the app).
if (identical(Sys.getenv("USER", "chen"), "root")) {
  ok(TRUE, "snapshot failure probe skipped for root (chmod is not enforced)")
} else {
  .ro_dir <- file.path(tempdir(), paste0("esv-ro-", Sys.getpid()))
  if (dir.exists(.ro_dir)) unlink(.ro_dir, recursive = TRUE)
  dir.create(.ro_dir, recursive = TRUE)
  Sys.chmod(.ro_dir, "0555")
  app_ro <- function(input, output, session) {
    omicsViewer:::app_module(
      "app", .dir = shiny::reactive(.ro_dir),
      ESVObj = shiny::reactive(.rt_dat), store = omicsViewer:::widget_store_new())
  }
  shiny::testServer(app_ro, {
    session$setInputs(`app-dataspace-eset` = "Feature")
    session$flushReact()
    session$setInputs(`app-snapshot_name` = "ro1")
    session$setInputs(`app-snapshot_save` = 1L)
    session$flushReact()
    ok(
      length(list.files(.ro_dir, pattern = "\\.ESS$", ignore.case = TRUE)) == 0L,
      "failed snapshot save writes no .ESS file"
    )
    ok(
      ut_cmp_identical(isFALSE(session$isClosed()), TRUE),
      "session survives a failed snapshot save"
    )
    # a second attempt must still work (the modal stays usable)
    session$setInputs(`app-snapshot_name` = "ro2")
    session$setInputs(`app-snapshot_save` = 2L)
    session$flushReact()
    ok(
      ut_cmp_identical(isFALSE(session$isClosed()), TRUE),
      "session survives repeated failed snapshot saves"
    )
  })
  Sys.chmod(.ro_dir, "0755")
  unlink(.ro_dir, recursive = TRUE)
}

ok(
  ut_cmp_identical(
    nchar(sanitize_snapshot_name(paste(replicate(300, "a", simplify = "character"), collapse = "")),
                    type = "bytes") <= 80L,
    TRUE
  ),
  "snapshot names are capped so file names stay under file system limits"
)
