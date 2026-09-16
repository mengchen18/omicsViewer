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
  ut_cmp_identical(state$schema_version, 1L),
  "migrated snapshot has schema version 1"
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
  ut_cmp_warning(
    validate_app_state(state2, dataset = es_changed, dataset_id = "changed.RDS"),
    c("fingerprint mismatch", "saved for dataset"),
    expected_count = 2L
  ),
  "dataset incompatibility warns instead of blocking"
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
