# Generator for the golden snapshot fixtures (todo 2.6/2.9).
#
# Run manually when the snapshot SCHEMA changes:
#   Rscript tests/fixtures/make_snapshot_fixtures.R
#
# The three files pin one representative of every readable generation:
#   snapshot_v0_legacy.ESS  pre-versioning flat list (eset_*/analyst_* names)
#   snapshot_v1.ESS         schema 1 (panel-carried axes/attr4/selectByCorner)
#   snapshot_v2.ESS         schema 2 (widget_store + selection records)
#
# Content is deterministic (fixed created_at/package_version, demo.RDS ids)
# and every fixture is restored through app_module in test_snapshotRoundTrip.

suppressMessages({
  library(Biobase)
  library(omicsViewer)
})
dir <- file.path("tests", "fixtures")
if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
demo_path <- system.file("extdata", "demo.RDS", package = "omicsViewer")
if (!nzchar(demo_path)) demo_path <- file.path("inst", "extdata", "demo.RDS")
dat <- readRDS(demo_path)
fd <- fData(dat)
pd <- pData(dat)
feat <- head(rownames(fd), 12)
samp <- head(colnames(pd), 6)
stamp <- "2020-01-01T00:00:00+0000"
pver <- "2.1.57"

ds_state <- omicsViewer:::dataset_state(dat, id = "ESVObj.RDS")

## ---- v0: legacy flat list --------------------------------------------------
v0 <- list(
  label = "golden legacy",
  created_at = stamp,
  package_version = pver,
  eset_active_tab = "Feature table",
  eset_selected_features = feat,
  eset_selected_samples = samp,
  analyst_active_tab = "ORA",
  active_feature = feat,
  active_sample = samp
)
saveRDS(v0, file.path(dir, "snapshot_v0_legacy.ESS"), compress = "xz")

## ---- v1: schema 1 (panel-carried axes/attr4, selectByCorner) ----------------
v1 <- omicsViewer:::new_app_state(
  dataset_id = "ESVObj.RDS",
  dataset = dat,
  selection = list(features = feat, samples = samp),
  app = list(data_active_tab = "Feature", analysis_active_tab = "Feature"),
  data_space = list(
    eset_active_tab = "Feature",
    eset_fdata_fig = list(
      axisMode = "custom",
      xax = list("ttest", "RE_vs_ME", "mean.diff"),
      yax = list("ttest", "RE_vs_ME", "log.fdr"),
      showRegLine = TRUE,
      selectByCorner = FALSE,
      selection_clicked = character(0),
      selection_selected = feat[1:3],
      attr4 = list(
        selectColor = list("General", "All", "Cell.line"),
        selectShape = NULL,
        xcut = "log10(2)",
        ycut = "-log10(0.05)",
        acorner = "None"
      )
    )
  ),
  label = "golden v1",
  package_version = pver,
  schema_version = 1L,
  created_at = stamp
)
saveRDS(v1, file.path(dir, "snapshot_v1.ESS"), compress = "xz")

## ---- v2: schema 2 (widget_store + selection records) -----------------------
v2 <- omicsViewer:::new_app_state(
  dataset_id = "ESVObj.RDS",
  dataset = dat,
  selection = list(
    features = feat[4:8],
    samples = samp[1:2],
    records = list(
      feature = list(ids = feat[4:8], clicked = character(0),
                     origin = "table", anchor = NULL,
                     mirror = feat[4:8], epoch = 0L),
      sample = list(ids = samp[1:2], clicked = character(0),
                    origin = "figure", anchor = NULL,
                    mirror = TRUE, epoch = 0L)
    )
  ),
  app = list(data_active_tab = "Sample", analysis_active_tab = "Feature"),
  data_space = list(
    eset_active_tab = "Sample",
    eset_pdata_fig = list(selectByCorner = FALSE)
  ),
  widget_store = list(
    values = list(
      `dataspace.feature_space.x_analysis` = "Cor",
      `dataspace.feature_space.x_subset` = "MDR",
      `dataspace.feature_space.x_variable` = "R",
      `dataspace.feature_space.y_analysis` = "Cor",
      `dataspace.feature_space.y_subset` = "MDR",
      `dataspace.feature_space.y_variable` = "logP",
      `dataspace.feature_space.axis_mode` = "quick",
      `dataspace.expr_heatmap.heatmap_colors` = "RdGy",
      `dataspace.feature_space.attr4.color_analysis` = NULL
    ),
    unset = "dataspace.feature_space.attr4.color_analysis"
  ),
  label = "golden v2",
  package_version = pver,
  schema_version = 2L,
  created_at = stamp
)
saveRDS(v2, file.path(dir, "snapshot_v2.ESS"), compress = "xz")

cat("fixtures written to", normalizePath(dir), ":\n")
print(list.files(dir, pattern = "\\.ESS$"))
