library(omicsViewer)
library(unittest, quietly = TRUE)

agent_environment_config <- omicsViewer:::agent_environment_config
agent_request_limit <- omicsViewer:::agent_request_limit
agent_validate_provider_config <- omicsViewer:::agent_validate_provider_config
agent_annotation_catalog <- omicsViewer:::agent_annotation_catalog
agent_search_annotations <- omicsViewer:::agent_search_annotations
agent_summarize_annotation <- omicsViewer:::agent_summarize_annotation
agent_normalize_state_update <- omicsViewer:::agent_normalize_state_update
agent_normalize_scatter_view <- omicsViewer:::agent_normalize_scatter_view

fd <- data.frame(
  `ttest|A_vs_B|log.fdr` = c(0.1, 0.4, 0.8),
  category = c("kinase", "phosphatase", "kinase"),
  check.names = FALSE,
  row.names = c("Gene1", "Gene2", "Gene3")
)
pd <- data.frame(
  group = c("WT", "WT", "KO", "KO"),
  batch = c(1, 2, 1, 2),
  row.names = c("S1", "S2", "S3", "S4")
)

config <- agent_validate_provider_config(
  provider = "openai", model = " test-model ", api_key = " secret ",
  base_url = " https://example.org/v1 "
)
ok(ut_cmp_identical(config$provider, "openai"), "provider config is normalized")
ok(ut_cmp_identical(config$model, "test-model"), "provider model is trimmed")
ok(ut_cmp_identical(config$base_url, "https://example.org/v1"), "provider base URL is trimmed")
ok(ut_cmp_identical(config$configured, TRUE), "provider config reports a key")
ok(
  ut_cmp_error(
    agent_validate_provider_config(provider = "openai", base_url = "ftp://example.org"),
    "API base URL must use HTTPS, or local HTTP for this machine only."
  ),
  "remote non-HTTPS provider base URLs are rejected"
)
ok(
  ut_cmp_identical(
    agent_validate_provider_config("openai", base_url = "http://127.0.0.1:8080/v1")$base_url,
    "http://127.0.0.1:8080/v1"
  ),
  "local HTTP provider endpoints are allowed"
)

old_request_limit <- Sys.getenv("OMICSVIEWER_LLM_MAX_REQUESTS")
Sys.setenv(OMICSVIEWER_LLM_MAX_REQUESTS = "1000")
on.exit(Sys.setenv(
  OMICSVIEWER_LLM_MAX_REQUESTS = if (identical(old_request_limit, "")) "" else old_request_limit
), add = TRUE)
ok(ut_cmp_identical(omicsViewer:::agent_request_limit(), 200L), "request limit is capped")
Sys.setenv(OMICSVIEWER_LLM_MAX_REQUESTS = "not-a-number")
ok(ut_cmp_identical(omicsViewer:::agent_request_limit(), 40L), "invalid request limit falls back to 40")
Sys.setenv(
  OMICSVIEWER_LLM_MAX_REQUESTS = if (identical(old_request_limit, "")) "" else old_request_limit
)

catalog <- agent_annotation_catalog(fd, pd)
ok(ut_cmp_identical(catalog$feature$rows, 3L), "feature catalog row count")
ok(
  ut_cmp_identical(catalog$feature$columns[["ttest|A_vs_B|log.fdr"]]$type, "numeric"),
  "feature catalog column type"
)
ok(
  ut_cmp_identical(catalog$sample$columns[["group"]]$unique_count, 2L),
  "sample catalog cardinality"
)

search <- agent_search_annotations(
  space = "feature", query = "gene1", feature_data = fd, sample_data = pd
)
ok(ut_cmp_identical(search$matching_ids, "Gene1"), "annotation search is case-insensitive")
ok(
  ut_cmp_identical(search$value_matches, list()),
  "annotation search does not invent value matches"
)
value_search <- agent_search_annotations(
  space = "sample", query = "ko", feature_data = fd, sample_data = pd
)
ok(
  ut_cmp_identical(value_search$value_matches[[1]]$matching_rows, 2L),
  "annotation value search counts matches"
)
ok(
  ut_cmp_identical(value_search$matching_id_count, 0L),
  "annotation value search does not confuse values with IDs"
)

summary <- agent_summarize_annotation(
  space = "feature", column = "ttest|A_vs_B|log.fdr", feature_data = fd, sample_data = pd
)
ok(ut_cmp_identical(summary$type, "numeric"), "numeric annotation summary type")
ok(ut_cmp_identical(summary$non_missing, 3L), "numeric annotation summary count")
ok(ut_cmp_identical(summary$summary[[1]], 0.1), "numeric annotation summary minimum")
ok(
  ut_cmp_error(
    agent_summarize_annotation("feature", "missing", fd, pd),
    "Unknown feature annotation column: missing"
  ),
  "unknown annotation columns are rejected"
)

update <- agent_normalize_state_update(
  update = list(
    data_space_tab = "Sample",
    features = c("Gene3", "Gene1", "Gene1"),
    samples = character()
  ),
  data_tabs = c("Feature", "Sample"),
  analysis_tabs = c("Feature", "ORA"),
  feature_ids = rownames(fd),
  sample_ids = rownames(pd)
)
ok(ut_cmp_identical(update$data_space_tab, "Sample"), "assistant tab update is validated")
ok(ut_cmp_identical(update$features, c("Gene3", "Gene1")), "assistant feature IDs are deduplicated")
ok(ut_cmp_identical(update$samples, character()), "assistant can explicitly clear samples")
unknown_id_error <- tryCatch(
  agent_normalize_state_update(
    list(features = "Missing"), c("Feature"), c("ORA"), rownames(fd), rownames(pd)
  ),
  error = function(e) conditionMessage(e)
)
ok(
  ut_cmp_identical(grepl("Unknown feature ID", unknown_id_error, fixed = TRUE), TRUE),
  "unknown assistant feature IDs are rejected"
)

views <- data.frame(
  id = "volcano", label = "Volcano", x = "ttest|A_vs_B|mean.diff",
  y = "ttest|A_vs_B|log.fdr", description = "volcano", source = "auto",
  stringsAsFactors = FALSE
)
view <- agent_normalize_scatter_view(
  space = "feature", quick_view_id = "volcano", quick_views = list(feature = views)
)
ok(ut_cmp_identical(view$x_axis, "ttest|A_vs_B|mean.diff"), "quick-view X axis is resolved")
ok(ut_cmp_identical(view$mode, "quick"), "quick-view mode is preserved")
ok(
  ut_cmp_error(
    agent_normalize_scatter_view(
      space = "feature", x_axis = "missing", y_axis = "ttest|A_vs_B|log.fdr",
      feature_columns = colnames(fd)
    ),
    "Unknown feature X-axis annotation: missing"
  ),
  "unknown custom scatter axes are rejected"
)
