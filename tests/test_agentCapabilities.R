# WP8/WP9 unit suite: capability registry + discovery helpers
# (auxi_agentCapabilities.R) and the mapping widget kind (auxi_widgetStore.R).

library(omicsViewer)
library(unittest, quietly = TRUE)

agent_capability_records <- omicsViewer:::agent_capability_records
agent_capability_search <- omicsViewer:::agent_capability_search
agent_capability_get <- omicsViewer:::agent_capability_get
agent_capability_summary <- omicsViewer:::agent_capability_summary
agent_tool_capabilities <- omicsViewer:::agent_tool_capabilities
agent_normalize_enrichment_update <- omicsViewer:::agent_normalize_enrichment_update
agent_normalize_table_view_update <- omicsViewer:::agent_normalize_table_view_update
agent_compact_state <- omicsViewer:::agent_compact_state
widget_store_new <- omicsViewer:::widget_store_new
widget_store_child <- omicsViewer:::widget_store_child
widget_binding <- omicsViewer:::widget_binding
store_register <- omicsViewer:::store_register
store_read <- omicsViewer:::store_read
store_apply <- omicsViewer:::store_apply

mk_store <- function() {
  s <- widget_store_new()
  ds <- widget_store_child(s, "dataspace")
  store_register(
    ds,
    widget_binding("feature_space.x_analysis", "select", label = "X analysis",
      help = "Feature scatter X-axis analysis",
      choices_provider = function(v) c("ttest", "PCA")),
    widget_binding("tab_pheno.multi_selection", "boolean",
      label = "Multiple selection",
      help = "Allow selecting more than one row in the sample table")
  )
  rs <- widget_store_child(s, "resultspace")
  store_register(
    rs,
    widget_binding("ora.xax_analysis", "select", label = "Collapse category",
      help = "Annotation category whose values the input features are collapsed on",
      choices_provider = function(v) c("ttest")),
    widget_binding("internal_counter", "numeric", internal = TRUE,
      label = "Internal", help = "never agent-visible")
  )
  s
}

## ------------------------------------------------- records generation ----
s <- mk_store()
records <- agent_capability_records(s)
ids <- vapply(records, function(r) r$id, character(1))
ok(
  all(c("dataspace.feature_space.x_analysis", "dataspace.tab_pheno.multi_selection",
        "resultspace.ora.xax_analysis") %in% ids) &&
    !any(grepl("internal_counter", ids)),
  "widget records are generated from agent-visible bindings only"
)
ok(
  identical(length(records), 3L + length(agent_tool_capabilities())),
  "records = widget records + every tool record"
)
ora_rec <- Filter(function(r) identical(r$id, "resultspace.ora.xax_analysis"),
                  records)[[1]]
ok(
  identical(ora_rec$panel, "Analysis space") &&
    identical(ora_rec$operation, "set_enrichment_parameters"),
  "widget records map to panels and their covering semantic tool"
)
scatter_rec <- Filter(
  function(r) identical(r$id, "dataspace.feature_space.x_analysis"), records)[[1]]
ok(
  identical(scatter_rec$panel, "Data space") &&
    identical(scatter_rec$operation, "set_scatter_view"),
  "scatter widgets map to set_scatter_view"
)
tab_rec <- Filter(
  function(r) identical(r$id, "dataspace.tab_pheno.multi_selection"), records)[[1]]
ok(
  identical(tab_rec$operation, "set_table_view"),
  "data-space table widgets map to set_table_view"
)
ok(
  identical(agent_capability_records(NULL) |> length(),
            length(agent_tool_capabilities())),
  "NULL store yields the tool records only"
)
ok(
  identical(setequal(
    vapply(Filter(function(r) startsWith(r$id, "tool:"), records),
           function(r) r$operation, character(1)),
    c("get_omics_viewer_state", "search_annotations", "summarize_annotation",
      "set_omics_viewer_state", "set_scatter_view", "set_enrichment_parameters",
      "set_table_view", "create_figure", "update_figure", "list_widgets",
      "get_widget", "set_widgets", "search_ui_capabilities",
      "get_ui_capability")),
    TRUE),
  "tool records cover the full assistant allowlist incl. WP8/WP9 tools"
)

## ---------------------------------------------------------- searching ----
res <- agent_capability_search(s, "enrichment")
ok(
  identical(res$capability_count, 3L + length(agent_tool_capabilities())) &&
    res$match_count >= 2L &&
    identical(vapply(res$capabilities, function(r) r$id, character(1))[1],
              "resultspace.ora.xax_analysis"),
  "search matches tool help text and widget help text by meaning"
)
ok(
  (function(x) length(x$capabilities) == 2L && isTRUE(x$truncated))(
    agent_capability_search(s, "s", max_results = 2L)),
  "search respects max_results and reports truncation"
)
ok(
  ut_cmp_error(agent_capability_search(s, "null"),
               "non-empty capability search query"),
  "sentinel query strings are rejected"
)
ok(
  ut_cmp_error(agent_capability_search(s, "x", max_results = 0L),
               "max_results must be an integer from 1 through 50"),
  "max_results bounds are enforced"
)
empty <- agent_capability_search(s, "zzzznope")
ok(
  identical(empty$match_count, 0L) && identical(empty$capabilities, list()),
  "zero-hit searches return an empty bounded result, not an error"
)

## -------------------------------------------------------- describing ----
d <- agent_capability_get(s, "resultspace.ora.xax_analysis")
ok(
  identical(d$operation, "set_enrichment_parameters") &&
    identical(d$kind, "select"),
  "get returns the widget record"
)
d2 <- agent_capability_get(s, "set_table_view")
ok(
  identical(d2$id, "tool:set_table_view") && identical(d2$writable, TRUE),
  "get accepts raw semantic tool ids"
)
d3 <- agent_capability_get(s, "tool:set_table_view")
ok(identical(d3$id, d2$id), "get accepts tool-prefixed ids")
ok(
  ut_cmp_error(agent_capability_get(s, "set_table_vie"),
               paste("Unknown capability id: set_table_vie.",
                     "Closest matches: set_table_view.")),
  "get suggests the closest id on typos"
)

## ---------------------------------------------------------- overview ----
sm <- agent_capability_summary(s)
ok(
  identical(sm$capability_count, 3L) &&
    setequal(sm$panels, c("Analysis space", "Data space")) &&
    "set_enrichment_parameters" %in% sm$semantic_tools,
  "summary carries counts, panels, and semantic tool ids only"
)
ok(identical(agent_capability_summary(NULL), NULL), "no store, no summary")

# through agent_compact_state (the get_omics_viewer_state overview)
state <- list(
  dataset = list(id = "x", class = "ExpressionSet",
                 dimensions = c(features = 3L, samples = 4L)),
  app = list(data_active_tab = "Feature", analysis_active_tab = "Feature"),
  selection = list(features = list(), samples = list()),
  panels = list(),
  policy = "test"
)
cs <- agent_compact_state(state = state, store = s)
ok(
  identical(cs$capabilities$capability_count, 3L) &&
    is.null(cs$capabilities$capabilities),
  "compact overview carries capability counts, never contents"
)
cs2 <- agent_compact_state(state = state, store = NULL)
ok(
  is.null(cs2$capabilities),
  "overview without a store omits the capabilities block"
)

## ---------------------------------------------- WP8 enrichment tool ----
fd <- data.frame(
  row.names = c("G1", "G2", "G3"),
  `ttest|A_vs_B|log.fdr` = c(1, 2, 3),
  `ttest|A_vs_B|md` = c(-1, 0, 1),
  check.names = FALSE
)
attr(fd, "GS") <- data.frame(
  featureId = c("G1", "G2", "G3"),
  gsId = c("gs1", "gs1", "gs2")
)

up <- agent_normalize_enrichment_update(
  list(method = "ora", collapse = "ttest|A_vs_B|log.fdr",
       selected_pathway = "gs1"), fd)
ok(
  identical(up$method, "ora") && identical(up$tab, "ORA") &&
    identical(up$patch$`resultspace.ora.xax_analysis`, "ttest") &&
    identical(up$patch$`resultspace.ora.xax_subset`, "A_vs_B") &&
    identical(up$patch$`resultspace.ora.xax_variable`, "log.fdr") &&
    identical(up$patch$`resultspace.ora.selected_row`, "gs1") &&
    identical(up$patch$`resultspace.analyst_tab`, "ORA"),
  "enrichment update builds the full store patch incl. tab switch"
)
ok(
  identical(up$method, "ora"),
  "method label mapping"
)
up2 <- agent_normalize_enrichment_update(
  list(method = "fgsea", collapse = "ttest|A_vs_B|md"), fd)
ok(
  identical(up2$patch$`resultspace.fgsea.xax_variable`, "md") &&
    identical(up2$patch$`resultspace.analyst_tab`, "fGSEA"),
  "fgsea method routes to the fgsea namespace"
)
ok(
  ut_cmp_error(
    agent_normalize_enrichment_update(list(method = "ora"), fd),
    "Enrichment update contains no changes."),
  "empty enrichment update is rejected"
)
ok(
  ut_cmp_error(
    agent_normalize_enrichment_update(
      list(method = "ora", collapse = "log.fdr"), fd),
    "collapse must be a full 'Category|Subcategory|Variable'"),
  "partial collapse triples are rejected with the naming convention"
)
ok(
  grepl("Unknown feature annotation",
        tryCatch(agent_normalize_enrichment_update(
          list(method = "ora", collapse = "ttest|A_vs_B|logg.fdrr"), fd),
          error = function(e) conditionMessage(e))) &&
    grepl("Closest matches: ttest\\|A_vs_B\\|log.fdr",
          tryCatch(agent_normalize_enrichment_update(
            list(method = "ora", collapse = "ttest|A_vs_B|logg.fdrr"), fd),
            error = function(e) conditionMessage(e))),
  "unknown collapse columns get closest-match suggestions"
)
ok(
  ut_cmp_error(
    agent_normalize_enrichment_update(
      list(method = "gsea", collapse = "ttest|A_vs_B|md"), fd),
    "Unknown enrichment method: gsea."),
  "unknown methods are rejected"
)
fd_no_gs <- fd
attr(fd_no_gs, "GS") <- NULL
ok(
  ut_cmp_error(
    agent_normalize_enrichment_update(
      list(method = "ora", collapse = "ttest|A_vs_B|md"), fd_no_gs),
    "no gene-set annotations"),
  "datasets without gene sets report the missing capability"
)
ok(
  ut_cmp_error(
    agent_normalize_enrichment_update(
      list(method = "ora", wrong = 1), fd),
    "Unknown enrichment field"),
  "unknown fields are rejected"
)

## ------------------------------------------------ WP8 table-view tool ----
tv <- agent_normalize_table_view_update(list(
  table = "feature_table",
  columns = c("General|Gene|symbol", "ttest|A_vs_B|log.fdr"),
  multi_selection = TRUE,
  column_filters = list("ttest|A_vs_B|log.fdr" = "2"),
  page = 3L
))
ok(
  identical(tv$tab, "Feature table") &&
    identical(tv$patch$`dataspace.tab_feature.columns`,
              c("General|Gene|symbol", "ttest|A_vs_B|log.fdr")) &&
    identical(tv$patch$`dataspace.tab_feature.multi_selection`, TRUE) &&
    identical(tv$patch$`dataspace.tab_feature.column_filters`,
              c("ttest|A_vs_B|log.fdr" = "2")) &&
    identical(tv$patch$`dataspace.tab_feature.page`, 3L) &&
    identical(tv$patch$`dataspace.active_tab`, "Feature table"),
  "table-view update builds the full store patch incl. tab switch"
)
tv2 <- agent_normalize_table_view_update(list(
  table = "sample_table", column_filters = '{"group": "KO"}'))
ok(
  identical(tv2$patch$`dataspace.tab_pheno.column_filters`, c(group = "KO")) &&
    identical(tv2$patch$`dataspace.active_tab`, "Sample table"),
  "column_filters accepts JSON-object strings"
)
ok(
  identical(tv2$patch$`dataspace.tab_pheno.column_filters`,
            agent_normalize_table_view_update(list(
              table = "sample_table",
              column_filters = list(group = "KO")))$patch$`dataspace.tab_pheno.column_filters`),
  "JSON-string and named-list column_filters normalize identically"
)
tv3 <- agent_normalize_table_view_update(list(
  table = "expression_table", clear_filters = "true"))
ok(
  identical(tv3$patch$`dataspace.tab_expr.column_filters`,
            setNames(character(0), character(0))),
  "clear_filters produces an explicit empty mapping"
)
ok(
  is.null(agent_normalize_table_view_update(list(
    table = "sample_table", page = 2L,
    column_filters = "{}"))$patch$`dataspace.tab_pheno.column_filters`),
  "sentinel column_filters strings are treated as omitted"
)
ok(
  ut_cmp_error(agent_normalize_table_view_update(list(table = "gene_table")),
               "Unknown table: gene_table."),
  "unknown tables are rejected with suggestions"
)
ok(
  ut_cmp_error(agent_normalize_table_view_update(list(
    table = "feature_table", columns = character())),
    "columns cannot be empty"),
  "empty column arrays are rejected (min 1 entry)"
)
ok(
  ut_cmp_error(agent_normalize_table_view_update(list(
    table = "feature_table", page = 0L)),
    "page must be an integer of at least 1"),
  "page bounds are enforced"
)
ok(
  ut_cmp_error(agent_normalize_table_view_update(list(table = "feature_table")),
               "Table-view update contains no changes."),
  "empty table-view update is rejected"
)
ok(
  ut_cmp_error(agent_normalize_table_view_update(list(
    table = "feature_table", bogus = 1)),
    "Unknown table-view field"),
  "unknown fields are rejected"
)

## --------------------------------------------- mapping widget kind ----
ms <- widget_store_new()
store_register(ms,
  widget_binding("t.filters", "mapping", label = "Column filters",
    help = "Per-column search patterns",
    choices_provider = function(v) c("a", "b", "c")))
r <- store_apply(ms, list(`t.filters` = list(a = "x", b = "y")), origin = "agent")
ok(
  identical(store_read(ms)$`t.filters`, c(a = "x", b = "y")),
  "mapping values store as named character vectors"
)
r <- store_apply(ms, list(`t.filters` = setNames(character(0), character(0))),
                origin = "agent")
ok(
  identical(store_read(ms)$`t.filters`, setNames(character(0), character(0))),
  "empty mapping clears all entries"
)
ok(
  ut_cmp_error(
    store_apply(ms, list(`t.filters` = list(z = "x")), origin = "agent"),
    "Unknown key\\(s\\) for t.filters: z."),
  "mapping keys validate against choices with suggestions"
)
ok(
  ut_cmp_error(
    store_apply(ms, list(`t.filters` = list(a = "")), origin = "agent"),
    "requires mapping values of 1-200"),
  "mapping values must be non-empty"
)
ok(
  ut_cmp_error(
    store_apply(ms, list(`t.filters` = list("x")), origin = "agent"),
    "requires a named key-to-value mapping"),
  "unnamed mapping lists are rejected"
)
ok(
  identical(store_apply(ms, list(`t.filters` = "null"), origin = "agent") |>
              (\(x) store_read(ms)$`t.filters`)(),
            setNames(character(0), character(0))),
  "literal sentinel strings for mappings are treated as omitted (kept)"
)
