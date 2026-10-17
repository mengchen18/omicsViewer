# WP-7 (F4): panel_results - readable enrichment tables.
#
# Pure helper tests for agent_enrichment_table_excerpt plus testServer
# integration of the result-space reader: computed vs uncomputed panels,
# method filtering, and the snapshot-safety of the attribute transport
# (result tables must never leak into .sanitize_widget_state output).
#
# Run: Rscript tests/test_agentPanelResults.R

library(omicsViewer)
library(unittest, quietly = TRUE)

agent_enrichment_table_excerpt <- omicsViewer:::agent_enrichment_table_excerpt
agent_compact_state <- omicsViewer:::agent_compact_state
build_app_state <- omicsViewer:::build_app_state
widget_store_new <- omicsViewer:::widget_store_new
store_apply <- omicsViewer:::store_apply
store_read <- omicsViewer:::store_read

## ----------------------------------------------- pure excerpt helper ----
ok(
  identical(agent_enrichment_table_excerpt(NULL, "ORA")$status, "not_run") &&
    grepl("set_enrichment_parameters",
          agent_enrichment_table_excerpt(NULL, "ORA")$message, fixed = TRUE),
  "excerpt: NULL table reports not_run with the run hint"
)
ok(
  identical(agent_enrichment_table_excerpt("Too few feature IDs", "ORA")$status,
            "message") &&
    identical(agent_enrichment_table_excerpt("Too few feature IDs", "ORA")$message,
              "Too few feature IDs"),
  "excerpt: character tables pass the panel message through"
)

# real ORA table from the demo dataset
dat <- readRDS(system.file(package = "omicsViewer", "extdata/demo.RDS"))
obj <- omicsViewer:::tallGS(dat)
fd <- Biobase::fData(obj)
fdgs <- attr(fd, "GS")
sel <- rownames(fd)[fd[["ttest|RE_vs_ME|pvalue"]] < 0.01]
ora_tab <- omicsViewer:::vectORATall(fdgs, i = sel, background = nrow(fd))
x <- agent_enrichment_table_excerpt(ora_tab, "ORA")
ok(identical(x$status, "ok"), "excerpt: real ORA table yields ok rows")
ok(
  identical(utils::head(x$columns, 3L), c("pathway", "size_overlap", "size_geneset")) &&
    !"overlap_ids" %in% x$columns,
  "excerpt: ORA key columns lead and the overlap_ids list column is dropped"
)
ok(
  x$row_count == nrow(ora_tab) &&
    x$rows[[1]]$pathway == ora_tab$pathway[[1]],
  "excerpt: ORA rows carry real pathway ids in table order"
)
ok(
  !any(vapply(x$rows, function(r)
    any(vapply(r, function(v) is.numeric(v) && is.finite(v) &&
                 nchar(format(v, scientific = FALSE)) > 8, logical(1))),
    logical(1))),
  "excerpt: numeric values are rounded"
)

# real fGSEA table with the module's renamed columns
stats <- fd[, "ttest|RE_vs_ME|mean.diff"]
names(stats) <- rownames(fd)
stats <- stats[!is.na(stats)]
fdgs2 <- fdgs[fdgs$featureId %in% names(stats), ]
fg_tab <- omicsViewer:::fgsea1(fdgs2, stats = stats, minSize = 3, maxSize = 500)
cn <- colnames(fg_tab)
cn[cn == "ES"] <- "enrichment score (ES)"
cn[cn == "NES"] <- "normalized ES"
colnames(fg_tab) <- cn
xf <- agent_enrichment_table_excerpt(fg_tab, "fGSEA")
ok(
  identical(utils::head(xf$columns, 2L), c("pathway", "normalized ES")) &&
    !"leadingEdge" %in% xf$columns,
  "excerpt: fGSEA key columns lead (NES first stat) and leadingEdge is dropped"
)
ok(
  xf$rows[[1]]$pathway == fg_tab$pathway[[1]] &&
    is.finite(xf$rows[[1]]$`normalized ES`),
  "excerpt: fGSEA rows carry real pathway ids and NES values"
)
ok(
  nchar(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null"), type = "bytes") <
    4096L,
  "excerpt: ORA excerpt stays well under the tool budget"
)
ok(
  nchar(jsonlite::toJSON(xf, auto_unbox = TRUE, null = "null"), type = "bytes") <
    4096L,
  "excerpt: fGSEA excerpt stays well under the tool budget"
)

# byte shaping under a tiny budget: columns drop, then rows (floor 3)
xt <- agent_enrichment_table_excerpt(ora_tab, "ORA", n = 15L, budget_bytes = 600L)
ok(
  length(xt$rows) >= 3L &&
    nchar(jsonlite::toJSON(
      list(columns = xt$columns, rows = xt$rows),
      auto_unbox = TRUE, null = "null"), type = "bytes") <= 620L,
  "excerpt: tiny budget keeps >=3 rows inside the cap (codec never truncates)"
)
ok(
  setequal(xt$columns_available,
           setdiff(setdiff(colnames(ora_tab), "overlap_ids"), xt$columns)),
  "excerpt: dropped columns are reported in columns_available"
)
ok(
  identical(agent_enrichment_table_excerpt(ora_tab, "ORA", n = 3L)$returned, 3L) &&
    isTRUE(agent_enrichment_table_excerpt(ora_tab, "ORA", n = 3L)$rows_truncated),
  "excerpt: n bounds the returned rows and flags truncation"
)

## ----------------------------------- state section: compose + guards ----
full_state <- build_app_state(
  dataset = NULL, dataset_id = "unit.RDS",
  data_status = list(eset_active_tab = "Feature"),
  result_status = list(analyst_active_tab = "Feature"),
  selected_features = "G1", selected_samples = character(),
  label = "unit"
)
overview <- agent_compact_state(state = full_state)
ok(
  identical(overview$available_sections,
            c("annotations", "quick_views", "panels", "panel_results",
              "figure_grammar")),
  "section menu now advertises panel_results"
)
ok(is.null(overview$panel_results),
   "overview omits panel_results unless requested")
with_pr <- agent_compact_state(
  state = full_state,
  panel_results = list(ora = x, fgsea = xf),
  sections = "panel_results"
)
ok(
  identical(with_pr$panel_results$ora, x) &&
    identical(with_pr$panel_results$fgsea, xf) &&
    is.null(with_pr$annotations),
  "requested panel_results section returns the reader payload only"
)
ok(
  ut_cmp_error(
    agent_compact_state(state = full_state, sections = "panel_result"),
    "Unknown state section"),
  "misspelled panel_results is still rejected with a suggestion"
)

## --------------------------------- result-space reader (testServer) -----
# Drive the REAL L1 result-space module: the enrichment modules expose
# their tables as attributes on their returns; the L1 status value carries
# the reader attribute; snapshots (.sanitize_widget_state) never see the
# tables themselves. Headless sessions report every output visible, so the
# visibility-gated fgsea mirror computes once the ranking settles.
.prd_fd <- fd[1:120, , drop = FALSE]
attr(.prd_fd, "GS") <- fdgs[fdgs$featureId %in% rownames(.prd_fd), ]
.prd_pd <- Biobase::pData(obj)[1:10, , drop = FALSE]
.prd_ex <- Biobase::exprs(dat)[1:120, 1:10, drop = FALSE]
.prd_sel <- rownames(.prd_fd)[.prd_fd[["ttest|RE_vs_ME|pvalue"]] < 0.05]
.prd_root <- widget_store_new()
.prd_env <- new.env(); .prd_env$ret <- NULL
.prd_host <- function() function(input, output, session) {
  .prd_env$ret <- omicsViewer:::L1_result_space_module(
    "rs",
    reactive_expr = shiny::reactive(.prd_ex),
    reactive_phenoData = shiny::reactive(.prd_pd),
    reactive_featureData = shiny::reactive(.prd_fd),
    reactive_i = shiny::reactive(.prd_sel),
    reactive_highlight = shiny::reactive(NULL),
    store = .prd_root)
}
.prd_flush <- function(session, n = 8)
  for (i in seq_len(n)) tryCatch(session$flushReact(), error = function(e) NULL)
shiny::testServer(.prd_host(), {
  .prd_flush(session, 5)
  val <- .prd_env$ret()
  reader <- attr(val, "read_enrichment_results")
  ok(is.function(reader), "L1 status value carries the reader attribute")
  ok(
    is.null(val$analyst_ora$table) && is.null(val$analyst_fgsea$table) &&
      is.null(attr(val$analyst_ora, "table")) &&
      is.null(attr(val$analyst_fgsea, "table")),
    "sanitized status never carries result tables (snapshot safety)"
  )
  # uncomputed fGSEA: the ranking req()s a settled variable, so nothing
  # ran yet. (ORA runs on the raw selection with an unset collapse triple
  # by historical contract - its not_run path is covered by the empty-
  # selection host below.)
  .prd_env$cold <- reader()
  # user-path cascade picks (testServer quirk: update*Input messages are
  # not applied, so drive the triselector inputs like test_agentWidgets)
  session$setInputs(`rs-fgsea-tris_fgsea-analysis` = "ttest")
  session$setInputs(`rs-fgsea-tris_fgsea-subset` = "RE_vs_ME")
  session$setInputs(`rs-fgsea-tris_fgsea-variable` = "mean.diff")
  session$setInputs(`rs-ora-tris_ora-analysis` = "General")
  session$setInputs(`rs-ora-tris_ora-subset` = "All")
  session$setInputs(`rs-ora-tris_ora-variable` = "Gene.name")
  .prd_flush(session)
  .prd_env$warm <- reader()
  .prd_env$warm_fg <- reader(method = "fgsea")
  .prd_env$bad <- tryCatch(reader(method = "gsea"),
                           error = function(e) conditionMessage(e))
})
ok(
  identical(.prd_env$cold$fgsea$status, "not_run") &&
    grepl("set_enrichment_parameters", .prd_env$cold$fgsea$message,
          fixed = TRUE),
  "reader: uncomputed fGSEA returns an explicit not_run note"
)
ok(
  identical(.prd_env$warm$fgsea$status, "ok") &&
    length(.prd_env$warm$fgsea$rows) > 0L &&
    identical(.prd_env$warm$fgsea$columns[[2]], "normalized ES"),
  "reader: fGSEA rows are readable after the ranking settles (NES column first)"
)
ok(
  identical(.prd_env$warm$ora$status, "ok") &&
    length(.prd_env$warm$ora$rows) > 0L &&
    identical(.prd_env$warm$ora$columns[[1]], "pathway"),
  "reader: ORA rows are readable after the collapse cascade commits"
)
# ORA not_run path: an empty selection leaves oraTab NULL (the observer
# req()s a selection); the excerpt must say so instead of erroring
.prd_env$cold_ora <- NULL
.prd_root_cold <- widget_store_new()
.prd_host_cold <- function() function(input, output, session) {
  .prd_env$ret <- omicsViewer:::L1_result_space_module(
    "rs",
    reactive_expr = shiny::reactive(.prd_ex),
    reactive_phenoData = shiny::reactive(.prd_pd),
    reactive_featureData = shiny::reactive(.prd_fd),
    reactive_i = shiny::reactive(NULL),
    reactive_highlight = shiny::reactive(NULL),
    store = .prd_root_cold)
}
shiny::testServer(.prd_host_cold(), {
  .prd_flush(session, 5)
  reader <- attr(.prd_env$ret(), "read_enrichment_results")
  .prd_env$cold_ora <- reader(method = "ora")
})
ok(
  identical(.prd_env$cold_ora$ora$status, "not_run") &&
    grepl("set_enrichment_parameters", .prd_env$cold_ora$ora$message,
          fixed = TRUE),
  "reader: ORA without a selection returns an explicit not_run note"
)
ok(
  identical(names(.prd_env$warm), c("ora", "fgsea")) &&
    identical(names(.prd_env$warm_fg), "fgsea"),
  "reader: omitted method returns both panels, exact method filters"
)
ok(
  grepl("Unknown enrichment method", .prd_env$bad, fixed = TRUE) &&
    grepl("'ora'", .prd_env$bad, fixed = TRUE),
  "reader: unknown methods are rejected with the allowed vocabulary"
)
# L0 fallback contract (read_agent_panel_results): a missing reader
# degrades to the same explicit not_run notes, never an error
ok(
  identical(agent_enrichment_table_excerpt(NULL, "fGSEA")$status, "not_run"),
  "L0 fallback contract: missing reader still yields not_run notes"
)
