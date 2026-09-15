library(omicsViewer)
library(Biobase)
library(unittest, quietly = TRUE)

prepare_quick_views <- omicsViewer:::prepare_quick_views
parse_quick_views <- omicsViewer:::parse_quick_views
detect_quick_views <- omicsViewer:::detect_quick_views
active_quick_view <- omicsViewer:::active_quick_view
quick_badges_module <- omicsViewer:::quick_badges_module

triset <- rbind(
  c("PCA", "All", "PC1(10%)"),
  c("PCA", "All", "PC2(7%)"),
  c("ttest", "A_B", "mean.diff"),
  c("ttest", "A_B", "log.fdr")
)

meta <- data.frame(
  "PCA|All|PC1(10%)" = 1:3,
  "PCA|All|PC2(7%)" = 3:1,
  "ttest|A_B|mean.diff" = c(-1, 0, 1),
  "ttest|A_B|log.fdr" = c(2, 1, 3),
  check.names = FALSE
)

views <- detect_quick_views(meta, triset)
ok(
  ut_cmp_equal(views$id, c("pca", "volcano_A_B")),
  "detect PCA and volcano quick views"
)
ok(
  ut_cmp_equal(views$x, c("PCA|All|PC1(10%)", "ttest|A_B|mean.diff")),
  "quick-view X axes"
)
ok(
  ut_cmp_equal(views$y, c("PCA|All|PC2(7%)", "ttest|A_B|log.fdr")),
  "quick-view Y axes"
)

meta2 <- meta
meta2$"ttest|A_B|log.fdr" <- NULL
meta2$"ttest|A_B|log.pvalue" <- 1:3
triset2 <- rbind(
  triset,
  c("ttest", "A_B", "log.pvalue")
)
views2 <- detect_quick_views(meta2, triset2)
ok(
  ut_cmp_equal(views2$y[views2$id == "volcano_A_B"], "ttest|A_B|log.pvalue"),
  "volcano falls back to log.pvalue"
)

custom <- list(
  PCA_all = c(
    "PCA|All|PC1(10%)",
    "PCA|All|PC2(7%)",
    name = "PCA all feature"
  )
)
attr(meta, "shortcut") <- custom
views3 <- prepare_quick_views(meta, triset)
ok(
  ut_cmp_equal(views3$id, c("PCA_all", "volcano_A_B")),
  "custom shortcut is used and auto views are retained"
)
ok(
  ut_cmp_equal(views3$label[1], "PCA all feature"),
  "compact shortcut name becomes the badge label"
)

invalid <- list(
  good = c("PCA|All|PC1(10%)", "PCA|All|PC2(7%)"),
  badColumn = c("PCA|All|PC1(10%)", "Not|Available|Here"),
  malformed = c("PCA|All|PC1(10%)", "PCA")
)
ok(
  ut_cmp_equal(nrow(parse_quick_views(invalid, triset)), 1),
  "invalid shortcut axes are dropped"
)

duplicate <- list(
  x = list(x = "PCA|All|PC1(10%)", y = "PCA|All|PC2(7%)"),
  x2 = list(x = "ttest|A_B|mean.diff", y = "ttest|A_B|log.fdr")
)
ok(
  ut_cmp_equal(parse_quick_views(duplicate, triset)$id, c("x", "x2")),
  "shortcut ids are preserved and made available"
)

ok(
  ut_cmp_equal(
    active_quick_view(
      views3,
      list(analysis = "PCA", subset = "All", variable = "PC1(10%)"),
      list(analysis = "PCA", subset = "All", variable = "PC2(7%)")
    ),
    "PCA_all"
  ),
  "active badge is derived from current triselector values"
)
ok(
  length(active_quick_view(
    views3,
    list(analysis = "PCA", subset = "All", variable = "PC1(10%)"),
    list(analysis = "ttest", subset = "A_B", variable = "log.fdr")
  )) == 0,
  "non-shortcut axes have no active badge"
)

metaCor <- data.frame(
  "Cor|MDR|R" = -1:1,
  "Cor|MDR|logP" = 1:3,
  check.names = FALSE
)
trisetCor <- rbind(
  c("Cor", "MDR", "R"),
  c("Cor", "MDR", "logP")
)
viewsCor <- detect_quick_views(metaCor, trisetCor)
ok(
  ut_cmp_equal(viewsCor$id, "cor_MDR"),
  "detect correlation quick view"
)
ok(
  ut_cmp_equal(viewsCor$x, "Cor|MDR|R"),
  "correlation quick-view X axis"
)
ok(
  ut_cmp_equal(viewsCor$y, "Cor|MDR|logP"),
  "correlation quick-view Y axis"
)

set.seed(1234)
exprTest <- matrix(rnorm(100), nrow = 5, ncol = 20)
rownames(exprTest) <- paste0("f", 1:5)
colnames(exprTest) <- paste0("s", 1:20)
pdTest <- data.frame(
  score = rnorm(20),
  group = rep(c("a", "b"), each = 10),
  row.names = colnames(exprTest)
)
fdTest <- data.frame(name = rownames(exprTest), row.names = rownames(exprTest))
prepared <- prepOmicsViewer(
  expr = exprTest,
  pData = pdTest,
  fData = fdTest,
  PCA = TRUE,
  pca.fillNA = FALSE,
  t.test = rbind(c("group", "a", "b")),
  ttest.fillNA = FALSE,
  SummarizedExperiment = FALSE
)
ok(
  ut_cmp_equal(
    c("cor_score", "volcano_a_vs_b") %in% attr(Biobase::fData(prepared), "quickViews")$id,
    c(TRUE, TRUE)
  ),
  "prepOmicsViewer adds correlation and t-test shortcuts"
)
ok(
  ut_cmp_equal("pca" %in% attr(Biobase::pData(prepared), "quickViews")$id, TRUE),
  "prepOmicsViewer adds PCA sample shortcut"
)

# Rendering only the active state must not recreate badge actionButtons: doing
# so resets their click counters and makes the first switch appear to fail.
quickClickViews <- data.frame(
  id = c("pca", "volcano"),
  label = c("PCA", "Volcano"),
  x = c("PCA|All|PC1", "ttest|A_B|mean.diff"),
  y = c("PCA|All|PC2", "ttest|A_B|log.fdr"),
  description = c("PCA", "Volcano"),
  source = "test",
  stringsAsFactors = FALSE
)
quickClickResult <- new.env(parent = emptyenv())

quickClickApp <- function(input, output, session) {
  active <- shiny::reactiveVal(NULL)
  selected <- quick_badges_module(
    "quick",
    views = shiny::reactive(quickClickViews),
    activeId = active
  )
  quickClickResult$selected <- selected
  shiny::observeEvent(selected()$trigger, active(selected()$view$id))
}

shiny::testServer(quickClickApp, {
  session$flushReact()
  session$flushOutput()
  session$setInputs(`quick-badge_pca` = 1)
  session$flushReact()
  session$flushOutput()
  ok(
    ut_cmp_equal(
      c(quickClickResult$selected()$trigger, quickClickResult$selected()$view$id),
      c(1, "pca")
    ),
    "first quick-view click switches immediately"
  )

  session$setInputs(`quick-badge_volcano` = 1)
  session$flushReact()
  session$flushOutput()
  ok(
    ut_cmp_equal(
      c(quickClickResult$selected()$trigger, quickClickResult$selected()$view$id),
      c(2, "volcano")
    ),
    "second quick-view click switches immediately"
  )
})
