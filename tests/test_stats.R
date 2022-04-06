library(unittest)
correlationAnalysis <- omicsViewer:::correlationAnalysis

ph <- data.frame(
  set  = 1:15
)
expr <- rbind(1:15, matrix(rnorm(150), 10, 15))
res <- correlationAnalysis(expr, ph, prefix = "test")
ok(ut_cmp_equal(
  colnames(res), 
  c("test|set|R", "test|set|N", "test|set|P", "test|set|logP", "test|set|range")),
  "correlationAnalysis - colnames"
  )
ok(
  ut_cmp_equal( nrow(res), 11 ),
  "correlationAnalysis - row numbers"
)
ok(
  ut_cmp_equal( res[1, 1], 1 ),
  "correlationAnalysis - correlation 1"
)


ph <- data.frame(
  var1  = rep(LETTERS[1:2], each = 6),
  var2  = rep(c("C", "D", "C", "D"), each = 3),
  stringsAsFactors = FALSE, 
  row.names = paste("S", 1:12, sep = "")
)
cmp <- rbind(
  c("var1", "A", "B"),
  c("var2", "C", "D")
)
expr <- cbind(matrix(rnorm(60), 10), matrix(rnorm(60, mean = 3), 10))
colnames(expr) <- rownames(ph)
rownames(expr) <- paste("P", 1:10, sep = "")

multi.t.test <- omicsViewer:::multi.t.test
res <- multi.t.test(x = expr, pheno = ph, compare = cmp)
ok(
  ut_cmp_equal(all(res$`ttest|A_vs_B|pvalue` < res$`ttest|C_vs_D|pvalue`), TRUE),
  "multi.t.test - significance"
)

exprspca <- omicsViewer:::exprspca
pcres <- exprspca(expr, n = 3)
ok(ut_cmp_equal(dim(pcres$samples), c(12, 3)), "exprspca sample")
ok(ut_cmp_equal(dim(pcres$features), c(10, 3)), "exprspca sample")

exprNA <- expr
exprNA[1, 1:11] <- NA
fillNA <- omicsViewer:::fillNA
res <- fillNA(exprNA)
ok(ut_cmp_equal(length(unique(res[1, 1:11])), 1), "fillNA")
