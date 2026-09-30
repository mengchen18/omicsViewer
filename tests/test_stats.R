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

# R-M12: var.equal must be a formal, not a ... entry - a documented
# var.equal = FALSE call used to collide with the hard-coded TRUE inside
# the per-row t.test call ("formal argument ... matched by multiple
# actual arguments"), the error was swallowed per row and every p-value
# came back NA.
res_ve <- multi.t.test(x = expr, pheno = ph, compare = cmp, var.equal = FALSE)
ok(
  ut_cmp_equal(all(!is.na(res_ve$`ttest|A_vs_B|pvalue`)), TRUE),
  "multi.t.test - var.equal = FALSE yields p-values (R-M12)"
)
ok(
  ut_cmp_equal(identical(
    res_ve$`ttest|A_vs_B|pvalue`,
    multi.t.test(x = expr, pheno = ph, compare = cmp, var.equal = TRUE)$`ttest|A_vs_B|pvalue`),
    FALSE),
  "multi.t.test - var.equal = FALSE differs from TRUE (R-M12, Welch)"
)

# R-L4: read.proteinGroups.tmt - with NO contaminator/reverse/site-only
# rows, ir was integer(0) and ab[-ir, ] subset every table to ZERO rows
# (the whole dataset silently disappeared).
tf <- tempfile(fileext = ".txt")
# full column families (summed numeric-suffix columns + plain individual
# channels), as in a real proteinGroups.txt
writeLines(c(
  paste("Majority.protein.IDs", "Only.identified.by.site",
        "Fraction.1", "Reporter.intensity.corrected.1", "Reporter.intensity.corrected.2",
        "Reporter.intensity.corrected", "Reporter.intensity.count.1",
        "Reporter.intensity.count", "Reporter.intensity.1", "Reporter.intensity",
        sep = "\t"),
  paste("P1", "-", "5", "10", "20", "10", "1", "1", "10", "10", sep = "\t"),
  paste("P2", "-", "6", "11", "21", "11", "1", "1", "11", "11", sep = "\t")), tf)
rp <- omicsViewer:::read.proteinGroups.tmt(tf)
ok(ut_cmp_equal(nrow(rp$Reporter.intensity.corrected), 2),
   "read.proteinGroups.tmt: no filtered rows keeps ALL rows (R-L4)")
ok(ut_cmp_equal(nrow(rp$annot), 2),
   "read.proteinGroups.tmt: annotation keeps ALL rows (R-L4)")
