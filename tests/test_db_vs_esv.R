library(unittest, quietly = TRUE)
library(RSQLite)
library(omicsViewer)


f <- system.file(package = 'omicsViewer', 'extdata/demo.RDS')
obj <- readRDS(f)
dd <- tools::R_user_dir("omicsViewer", which="cache")
dir.create(dd)
db <- tempfile(tmpdir = dd, fileext = ".db")
savedPath <- saveOmicsViewerDb(obj, db)
ok(ut_cmp_identical(is.character(savedPath), TRUE), "save sqlite database")

esv1 <- readESVObj(f)
ok(ut_cmp_identical(inherits(esv1, "ExpressionSet"), TRUE), "read expressionset")

esv2 <- readESVObj(db)
ok(ut_cmp_identical(inherits(esv2, "SQLiteConnection"), TRUE), "connect to database")

getExprs <- omicsViewer:::getExprs
expr1 <- getExprs(esv1)
expr2 <- getExprs(esv2)
ok(ut_cmp_identical(expr1, expr2), "db vs esv - expression")

getExprsImpute <- omicsViewer:::getExprsImpute
expr1 <- getExprsImpute(esv1)
expr2 <- getExprsImpute(esv2)
ok(ut_cmp_identical(expr1, expr2), "db vs esv - expression imputed")


getPData <- omicsViewer:::getPData
getFData <- omicsViewer:::getFData
getAx <- omicsViewer:::getAx
getDend <- omicsViewer:::getDend

pd1 <- getPData(esv1)
pd2 <- getPData(esv2)
ok(ut_cmp_identical(colnames(pd1), colnames(pd2)), "db vs esv - phenotype 1")
ok(ut_cmp_identical(rownames(pd1), rownames(pd2)), "db vs esv - phenotype 2")

fd1 <- getFData(esv1)
fd2 <- getFData(esv2)
ok(ut_cmp_identical(rownames(fd1), rownames(fd2)), "db vs esv - feature data")

for (i in c("sx", "sy", "fx", "fy")) {
  ax1 <- getAx(esv1, i)
  ax2 <- getAx(esv2, i)
  ok(ut_cmp_identical(ax1, ax2), sprintf("db vs esv - get axis - %s", i))
}

dd2 <- getDend(esv2)
ok(ut_cmp_identical(dd2, NULL), "db get dend")


# ---------------- M4: ESVObj = SummarizedExperiment ----------------
# Documented `omicsViewer(dir, ESVObj = se)` failed: the L0 path called
# tallGS directly on the SE (fData/pData don't exist) and the getters
# silently returned NULL for unsupported classes. The fix converts via
# asEsetWithAttr first and the getters raise an explicit error.
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(exprs = Biobase::exprs(esv1)),
  colData = S4Vectors::DataFrame(Biobase::pData(esv1)),
  rowData = S4Vectors::DataFrame(Biobase::fData(esv1)))
# gene-set matrix attribute must survive the conversion chain
rd <- SummarizedExperiment::rowData(se)
attr(rd, "GS") <- attr(Biobase::fData(esv1), "GS")
SummarizedExperiment::rowData(se) <- rd
esv3 <- omicsViewer:::tallGS(omicsViewer:::asEsetWithAttr(se))
ok(ut_cmp_identical(inherits(esv3, "ExpressionSet"), TRUE),
   "M4: ESVObj SummarizedExperiment converts to ExpressionSet")
ok(ut_cmp_identical(
   is.data.frame(attr(Biobase::fData(esv3), "GS")), TRUE),
   "M4: GS attribute survives the conversion chain")
ok(ut_cmp_identical(
   identical(getExprs(esv3), Biobase::exprs(obj)), TRUE),
   "M4: converted SE reads through getExprs")
ok(ut_cmp_identical(
   inherits(tryCatch(getExprs(se), error = function(e) e), "error"), TRUE),
   "M4: getExprs raises an explicit error on a raw SE")
