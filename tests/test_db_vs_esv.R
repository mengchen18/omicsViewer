library(unittest, quietly = TRUE)
library(RSQLite)
library(omicsViewer)


f <- system.file(package = 'omicsViewer', 'extdata/demo.RDS')
obj <- readRDS(f)
dd <- tools::R_user_dir("omicsViewer", which="cache")
dir.create(dd)
db <- file.path(dd, "temp.db")
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

