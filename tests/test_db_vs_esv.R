library(unittest, quietly = TRUE)
library(RSQLite)
library(DBI)
library(shiny)
library(omicsViewer)


f <- system.file(package = 'omicsViewer', 'extdata/demo.RDS')
obj <- readRDS(f)
dd <- tools::R_user_dir("omicsViewer", which="cache")
dir.create(dd)
db <- tempfile(tmpdir = dd, fileext = ".db")
savedPath <- saveOmicsViewerDb(obj, db)
ok(ut_cmp_identical(is.character(savedPath), TRUE), "save sqlite database")
ok(ut_cmp_identical("meta" %in% RSQLite::dbListTables(RSQLite::dbConnect(SQLite(), db)), TRUE),
   "save writes the meta schema marker")

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

# =====================================================================
# sqlite-revival WP-0/1/2: content-level round-trip parity. The database
# must serve the SAME app-ready view as the .RDS load path: gene sets
# attached (not an empty GS table + 568 wide GS|* columns polluting the
# feature space), content-identical annotations.
# =====================================================================

ok(ut_cmp_identical(ncol(fd2), ncol(fd1)), "feature table is the app-ready view (no wide GS|* columns)")
ok(ut_cmp_identical(length(grep("^GS\\|", colnames(fd2))), 0L),
   "no GS|* membership columns leak into the db feature table")

g1 <- attr(fd1, "GS")
g2 <- attr(fd2, "GS")
ok(ut_cmp_identical(is.null(g1), FALSE), "RDS path carries the gene-set annotation")
ok(ut_cmp_identical(is.null(g2), FALSE), "db path carries the gene-set annotation (was lost)")
ok(ut_cmp_identical(dim(g1), dim(g2)), "db vs esv - GS dimensions")
ok(ut_cmp_identical(as.character(g1$featureId), as.character(g2$featureId)), "db vs esv - GS featureId")
ok(ut_cmp_identical(as.character(g1$gsId), as.character(g2$gsId)), "db vs esv - GS gsId")
ok(ut_cmp_identical(as.numeric(g1$weight), as.numeric(g2$weight)), "db vs esv - GS weight")
ok(ut_cmp_identical(setequal(levels(g1$gsId), levels(g2$gsId)), TRUE), "db vs esv - GS gsId level set")
# melt's permuted row names carry no semantics; compare with them
# normalized away on both sides
strip_gs_rownames <- function(g) { if (!is.null(g)) rownames(g) <- NULL; g }
ok(ut_cmp_identical(strip_gs_rownames(g1), strip_gs_rownames(g2)),
   "db vs esv - GS annotation identical (factors, levels, weights)")

# NaN normalizes to NA at save (SQLite cannot represent NaN); GS row
# names (melt artifact) normalized away on both sides
norm_fd <- function(fd) {
  g <- attr(fd, "GS")
  if (!is.null(g)) { rownames(g) <- NULL; attr(fd, "GS") <- g }
  fd
}
nan2na_df <- function(df) {
  for (i in which(vapply(df, is.numeric, logical(1)))) {
    v <- df[[i]]
    if (any(idx <- is.nan(v))) { v[idx] <- NA; df[[i]] <- v }
  }
  df
}
ok(ut_cmp_identical(nan2na_df(norm_fd(fd1)), norm_fd(fd2)), "db vs esv - feature data content identical")
ok(ut_cmp_identical(nan2na_df(pd1), pd2), "db vs esv - phenotype data content identical")

# =====================================================================
# sqlite-revival WP-5: type fidelity - factors rebuilt from recorded
# levels, logicals and Dates round trip
# =====================================================================
ok(ut_cmp_identical(
  vapply(pd1, function(x) class(x)[1], character(1)),
  vapply(pd2, function(x) class(x)[1], character(1))),
  "db vs esv - phenotype column classes identical")

obj_tf <- readRDS(f)
pd_tf <- Biobase::pData(obj_tf)
qv <- attr(pd_tf, "quickViews")
pd_tf$tf_logical <- rep(c(TRUE, FALSE), 30)
pd_tf$tf_date <- as.Date("2020-01-01") + (0:59)
pd_tf$tf_factor <- factor(rep(c("low", "mid", "high"), 20), levels = c("low", "mid", "high"))
Biobase::pData(obj_tf) <- pd_tf
attr(Biobase::pData(obj_tf), "quickViews") <- qv
db_tf <- tempfile(tmpdir = dd, fileext = ".db")
saveOmicsViewerDb(obj_tf, db_tf)
pd_tfr <- getPData(readESVObj(db_tf))
ok(ut_cmp_identical(class(pd_tfr$tf_logical)[1], "logical"), "logical column round trips")
ok(ut_cmp_identical(class(pd_tfr$tf_date)[1], "Date"), "Date column round trips")
ok(ut_cmp_identical(class(pd_tfr$tf_factor)[1], "factor"), "factor column round trips")
ok(ut_cmp_identical(levels(pd_tfr$tf_factor), c("low", "mid", "high")),
   "factor levels round trip in order (incl. non-alphabetical)")
ok(ut_cmp_identical(pd_tfr$tf_factor, pd_tf$tf_factor), "factor values round trip")

# =====================================================================
# sqlite-revival WP-3: attribute round-trip through the attrs table
# =====================================================================
ok(ut_cmp_identical(attr(fd2, "quickViews"), attr(fd1, "quickViews")),
   "user-defined feature quickViews round trip")
ok(ut_cmp_identical(attr(pd2, "quickViews"), attr(pd1, "quickViews")),
   "user-defined sample quickViews round trip")

set.seed(42)
obj_at <- readRDS(f)
# precomputed expression dendrogram rides on the exprs matrix (the
# heatmap HCL choices surface reads exactly this attribute). NOTE:
# assayDataElement<- strips matrix attributes, so the assayData slot is
# rebuilt via assayDataNew (the same way prepOmicsViewer stores attrs)
ex_at <- Biobase::exprs(obj_at)
hc <- hclust(dist(t(matrix(rnorm(2702 * 4), 2702))))
attr(ex_at, "rowDendrogram") <- list(ward = list(ord = hc$order, hcl = as.dendrogram(hc)))
Biobase::assayData(obj_at) <- Biobase::assayDataNew("exprs" = ex_at)
attr(obj_at, "cormat") <- matrix(runif(3600), 60,
  dimnames = list(colnames(ex_at), colnames(ex_at)))
attr(obj_at, "S6.6_drc") <- c("drc1" = "AUC|dose", "drc2" = "Hill|dose")
db_at <- tempfile(tmpdir = dd, fileext = ".db")
saveOmicsViewerDb(obj_at, db_at)
con_at <- readESVObj(db_at)
ok(ut_cmp_identical(getAx(con_at, "cormat"), attr(obj_at, "cormat")),
   "cormat round trips through the getAx seam")
ok(ut_cmp_identical(getAx(con_at, "S6.6_drc"), attr(obj_at, "S6.6_drc")),
   "S6.6_drc round trips through the getAx seam")
ok(ut_cmp_identical(
  attr(getExprs(con_at), "rowDendrogram"), attr(ex_at, "rowDendrogram")),
  "exprs rowDendrogram round trips onto the expression matrix")
# the ExpressionSet side of the seam still serves plain attributes
ok(ut_cmp_identical(getAx(esv1, "S6.6_drc"), attr(esv1, "S6.6_drc")),
   "getAx ExpressionSet method serves arbitrary attributes")

# =====================================================================
# sqlite-revival WP-1: stale exprsimpute must not survive a re-save
# without imputed values
# =====================================================================
db_st <- tempfile(tmpdir = dd, fileext = ".db")
obj_wi <- readRDS(f)
set.seed(1)
Biobase::assayData(obj_wi) <- Biobase::assayDataNew(
  "exprs" = Biobase::exprs(obj_wi),
  "exprs_impute" = jitter(Biobase::exprs(obj_wi)))
saveOmicsViewerDb(obj_wi, db_st)
ok(ut_cmp_identical(is.null(getExprsImpute(readESVObj(db_st))), FALSE),
   "imputed matrix written when present")
obj_ni <- readRDS(f)
Biobase::assayData(obj_ni) <- Biobase::assayDataNew("exprs" = Biobase::exprs(obj_ni))
saveOmicsViewerDb(obj_ni, db_st)
ok(ut_cmp_identical(is.null(getExprsImpute(readESVObj(db_st))), TRUE),
   "stale exprsimpute table dropped on re-save without imputed values")

# =====================================================================
# sqlite-revival WP-1: guards, overwrite=FALSE, failure atomicity
# =====================================================================
db_g <- tempfile(tmpdir = dd, fileext = ".db")
saveOmicsViewerDb(obj, db_g)

obj_rn <- readRDS(f)
pd_rn <- Biobase::pData(obj_rn)
pd_rn$rowname <- "x"
Biobase::pData(obj_rn) <- pd_rn
err <- tryCatch(saveOmicsViewerDb(obj_rn, db_g), error = function(e) e)
ok(ut_cmp_identical(inherits(err, "error"), TRUE), "rowname column rejected")
ok(ut_cmp_identical(grepl("rowname", conditionMessage(err)), TRUE),
   "rowname rejection names the column")

# duplicated column names: the container setters reject them long before
# the database layer, so the guard is probed directly (its wiring into
# saveOmicsViewerDb is covered by the rowname/unsupported-class cases)
df_dup <- data.frame(a = 1:3, b = 4:6, check.names = FALSE)
colnames(df_dup)[2] <- "a"
err <- tryCatch(omicsViewer:::`.db_guard_table`(df_dup, "pData"), error = function(e) e)
ok(ut_cmp_identical(inherits(err, "error"), TRUE), "duplicate column names rejected")
ok(ut_cmp_identical(grepl("duplicated column", conditionMessage(err)), TRUE),
   "duplicate rejection names the problem")

obj_lc <- readRDS(f)
pd_lc <- Biobase::pData(obj_lc)
pd_lc$lst <- I(as.list(rnorm(60)))
Biobase::pData(obj_lc) <- pd_lc
err <- tryCatch(saveOmicsViewerDb(obj_lc, db_g), error = function(e) e)
ok(ut_cmp_identical(inherits(err, "error"), TRUE), "unsupported column class rejected")
ok(ut_cmp_identical(grepl("unsupported", conditionMessage(err)), TRUE),
   "unsupported class rejection names the class")

# guards fire before the connection: an existing database is untouched
ok(ut_cmp_identical(ncol(getFData(readESVObj(db_g))), ncol(fd1)),
   "guard failures leave an existing database unchanged")

err <- tryCatch(saveOmicsViewerDb(obj, db_g, overwrite = FALSE), error = function(e) e)
ok(ut_cmp_identical(inherits(err, "error"), TRUE), "overwrite=FALSE refuses existing tables")
ok(ut_cmp_identical(ncol(getFData(readESVObj(db_g))), ncol(fd1)),
   "refused overwrite leaves the database intact")

# =====================================================================
# sqlite-revival WP-2: legacy (pre-2.1.81) database recovery - wide
# GS|* columns in the feature table with an empty GS table
# =====================================================================
ldb <- tempfile(tmpdir = dd, fileext = ".db")
lcon <- dbConnect(RSQLite::SQLite(), ldb)
fdw <- Biobase::fData(obj)
invisible(dbWriteTable(lcon, "exprs",
  data.frame(Biobase::exprs(obj), rowname = rownames(Biobase::exprs(obj)), check.names = FALSE)))
invisible(dbWriteTable(lcon, "feature",
  data.frame(fdw, rowname = rownames(fdw), check.names = FALSE)))
invisible(dbWriteTable(lcon, "sample",
  data.frame(Biobase::pData(obj), rowname = rownames(Biobase::pData(obj)), check.names = FALSE)))
invisible(dbWriteTable(lcon, "GS",
  data.frame(featureId = character(0), gsId = character(0), weight = numeric(0))))
invisible(dbWriteTable(lcon, "axes", data.frame(axis = character(0), value = character(0))))
invisible(dbWriteTable(lcon, "dendrogram",
  data.frame(name = character(0), object = character(0), dim = integer(0))))
dbDisconnect(lcon)

lread <- readESVObj(ldb)
fdl <- getFData(lread)
ok(ut_cmp_identical(ncol(fdl), ncol(fd1)),
   "legacy db: wide GS|* columns stripped from the feature table")
gl <- attr(fdl, "GS")
ok(ut_cmp_identical(is.null(gl), FALSE), "legacy db: gene sets recovered from wide columns")
ok(ut_cmp_identical(dim(gl), dim(g1)), "legacy db: recovered GS dimensions match")
ok(ut_cmp_identical(as.character(gl$gsId), as.character(g1$gsId)),
   "legacy db: recovered GS content matches")

# =====================================================================
# sqlite-revival WP-4: readESVObj validation and read-only connections
# =====================================================================
miss <- tempfile(tmpdir = dd, fileext = ".db")
err <- tryCatch(readESVObj(miss), error = function(e) e)
ok(ut_cmp_identical(inherits(err, "error"), TRUE), "missing database file rejected")
ok(ut_cmp_identical(file.exists(miss), FALSE),
   "missing database file is not silently created")

junk <- tempfile(tmpdir = dd, fileext = ".sqlite")
writeLines("this is not a database", junk)
err <- tryCatch(readESVObj(junk), error = function(e) e)
ok(ut_cmp_identical(inherits(err, "error"), TRUE), "non-sqlite file rejected")
ok(ut_cmp_identical(grepl("not a valid SQLite database", conditionMessage(err)), TRUE),
   "non-sqlite rejection explains the problem")

empty <- tempfile(tmpdir = dd, fileext = ".db")
econ <- dbConnect(RSQLite::SQLite(), empty)
dbDisconnect(econ)
err <- tryCatch(readESVObj(empty), error = function(e) e)
ok(ut_cmp_identical(inherits(err, "error"), TRUE), "database without omicsViewer tables rejected")
ok(ut_cmp_identical(grepl("missing required table", conditionMessage(err)), TRUE),
   "missing-table rejection names the tables")

ok(ut_cmp_identical(isTRUE(DBI::dbIsValid(esv2)), TRUE), "dataset connection is valid")
ro_err <- tryCatch(dbExecute(esv2, "CREATE TABLE nope (x)"), error = function(e) e)
ok(ut_cmp_identical(inherits(ro_err, "error"), TRUE), "dataset connection is read-only")

con_close <- readESVObj(db)
omicsViewer:::esv_close_connection(con_close)
ok(ut_cmp_identical(isTRUE(DBI::dbIsValid(con_close)), FALSE), "esv_close_connection disconnects")
omicsViewer:::esv_close_connection(con_close)  # double close is quiet
omicsViewer:::esv_close_connection(NULL)       # NULL is quiet

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

# =====================================================================
# sqlite-revival WP-6: app-level boot on a database dataset and dataset
# switch with connection cleanup, through the real app_module
# =====================================================================
if (identical(Sys.getenv("NOT_CRAN"), "true")) {

dbdir <- file.path(tempdir(), paste0("db-boot-", Sys.getpid()))
if (dir.exists(dbdir)) unlink(dbdir, recursive = TRUE)
dir.create(dbdir, recursive = TRUE, showWarnings = FALSE)
file.copy(f, file.path(dbdir, "demo.RDS"))
saveOmicsViewerDb(readRDS(f), file.path(dbdir, "demo.db"))

.db_recorded <- list()
app_db <- function(input, output, session) {
  omicsViewer:::app_module(
    "app",
    .dir = shiny::reactive(dbdir),
    esetLoader = function(path) {
      v <- readESVObj(path)
      if (inherits(v, "SQLiteConnection"))
        .db_recorded[[length(.db_recorded) + 1L]] <<- v
      v
    })
}

shiny::testServer(app_db, {
  # the navbar input must exist before the first flush (test_appState note)
  session$setInputs(`app-dataspace-eset` = "Feature")
  session$flushReact()
  session$setInputs(`app-selectFile` = "demo.db")
  for (i in 1:10) session$flushReact()
  session$flushOutput()

  sm <- as.character(output[["app-summary"]])
  ok(ut_cmp_identical(any(grepl("features and 60 samples", sm, fixed = TRUE)), TRUE),
     "app boots on a database dataset (2702 features, 60 samples)")
  ok(ut_cmp_identical(
    any(grepl("Dataset loaded successfully", as.character(output[["app-loadingStatus"]]), fixed = TRUE)),
    TRUE), "app reports the database dataset as loaded")
  ok(ut_cmp_identical(length(.db_recorded) == 1L, TRUE),
     "exactly one dataset connection was opened")

  # dataset switch db -> RDS: the connection must be closed, the new
  # dataset must load
  session$setInputs(`app-selectFile` = "demo.RDS")
  for (i in 1:10) session$flushReact()
  session$flushOutput()
  sm <- as.character(output[["app-summary"]])
  ok(ut_cmp_identical(any(grepl("features and 60 samples", sm, fixed = TRUE)), TRUE),
     "dataset switch db -> RDS loads the new dataset")
  ok(ut_cmp_identical(
    all(!vapply(.db_recorded, DBI::dbIsValid, logical(1))), TRUE),
    "database connection closed on dataset switch (no handle leak)")

  # and back onto the database
  n0 <- length(.db_recorded)
  session$setInputs(`app-selectFile` = "demo.db")
  for (i in 1:10) session$flushReact()
  session$flushOutput()
  ok(ut_cmp_identical(any(grepl("features and 60 samples",
    as.character(output[["app-summary"]]), fixed = TRUE)), TRUE),
    "dataset switch RDS -> db loads the database dataset again")
  ok(ut_cmp_identical(
    all(!vapply(.db_recorded[seq_len(n0)], DBI::dbIsValid, logical(1))), TRUE),
    "all replaced connections stay closed")
})

}
