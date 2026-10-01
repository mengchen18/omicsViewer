library(unittest, quietly = TRUE)
library(methods)
library(omicsViewer)
library(SummarizedExperiment)

f <- system.file(package = 'omicsViewer', 'extdata/demo.RDS')
obj <- readRDS(f)

err <- function(expr) tryCatch(expr, error = function(e) conditionMessage(e))

# ---------------- dispatch parity on ExpressionSet ----------------
# every method body is a verbatim port; results must stay identical

ok(ut_cmp_identical(getExprs(obj), Biobase::exprs(obj)),
   "S4: getExprs ExpressionSet dispatch identical to exprs()")
ok(ut_cmp_identical(getPData(obj), Biobase::pData(obj)),
   "S4: getPData ExpressionSet dispatch identical to pData()")
ok(ut_cmp_identical(getFData(obj), Biobase::fData(obj)),
   "S4: getFData ExpressionSet dispatch identical to fData()")
ok(ut_cmp_identical(omicsViewer:::getExprsImpute(obj), omicsViewer:::exprsImpute(obj)),
   "S4: getExprsImpute ExpressionSet dispatch identical to exprsImpute()")
ok(ut_cmp_identical(getAx(obj, "sx"), attr(obj, "sx")),
   "S4: getAx ExpressionSet dispatch reads attribute")
ok(ut_cmp_identical(getAx(obj, "nosuchaxis"), NULL),
   "S4: getAx unset axis returns NULL")

# ---------------- the six generics are exported S4 generics ----------------
for (g in c("asEsetWithAttr", "getExprs", "getExprsImpute", "getPData", "getFData", "getAx")) {
  ok(ut_cmp_identical(
     is(getExportedValue("omicsViewer", g), "standardGeneric"), TRUE),
     sprintf("S4: %s is an exported S4 generic", g))
}

# ---------------- ANY fallbacks keep the historical errors ----------------
for (entry in list(
     list(getExprs, "getExprs"), list(omicsViewer:::getExprsImpute, "getExprsImpute"),
     list(getPData, "getPData"), list(getFData, "getFData"))) {
  m <- err(entry[[1]](42))
  ok(ut_cmp_identical(grepl(paste0(entry[[2]], ": unsupported class"), m), TRUE),
     sprintf("S4: %s ANY method raises explicit error", entry[[2]]))
}
ok(ut_cmp_identical(grepl("getAx: unsupported class", err(getAx(42, "sx"))), TRUE),
   "S4: getAx ANY method raises explicit error")
ok(ut_cmp_identical(grepl("x should be either", err(asEsetWithAttr(42))), TRUE),
   "S4: asEsetWithAttr ANY method raises explicit error")

# ---------------- asEsetWithAttr conversion ----------------
ok(ut_cmp_identical(asEsetWithAttr(obj), obj),
   "S4: asEsetWithAttr identity on ExpressionSet")

se <- SummarizedExperiment(
  assays = list(exprs = Biobase::exprs(obj)),
  colData = S4Vectors::DataFrame(Biobase::pData(obj)),
  rowData = S4Vectors::DataFrame(Biobase::fData(obj)))
e2 <- asEsetWithAttr(se)
ok(ut_cmp_identical(inherits(e2, "ExpressionSet"), TRUE),
   "S4: asEsetWithAttr converts SummarizedExperiment")
ok(ut_cmp_identical(getExprs(e2), Biobase::exprs(obj)),
   "S4: converted SE reads identically through getExprs")
ok(ut_cmp_identical(inherits(tryCatch(getExprs(se), error = function(e) e), "error"), TRUE),
   "S4: raw SE still has no getExprs method (must convert first)")

# ---------------- tallGS ----------------
ok(ut_cmp_identical(omicsViewer:::tallGS(42), 42),
   "S4: tallGS ANY identity default")
ok(ut_cmp_identical(omicsViewer:::tallGS(letters), letters),
   "S4: tallGS ANY identity on non-Eset container")
t1 <- omicsViewer:::tallGS(obj)
ok(ut_cmp_identical(dim(Biobase::fData(t1)), c(2702L, 145L)),
   "S4: tallGS ExpressionSet method collapses wide GS columns (713 -> 145)")
ok(ut_cmp_identical(
   identical(names(attr(Biobase::fData(t1), "GS")), c("featureId", "gsId", "weight")), TRUE),
   "S4: tallGS ExpressionSet method produces long GS data.frame")
t2 <- omicsViewer:::tallGS(t1)
ok(ut_cmp_identical(dim(Biobase::fData(t2)), dim(Biobase::fData(t1))),
   "S4: tallGS is idempotent on its own output")
ok(ut_cmp_identical(
   identical(attr(Biobase::fData(t2), "GS"), attr(Biobase::fData(t1), "GS")), TRUE),
   "S4: tallGS output GS attribute stable")

# ---------------- third-party class simulation ----------------
# The extension pattern the conversion enables: a package defines its own
# class (here: an SE subclass) and registers methods for omicsViewer's
# generics - omicsViewer itself never references the class.
setClass("DummyThirdPartySE", contains = "SummarizedExperiment")
setMethod("asEsetWithAttr", "DummyThirdPartySE", function(x) x)
setMethod("getExprs", "DummyThirdPartySE", function(x) SummarizedExperiment::assay(x))
setMethod("getFData", "DummyThirdPartySE", function(x) as.data.frame(SummarizedExperiment::rowData(x)))
setMethod("getPData", "DummyThirdPartySE", function(x) as.data.frame(SummarizedExperiment::colData(x)))
setMethod("getAx", signature(x = "DummyThirdPartySE", what = "ANY"),
          function(x, what) attr(x, what))

d <- new("DummyThirdPartySE", se)
attr(d, "sx") <- "sample_group"

ok(ut_cmp_identical(asEsetWithAttr(d), d),
   "S4: third-party identity method wins over SE coercion")
ok(ut_cmp_identical(getExprs(d), SummarizedExperiment::assay(se)),
   "S4: third-party getExprs dispatch")
ok(ut_cmp_identical(getFData(d), as.data.frame(SummarizedExperiment::rowData(se))),
   "S4: third-party getFData dispatch")
ok(ut_cmp_identical(getPData(d), as.data.frame(SummarizedExperiment::colData(se))),
   "S4: third-party getPData dispatch")
ok(ut_cmp_identical(getAx(d, "sx"), "sample_group"),
   "S4: third-party getAx dispatch")
ok(ut_cmp_identical(omicsViewer:::tallGS(d), d),
   "S4: tallGS identity covers third-party class")
# the exact reactive_eset() composition in L0_module_app passes through
ok(ut_cmp_identical(omicsViewer:::tallGS(asEsetWithAttr(d)), d),
   "S4: L0 hot path composition (tallGS(asEsetWithAttr(x))) passes third-party objects through")

# ---------------- iheatmap extraction goes through the generics ----------------
app <- omicsViewer:::iheatmap(obj)
ok(ut_cmp_identical(inherits(app, "shiny.appobj"), TRUE),
   "S4: iheatmap extracts ExpressionSet through the generics")
mm <- matrix(rnorm(30), 6)
app2 <- omicsViewer:::iheatmap(mm, fData = data.frame(a = 1:5), pData = data.frame(b = 1:6))
ok(ut_cmp_identical(inherits(app2, "shiny.appobj"), TRUE),
   "S4: iheatmap matrix path unchanged")
ok(ut_cmp_identical(
   grepl("unsupported class", err(omicsViewer:::iheatmap(data.frame(a = 1:3)))), TRUE),
   "S4: iheatmap non-matrix without methods raises the generic's error")
