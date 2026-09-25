# Mock up SeqLogo (PTM motif) and ResponseCurve (dose-response) content
# in the demo ExpressionSet so the corresponding right-panel tabs can be
# exercised and debugged without a dedicated PTM / dose-response dataset.
#
# Run from the package root:
#   Rscript inst/script/mock_demo_tabs.R
#
# Requirements implemented here (see R/module_PTMotif.R, R/module_doseResponse.R,
# R/proc_doseCurve.R):
#  - SeqLogo tab: fData column `SeqLogo|<category>|<subcategory>` holding
#    peptide sequence windows; `;`-separated multi-windows allowed; every
#    unique sequence must have the SAME length (odd, e.g. 15 aa, so the
#    dashed center line marks the modification site); only the 20 standard
#    uppercase amino acids are counted by aaFreq().
#  - Response tab: fData columns `ResponseCurve|<curveid>|<param>` (canonical
#    source: drmMat() + extractParamDCList()), pData columns
#    `General|All|<dose_col>` + `General|All|<curveid_col>`, referenced via
#    attr(eset, "S6.6_drc") = c(dose_col = ..., curveid_col = ...).
#    Curve ids in the fData column names must match the pdata curve values.

suppressMessages({
  library(Biobase)
  library(drc)
})

set.seed(42)

f_demo <- file.path("inst", "extdata", "demo.RDS")
eset <- readRDS(f_demo)

fd <- fData(eset)
pd <- pData(eset)
fdAttrs <- attributes(fd)[setdiff(names(attributes(fd)),
  c("names", "class", "row.names"))]
pdAttrs <- attributes(pd)[setdiff(names(attributes(pd)),
  c("names", "class", "row.names"))]

## ------------------------------------------------------------------
## 1. SeqLogo mock: 15-aa windows, SP-proline-directed motif enriched
##    among features up-regulated in RE_vs_ME (so selecting those in the
##    Feature table yields a visible foreground motif enrichment).
## ------------------------------------------------------------------
AA <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
        "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

# proline-directed kinase window: ...RxxSP.... (S at center position 8)
motifWindow <- function() {
  s <- sample(AA, 15, replace = TRUE)
  s[5] <- "R"
  s[8] <- "S"
  s[9] <- "P"
  paste(s, collapse = "")
}
randomWindow <- function() paste(sample(AA, 15, replace = TRUE), collapse = "")

enriched <- which(fd[["ttest|RE_vs_ME|fdr"]] < 0.05 &
                  fd[["ttest|RE_vs_ME|mean.diff"]] > 0)

wins <- character(nrow(fd))
for (i in seq_len(nrow(fd))) {
  if (i %in% enriched) {
    # enriched features: motif window (some carry a second, random window
    # to exercise the `;`-separated multi-window path; a few are empty to
    # exercise the zero-length filter in cleanSeqs())
    w <- motifWindow()
    if (i %% 5L == 0L) w <- paste(w, randomWindow(), sep = ";")
    if (i %% 37L == 0L) w <- ""
    wins[i] <- w
  } else {
    wins[i] <- randomWindow()
  }
}
stopifnot(all(nchar(gsub(";", "", wins)) %% 15 == 0 | wins == ""))
fd[["SeqLogo|All|seq.window"]] <- wins

## ------------------------------------------------------------------
## 2. ResponseCurve mock: 4 curves x 15 samples (5 log-spaced doses,
##    3 replicates each); LL.4 fits via the canonical drmMat pipeline.
##    Fits that do not converge become NA parameter rows (valid case:
##    empty curve, working UI).
## ------------------------------------------------------------------
curves <- paste0("curve", 1:4)
doses <- c(0.1, 0.4, 1.6, 6.4, 25.6)
curveCol <- rep(curves, each = 15)
doseCol <- rep(rep(doses, each = 3), 4)

pd[["General|All|Dose"]] <- doseCol
pd[["General|All|Curve"]] <- curveCol
attr(eset, "S6.6_drc") <- c(dose_col = "Dose", curveid_col = "Curve")

mods <- omicsViewer:::drmMat(
  exprs(eset), fitvar = doseCol, fitvar.name = "Dose",
  curveid = curveCol, fct.name = "LL.4()")
ok <- sapply(mods, inherits, "drc")
message("dose-response fits: ", sum(ok), "/", length(ok), " converged")
drcPar <- omicsViewer:::extractParamDCList(mods)
fd <- cbind(fd, drcPar)

## ------------------------------------------------------------------
## re-assemble, preserving dropped data.frame attributes (quickViews)
## ------------------------------------------------------------------
for (a in names(fdAttrs)) attr(fd, a) <- fdAttrs[[a]]
for (a in names(pdAttrs)) attr(pd, a) <- pdAttrs[[a]]
fData(eset) <- fd
pData(eset) <- pd

stopifnot(
  any(grepl("^SeqLogo\\|", colnames(fData(eset)))),
  any(grepl("^ResponseCurve\\|", colnames(fData(eset)))),
  identical(attr(eset, "S6.6_drc")["dose_col"], c(dose_col = "Dose")),
  identical(attr(eset, "S6.6_drc")["curveid_col"], c(curveid_col = "Curve")),
  all(nchar(unique(unlist(strsplit(
    fData(eset)[["SeqLogo|All|seq.window"]], ";")))) == 15)
)

saveRDS(eset, f_demo)
message("updated ", f_demo, ": +", 1L, " SeqLogo column, +",
        ncol(drcPar), " ResponseCurve columns, +2 pData columns")
