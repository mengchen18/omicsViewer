library(omicsViewer)
library(unittest, quietly = TRUE)
library(Matrix)
library(fastmatch)

# =========================== conversion ========================
gs <- cbind(
  s1 = c(1, 1, 1, 1, 0, 0, 0, 0),
  s2 = c(0, 0, 0, 0, 1, 1, 1, 1),
  s3 = c(0, 0, 1, 1, 1, 1, 0, 0))
stats <- 1:8
rownames(gs) <- names(stats) <- paste0("g", 1:8)
gss <- as(gs, "dgCMatrix")
ok(ut_cmp_equal(omicsViewer:::csc2list(gs), omicsViewer:::csc2list(gss)), 
  "convert matrix to data.frame for fgsea")
df <- omicsViewer:::csc2list(gs)
ok(
  ut_cmp_equal(
    gss, omicsViewer:::list2csc(df, dimnames = dimnames(gs)), check.attributes = FALSE
    ), "convert data.frame to matrix"
   )

vectORA <- omicsViewer:::vectORA
r1 <- vectORA( gs, i = 1:4, minOverlap = 1, minSize = 1 )
r2 <- vectORA( gss, i = 1:4, minOverlap = 1, minSize = 1 )
vectORATall <- omicsViewer:::vectORATall
r3 <- vectORATall( df, i=rownames(gs)[1:4], background=8 )

ok(
  ut_cmp_equal(r1$pathway, c("s1", "s3"), check.attributes = FALSE),
  "vectORA 1"
  )
ok(ut_cmp_equal(r1, r2, check.attributes = FALSE), "vectORA 2")
ok(ut_cmp_equal(r1, r3), "vectORATall")



res <- omicsViewer:::fgsea1(
  gs = gss, stats = stats, minSize = 1, maxSize = 500, sampleSize = 3)
ut_cmp_identical(res$pathway, colnames(gss))
ok(ut_cmp_equal(res$pval[1] < 0.2, TRUE), "fgsea significant pathway 1")
ok(ut_cmp_equal(res$pval[2] < 0.2, TRUE), "fgsea significant pathway 2")
ok(ut_cmp_equal(res$pval[3] > 0.6, TRUE), "fgsea insignificant pathway")

totall <- omicsViewer:::totall
tgs <- totall(gs)
ok(ut_cmp_identical(colnames(tgs), c("featureId", "gsId", "weight")), 
  "totall 1")
ok(ut_cmp_equal(nrow(tgs), 12), "totall 2")


# ====================== trisetters ==========================

expr <- matrix(1:6, 3, 2)
rownames(expr) <- c("g1", "g2", "g3")
colnames(expr) <- c("s1", "s2")

m1 <- data.frame(
  "F1|pos2|pos3" = 1:3,
  "F2|pos2|pos3" = 1:3,
  check.names = FALSE
)
rownames(m1) <- rownames(expr)

m2 <- data.frame(
  "var1|pos2|pos3" = 1:2,
  "var2|pos2|pos3" = 1:2,
  check.names = FALSE
)
rownames(m2) <- colnames(expr)

trisetter <- omicsViewer:::trisetter

f1 <- rbind(c("F1", "pos2", "pos3"),
            c("F2", "pos2", "pos3"))
ok(
  ut_cmp_equal(
    trisetter(meta = m1, expr=NULL, combine="none"), f1,
    check.attributes = FALSE
  ),
  "trisetter feature meta"
)

var1 <- rbind(c("var1", "pos2", "pos3"),
              c("var2", "pos2", "pos3"))
ok(
  ut_cmp_equal(
    trisetter(meta = m2, expr=NULL, combine="none"), var1,
    check.attributes = FALSE
  ),
  "trisetter pheno meta"
)

cf1 <- rbind(f1, c("Sample", "Auto", "s1"),
             c("Sample", "Auto", "s2"))
ok(
  ut_cmp_equal(
    trisetter(meta = m1, expr=expr, combine="feature"), cf1,
    check.attributes = FALSE
  ),
  "trisetter feature combined with expr"
)
cvar1 <- rbind(var1, 
               c("Feature", "Auto", "g1"),
               c("Feature", "Auto", "g2"),
               c("Feature", "Auto", "g3"))
ok(
  ut_cmp_equal(
    trisetter(meta = m2, expr=expr, combine="pheno"), cvar1,
    check.attributes = FALSE
  ),
  "trisetter pheno combined with expr"
)

varSelector <- omicsViewer:::varSelector
l1 <- list(analysis = "Feature", subset= "Auto", variable = "g1")
ok(
ut_cmp_equal(
  varSelector(x = l1, expr = expr, meta = m2),
  c(1, 4),
  check.attributes = FALSE
  ),
"select from triselector - feature"
)
  
l1 <- list(analysis = "Sample", subset= "Auto", variable = "s1")
ok(
ut_cmp_equal(
  varSelector(x = l1, expr = expr, meta = m2), 
  1:3, 
  check.attributes = FALSE),
"select from triselector - sample"
)
text2num <- omicsViewer:::text2num
ok(ut_cmp_equal(text2num("-log10(0.01)"), -log10(0.01)), "text2num - 1")
ok(ut_cmp_equal(text2num("0.05"), 0.05), "text2num - 2")

######################################
terms <- data.frame(
  id = c("ID1", "ID2", "ID1", "ID2", "ID8", "ID10"),
 term = c("T1", "T1", "T2", "T2", "T2", "T2"),
  stringsAsFactors = FALSE
)
features <- list(c("ID1", "ID2"), c("ID13"), c("ID4", "ID8", "ID10"))
gsAnnotIdList(idList = features, gsIdMap = terms, minSize = 1, maxSize = 500)

terms <- data.frame(
id = c("ID1", "ID2", "ID1", "ID2", "ID8", "ID10", "ID4", "ID4"),
term = c("T1", "T1", "T2", "T2", "T2", "T2", "T1", "T2"),
stringsAsFactors = FALSE
)
features <- list(F1 = c("ID1", "ID2", "ID4"), F2 = c("ID13"), F3 = c("ID4", "ID8", "ID10"))

res <- data.frame(
  featureId = c(1, 1, 3, 3),
  gsId = c("T1", "T2", "T2", "T1"),
  weight = rep(1, 4),
  stringsAsFactors = FALSE
)
r1 <- gsAnnotIdList(features, gsIdMap = terms, data.frame = TRUE, minSize = 1)
ok(ut_cmp_equal(r1, res, check.attributes = FALSE), 
   "gsAnnotIdList data.frame"
   )

res <- sparseMatrix(i = c(1, 1, 3, 3), j = c(1, 2, 1, 2), x = 1)
colnames(res) <- c("T1", "T2")
r2 <- gsAnnotIdList(features, gsIdMap = terms, data.frame = FALSE, minSize = 1)
ok(ut_cmp_equal(r2, res), "gsAnnotIdList sparseMatrix")

# ==============================================
xq <- rbind(c(4, 2, 4),
            c(20, 40, 10),
            c(11, 234, 10),
            c(200, 1000, 100))

vectORA.core <- omicsViewer:::vectORA.core
r1 <- vectORA.core(xq[1, ], xq[2, ], xq[3, ], xq[4, ])
r2 <- vectORA.core(xq[1, ], xq[2, ], xq[3, ], xq[4, ], unconditional.or = FALSE)

ok(ut_cmp_equal(r1$p.value, r2$p.value), "vecORA.core 1")

# fisher's test
pv <- t(apply(xq, 2, function(x1) {
  m <- rbind(c(x1[1], x1[2]-x1[1]),
             c(x1[3]-x1[1], x1[4] - x1[2] - x1[3] + x1[1]))
  v <- fisher.test(m, alternative = "greater")
  c(p.value = v$p.value, v$estimate)
}))

ok(ut_cmp_equal(r1$p.value, pv[, "p.value"]), "vectORA.core - conditional OR")
ok(ut_cmp_equal(r2$OR, pv[, "odds ratio"]), "vectORA - unconditional OR")

