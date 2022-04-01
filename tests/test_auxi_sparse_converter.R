library(omicsViewer)
library(Matrix)
library(unittest, quietly = TRUE)

# source("R/auxi_sparse_converter.R")

m <- matrix(0, 10, 5, dimnames = list(paste0("R", 1:10), paste0("C", 1:5)))
for (i in 2:ncol(m)) {
   m[i*2, i] <- 1
}

csc <- as(m, "dgCMatrix")
ok(ut_cmp_equal(colSums(csc), colSums(m)), "as sparse")

scc <- omicsViewer:::csc2list(csc)
csc2 <- omicsViewer:::list2csc(scc, dimnames = list(rownames(m)))
ok(ut_cmp_equal(colSums(csc2), colSums(m)[-1], check.attributes = FALSE), 
   "only rownames, missing 0 column")

csc3 <- omicsViewer:::list2csc(scc, dimnames = list(rownames(m), colnames(m)))
ok(ut_cmp_equal(colSums(csc3), colSums(m), check.attributes = FALSE), 
   "both rownames and columns given")
ok(ut_cmp_equal(as(csc3, "matrix"), m, check.attributes = FALSE), 
   "regular matrix")

scc2 <- scc
scc2$weight <- NULL
csc4 <- omicsViewer:::list2csc(scc2, dimnames = list(rownames(m), colnames(m)))
ok(ut_cmp_equal(colSums(csc4), colSums(m), check.attributes = FALSE), 
   "as binary matrix")


# ================= conversion of hclust object =================
m <- matrix(rnorm(50), 25, dimnames = list(paste0("g", 1:25), c('c1', 'c2')))
hc <- hclust(dist(m))
plot(hc)
te <- omicsViewer:::hclust2str(hc)
hc2 <- omicsViewer:::str2hclust(te)
plot(hc2)
hc_elements <- names(hc)
ok(ut_cmp_equal(all(names(hc2) %in% hc_elements), TRUE), "Element name test")
j <- sapply(setdiff(hc_elements, "call"), function(x) 
  all.equal(hc[[x]], hc2[[x]], tol = 5e-5, check.attributes = FALSE))
ok(ut_cmp_equal(all(j), TRUE), "HCL element conversion test.")


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


res <- omicsViewer:::fgsea1(
  gs = gss, stats = stats, minSize = 1, maxSize = 500, sampleSize = 3)
ut_cmp_identical(res$pathway, colnames(gss))
ok(ut_cmp_equal(res$pval[1] < 0.2, TRUE), "fgsea significant pathway 1")
ok(ut_cmp_equal(res$pval[2] < 0.2, TRUE), "fgsea significant pathway 2")
ok(ut_cmp_equal(res$pval[3] > 0.6, TRUE), "fgsea insignificant pathway")

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

# ===================== select highlight block ===================
cl <- list(x = -5:5, y = c(5:0, 1:5))
line_rect <- omicsViewer:::line_rect
r <- line_rect(coord = cl, l = list(x = 2.5, y = 2.5, corner = "topright"))
ok(
  ut_cmp_equal(
    r$rect[[1]][c("x0", "y0")], c(2.5, 2.5), check.attributes = FALSE
    ), "line_rect - topright"
)

r <- line_rect(coord = cl, l = list(x = 2.5, y = 2.5, corner = "volcano"))
ok(
  ut_cmp_identical(length(r$rect), as.integer(2)),
  "line_rect - volcano n block"
)
ok(
  ut_cmp_equal(
    r$rect[[1]][c("x1", "y0")], c(-2.5, 2.5), check.attributes = FALSE
    ), "line_rect - volcano block 1"
)
ok(
  ut_cmp_equal(
    r$rect[[2]][c("x0", "y0")], c(2.5, 2.5), check.attributes = FALSE
    ), "line_rect - volcano block 2"
)

ok(
  ut_cmp_identical(
    line_rect(coord = cl, l = list(y = 2.5, corner = "volcano")), NULL
    ), "line_rect - NULL 1"
)

ok(
  ut_cmp_identical(
    line_rect(coord = cl, l = list(x = 2.5, corner = "volcano")), NULL
  ), "line_rect - NULL 2"
)

ok(
  ut_cmp_identical(
    line_rect(coord = cl, l = list(x = 2.5, y = 2.5, corner = "None")), NULL
  ), "line_rect - NULL 3"
)

