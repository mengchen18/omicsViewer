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

