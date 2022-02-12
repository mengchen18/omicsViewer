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






