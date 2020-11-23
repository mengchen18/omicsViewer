# library(fgsea)
# library(Biobase)
# dat <- readRDS("Dat/exampleEset.RDS")
# fd <- fData(dat)
# fdgs <- fd[, grep("^GS\\|", colnames(fd))]
# res <- fgsea0(fdgs, stats = fd$`t-test|OV_BR|md`, nperm = 1000, minSize = 5, maxSize = 500)
# res <- fgsea0(fdgs, stats = fd$`t-test|OV_BR|md`, nperm = 1000, minSize = 5, maxSize = 500, gs_desc = colnames(fdgs))

#' @param gs a matrix where rows are features and columns are gene sets. The zeros in the matrix 
#'   indicate a gene is not a memeber of the gene sets in the column. 
#' @param stats other parameter passed to fgsea
#' @param ... other parameters passed to fgsea
#' 
fgsea0 <- function(gs, stats, nperm = 100, gs_desc = NULL, feature_desc = NULL, ...) {
  
  rn <- rownames(gs)
  pw <- lapply(1:ncol(gs), function(i) {
    rn[gs[, i] != 0]
  })
  names(pw) <- colnames(gs)
  
  if (is.null(names(stats)))
    names(stats) <- rownames(gs)
  
  params <- list(pathways = pw, stats = stats, nperm = nperm, ...)
  res <- do.call(fgsea, params)
  if (!is.null(gs_desc)) {
    if (length(gs_desc) != ncol(gs))
      stop("Length of gs_desc != ncol(gs) !!")
    res$desc <- gs_desc[match(res$pathway, colnames(gs))]
  }
    
  as.data.frame(res)
}

