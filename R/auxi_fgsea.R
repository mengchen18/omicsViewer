#' Wrapper of fgsea function to take binary gene set matrix as input
#' @param gs a data.frame object with at least three columns named as "featureId", "gsId" and "weight".
#' @param stats other parameter passed to fgsea
#' @param nperm number of permutation
#' @param gs_desc description of gene sets, it should be a named vector and the names
#'  should be the same as colnames(gs)
#' @param ... other parameters passed to fgsea
#' @importFrom fgsea fgsea
#' @examples
#' ## not for users
#' # library(fgsea)
#' # library(Biobase)
#' # dat <- readRDS("Dat/exampleEset.RDS")
#' # fd <- fData(dat)
#' # fdgs <- fd[, grep("^GS\\|", colnames(fd))]
#' # res <- fgsea0(fdgs, stats = fd$`t-test|OV_BR|md`, nperm = 1000, minSize = 5, maxSize = 500)
#' # res <- fgsea0(
#' #   fdgs, stats = fd$`t-test|OV_BR|md`, nperm = 1000, 
#' #   minSize = 5, maxSize = 500, gs_desc = colnames(fdgs))
#' @return a \code{data.frame} of fgsea results

fgsea1 <- function(gs, stats, nperm = 100, gs_desc = NULL, ...) {
  
  pw <- split(gs$featureId, gs$gsId)
  
  params <- list(pathways = pw, stats = stats, nperm = nperm, ...)
  res <- do.call(fgsea, params)
  if (!is.null(gs_desc)) {
    res$desc <- as.character(gs_desc[res$pathway])
  }
  as.data.frame(res)
}
