#' Wrapper of fgseaMultilevel function to take binary gene set matrix as input
#' @param gs either a data.frame or a (sparse) matrix input. If a data.frame object 
#'   is given, it should have at least three columns named as "featureId", "gsId" and "weight".
#'   If a matrix is given, the matrix is binary matrix where rows are features and 
#'   columns are gene sets. The values in the matrix should be either 1 or 0 representing
#'   the presence and absence of a feature in the genesets, repectively. 
#' @param stats other parameter passed to fgseaMultilevel
#' @param gs_desc description of gene sets, it should be a named vector and the names
#'  should be the same as colnames(gs)
#' @param ... other parameters passed to fgseaMultilevel
#' @importFrom fgsea fgsea fgseaMultilevel
#' @examples
#' ## not for users
#' # library(fgsea)
#' # library(Biobase)
#' # dat <- readRDS("Dat/exampleEset.RDS")
#' # fd <- fData(dat)
#' # fdgs <- fd[, grep("^GS\\|", colnames(fd))]
#' # res <- fgsea1(fdgs, stats = fd$`t-test|OV_BR|md`, minSize = 5, maxSize = 500)
#' # res <- fgsea1(
#' #   fdgs, stats = fd$`t-test|OV_BR|md`,  
#' #   minSize = 5, maxSize = 500, gs_desc = colnames(fdgs))
#' @return a \code{data.frame} of fgsea results

fgsea1 <- function(gs, stats, gs_desc = NULL, ...) {
  if (inherits(gs, c("matrix", "dgCMatrix"))) {
    gs <- csc2list(gs)
  } 
  pw <- split(gs$featureId, gs$gsId)
  params <- list(pathways = pw, stats = stats, ...)
  res <- do.call(fgseaMultilevel, params)
  if (!is.null(gs_desc)) {
    res$desc <- as.character(gs_desc[res$pathway])
  }
  as.data.frame(res)
}
