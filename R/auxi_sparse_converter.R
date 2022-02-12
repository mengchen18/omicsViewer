#' convert a column compressed sparse matrix to a list
#' @param x a CsparseMatrix object
#' @importFrom fastmatch fmatch
#' @rawNamespace import(Matrix, except = image)
#' @return a sparse frame in data.frame 
#'  
csc2list <- function(x) {
  if (is.null(x@Dimnames[[1]]) || is.null(x@Dimnames[[2]]))
    stop("csmlist: x need to have dimnames!")
  df <- data.frame(
    featureId = x@Dimnames[[1]][x@i+1],
    gsId = rep(x@Dimnames[[2]], x@p[-1] - x@p[-length(x@p)]),
    stringsAsFactors = TRUE
  )
  if (hasAttr(x, "x"))
    df$weight <- x@x
  df
}

#' convert a list to column compressed sparse matrix
#' @param l a data.frame with at least two columns - featureId, gsId; optionally
#'   a "weight" column. 
#' @param dimnames a list of dimnames, should contain at least one element for the
#'   row names. 
#' @importFrom fastmatch fmatch
#' @return a sparse matrix, CsparseMatrix, column compressed
#' 
list2csc <- function(l, dimnames) {
  if (!is.list(dimnames))
    stop("list2csc: dimnames should be a list!")
  if (is.null(dimnames[[1]]))
    stop("list2csc: dimnames should contain as least one element for the rownames of the matrix")
  if (!all(colnames(l) %in% c("featureId", "gsId", "weight")))
    stop('list2csc: colnames of l should be "featureId", "gsId", "weight"')
  
  names <- dimnames[[1]]
  if (length(dimnames) > 1)
    cn <- dimnames[[2]] else
      cn <- unique(l$gsId)
  args <- list(
    i = c(fmatch( l$featureId, names )), 
    j = c(fmatch(l$gsId, cn)),
    dims = c(length(names), length(cn)),
    dimnames = list(names, cn)
  )
  if ("weight" %in% colnames(l))
    args$x <- l$weight
  do.call(sparseMatrix, args)
}
