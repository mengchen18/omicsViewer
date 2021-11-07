getSearchCols <- function(x) {
  if (is.null(x))
    return(NULL)
  lapply(x$columns, "[[", "search")
}

getOrderCols <- function(x) {
  if (is.null(x))
    return(NULL)
  x[["order"]]
}

null2empty <- function(x) {
  if (is.null(x))
    return("")
  x
}

exprsImpute <- function(x) { 
	v <- try(x@assayData$exprs_impute, silent = TRUE)
    if (inherits(x, "try-error")) 
      v <- NULL
    v
    }