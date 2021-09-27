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