#' Create a nx3 matrix that can be use for triselector given a meta and expression table
#' @description only used inside reactive
#' @param meta a meta data, usually either phenotype data or feature data
#' @param expr expression matrix, optional.
#' @param combine how the meta and expression to be combined. Should be either 
#'  "pheno" or "feature' or "none".
#' @return a nx3 matrix
#' @importFrom shiny req
#' @importFrom stringr str_split_fixed
# 
trisetter <- function(meta, expr=NULL, combine) {
  req(meta)
  nm <- colnames(meta)
  if (!is.null(expr)) {
    combine <- combine[1]
    combine <- match.arg(combine, choices = c("pheno", "feature", "none"))
    if (combine == "feature") {
      if (nrow(expr) != nrow(meta))
        stop("nrow(meta) needs to equal nrow(expr) when combined by 'feature'!")
      rt <- paste0('Sample|Auto|', colnames(expr))
    } else if (combine == "pheno") {
      if (ncol(expr) != nrow(meta))
        stop("nrow(meta) needs to equal ncol(expr) when combined by 'pheno'!")
      rt <- paste0('Feature|Auto|', rownames(expr))
    } else {
      rt <- NULL
    }
    nm <- c(nm, rt)
  }
  str_split_fixed(nm, "\\|", n = 3)
}

#' variable selector 
#' @param x variable return by triselector, a list of length three named as "analysis",
#'  "subset" and "variable'
#' @param expr the expression matrix
#' @param meta a meta matrix
#' @param alternative alternative value to be returned when nothing to select
#' 
varSelector <- function(x, expr, meta, alternative = NULL) {
  if (x$variable %in% c("Select a variable!", "")) {
    if (is.null(alternative))
      return(NULL)
    if (alternative %in% c("", "Select a variable!"))
      return(NULL)
    # alternative <- gsub(" ", "<br>", alternative)
    return(alternative)
  }
  lab <- paste(x, collapse = "|")

  selectIfCan <- function(x, i, dim) {
    if (dim == 1) {
      if (!all(i %in% rownames(x)))
        return(NULL)
      r <- x[i, ]
    } else if (dim == 2) {
      if (!all(i %in% colnames(x))) 
        return(NULL)
      r <- x[, i]
    } else
      stop("dim should be either 1 or 2!")
    r
  }

  if (x$analysis == "Feature" && x$subset == "Auto") {
    x <- selectIfCan(expr, x$variable, dim = 1)
  } else if (x$analysis == "Sample" && x$subset == "Auto") {
    x <- selectIfCan(expr, x$variable, dim = 2)
  } else {
    x <- selectIfCan(meta, lab, dim = 2)
  }
  if (!is.null(x)) attr(x, "label") <- lab
  x
}

#' convert text to number
#' @param x a number or a string can be calculated
#' 
text2num <- function(x) {
  x0 <- suppressWarnings(as.numeric(x))
  if (is.na(x0))
    x0 <-  try(eval(parse(text = x)), silent = TRUE)
  if (!is.numeric(x0)) {
    warning("text2num: cannot convert x to num!")
    return(NULL)
  }
  x0
}

###
line_rect <- function(l, coord) {
  
  x <- l$x
  y <- l$y
  if (is.null(x) && is.null(y))
    return(NULL)
  
  minx <- min(coord$x, na.rm = TRUE)
  maxx <- max(coord$x, na.rm = TRUE)
  miny <- min(coord$y, na.rm = TRUE)
  maxy <- max(coord$y, na.rm = TRUE)
  x.exp <- (maxx-minx)*0.05
  y.exp <- (maxy-miny)*0.05
  if (x.exp == 0 || y.exp == 0)
    return(NULL)
  
  if (grepl("top", l$corner)) {
    y0 <- max(l$y, miny - y.exp)
    y1 <- maxy + y.exp
  } else if (grepl("bottom", l$corner)) {
    y0 <- miny - y.exp
    y1 <- min(l$y, maxy + y.exp)
  } else {
    y0 <- miny - y.exp
    y1 <- maxy + y.exp
  }
  
  if (grepl("left", l$corner)) {
    x0 <- minx - x.exp
    x1 <- min(l$x, maxx + x.exp)
  } else if (grepl("right", l$corner)) {
    x0 <- max(l$x, minx - x.exp)
    x1 <- maxx + x.exp
  } else {
    x0 <- minx - x.exp
    x1 <- maxx + x.exp
  }
  if ((x0 > x1) || (y0 > y1))
    x <- y <- x0 <- x1 <- y0 <- y1 <- NULL
  
  list(
    x = x, y = y, rect = c(x0 = x0, x1 = x1, y0 = y0, y1 = y1)
  )
}

