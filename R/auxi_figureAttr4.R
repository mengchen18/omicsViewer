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
  stringr::str_split_fixed(nm, "\\|", n = 3)
}

varSelector <- function(x, expr, meta, alternative = NULL) {
  if (x$variable %in% c("Select a variable!", "")) {
    if (is.null(alternative))
      return(NULL)
    if (alternative %in% c("", "Select a variable!"))
      return(NULL)
    alternative <- gsub(" ", "<br>", alternative)
    return(alternative)
  }
  lab <- paste(x, collapse = "|")
  if (x$analysis == "Feature" && x$subset == "Auto") {
    x <- expr[x$variable, ]
  } else if (x$analysis == "Sample" && x$subset == "Auto") {
    x <- expr[, x$variable]
  } else
    x <- meta[, lab]
  attr(x, "label") <- lab
  x
}