#' Function to perform pairwse t-test
#' @param x an expression matrix, log10 transformed
#' @param pheno phenotype data of x, nrow(pheno) must equal ncol(x)
#' @param compare comparison want to do. If not null, it a nx3 matrix. 
#'   The first column should be column headers in pData, the second and 
#'   third columns should be two values in the columns of pData selected by the header 
#'   (first column). The samples mapped to the two values are compared. If paired comparisons 
#'   to be done, the orders of samples should be mapped
#' @param fillNA Whether NA should be filled. If TRUE, the missing value will
#'   be replaced by a constant a slightly lower than the min detected value (-log10(2))
#' @param ... other parameters passed to t.test
#' @importFrom methods is
#' @export
#' @examples 
#' # reading expression
#' packdir <- system.file("extdata", package = "ExpressionSetViewer")
#' expr <- read.delim(file.path(packdir, "expressionMatrix.tsv"), stringsAsFactors = FALSE)
#' # reading phenotype data
#' pd <- read.delim(file.path(packdir, "sampleGeneral.tsv"), stringsAsFactors = FALSE)
#' 
#' ## Single t-test
#' head(pd)
#' # define comparisons
#' tests <- c("Origin", "RE", "ME")
#' tres <- multi.t.test(x = expr, pheno = pd, compare = tests)
#' 
#' ## multiple t-test
#' head(pd)
#' # define comparisons
#' tests <- rbind(
#' c("Origin", "RE", "ME"),
#' c("Origin", "RE", "LE"),
#' c('TP53.Status', "MT", "WT")
#' )
#' tres <- multi.t.test(x = expr, pheno = pd, compare = tests)
#' 
multi.t.test <- function(x, pheno, compare = NULL, fillNA = FALSE, ...) {
  
  x0 <- x
  if( is.vector(compare) || length(compare) == 3)
    compare <- matrix(compare, nrow = 1)    
  if (fillNA) x <- fillNA(x)  
  if (is.null(compare)) return(NULL)

  tl <- lapply(unique(compare[, 1]), function(x) {
    x <- compare[x == compare[, 1], -1, drop = FALSE]
    unique(unlist(split(x, row(x))))
  })
  names(tl) <- unique(compare[, 1])
  
  df <- data.frame(row.names = rownames(x))
  for ( i in names(tl) ) {
    for ( j in tl[[i]] ) {
      m <- x[, which(pheno[[i]] == j), drop = FALSE]
      m0 <- x0[, which(pheno[[i]] == j), drop = FALSE]
      rv <- rowSums(!is.na(m0))
      rm <- rowMeans(m, na.rm = TRUE)
      rq <- rank(rm, na.last = TRUE)/sum(!is.na(rm))
      rq[is.na(rm)] <- NA
      df[[paste("mean", i, j, sep = "|")]] <- rm
      df[[paste("n value", i, j, sep = "|")]] <- rv
      df[[paste("quantile", i, j, sep = "|")]] <- rq
    }
  } 
  for ( i in 1:nrow(compare) ) {
    v <- compare[i, ]
    i1 <- which(pheno[[v[1]]] == v[2])
    i2 <- which(pheno[[v[1]]] == v[3])
    
    if (length(i1) == 0 || length(i2) == 0)
      stop(paste("Didn't find var:", paste(v, collapse = "-")))
    
    tv <- apply(x, 1, function(xx) {
      t <- try(t.test(xx[i1], xx[i2], var.equal = TRUE, ...), silent = TRUE) # 
      if (!is(t, "htest"))
        return(c(pvalue = NA, mean.diff = NA))
      if (length(t$estimate) == 1)
        md <- t$estimate[[1]] else
          md <- t$estimate[1] - t$estimate[2]
        c(pvalue = t$p.value, mean.diff = md)
    })
    
    pv <- tv[1, ]
    fdr <- p.adjust(pv, method = "fdr")
    
    df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "pvalue", sep = "|")]] <- pv
    df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "log.pvalue", sep = "|")]] <- -log10(pv)
    df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "fdr", sep = "|")]] <- fdr
    df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "log.fdr", sep = "|")]] <- -log10(fdr)
    df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "mean.diff", sep = "|")]] <- tv[2, ]
  }
  df
}

#' perform PCA and prepare for ExpressionSetViewer
#' @param x expression matrix, where rows are feature and samples are columns. Don't need to 
#'   transpose the matrix. 
#' @param n number of components to keep
#' @param prefix prefix of the names
#' @param fillNA Whether NA should be filled. If TRUE, the missing value will
#'   be replaced by a constant a slightly lower than the min detected value (-log10(2)).
#'   If set to false, na.omit will be called before pca
#' @export
#' @importFrom stats prcomp
#' @examples 
#' # reading expression
#' packdir <- system.file("extdata", package = "ExpressionSetViewer")
#' expr <- read.delim(file.path(packdir, "expressionMatrix.tsv"), stringsAsFactors = FALSE)
#' # call PCA
#' pc <- exprspca(expr)
#' head(pc$samples)
#' head(pc$features)
#'  
exprspca <- function(x, n = min(8, ncol(x)-1), prefix = "PCA|All", fillNA = FALSE) {
  
  writePC <- function(x, n) {
    n <- min(n, length(x$sdev))
    var <- round(x$sdev[1:n]^2/(sum(x$sdev^2))*100, digits = 1)
    xx <- x$x[, 1:min(n, ncol(x$x))]
    colnames(xx) <- paste0(prefix, "|", colnames(xx), "(", var, "%", ")")
    pp <- x$rotation[, 1:min(n, ncol(x$x))]
    colnames(pp) <- paste0(prefix, "|", colnames(pp), "(", var, "%", ")")
    list(samples = xx, features = pp)
  }
  
  if (fillNA) {
    x <- fillNA(x)
    pc <- prcomp(t(x))
  } else {
    nr <- nrow(x)
    
    x <- na.omit(x)
    pc <- prcomp(t(x))
    pos <- setdiff(1:nr, attr(x, "na.action"))
    rotation <- matrix(NA, nrow = nr, ncol = ncol(pc$rotation))
    rotation[pos, ] <- pc$rotation
    colnames(rotation) <- colnames(pc$rotation)
    pc$rotation <- rotation
  }
  
  writePC(pc, n = n)
}

#' Filling NA using a constant or half lowest detected value (for each protein) method
#' @param x matrix with NA values
#' @param maxfill the maximum of filled value
#' @param fillingFun function to calculate half value
#'   e.g.
#'   function(x) x - log10(2) # default, when x is in log10 scale
#'   function(x) x - 1 # half lowest detected value when x is in log2 scale
#'   function(x) 0 # replace NA by 0
#' @export

fillNA <- function(x, maxfill = quantile(x, probs = 0.15, na.rm = TRUE), fillingFun = function(x) x - log10(2)) {
  x <- apply(x, 1, function(xx) {
    x3 <- xx
    x3[is.na(x3)] <- min(maxfill, fillingFun(min(x3, na.rm = TRUE)))
    x3
  })
  t(x)
}

