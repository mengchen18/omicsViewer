#' Correlating a expression matrix with phenotypical variables
#' @description This is a convenience function to perform correlation analysis, 
#'   the output is in a format ready to be incorporated into object to be 
#'   visualized by \code{omicsViewer}. 
#' @param x an expression matrix, rows are the features (e.g. proteins), 
#'   columns are the samples
#' @param pheno  a \code{data.frame} storing the numerical phenotypical variable to 
#'   be correlated with the rows (features) in expression matrix.
#' @param min.value the minimum number of samples required in the correlation 
#'   analysis, if lower than this number, NA will be returned. 
#' @param prefix prefix of the names. Usually don't need to be changed by the user. 
#'   When changes are needed, the prefix should be in a format like 
#'   [analysis name]|[subset] so the "analysis name" and "subset" can be selected
#'   in the \code{omicsViewer}.
#' @importFrom psych corr.test
#' @importFrom matrixStats rowMaxs rowMins
#' @export
#' @examples 
#' e1 <- matrix(rnorm(500), 50, 10)
#' rownames(e1) <- paste0("FT", 1:50)
#' p1 <- matrix(rnorm(50), 10, 5)
#' colnames(p1) <- paste0("PH", 1:5)
#' colnames(e1) <- rownames(p1) <- paste0("S", 1:10)
#' correlationAnalysis(x = e1, pheno = p1, min.value = 8)
#' @return Every correlation analysis returns a \code{data.frame} with five 
#'   columns:
#'   \code{R} - pearson correlation coefficient
#'   \code{N} - number of values used in the analysis
#'   \code{P} - p-values returned by pearson correlation analysis
#'   \code{logP} - log transformed p-values
#'   \code{range} - the range of values in expression matrix used in the analysis

correlationAnalysis <- function(x, pheno, min.value = 12, prefix = "Cor") {
  
  if (is.data.frame(pheno)) {
    cn <- colnames(pheno)
    rn <- rownames(pheno)
    pheno <- apply(pheno, 2, function(x) as.numeric(as.character(x)))
    colnames(pheno) <- cn
    rownames(pheno) <- rn
  }
  
  if (is.null(colnames(pheno)))
    colnames(pheno) <- make.names(seq_len(ncol(pheno)))
  
  if (ncol(x) != nrow(pheno))
    stop("ncol x != nrow pheno")
  
  r <- psych::corr.test(t(x), pheno, use = "pair", adjust = "none", ci = FALSE)
  
  rs <- lapply(seq_len(ncol(pheno)), function(i) {
    if (length(r$n) == 1)
      nv <- r$n else
        nv <- r$n[, i]
      df <- data.frame(
        R = r$r[, i],
        N = nv,
        P = r$p[, i],
        logP = -log10(r$p[, i])
      )
      texp <- x[, !is.na(pheno[, i])]
      df$range <- rowMaxs(texp, na.rm = TRUE) - rowMins(texp, na.rm = TRUE)
      
      ii <- which(df$N < min.value)
      df[c("R", "P", "logP")] <- lapply(df[c("R", "P", "logP")], function(x) {
        x[ii] <- NA
        x
      })
      
      if (all(is.na(df$R)))
        return(df[, character(0)])
      
      colnames(df) <- paste(colnames(pheno)[i], colnames(df), sep = "|")
      df
  })
  rs <- do.call(cbind, rs)
  if (ncol(rs) > 0)
    colnames(rs) <- paste(prefix, colnames(rs), sep = "|")
  rs
}


#' Function to perform multiple t-tests on an expression matrix
#' @description This is a convenience function to perform multiple student's t-test. 
#'   The output is in a format ready to be incorporated into object to be visualized by 
#'   \code{omicsViewer}. This function use \code{\link{t.test}}. 
#' @param x an expression matrix, usually log10 transformed.
#' @param pheno phenotype data of x, the number of rows in \code{pheno} must equal
#'   the number of columns of \code{x}. Please refer to examples for more details. 
#' @param compare NULL or a matrix with three columns to define the comparisons to do. 
#'   When a matrix is given, the first column should be one of the column headers 
#'   in \code{pheno}; then the second and third columns should be two values 
#'   presented (more than once) in the columns of \code{pheno} selected by the
#'   values in the first column. The samples mapped to the two values are compared. 
#'   If paired comparisons to be done, the orders of samples should be mapped
#' @param fillNA logical; whether NA should be filled? If FALSE (default), t test 
#'   will be performed whenever possible. If not possible, then NA will be returned. 
#'   If TRUE, the missing value will be replaced using \code{\link{fillNA}}. 
#' @param ... other parameters passed to \code{\link{t.test}}
#' @importFrom methods is
#' @return a \code{data.frame} stores the t-test results with the follow columns:
#' \code{mean|[selected header in pheno]|[group 1 in test]} - The mean value of group 1
#' \code{n value|[selected header in pheno]|[group 1 in test]} - The number of value used in the test for group 1
#' \code{quantile|[selected header in pheno]|[group 1 in test]} - The quantile of means values in group 1
#' \code{mean|[selected header in pheno]|[group 2 in test]} - The mean value of group 2
#' \code{n value|[selected header in pheno]|[group 2 in test]} - The number of value used in the test for group 2
#' \code{quantile|[selected header in pheno]|[group 2 in test]} - The quantile of means values in group 2
#' \code{ttest|[group 1 in test]_vs_[group 2 in test]|pvalue} - The p-value return by \code{\link{t.test}}
#' \code{ttest|[group 1 in test]_vs_[group 2 in test]|log.pvalue} - The -log10 transformed p-value
#' \code{ttest|[group 1 in test]_vs_[group 2 in test]|fdr} - The BH method corrected p-values, e.g. FDR
#' \code{ttest|[group 1 in test]_vs_[group 2 in test]|log.fdr} - The -log10 transformed FDR
#' \code{ttest|[group 1 in test]_vs_[group 2 in test]|mean.diff} - The difference between the means of the two groups, e.g. fold change
#' @export 
#' @examples 
#' # reading expression
#' packdir <- system.file("extdata", package = "omicsViewer")
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
  for ( i in seq_len(nrow(compare)) ) {
    v <- compare[i, ]
    i1 <- which(pheno[[v[1]]] == v[2])
    i2 <- which(pheno[[v[1]]] == v[3])
    
    if (length(i1) == 0 || length(i2) == 0)
      stop("Didn't find var: ", paste(v, collapse = "-"))
    
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

#' Perform PCA and prepare results for \code{omicsViewer}
#' @description This is a convenience function to perform PCA on expression matrix, 
#'   the output of PCA
#'   will be in a format ready to be incorporated into object to be visualized by 
#'   \code{omicsViewer}.
#' @param x an expression matrix, where rows are features and samples are on columns. 
#' @param n number of components to keep
#' @param prefix prefix of the names. Usually don't need to be changed by the user. 
#'   When changes are needed, the prefix should be in a format like 
#'   [analysis name]|[subset] so the "analysis name" and "subset" can be selected
#'   in the \code{omicsViewer}.
#' @param fillNA logical; whether NA should be filled? If FALSE (default), 
#'   \code{na.omit} will be called before PCA. If TRUE, the missing value 
#'   will be replaced using \code{\link{fillNA}}. 
#' @param ... other parameters passed to \code{\link{prcomp}}
#' @export
#' @importFrom stats prcomp
#' @examples 
#' # reading expression
#' packdir <- system.file("extdata", package = "omicsViewer")
#' expr <- read.delim(file.path(packdir, "expressionMatrix.tsv"), stringsAsFactors = FALSE)
#' # call PCA
#' pc <- exprspca(expr)
#' head(pc$samples)
#' head(pc$features)
#' @return a \code{data.frame} storing the PCA results

exprspca <- function(x, n = min(8, ncol(x)-1), prefix = "PCA|All", fillNA = FALSE, ...) {
  
  writePC <- function(x, n) {
    n <- min(n, length(x$sdev))
    var <- round(x$sdev[seq_len(n)]^2/(sum(x$sdev^2))*100, digits = 1)
    xx <- x$x[, seq_len(min(n, ncol(x$x)))]
    colnames(xx) <- paste0(prefix, "|", colnames(xx), "(", var, "%", ")")
    pp <- x$rotation[, seq_len(min(n, ncol(x$x)))]
    colnames(pp) <- paste0(prefix, "|", colnames(pp), "(", var, "%", ")")
    list(samples = xx, features = pp)
  }
  
  if (fillNA) {
    x <- fillNA(x)
    pc <- prcomp(t(x))
  } else {
    nr <- nrow(x)
    x <- na.omit(x)
    pc <- prcomp(t(x), ...)
    pos <- setdiff(seq_len(nr), attr(x, "na.action"))
    rotation <- matrix(NA, nrow = nr, ncol = ncol(pc$rotation))
    rotation[pos, ] <- pc$rotation
    colnames(rotation) <- colnames(pc$rotation)
    pc$rotation <- rotation
  }
  
  writePC(pc, n = n)
}

#' Filling NAs in a matrix using constants calculated from user the defined function
#' @description This function is usually use to impute missing values in expression matrix, where the rows are
#'   feature and columns are samples. This function impute the missing values on the row-wise, that is, every 
#'   row will be imputed using different constant. 
#' @param x a matrix with NA values
#' @param maxfill the maximum filled value, if the value calculated by \code{fillingFun} is greater than 
#'   \code{maxfill}, then \code{maxfill} will the used to replace NAs. 
#' @param fillingFun function to calculate the filling values. It should be a function accept at least one
#'   argument "x", which is a row of input expression matrix. The default is 
#'   \code{function(x) min(x, na.rm = TRUE) - log10(2)} 
#'   corresponds to the "half of lowest detected values" if the expression matrix is log10 transformed.
#'   More examples:#'   
#'   \code{function(x) min(x, na.rm = TRUE) - 1 }# half of lowest detected value when expression matrix
#'     is in log2 scale
#'   \code{function(x) 0} # replace NA by 0
#' @note The returned matrix may have -Inf, which may need to be filtered/replaced additionally
#' @export
#' @examples
#' m <- matrix(rnorm(200), 20, 10)
#' m[sample(1:200, size = 20)] <- NA
#' mf <- fillNA(m)
#' @return a matrix without NAs
fillNA <- function(x, maxfill = quantile(x, probs = 0.15, na.rm = TRUE), fillingFun = function(x) min(x, na.rm = TRUE) - log10(2)) {
  xf <- apply(x, 1, function(xx) {
    x3 <- xx
    x3[is.na(x3)] <- min(maxfill, fillingFun(xx))
    x3
  })
  xf <- t(xf)
  rownames(xf) <- rownames(x)
  colnames(xf) <- colnames(xf)
  xf
}

