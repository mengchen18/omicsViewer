#' Prepare data that can be viewer by ExpressionSetViewer
#' @param expr expression matrix where the rows are feature and columns are samples, matrix should be log10 transformed
#' @param pData phenotype data
#' @param fData feature data
#' @param PCA pca
#' @param ncomp number of components to keep 
#' @param pca.fillNA logical, whether the NA should be filled with a constant in PCA. 
#' @param t.test How t-test should be done. it should be null (no t-test will be done)
#'   or a nx3 matrix. The first column should be column headers in pData, the second and 
#'   third columns should be two values in the columns of pData selected by the header 
#'   (first column). The samples mapped to the two values are compared. If paired comparisons 
#'   to be done, the orders of samples should be mapped. 
#' @param ttest.fillNA logical, whether the NA should be filled with a constant in t-test. 
#' @param gs gene set data
#' @param stringDB string database
#' @param surv survival data
#' @param ... arguments passed to t.test, such as paired.
#' @importFrom Biobase ExpressionSet AnnotatedDataFrame
#' @details expr should be a matrix with unique rownames and colnames. 
#'   pData and fData don't need to hav row names, row names will be assigned
#'   according to row names and columns of expr. 
#'   PCA and pairwise t-test could be done. The t.test argument should be a 
#'   nx3 matrix, the i
#' @export
#' @examples 
#' packdir <- system.file("extdata", package = "ExpressionSetViewer")
#' # reading expression
#' expr <- read.delim(file.path(packdir, "expressionMatrix.tsv"), stringsAsFactors = FALSE)
#' colnames(expr) <- make.names(colnames(expr))
#' rownames(expr) <- make.names(rownames(expr))
#' # reading feature data
#' fd <- read.delim(file.path(packdir, "featureGeneral.tsv"), stringsAsFactors = FALSE)
#' # reading phenotype data
#' pd <- read.delim(file.path(packdir, "sampleGeneral.tsv"), stringsAsFactors = FALSE)
#' 
#' #  reading other datasets
#' drugData <- read.delim(file.path(packdir, "sampleDrug.tsv"))
#' # survival data
#' # this data is from cell line, the survival data are fake data to 
#' # show how to use the survival data in #' ExpressionSetViewer
#' surv <- read.delim(file.path(packdir, "sampleSurv.tsv"))
#' # gene set information
#' genesets <- read.delim(file.path(packdir, "featureGS.tsv"))
#' 
#' # Define t-test to be done
#' tests <- rbind(
#' c("Origin", "RE", "ME"),
#' c("Origin", "RE", "LE"),
#' c('TP53.Status', "MT", "WT")
#' )
#' 
#' # prepare column for stringDB query
#' strid <- sapply(strsplit(fd$Protein.ID, ";|-"), "[", 1)
#' 
#' ###
#' d <- prepEsetViewer(
#' expr = expr, pData = pd, fData = fd, 
#' PCA = TRUE, pca.fillNA = TRUE,
#' t.test = tests, ttest.fillNA = FALSE, 
#' gs = genesets, stringDB = strid, surv = surv)
#' 
#' saveRDS(d, file = "dtest.RDS")
#' ##  to open the viewer
#' # ExpressionSetViewer("./")
#' @return an object of \code{ExpressionSet} that can be visualized using
#' \code{ExpressionSetViewer}

prepEsetViewer <- function(
  expr, pData, fData, 
  PCA = TRUE, ncomp = min(8, ncol(expr)), pca.fillNA = TRUE,
  t.test = NULL, ttest.fillNA = FALSE, ..., 
  gs = NULL, stringDB = NULL, surv = NULL) {
  
  
  p0 <- pData
  
  ## ======================= check dimension and names  ============================
  de <- dim(expr)
  if (nrow(pData) != de[2])
    stop("nrow pData != ncol(exprs)")
  if (nrow(fData) != de[1])
    stop("nrow pData != nrow(exprs)")
  
  expr_rn <- rownames(expr)
  if (is.null(expr_rn) && !is.null(rownames(fData)))
    expr_rn <- rownames(fData)
  if (is.null(expr_rn))
    expr_rn <- str_pad(1:nrow(expr), width = nchar(nrow(expr)), pad = "0")
  expr_rn <- make.names(expr_rn, unique = TRUE)
  expr_cn <- colnames(expr)
  if (is.null(expr_cn) && !is.null(rownames(pData)))
    expr_cn <- rownames(pData)
  if (is.null(expr_cn))
    expr_cn <- str_pad(1:ncol(expr), width = nchar(ncol(expr)), pad = "0")
  expr_cn <- make.names(expr_cn, unique = TRUE)
  
  pData_rn <- rownames(pData)
  fData_rn <- rownames(fData)
  
  if (is.null(expr_rn) || is.null(expr_cn))
    stop("expr needs to be a matrix with dimnames")
  
  if (!identical(expr_cn, pData_rn)) {
    message("rownames of pData reset by colnames of exprs")
    rownames(pData) <- expr_cn
  }
  if (!identical(expr_rn, fData_rn)) {
    message("rownames of fData reset by rownames of exprs")
    rownames(fData) <- expr_rn
  }
  rownames(expr) <- expr_rn
  colnames(expr) <- expr_cn
  pData <- cbind(pData, numberOfFeatures = colSums(!is.na(expr)))
  
  colnames(pData) <- paste0('General|All|', trimws(colnames(pData)))
  colnames(fData) <- paste0('General|All|', trimws(colnames(fData)))
  
  
  ## ======================= PCA ============================
  if (PCA) {
    pc <- exprspca(expr, n = ncomp, fillNA =  pca.fillNA)
    pData <- cbind(pData, pc$samples)
    fData <- cbind(fData, pc$features)
  }
  
  
  ## ======================= perform t-test ============================
  if (!is.null(t.test)) {
    tres <- multi.t.test(x = expr, pheno = p0, compare = t.test, fillNA = ttest.fillNA, ...)
    fData <- cbind(fData, tres)
  }
  
  ## ======================= ranking ==========================
  rk <- data.frame(apply(expr, 2, rank), stringsAsFactors = FALSE)
  colnames(rk) <- paste0("Rank|All|", colnames(rk))
  fData <- cbind(fData, rk)
  
  ##  ================== make headers for others =========================
  # gs
  if (!is.null(gs)) {
    if (is.data.frame(gs))
      gs <- as.matrix(gs)
    if (is.vector(gs) && length(gs) == nrow(expr)) {
      gs <- matrix(gs, ncol = 1)
      colnames(gs) <- "geneset"
    } else if ( !(is.matrix(gs) && nrow(gs) == nrow(expr)) )
      stop("incompatible 'gs'")
    if (is.null(colnames(gs)))
      stop("colnames of gs should not bu null!")
    if (any(duplicated(colnames(gs))))
      stop('colnames of gs should be unique!')
    colnames(gs) <- paste0('GS|All|', colnames(gs))
    gs <- data.frame(gs, check.names = FALSE, stringsAsFactors = FALSE)
    
    fData <- cbind(fData, gs)
  }
  
  # stringDB
  if (!is.null(stringDB)) {
    if (length(stringDB) == nrow(expr))
      strdb <- data.frame(
        "StringDB|All|ID" = stringDB, 
        stringsAsFactors = FALSE, 
        check.names = FALSE) else
          stop('stringDB should be a vector same length as nrow(expr), containing the IDs can be used to query stringDB')
    fData <- cbind(fData, strdb)
  }
  
  # surv
  if (!is.null(surv)) {
    if (is.vector(surv) && length(surv) == nrow(expr)) {
      surv <- data.frame( "Surv|all|surv" = surv, 
                          stringsAsFactors = FALSE, 
                          check.names = FALSE)
    } else if (is.matrix(surv) || is.data.frame(surv)) {
      if (nrow(surv) != ncol(expr))
        stop( "nrow(surv) should equals ncol(expr)" )
      if (is.null(colnames(surv)))
        colnames(surv) <- paste0("Surv", 1:ncol(surv))
      surv <- data.frame(surv, stringsAsFactors = FALSE, check.names = TRUE)
      colnames(surv) <- paste0("Surv|all|", trimws(colnames(surv)))
    } else
      stop("incompatible 'surv'")
    sv <- sapply(surv, function(x) {
      all(grepl("\\D|+", sub("+$", "", surv)))
    })
    if (any(!sv))
      stop("only numbers and + are allowed in surv")
    
    pData <- cbind(pData, surv)
  }
  
  ExpressionSet(assayData = as.matrix(expr), 
                phenoData = AnnotatedDataFrame(pData), 
                featureData = AnnotatedDataFrame(fData))
}


#' Correlating a expression matrix with phenotypical variables
#' @param expr expression matrix, rows are the features (e.g. proteins), colnames are the samples
#' @param pheno  a matrix stores the phenotypical variable to be correlated with expr. The rows 
#'   should be samples, columns are phenotypical variables.
#' @param min.value the minimum number of samples required in the correlation analysis
#' @param prefix prefix used for columns headers
#' @importFrom psych corr.test
#' @importFrom matrixStats rowMaxs rowMins
#' @export
#' @examples 
#' e1 <- matrix(rnorm(500), 50, 10)
#' rownames(e1) <- paste0("FT", 1:50)
#' p1 <- matrix(rnorm(50), 10, 5)
#' colnames(p1) <- paste0("PH", 1:5)
#' colnames(e1) <- rownames(p1) <- paste0("S", 1:10)
#' correlationAnalysis(expr = e1, pheno = p1, min.value = 8)
#' @return a object of \code{data.frame} containing the results of correlation analysis

correlationAnalysis <- function(expr, pheno, min.value = 12, prefix = "Cor") {
  
  if (is.null(colnames(pheno)))
    colnames(pheno) <- make.names(1:ncol(pheno))
  
  if (ncol(expr) != nrow(pheno))
    stop("ncol expr != nrow pheno")
  
  r <- psych::corr.test(t(expr), pheno, use = "pair", adjust = "none", ci = FALSE)
  
  rs <- lapply(1:ncol(pheno), function(i) {
    if (length(r$n) == 1)
      nv <- r$n else
        nv <- r$n[, i]
    df <- data.frame(
      R = r$r[, i],
      N = nv,
      P = r$p[, i],
      logP = -log10(r$p[, i])
    )
    texp <- expr[, !is.na(pheno[, i])]
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

#' Annotation of gene/protein function using multiple IDs.
#' @param idList list of protein IDs, e.g. list(c("ID1", "ID2"), c("ID13"), c("ID4", "ID8", "ID10"))
#' @param gsIdMap a data frame for geneset to id map, it has two columns
#'   - id: the ID column
#'   - term: annotation terms
#'   e.g. 
#'   gsIdMap <- data.frame(
#'     id = c("ID1", "ID2", "ID1", "ID2", "ID8", "ID10"),
#'     term = c("T1", "T1", "T2", "T2", "T2", "T2"),
#'     stringsAsFactors = FALSE
#'     )
#' @param minSize minimum size of gene sets
#' @param maxSize maximum size of gene sets
#' @importFrom matrixStats colSums2
#' @export
#' @examples
#' terms <- data.frame(
#'   id = c("ID1", "ID2", "ID1", "ID2", "ID8", "ID10"),
#'  term = c("T1", "T1", "T2", "T2", "T2", "T2"),
#'   stringsAsFactors = FALSE
#' )
#' features <- list(c("ID1", "ID2"), c("ID13"), c("ID4", "ID8", "ID10"))
#' gsAnnotIdList(idList = features, gsIdMap = terms, minSize = 1, maxSize = 500)
#' @return a binary matrix. The number of rows is the same with length of idList, the columns 
#' are the annotated gene set. 

gsAnnotIdList <- function(idList, gsIdMap, minSize = 5, maxSize = 500) {
  pid <- data.frame(
    id = unlist(idList),
    index = rep(1:length(idList), times = sapply(idList, length)),
    stringsAsFactors = FALSE
  )
  
  gsIdMap <- gsIdMap[gsIdMap$id %in% pid$id, ]
  gsIdMap <- split(gsIdMap$id, gsIdMap$term)
  v <- sapply(gsIdMap, function(x, s) {
    s[unique(pid$index[fastmatch::"%fin%"(pid$id, x)])] <- 1
    s
  }, s = rep(0, length(idList)))
  colnames(v) <- names(gsIdMap)
  
  nm <- matrixStats::colSums2(v)
  v <- v[, which(nm >= minSize & nm <= maxSize)]
  v
}
