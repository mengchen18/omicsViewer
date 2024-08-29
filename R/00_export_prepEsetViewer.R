#' Prepare object to be viewed by omicsViewer
#' @description This is a convenience function to prepare the data to be visualized using \code{\link{omicsViewer}}.
#'   The result of PCA and t-test could be included directly. 
#' @param expr expression matrix where the rows are feature and columns are samples, 
#'   matrix should be log10 transformed and have unique row and column names
#' @param pData phenotype data
#' @param fData feature data
#' @param PCA pca
#' @param ncomp number of components to keep 
#' @param pca.fillNA logical, whether the NA should be filled with a constant in PCA. 
#' @param t.test will be passed to the \code{compare} argument in \code{\link{multi.t.test}}
#' @param ttest.fillNA logical, whether the NA should be filled with a constant in t-test. 
#' @param gs gene-set data, please refer to examples for more details about the format
#' @param stringDB the IDs that can be used in the STRING database (https://string-db.org/) query. 
#' @param surv survival data, please refer to examples for more details about the format
#' @param SummarizedExperiment logical; whether to return an object of class \code{SummarizedExperiment}. 
#'   If set to FALSE, the function will return an \code{ExpressionSet} object.
#' @param ... arguments passed to \code{\link{t.test}}, such as \code{paired}.
#' @importFrom Biobase ExpressionSet AnnotatedDataFrame
#' @export
#' @examples 
#' packdir <- system.file("extdata", package = "omicsViewer")
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
#' # show how to use the survival data in #' omicsViewer
#' surv <- read.delim(file.path(packdir, "sampleSurv.tsv"))
#' # gene set information
#' genesets <- read_gmt(file.path(packdir, "geneset.gmt"), data.frame = TRUE)
#' gsannot <- gsAnnotIdList(idList = rownames(fd), gsIdMap = genesets, data.frame = TRUE)
#' 
#' # Define t-test to be done, a matrix nx3
#' # every row define a t-test, the format
#' # [column header] [group 1 in the test] [group 2 in the test]
#' tests <- rbind(
#'  c("Origin", "RE", "ME"),
#'  c("Origin", "RE", "LE"),
#'  c('TP53.Status', "MT", "WT")
#'  )
#' # prepare column for stringDB query
#' strid <- sapply(strsplit(fd$Protein.ID, ";|-"), "[", 1)
#' ###
#' d <- prepOmicsViewer(
#'   expr = expr, pData = pd, fData = fd, 
#'   PCA = TRUE, pca.fillNA = TRUE,
#'   t.test = tests, ttest.fillNA = FALSE, 
#'   gs = gsannot, stringDB = strid, surv = surv)
#' # feature space - default x axis
#' attr(d, "fx") <- "ttest|RE_vs_ME|mean.diff"
#' # feature space - default y axis
#' attr(d, "fy") <- "ttest|RE_vs_ME|log.fdr"
#' # sample space - default x axis
#' attr(d, "sx") <- "PCA|All|PC1("
#' # sample space - default y axis
#' attr(d, "sy") <- "PCA|All|PC2("
#' # Save object and view
#' # saveRDS(d, file = "dtest.RDS")
#' ##  to open the viewer
#' # omicsViewer("./")
#' @return an object of \code{ExpressionSet} or \code{SummarizedExperiment} that can be visualized using
#' \code{omicsViewer}

prepOmicsViewer <- function(
  expr, pData, fData, 
  PCA = TRUE, ncomp = min(8, ncol(expr)), pca.fillNA = TRUE,
  t.test = NULL, ttest.fillNA = FALSE, ..., 
  gs = NULL, stringDB = NULL, surv = NULL, 
  SummarizedExperiment = TRUE) {
  
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
    expr_rn <- str_pad(seq_len(nrow(expr)), width = nchar(nrow(expr)), pad = "0")
  expr_rn <- make.names(expr_rn, unique = TRUE)
  expr_cn <- colnames(expr)
  if (is.null(expr_cn) && !is.null(rownames(pData)))
    expr_cn <- rownames(pData)
  if (is.null(expr_cn))
    expr_cn <- str_pad(seq_len(ncol(expr)), width = nchar(ncol(expr)), pad = "0")
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
    if (pca.fillNA) {
      pc <- exprspca(expr, n = ncomp, fillNA =  pca.fillNA)
      pData <- cbind(pData, pc$samples)
      fData <- cbind(fData, pc$features)
    }
    pc <- try( exprspca(expr, n = ncomp, fillNA = FALSE, prefix = "PCA|removeMissing") )
    if (!inherits(pc, "try-error")) {
      pData <- cbind(pData, pc$samples)
      fData <- cbind(fData, pc$features)
    } else 
      warning("PCA on missing value removed matrix cannot be done. It's like because 
              no feature remains after missing value filtering. ")
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
  
  # ==================  surv ================== 
  if (!is.null(surv)) {
    if (is.vector(surv) && length(surv) == ncol(expr)) {
      surv <- data.frame( "Surv|all|surv" = surv, 
                          stringsAsFactors = FALSE, 
                          check.names = FALSE)
    } else if (is.matrix(surv) || is.data.frame(surv)) {
      if (nrow(surv) != ncol(expr))
        stop( "nrow(surv) should equals ncol(expr)" )
      if (is.null(colnames(surv)))
        colnames(surv) <- paste0("Surv", seq_len(ncol(surv)))
      surv <- data.frame(surv, stringsAsFactors = FALSE, check.names = TRUE)
      colnames(surv) <- paste0("Surv|all|", trimws(colnames(surv)))
    } else
      stop("incompatible 'surv'")
    sv <- vapply(surv, function(x) {
      all(grepl("\\D|+", sub("+$", "", surv)))
    }, logical(1))
    if (any(!sv))
      stop("only numbers and + are allowed in surv")
    
    pData <- cbind(pData, surv)
  }
  
  # ==================  gs ================== 
  if (!is.null(gs)) {
    if (all(colnames(gs) %in% c("featureId", "gsId", "weight"))) {
      gs$featureId <- factor(rownames(fData)[gs$featureId])
      gs$gsId <- factor(gs$gsId)
      gs$weight <- as.integer(gs$weight)
    } else {
      if (is.vector(gs) && length(gs) == nrow(expr)) {
        gs <- matrix(gs, ncol = 1)
        colnames(gs) <- "geneset"
      }
      if (!inherits(gs, c("matrix", "dgCMatrix")))
        stop("gs should either be a (sparse) matrix or data.frame with three columns: featureId, gsId, weight!")
      if (is.null(rownames(gs)))
        rownames(gs) <- rownames(fData)
      if ( !(nrow(gs) == nrow(expr)) )
        stop("incompatible 'gs'")
      if (is.null(colnames(gs)))
        stop("colnames of gs should not be null!")
      if (any(duplicated(colnames(gs))))
        stop('colnames of gs should be unique!')
    }
    attr(fData, "GS") <- gs
  }
  
  # options to set default axis
  fx1 <- grep("ttest\\|(.*?)_vs_(.*?)\\|mean.diff", colnames(fData), value = TRUE)
  fy1 <- intersect(colnames(fData), sub("mean.diff$", "log.fdr", fx1))
  fx2 <- grep("PCA\\|All\\|PC1\\(", colnames(fData), value = TRUE)
  fy2 <- grep("PCA\\|All\\|PC2\\(", colnames(fData), value = TRUE)
  px <- grep("PCA\\|All\\|PC1\\(", colnames(pData), value = TRUE)
  py <- grep("PCA\\|All\\|PC2\\(", colnames(pData), value = TRUE)

  exprsWithAttr <- function(x, fillNA = FALSE, environment = FALSE, attrs = c("rowDendrogram", "colDendrogram")) {
    if (environment)
     aenv <- new.env() else 
       aenv <- list()
    mx <- as.matrix(x)    
    for (i in attrs) attr(mx, i) <- attr(x, i)
    aenv$exprs <- mx
    if (fillNA) {      
      mxf <- fillNA(mx)
      for (i in attrs) attr(mxf, i) <- attr(x, i)
      aenv$exprs_impute <- mxf
    }
    aenv
  }
  
  # prep object
  if (!SummarizedExperiment) {
    aenv <- exprsWithAttr(expr, fillNA = pca.fillNA || ttest.fillNA, environment = TRUE)
    res <- ExpressionSet(
      assayData = aenv, 
      phenoData = AnnotatedDataFrame(pData), 
      featureData = AnnotatedDataFrame(fData))
  } else {
    
    DataFrameWithAttr <- function(x) {
      attrs <- setdiff(names(attributes(x)), c("names", "class", "row.names"))
      attr_list <- lapply(attrs, function(attr_name) attr(x, attr_name))
      names(attr_list) <- attrs
      x <- DataFrame(x, check.names = FALSE)
      for (i in names(attr_list)) 
        attr(x, i) <- attr_list[[i]]
      x
    }
    
    aenv <- exprsWithAttr(expr, fillNA = pca.fillNA || ttest.fillNA, environment = FALSE)
    pd <- DataFrameWithAttr(pData)
    fd <- DataFrameWithAttr(fData)
    
    res <- SummarizedExperiment(
      assays = aenv,
      rowData = fd, 
      colData = pd
    )
  }
  
  # set default axes
  if (length(fx1) >= 1 && length(fy1) >= 1) {
    attr(res, "fx") <- fx1[1]
    attr(res, "fy") <- fy1[1]
  } else if (length(fx2) >= 1 && length(fy2) >= 1) {
    attr(res, "fx") <- fx2[1]
    attr(res, "fy") <- fy2[1]
  } 
  
  if (length(px) >= 1 && length(py) >= 1) {
    attr(res, "sx") <- px[1]
    attr(res, "sy") <- py[1]
  }
  
  res
}
