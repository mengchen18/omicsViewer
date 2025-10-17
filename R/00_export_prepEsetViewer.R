#' Prepare Omics Data for Visualization with omicsViewer
#'
#' @description
#' A comprehensive data preparation function that processes expression matrices and associated
#' metadata for interactive visualization with \code{\link{omicsViewer}}. Automatically performs
#' dimensionality reduction (PCA), statistical testing (t-tests), and integrates gene set
#' annotations, STRING database IDs, and survival data.
#'
#' @param expr Numeric matrix. Expression data with features in rows and samples in columns.
#'   Should be log-transformed (e.g., log2 or log10). Row and column names must be unique.
#'   Missing values (NA) are permitted if \code{pca.fillNA} or \code{ttest.fillNA} are TRUE.
#' @param pData Data.frame. Sample/phenotype metadata with one row per sample. Row names must
#'   match column names of \code{expr}. Should contain grouping variables for statistical tests.
#' @param fData Data.frame. Feature metadata with one row per feature. Row names must match
#'   row names of \code{expr}. Can include gene symbols, descriptions, database IDs, etc.
#' @param PCA Logical. Whether to perform Principal Component Analysis. Default: TRUE.
#'   Results are added to both sample and feature metadata.
#' @param ncomp Integer. Number of principal components to compute. Default: minimum of 8
#'   or the number of samples. Ignored if \code{PCA = FALSE}.
#' @param pca.fillNA Logical. If TRUE, missing values in \code{expr} are imputed before PCA
#'   by replacing with minimum value * 0.9. Default: TRUE. Two PCAs are performed:
#'   one with imputation and one without (if possible).
#' @param t.test Matrix or NULL. Definition of t-tests to perform. Should be an n×3 matrix where
#'   each row specifies: [column_name, group1, group2]. The column should exist in \code{pData}.
#'   Example: \code{rbind(c("Treatment", "Drug", "Control"), c("Genotype", "WT", "KO"))}.
#'   Results are added as columns to \code{fData}. NULL = no t-tests.
#' @param ttest.fillNA Logical. Whether to impute missing values before t-tests.
#'   Default: FALSE (features with NAs are excluded from testing).
#' @param gs Gene set annotations in one of two formats:
#'   \itemize{
#'     \item Data.frame with columns: \code{featureId} (indices), \code{gsId} (gene set IDs),
#'           \code{weight} (optional weights). See \code{\link{gsAnnotIdList}}.
#'     \item Matrix or sparse matrix (dgCMatrix) with features in rows and gene sets in columns.
#'           Values indicate membership (0/1 or weights).
#'   }
#'   NULL = no gene set annotations. Enables ORA and GSEA analyses in viewer.
#' @param stringDB Character vector of length \code{nrow(expr)}. Protein/gene identifiers
#'   compatible with STRING database queries (e.g., Ensembl protein IDs, gene names).
#'   NULL = STRING network analysis disabled.
#' @param surv Survival data in one of three formats:
#'   \itemize{
#'     \item Vector of length \code{ncol(expr)}: single survival time with censoring indicated
#'           by "+" suffix (e.g., "120+", "45").
#'     \item Matrix/data.frame: multiple survival endpoints with samples in rows. Column names
#'           will be prefixed with "Surv|all|". Values must be numeric with optional "+" suffix.
#'   }
#'   NULL = no survival analysis.
#' @param SummarizedExperiment Logical. If TRUE, returns a \code{SummarizedExperiment} object;
#'   if FALSE, returns an \code{ExpressionSet}. Default: TRUE.
#' @param ... Additional arguments passed to \code{\link{t.test}}, such as \code{paired = TRUE}
#'   for paired t-tests or \code{var.equal = TRUE} for equal variance assumption.
#'
#' @return
#' A \code{SummarizedExperiment} or \code{ExpressionSet} object ready for visualization with
#' \code{\link{omicsViewer}}. The object includes:
#' \itemize{
#'   \item Expression matrix (and optionally imputed matrix)
#'   \item Enhanced metadata with PCA results, t-test statistics, rankings
#'   \item Gene set annotations (as attributes)
#'   \item Default axis selections (as attributes: "sx", "sy", "fx", "fy")
#' }
#'
#' @details
#' The function performs the following processing steps:
#' \enumerate{
#'   \item Validates dimensions and ensures unique row/column names
#'   \item Standardizes column names by prefixing with data type (e.g., "General|All|")
#'   \item Performs PCA on expression data (with and without imputation)
#'   \item Conducts statistical tests (t-tests) between specified groups
#'   \item Computes feature rankings across samples
#'   \item Integrates gene set, STRING, and survival annotations
#'   \item Sets sensible default axes for visualization
#' }
#'
#' All metadata columns are prefixed with standardized headers following the pattern
#' "Category|Subcategory|Variable" to organize variables in the viewer interface.
#'
#' @importFrom Biobase ExpressionSet AnnotatedDataFrame
#' @export
#'
#' @seealso
#' \code{\link{omicsViewer}} for launching the viewer.
#' \code{\link{multi.t.test}} for details on t-test implementation.
#' \code{\link{gsAnnotIdList}} for gene set annotation formatting.
#'
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
