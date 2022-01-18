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
#' @param data.frame logical; whether to organize the result into \code{data.frame} format, 
#'   see "Value" section. 
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
#' 
#' terms <- data.frame(
#' id = c("ID1", "ID2", "ID1", "ID2", "ID8", "ID10", "ID4", "ID4"),
#' term = c("T1", "T1", "T2", "T2", "T2", "T2", "T1", "T2"),
#' stringsAsFactors = FALSE
#' )
#' features <- list(F1 = c("ID1", "ID2", "ID4"), F2 = c("ID13"), F3 = c("ID4", "ID8", "ID10"))
#' gsAnnotIdList(features, gsIdMap = terms, data.frame = TRUE, minSize = 1)
#' gsAnnotIdList(features, gsIdMap = terms, data.frame = FALSE, minSize = 1)
#' 
#' @return A binary matrix (if \code{data.frame} = FALSE),
#' the number of rows is the same with length of idList, the columns 
#' are the annotated gene set; or a \code{data.frame} (if \code{data.frame} = 
#' TRUE) with three columns: featureId, gsId, weight.
#' 
#' 
gsAnnotIdList <- function(idList, gsIdMap, minSize = 5, maxSize = 500, data.frame = FALSE) {
  pid <- data.frame(
    id = unlist(idList),
    index = rep(seq_along(idList), times = vapply(idList, length, integer(1))),
    stringsAsFactors = FALSE
  )
  gsIdMap <- gsIdMap[gsIdMap$id %in% pid$id, ]
  
  if (data.frame) {
    v <- split(pid$index, pid$id)[gsIdMap$id]
    v <- data.frame(
      featureId = seq_along(idList)[unlist(v)],
      gsId = rep(gsIdMap$term, vapply(v, length, integer(1))),
      weight = 1,
      stringsAsFactors = FALSE
    )
    v <- unique(v)
    n <- table(v$gsId)
    exc <- names(n[n > maxSize | n < minSize])
    v <- v[!v$gsId %fin% exc, ]
  } else {
    gsIdMap <- split(gsIdMap$id, gsIdMap$term)
    v <- vapply(gsIdMap, function(x, s) {
      s[unique(pid$index[fastmatch::"%fin%"(pid$id, x)])] <- 1
      s
    }, s = rep(0, length(idList)), FUN.VALUE = numeric(length(idList)))
    colnames(v) <- names(gsIdMap)
    nm <- matrixStats::colSums2(v)
    v <- v[, which(nm >= minSize & nm <= maxSize)]
  }
  v
}

totall <- function(gsmat) {
  gs <- as.matrix(gsmat)
  gs[gs == 0] <- NA
  gs <- reshape2::melt(gs, na.rm = TRUE)
  colnames(gs) <- c("featureId", "gsId", "weight")
  gs
}

#' @importFrom Biobase fData
tallGS <- function(obj) {
  fd <- Biobase::fData(obj)
  scn <- str_split_fixed(colnames(fd), "\\|", n = 3)
  ir <- which(scn[, 1] == "GS")
  if (length(ir) == 0)
    return(obj)
  
  igs <- fd[, ir, drop = FALSE]
  colnames(igs) <- make.names(scn[ir, 3])
  gs <- totall(igs)
  fd <- fd[, -ir]
  attr(fd, "GS") <- gs
  Biobase::fData(obj) <- fd
  obj
}

#' Read the object of SummarizedExperiment or ExpressetSet to be visualized using omicsViewer
#' @description Similar to \code{readRDS}. It reads the object to R working environment and perform extra two things. 
#'   1. If the loaded data an class of \code{SummarizedExperiment}, it will be converted to \code{ExpressionSet};
#'   2. If the gene set annotatio is in matrix format, the gene set annotation is converted to \code{data.frame} format.
#' @return an object of class \code{ExpressionSet} or \code{SummarizedExperiment} to be visualzied.
#' @param x the path of an object of \code{SummarizedExperiment} or \code{ExpressionSet}, passed to \link{readRDS}
#' @importFrom Biobase pData fData fData<- pData<-
#' @importFrom methods as
#' @export
#' @examples 
#' file <- system.file("extdata/demo.RDS", package = "omicsViewer")
#' obj <- readESVObj(file)
#' 
readESVObj <- function(x) {
  
  asEsetWithAttr <- function(x) {
    if (inherits(x, "SummarizedExperiment")) {
      eset <- as(x, "ExpressionSet")
      colnames(pData(eset)) <- colnames(colData(x))
      colnames(fData(eset)) <- colnames(rowData(x))
      
      DFattrs <- c("rownames", "nrows", "listData", "elementType", "elementMetadata", "metadata", "class")
      for (i in setdiff(names(attributes(colData(x))), DFattrs)) 
        attr(pData(eset), i) <- attr(colData(x), i)
      for (i in setdiff(names(attributes(rowData(x))), DFattrs))
        attr(fData(eset), i) <- attr(rowData(x), i)
      SEattrs <- c("assays", "colData", "NAMES", "elementMetadata", "metadata", "class")
      for (i in setdiff(names(attributes(x)), SEattrs))
        attr(eset, i) <- attr(x, i)
    } else if (inherits(x, "ExpressionSet")) {
      eset <- x
    } else
      stop("x should be either an SummarizedExperiment or ExpressionSet")
    eset
  }
  
  x <- asEsetWithAttr( readRDS(x) )
  tallGS(x)
}

vectORATall <- function(
  gs, i, background, minOverlap = 2, minSize=2, maxSize=Inf, gs_desc = NULL, feature_desc = NULL,
  unconditional.or = TRUE, mtc.method = "fdr", sort = c("none", "p.value", "OR")[1]
) {
  
  if (!is.factor(gs$featureId) || !is.factor(gs$gsId))
    stop("gs should have two columns of FACTOR named as featureId and gsId!")
  
  cnt_gs <- table(gs$gsId)
  cnt_gs <- cnt_gs[cnt_gs >= minSize & cnt_gs <= maxSize]
  if (length(cnt_gs) == 0) {
    message("No gene set in range c(minSize, maxSize), return NULL!")
    return(NULL)
  }
  gs <- gs[gs$gsId %in% names(cnt_gs), ]
  
  gsi <- gs[gs$featureId %in% i, ]
  cnt_gsi <- table(gsi$gsId)
  cnt_gsi <- cnt_gsi[cnt_gsi >= minOverlap]
  if (length(cnt_gsi) == 0) {
    message("No gene set have overlap >= minOverlap, return NULL!")
    return(NULL)
  }
  gsi <- gsi[gsi$gsId %in% names(cnt_gsi), ]
  
  bdf <- vectORA.core(
    n.overlap = c(cnt_gsi), 
    n.de = rep(length(i), length(cnt_gsi)), 
    n.gs = c(cnt_gs[names(cnt_gsi)]), 
    n.bkg = background, 
    unconditional.or = unconditional.or, mtc.method = mtc.method)
  
  gs_annot_x <- ""
  if (!is.null(gs_desc))
    gs_annot_x <- gs_desc[names(cnt_gsi)]
  if (is.null(feature_desc))
    feature_desc <- as.character(gsi$featureId)
  rs <- cbind(
    pathway = names(cnt_gsi), 
    desc = gs_annot_x, 
    bdf
  )
  overlap <- split(feature_desc, gsi$gsId)[as.character(rs$pathway)]
  rs$overlap_ids <- unname(overlap)
  sort <- sort[1]
  if (sort == "p.value") {
    rs <- rs[order(rs$p.value, decreasing = FALSE), ]
  } else if (sort == "OR") {
    rs <- rs[order(rs$OR, decreasing = TRUE), ]
  } else if (sort != "none") 
    warning("Unknown sort method, the results are not sorted!")
  rs
}