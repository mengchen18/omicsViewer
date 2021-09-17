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

#' Read the object of ExpressetSet to be visualized using ExpressionSetViewer
#' @description This function is similar to readRDS. The only difference is
#' this function will convert gene set data to data.frame format, if they are 
#' given in a matrix format. 
#' @return An object of ExpressionSet
#' @param x a name of the file to be read, passed to \link{readRDS}
#' @export
#' @examples 
#' file <- system.file("extdata/demo.RDS", package = "ExpressionSetViewer")
#' obj <- readESVObj(file)
readESVObj <- function(x) {
  x <- readRDS(x) 
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