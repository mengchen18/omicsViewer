#' @description Vectorized over-representation analysis using fisher's exact test
#' @param gs a matrix where rows are features and columns are gene sets. The zeros in the matrix 
#'   indicate a gene is not a memeber of the gene sets in the column. 
#' @param i the row number selected from the gs
#' @param background An integer to indicate the size of background; if NA, background = nrow(gs)
#' @param minOverlap the minimum required overlap between gene set and gene list, if the overlap is lower
#'  than this value, no test would be done on this gene set
#' @param minSize the minimum size of gene sets should be tested
#' @param maxSize the maximum size of genesets should be tested
#' @param gs_desc description of gene sets, a name character vector storing the description of
#'  of the gene set, the names should be the same as colnames of gs
#' @param feature_desc feature description, e.g. gene symbols, names
#' @param unconditional.or logical. Whether calculate odds ratio using Maximum Likelihood Estimate (the sample odds ratio). 
#'  Note that the conditional Maximum Likelihood Estimate (cMLE, set this parameter to FALSE) is used in fisher.test. 
#' @param mtc.method multiple test correction methods, passed to p.adjust function
#' @param sort could be one of c("none", "p.value", "OR") to indicate how the result should be sorted. 
#' @importFrom fastmatch '%fin%'
#' @examples 
#' 1
#' # library(Biobase)
#' # dat <- readRDS("Dat/exampleEset.RDS")
#' # fd <- fData(dat)
#' # fdgs <- fd[, grep("^GS\\|", colnames(fd))]
#' # ir <- which(fd$`t-test|OV_BR|pval` < 0.05)
#' # res <- vectORA(gs = fdgs, i = ir, unconditional.or = FALSE)
#' # res2 <- vectORA(gs = fdgs, i = ir, unconditional.or = TRUE)

vectORA <- function(
  gs, i, background=NA, minOverlap = 3, minSize=5, maxSize=Inf, gs_desc = NULL, feature_desc = NULL,
  unconditional.or = TRUE, mtc.method = "fdr", sort = c("none", "p.value", "OR")[1]
  ) {
  
  if (is.na(background))
    background <- nrow(gs)
  mat <- gs != 0
  mat[is.na(mat)] <- FALSE
  ncand <- colSums(mat != 0)
  ic <- which(ncand >= minSize & ncand <= maxSize)
  if (length(ic) == 0) {
    message("No gene set in range c(minSize, maxSize), return NULL!")
    return(NULL)
  }
  mat <- mat[, ic, drop = FALSE]
  
  nol <- colSums(mat[i, ])
  ic <- which(nol >=  minOverlap)
  if (length(ic) == 0) {
    message("No gene set have overlap >= minOverlap, return NULL!")
    return(NULL)
  }
  mat <- mat[, ic, drop = FALSE]
  nol <- nol[ic]
  
  bdf <- vectORA.core(
    n.overlap = nol, 
    n.de = rep(length(i), ncol(mat)), 
    n.gs = colSums(mat), 
    n.bkg = background, 
    unconditional.or = unconditional.or, mtc.method = mtc.method)
  
  gs_annot_x <- ""
  if (!is.null(gs_desc))
    gs_annot_x <- gs_desc[colnames(mat)]
  if (is.null(feature_desc))
    feature_desc <- rownames(gs)
  fi <- feature_desc[i]
  overlap <- lapply(seq_len(ncol(mat)), function(j) {
    x <- mat[i, j, drop = FALSE]
    fi[x != 0]
  })
  
  rs <- cbind(
    pathway = colnames(mat), desc = gs_annot_x, bdf
  )
  rs$overlap_ids <- overlap
  sort <- sort[1]
  if (sort == "p.value") {
    rs <- rs[order(rs$p.value, decreasing = FALSE), ]
  } else if (sort == "OR") {
    rs <- rs[order(rs$OR, decreasing = TRUE), ]
  } else if (sort != "none") 
    warning("Unknown sort method, the results are not sorted!")
  rs
}

#' @description Internal function used by vectORA 
#' @param n.overlap the number of overlap between de and gs. The number of white balls drawn 
#'  without replacement from an urn which contains both black and white balls (compared to hyper).
#' @param n.de the number of DE gene. The number of balls drawn from the urn (compared to hyper).
#' @param n.gs the size of gene set. The number of white balls in the urn (compared to hyper).
#' @param n.bkg the background size.
#' @param unconditional.or calculate odds ratio using Maximum Likelihood Estimate (the sample odds ratio). 
#'  Note that the conditional Maximum Likelihood Estimate (MLE) is used in fisher.test. 
#' @param mtc.method multiple test correction methods, passed to p.adjust function
#' @importFrom stats phyper dhyper
#' @examples
#' 1
#' # xq <- rbind(c(4, 2, 4),
#' #             c(20, 40, 10),
#' #             c(11, 234, 10),
#' #             c(200, 1000, 100))
#' # 
#' # vectORA.core(xq[1, ], xq[2, ], xq[3, ], xq[4, ])
#' # vectORA.core(xq[1, ], xq[2, ], xq[3, ], xq[4, ], unconditional.or = TRUE)
#' # 
#' ## fisher's test
#' # t(apply(xq, 2, function(x1) {
#' #   v <- fisher.test(rbind(c(x1[1], x1[2]-x1[1]), c(x1[3]-x1[1], x1[4] - x1[2] - x1[3] + x1[1])), alternative = "greater")
#' #   c(p.value = v$p.value, v$estimate)
#' # }))


vectORA.core <- function(n.overlap, n.de, n.gs, n.bkg, unconditional.or = TRUE, mtc.method = "fdr") {
  
  pval <- phyper(q = n.overlap-1, m = n.gs, n = n.bkg - n.gs, k = n.de, lower.tail = FALSE)
  if(length(pval) == 0)
    return(
      data.frame(
        p.value = numeric(0),
        p.adjusted = numeric(0),
        OR = numeric(0),
        size_overlap = numeric(0),
        size_geneset = numeric(0),
        size_input = numeric(0),
        size_backgroung = numeric(0),
        stringsAsFactors = FALSE
      )
    )
  
  if (unconditional.or)
    or <- (n.overlap/(n.de - n.overlap))/((n.gs-n.overlap)/(n.bkg - n.gs - n.de + n.overlap)) else {
      
      or <- function(n.overlap, n.gs, n.de, n.bkg) {        
        m <- n.gs
        n <- n.bkg - n.gs
        k <- n.de
        x <- n.overlap
        lo <- pmax(0L, k - n)
        hi <- pmin(k, m)
        
        supportl <- mapply(":", lo, hi, SIMPLIFY = FALSE)
        
        vapply(seq_len(length(x)), function(i) {
          support <- supportl[[i]]
          logdc <- dhyper(support, m[i], n[i], k[i], log = TRUE)
          
          dnhyper <- function(ncp) {
            d <- logdc + log(ncp) * support
            d <- exp(d - max(d))
            d/sum(d)
          }
          mnhyper <- function(ncp) {
            if (ncp == 0) 
              return(lo[i])
            if (ncp == Inf) 
              return(hi[i])
            sum(support * dnhyper(ncp))
          }
          mle <- function(x) {
            if (x == lo[i]) 
              return(0)
            if (x == hi[i]) 
              return(Inf)
            
            mu <- mnhyper(1)
            if (mu > x) 
              uniroot(function(t) mnhyper(t) - x, c(0, 1))$root
            else if (mu < x) 
              1/uniroot(function(t) mnhyper(1/t) - x, c(.Machine$double.eps, 1))$root
            else 1
          }
          mle(x[i])
        }, numeric(1))  
      }
      or <- or(n.overlap=n.overlap, n.gs=n.gs, n.de=n.de, n.bkg=n.bkg)
    }
  
  data.frame(
    p.value = pval,
    p.adjusted = p.adjust(pval, method = mtc.method),
    OR = or,
    size_overlap = n.overlap,
    size_geneset = n.gs,
    size_input = n.de,
    size_backgroung = n.bkg,
    stringsAsFactors = FALSE
  )
}





