#' Extract function annotation from uniprot .dat file
#' @param file the .dat or .dat.gz file
#' @param outputDir dir of output file
#' @param ... other parameters passed to readLines
#' @import stringr
#' @return a data.frame parse from .dat file
#' 
parseDatTerm <- function(file, outputDir = NULL, ...) {
  cat('Reading dat file ...\n')
  d0 <- readLines(file, ...)
  d0 <- split(d0, cumsum(d0 == "//"))
  
  org <- make.names(trimws(sub("OS", "", grep("^OS", d0[[1]], value = TRUE))))
  while(grepl("\\.\\.", org)) org <- gsub("\\.\\.", ".", org)
  
  fn <- basename(file)
  vn <- paste0("_", gsub("-", "", Sys.Date()), "_", org, "annot", sep = "")
  
  if (is.null(outputDir))
    outputDir <- dirname(file)
  outputFile <- file.path(outputDir, sub(".dat$", vn, fn))
  
  cat('Processing ...\n')
  dd <- lapply(d0, function(x) {
    
    dr <- grep("^DR", x, value = TRUE)
    an <- stringr::str_split_fixed(trimws(sub("^DR", "", dr)), ";", 4)
    if (nrow(an) == 0)
      return(NULL)
    
    ac <- stringr::str_split_fixed(trimws(sub("^AC", "", grep("^AC", x, value = TRUE))), ";", n = 2)[1]
    name <- stringr::str_split_fixed(trimws(sub("^ID", "", grep("^ID", x, value = TRUE))), " ", 2)[1]
    gn <- trimws(sub("^GN", "", grep("^GN", x, value = TRUE)))
    gn <- grep("Name=", gn, value = TRUE)
    if (length(gn) == 0)
      gn <- NA else {
        gn <- gsub("Name=|;$", "", strsplit(gn[1], " ")[[1]][1])
      }
    data.frame(
      ID = name,
      ACC = ac,
      geneName = gn,
      source = an[, 1],
      term = an[, 2],
      desc = an[, 3],
      stringsAsFactors = FALSE
    )
  })
  dd <- do.call(rbind, dd)
  uid <- paste(dd$source, dd$term)
  tb <- table(uid)
  i <- which(uid %in% names(tb[tb >= 5]) & uid %in% names(tb[tb < 0.25*length(d0)]))
  dd <- dd[i, ]
  for(i in colnames(dd)) 
    dd[[i]][is.na(dd[[i]])] <- "_NA_"
  if (!is.null(dd)) {
    cat('Writing table ...\n')
    write.table(dd, file = outputFile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  }
  invisible(dd)
}

#' get uniprot reference proteome IDs
#' @return a character vector of UP ids
#' @importFrom stringr str_match
#' @importFrom curl curl
#' @describeIn downloadUPRefProteome get uniprot reference protein IDs
getUPRefProteomeID <- function(domain = c("Eukaryota", "Archaea", "Bacteria", "Viruses")[1]) {
  url <- curl(sprintf("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/%s/", domain))
  on.exit(close(url))
  con <- readLines(url)
  na.omit(stringr::str_match(con, "\\>\\s*(.*?)\\s*/\\<")[, 2])
}

#' get uniprot reference proteome IDs
#' @param id the UP id to download
#' @param destdir destination directory
#' @param domain the domain, one of "Eukaryota", "Archaea", "Bacteria" or "Viruses"
#' @return a character vector of UP ids
#' @importFrom stringr str_match
#' @importFrom curl curl
#' @importFrom utils download.file

downloadUPRefProteome <- function(
  id, domain = c("Eukaryota", "Archaea", "Bacteria", "Viruses")[1], destdir = "./"
) {
  url <- sprintf(
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/%s/%s", 
    domain, id)
  url2 <- curl(url)
  on.exit(close(url2))
  con <- readLines(url2)
  con <- na.omit(stringr::str_match(con, paste0(id, "_s*(.*?)\\s*.dat.gz"))[, 1])
  con <- con[!grepl("additional", con)]
  durl <- paste(url, con, sep = "/")
  download.file(durl, destfile = file.path(destdir, con))
  con
}

# 
getSteps <- function(fdata, stepList) {
  sts <- lapply(stepList, function(cc) {
    cc <- unlist(cc)
    pname <- cc[1]
    gname <- cc[-1]
    n <- length(gname)
    
    m1 <- fdata[[paste("mean", pname, gname[1], sep = "|")]]
    mn <- fdata[[paste("mean", pname, gname[n], sep = "|")]]
    mm <- abs(m1 - mn)
    
    fillFalse <- function(x) {
      x[is.na(x)] <- FALSE
      x
    }
    s <- lapply(seq_len(length(gname) - 1), function(i) {
      x <- gname[c(i, i+1)]
      ifdr <- sprintf("ttest|%s_vs_%s|fdr", x[1], x[2])
      imd <- sprintf("ttest|%s_vs_%s|mean.diff", x[1], x[2])
      list(
        signif = cbind(
          fdata[, ifdr] < 0.01,
          fdata[, ifdr] < 0.05,
          fdata[, ifdr] < 0.1),
        down = fillFalse(fdata[, imd] > 0),
        up = fillFalse(fdata[, imd] < 0)
      )
    })
    
    # careful with NA
    down <- lapply(s, "[[", "down")
    down <- Reduce("&", down)
    
    up <- lapply(s, "[[", "up")
    up <- Reduce("&", up)
    
    sig <- lapply(s, "[[", "signif")
    sig <- Reduce("+", sig)
    sig[!(up | down), ] <- 0
    sig[down, ] <- -sig[down, ]
    
    colnames(sig) <- paste("x.fdr", c("01", "05", "1"), sep = ".")
    res <- cbind(sig, "y.mean.diff" = mm)
    colnames(res) <- paste(paste(cc, collapse = ":"), colnames(res), sep = "|")
    res
  })
  r <- do.call(cbind, sts)
  colnames(r) <- paste("step", colnames(r), sep = "|")
  r
}

getOutliersUp <- function(x, na.replace = quantile(x, 0.25, na.rm = TRUE), threshold = 0.1) {
  sec <- na.replace
  f <- apply(x, 1, function(x, sec) {
    x <- sort(x, decreasing = TRUE)
    if (!is.na(x[2]))
      sec <- x[2]
    c(max(x[1] - sec, 0), names(x[1]))
  }, sec = sec)
  df <- data.frame(
    fold.change.log10 = as.numeric(f[1, ]),
    sample = f[2, ],
    stringsAsFactors = FALSE
  )
  df[which(abs(df$fold.change.log10) < threshold), ] <- NA
  df
}

getOutliersDown <- function(x, threshold = 0.1) {
  f <- apply(x, 1, function(x) {
    x <- sort(x, decreasing = FALSE)
    if (is.na(x[2]))
      return(c(NA, NA))
    c(min(x[1] - x[2], 0), names(x[1]))
  })
  df <- data.frame(
    "fold.change.log10" = as.numeric(f[1, ]),
    sample = f[2, ],
    stringsAsFactors = FALSE
  )
  df[which(abs(df$fold.change.log10) < threshold), ] <- NA
  df
}

#' Row-wise normalization of expression matrix with or without reference sample
#' @param x an expression matrix where rows are features, e.g. genes, proteins 
#'   and columns are samples. The values in the matrix are usually log transformed.
#' @param batch a factor or vector has the same length as \code{ncol(x)} to indicate
#'   the batch assignment of samples.
#' @param ref a logical vector has the same length as \code{ncol(x)} to indicated which
#'   columns are the common references among batches. If it is NULL (by default), 
#'   the mean of all channels will be used as batch reference. When NA present in the 
#'   reference channels, the mean values will be used in correction.
#' @param useMean logical; whether to use means of batches, usually set to TRUE 
#'   when no reference available
#' @import matrixStats
#' @return a matrix (hopefully without/with less batch effect)
#' @export
#' @examples 
#' e1 <- matrix(rnorm(5000), 500, 10)
#' e1[, 6:10] <- e1[, 6:10] + 3
#' boxplot(e1)
#' f <- rep(c("a", "b"), each = 5)
#' e2 <- rowshift(x = e1, batch = f)
#' boxplot(e2)

rowshift <- function(x, batch, ref=NULL, useMean = FALSE) {
  
  if (is.data.frame(x))
    x <- apply(x, 2, as.numeric)
  if (is.null(ref))
    ref <- seq_len(ncol(x))
  
  b_ref <- batch[ref]
  expr_ref <- x[, ref, drop=FALSE]
  grandmeans <- rowMeans(expr_ref, na.rm = TRUE)
  grandmeans2 <- rowMeans(x, na.rm = TRUE)
  
  for (i in unique(batch)) {
    off <- rowMeans(expr_ref[, b_ref == i, drop = FALSE], na.rm = TRUE) - grandmeans
    if (useMean) {
      off2 <- rowMeans(x[, batch == i, drop = FALSE], na.rm = TRUE) - grandmeans2
      ona <- is.na(off)
      off[ona] <- off2[ona]
    }
    x[, batch == i] <- x[, batch == i] - off
  }
  x
}

#' Normalization using n quantiles
#' @param x an expression matrix, usually log transformed
#' @param probs the quantiles to be aligned across samples. If \code{probs} is a length 1 numerical 
#'   vector, the quantiles will aligned. As a special case, \code{probs = 0.5}
#'   equals the median centering. If \code{probs}' length is > 1, a shift and scaling factor
#'   of samples will be calculating by fitting linear models using quantiles of samples, 
#'   the median and variance of samples will be corrected using the intersect and
#'   slope of the fitted model. 
#' @param shareFeature logocal; if TRUE, the normalization will be based on the 
#'   shared features between samples
#' @param ref the columns name or index to specify the reference sample, only used 
#'   when \code{shareFeature = TRUE}
#' @importFrom matrixStats colQuantiles
#' @export
#' @return a normalized matrix
#' @examples 
#' e1 <- matrix(rnorm(5000), 500, 10)
#' e1[, 6:10] <- 0.3 *e1[, 6:10] + 3
#' boxplot(e1)
#' # median centering, no variance correction
#' e2 <- normalize.nQuantiles(x = e1, probs = 0.5)
#' boxplot(e2)
#' # median centering + variance stablization
#' e3 <- normalize.nQuantiles(x = e1, probs = seq(0.25, 0.75, by = 0.1))
#' boxplot(e3)

normalize.nQuantiles <- function(x, probs = 0.5, shareFeature = FALSE, ref = 1) {
  
  if (is.data.frame(x))
    x <- apply(x, 2, as.numeric)
  
  if (shareFeature) {
    if (length(probs) == 1) {
      fac <- vapply(seq_len(ncol(x)), function(i) {
        ir <- !is.na(x[, ref]) & !is.na(x[, i])
        quantile(x[ir, ref], probs = probs) - quantile(x[ir, i], probs = probs)
      }, numeric(1))
      x <- sweep(x, 2, fac, "+")
    } else {
      for (i in seq_len(ncol(x))) {
        ir <- !is.na(x[, ref]) & !is.na(x[, i])
        refquant <- quantile(x[ir, ref], probs = probs, na.rm = TRUE)
        colquant <- quantile(x[ir, i], probs = probs, na.rm = TRUE)
        mod <- lm(refquant ~ colquant)
        x[, i] <- mod$coefficients[[2]] * x[, i] + mod$coefficients[[1]]
      }
    }
  } else {
    if (length(probs) == 1) {
      colquant <- colQuantiles(x, probs = probs, na.rm = TRUE)
      grandquant <- quantile(x, probs = probs, na.rm = TRUE)
      off <- colquant - grandquant
      x <- sweep(x, 2, off, "-")
    } else {
      colquant <- t(colQuantiles(x, probs = probs, na.rm = TRUE))
      grandquant <- quantile(x, probs = probs, na.rm = TRUE)
      for (i in seq_len(ncol(x))) {
        mod <- lm(grandquant ~ colquant[, i])
        x[, i] <- mod$coefficients[[2]] * x[, i] + mod$coefficients[[1]]
      }
    }
  }
  x
}

#' Normalize total sum
#' @param x a log10 transformed expression matrix
#' @importFrom matrixStats colQuantiles
#' @importFrom stats median
#' @return a normalized matrix
#' @export
#' @examples 
#' e1 <- matrix(rnorm(5000), 500, 10)
#' e1[, 6:10] <- e1[, 6:10]+3
#' boxplot(e1)
#' e2 <- normalize.totsum(x = e1)
#' boxplot(e2)

normalize.totsum <- function(x) { 
  x <- 10^x
  xs <- colSums(x, na.rm = TRUE)
  x <- sweep(x, 2, xs, '/')
  x <- x*median(xs)
  log10(x)
}



#' Removing variance of reference samples
#' @description This normalization removes the variance in reference samples. The 
#'   method do not need to specific the batch assignment but cannot work with data
#'   contains less than five common reference samples. A typical use of this normalization
#'   is to correct some drifting effect in mass spec based label free proteomics or 
#'   untargeted metabolomics experiment. Usually, this is a very strong normalization 
#'   should only be used with good reasons. 
#' @param x an expression matrix
#' @param ref the index of reference samples
#' @param positive logical; force only positive values in the resulted matrix
#' @param ... if given, \code{\link{normalize.nQuantiles}} will be called first, 
#'   the arguments here will be passed to \code{\link{normalize.nQuantiles}}
#' @export
#' @return a normalized matrix
#' @examples 
#' e1 <- matrix(rnorm(5000), 100, 50)+10
#' e2 <- removeVarQC(x = e1, ref = seq(5, 45, by = 10))
#' boxplot(e2)

removeVarQC <- function(x, ref, positive = TRUE, ...) {
  ls <- list(...)
  if (length(ls) > 0)
    x <- normalize.nQuantiles(x, ...)
  x0 <- x[, ref]
  decomp0 <- svd(x0)
  m0 <- decomp0$u %*% t(t(x) %*% decomp0$u)
  mm <- x - m0
  mm <- mm - rowMedians(mm)
  mm <- mm + rowMedians(x)
  if (positive)
    mm[which(mm<0)] <- 0
  mm
}

#' Normalized expression matrix
#' @description A wrapper function of all normalization methods, including row-wise or column-wise normalization. 
#' @param x an expression matrix where rows are features and columns are samples, usually log transformed. 
#' @param ref index of reference samples
#' @param batch batch factor
#' @param colWise column-wise normalization method to use, see \code{\link{normalizeColWise}}
#' @param rowWise row-wise normalization method to used{}
#'   \code{Reference} - using \code{\link{removeVarQC}} method
#'   \code{Batch mean} - using \code{\link{rowshift}} method without reference samples
#'   \code{Batch reference} - using \code{\link{rowshift}} method with reference samples
#' @export
#' @return a normalized matrix
#' @examples 
#' e1 <- matrix(rnorm(5000), 100, 50)+10
#' boxplot(e1)
#' e2 <- normalizeData(x = e1, ref = seq(5, 45, by = 10), rowWise = "Reference")
#' boxplot(e2)
#' 
normalizeData <- function(
  x,
  colWise = c("None", "Median centering", "Median centering (shared ID)", "Total sum", "median centering + variance stablization")[1], 
  rowWise = c("None", "Reference", "Batch mean", "Batch reference")[1],
  ref = NULL, batch = NULL
) {
  d <- normalizeColWise( x, method = colWise )
  if (rowWise == "Reference") {
    if (is.null(ref))
      stop("Reference not given!")
    ina <- is.na(d)
    d[ina] <- min(d, na.rm = TRUE)
    d <- removeVarQC(d, ref)
    d[ina] <- NA
  } else if (grepl("Batch", rowWise)) {
    if (is.null(batch))
      stop("Batch not given!")
    if (is.null(ref) && rowWise == "Batch reference")
      stop("Reference not given!")
    d <- rowshift(d, batch = batch, ref=ref, useMean = rowWise == "Batch mean")
  } 
  d
}

#' Column-wise normalization of expression matrix
#' @description A wrapper function of all column-wise normalization methods
#' @param x an expression matrix where rows are features and columns are samples, usually log transformed. 
#' @param method normalization method to use
#'  "Median centering" - median centering, see \code{\link{normalize.nQuantiles}}
#'  "Median centering (shared ID)" - median centering using shared features, see \code{\link{normalize.nQuantiles}}
#'  "Total sum" - total sum normalization
#'  "median centering + variance stablization" - 10 quantile normalization using 0.25, 0.3, ..., 0.75, 
#'  see \code{\link{normalize.nQuantiles}}
#' @export
#' @return a normalized matrix
#' @examples 
#' e1 <- matrix(rnorm(5000), 100, 50)+10
#' boxplot(e1)
#' e2 <- normalizeColWise(x = e1, method = "Median centering")
#' boxplot(e2)
#' 
normalizeColWise <- function(
  x, method  = c("Median centering", "Median centering (shared ID)", "Total sum", "median centering + variance stablization")[1]
) {
  
  nr <- which.max(colSums(!is.na(x)))
  
  normfun <- switch (
    method,
    "None" = function(x) x,
    "Median centering" = function(x) normalize.nQuantiles(x, probs = 0.5, shareFeature = FALSE, ref = nr), 
    "Median centering (shared ID)" = function(x) normalize.nQuantiles(x, probs = 0.5, shareFeature = TRUE, ref = nr), 
    "Total sum" = normalize.totsum, 
    "median centering + variance stablization" = function(x) 
      normalize.nQuantiles(x, probs = seq(0.25, 0.75, by = 0.05), shareFeature = FALSE, ref = nr)
  )
  
  normfun(x)
}

#' Reading proteinGroup table of MaxQuant output
#' @description A convenience function to read the proteinGroups table of MaxQuant 
#'   output. The function organize the result into different tables, e.g. iBAQ. 
#' @param x the proteinGroup.txt file returned by MaxQuant search
#' @param quant the quantification method, LF or TMT
#' @return a list of tables extracted from proteinGroups.txt file
read.proteinGroups <- function(x, quant = c("LF", "TMT")[1]) {
  func <- read.proteinGroups.lf
  getName <- function(x) {    
    nm <- colnames(x$iBAQ)
    if (is.null(nm))
      nm <- colnames(x[[grep("LFQ", names(x))]])
    if (is.null(nm))
      nm <- colnames(x$Intensity)
    nm
  }
  if (quant == "TMT") {
    getName <- function(x) {
      colnames(v$Reporter.intensity.corrected)
    }
    func <- read.proteinGroups.tmt
  }
  
  v <- func(x)
  attr(v, "label") <- getName(v)
  
  id <- str_split_fixed(v$annot$Protein.IDs, pattern = ";", 2)[, 1]
  gn <- str_split_fixed(v$annot$Gene.names, pattern = ";", 2)[, 1]
  rn <- make.names(paste(gn, id, sep = "_"))
  ss <- names(which(unlist(vapply(v, nrow, FUN.VALUE = integer(1))) == nrow(v$annot)))
  for (i in ss)
    rownames(v[[i]]) <- rn
  v
}

#' Read protein groups output of maxquant output and split
#'   it to columns
#' @param file Maxquant proteinGroup.txt file path
#' @importFrom utils read.delim
#' @return a list of tables extracted from proteinGroups.txt file

read.proteinGroups.lf <- function(file) {
  pg <- read.delim(file, stringsAsFactors = FALSE)
  df <- data.frame(val = c("iBAQ.",
                           "LFQ.intensity.", 
                           "Peptides.", 
                           "Razor...unique.peptides.",
                           "Unique.peptides.",
                           "Sequence.coverage.",
                           "Intensity.",
                           "MS.MS.Count.",
                           "MS.MS.count.",
                           "Identification.type."
  ), 
  log = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), 
  stringsAsFactors = FALSE)
  
  vi <- vapply(df$val, function(x) length(grep(x, colnames(pg))) > 0, FUN.VALUE = logical(1))
  df <- df[vi, ]             
  
  i <- ! (grepl("^REV_", pg$Majority.protein.IDs) |
            grepl("^CON_", pg$Majority.protein.IDs) |
            pg$Only.identified.by.site == "+")
  
  annot <- pg[i, -grep(paste(df$val, collapse = "|"), colnames(pg))]
  
  getExpr <- function(x, type = "iBAQ.", log = TRUE, keep.row = NULL) {
    ic <- grep(type, colnames(x), ignore.case = FALSE, value = TRUE)
    ic <- setdiff(ic, "iBAQ.peptides")
    val <- apply(pg[, ic], 2, as.numeric)
    if (log) {
      val <- log10(val)
      val[is.infinite(val)] <- NA
    }
    if (!is.null(keep.row))
      val <- val[keep.row, ]
    colnames(val) <- gsub(type, "", colnames(val))
    val
  }
  
  ml <- mapply(function(val, log) getExpr(pg, type = val, log = log, keep.row = i), 
               val = df$val, log = df$log)
  names(ml) <- gsub("\\.$", "", df$val)  
  ml$annot <- annot
  
  if (!is.null(ml$iBAQ)) {
    i <- which(rowSums(ml$iBAQ, na.rm = TRUE) == 0)
    if (length(i) > 0)
      ml <- lapply(ml, function(x) x[-i, ])
    ml$iBAQ_mc <- sweep(ml$iBAQ, 2, matrixStats::colMedians(ml$iBAQ, na.rm = TRUE), "-") + median(ml$iBAQ, na.rm = TRUE) 
  }  
  ml
}


read.proteinGroups.tmt <- function(file, xref=NULL) {
  ab <- read.delim(file, stringsAsFactors = FALSE)
  
  ir <- c(grep("^CON_", ab$Majority.protein.IDs), 
          grep("^REV_", ab$Majority.protein.IDs), 
          which(ab$Only.identified.by.site == "+"))
  
  eSum <- c("Fraction", "Reporter.intensity.corrected", "Reporter.intensity", "Reporter.intensity.count")
  ls <- list()
  for (i in eSum) {
    gb <- grep(paste0(i, ".[0-9]*$"), colnames(ab), value = TRUE)
    ls[[i]] <- apply(ab[-ir, gb, drop = FALSE], 2, as.numeric)
    ab[gb] <- NULL
  }
  
  lsind <- list()
  eInd <- c("Reporter.intensity.corrected", "Reporter.intensity.count", "Reporter.intensity")
  for (i in eInd) {
    gb <- grep(i, colnames(ab), value = TRUE)
    lsind[[i]] <- apply(ab[-ir, gb, drop = FALSE], 2, as.numeric)
    ab[gb] <- NULL
  }
  
  lsind$Reporter.intensity.corrected.log10 <- log10(lsind$Reporter.intensity.corrected)
  lsind$Reporter.intensity.corrected.log10[is.infinite(lsind$Reporter.intensity.corrected.log10)] <- NA
  lsind$annot <- ab[-ir, ]
  lsind$Summed <- ls
  colnames(lsind$Reporter.intensity.corrected.log10) <- make.names(
    sub("Reporter.intensity.corrected.", "", colnames(lsind$Reporter.intensity.corrected.log10)) 
  )
  
  ec <- intersect(eSum, names(lsind))
  for (i in ec) {
    nn <- sub(i, "", colnames(lsind[[i]]))
    nn <- make.names(sub("^.", "", nn))
    colnames(lsind[[i]]) <- nn
  }
  
  fn <- make.names(xref$label)
  if (!is.null(xref)) {
    lab <- make.names(paste(xref$channel, xref$mix))
    if (!identical(lab, colnames(lsind[[ec[[1]]]])))
      stop("columne does not match!")
    for (ii in ec)
      colnames(lsind[[ii]]) <- fn
    if (!is.null(lsind$Reporter.intensity.corrected.log10))
      colnames(lsind$Reporter.intensity.corrected.log10) <- fn
    lsind$xref <- xref
  }
  lsind
}


#' @importFrom stringr str_split_fixed str_count
phenoTemplate <- function(label, quant = c("LF", "TMT", "undefined")[1]) {
  n <- min(str_count(label, "_"))+1
  tab <- data.frame(
    Label = label,
    stringsAsFactors = FALSE
  )  
  if (quant %in% c("TMT", "undefined")) {
    cn <- str_split_fixed(label, "\\.", 2)    
    tab$Batch <- cn[, 2]
    tab$Channel <- str_extract(cn[, 1], "(\\d)+")
  } else {
    tab$Batch <- "1"    
  }
  tab$Reference <- FALSE
  if (n > 1) {
    m <- str_split_fixed(label, pattern = "_", n = n)
    colnames(m) <- paste0("Var", seq_len(ncol(m)) )
    tab <- cbind(tab, m)
  }
  it <- which(vapply(tab, is.factor, logical(1)))
  if (length(it) > 0)
    tab[it] <- lapply(tab[it], as.character)
  tab
}

#' Filter out rows of expression matrix
#' @description The function is used to filter rows with values of low intensities or
#'   do not reproducible presented in replicates. 
#' @param x an expression matrix
#' @param max.quantile a single numerical value between (0, 1), if the row maximum
#'   is smaller than this quantile (calculated from the whole matrix), 
#'   the row will be removed.
#' @param max.value a single numerical value, if the the maximum value of a rwo is
#'   smaller than this value, the row will be removed. Only used if \code{max.quantile}
#'   is set to "NULL". 
#' @param var variables has the same length as the column number in \code{x} to 
#'   indicate which sample is from which group
#' @param min.rep the minimum number of replicate in at least one of the groups,
#'   if less than this value, the row will be removed.
#' @return a logical vector where the TRUE means row to keep 
#' @export
#' @examples 
#' e1 <- matrix(rnorm(5000, sd = 0.3), 500, 10) + rnorm(500)
#' f <- filterRow(x = e1, max.quantile = 0.25)
#' table(f)
#' 
filterRow <- function(x, max.quantile = NULL, max.value = NULL, var = NULL, min.rep = 2) {
  rmi <- rowMaxs(x, na.rm = TRUE)
  if (!is.null(max.value)) {
    qt <- max.value
  } else if (!is.null(max.quantile)) {
    qt <- quantile(x, na.rm = TRUE, probs = max.quantile)
  } else
    qt <- NULL
  
  f1 <- FALSE
  if (!is.null(qt))
    f1 <- rmi < qt
  
  f2 <- FALSE
  if (!is.null(var)) {
    sn <- vapply(unique(var), function(xx) {
      rowSums(!is.na(x[, var == xx, drop = FALSE]))
    }, double(nrow(x)))
    rsn <- rowMaxs(sn)
    f2 <- rsn < min.rep
  }
  !(f1 | f2)
}

#' MQ folder validator
#' Validate whether a folder is a MQ output folder
#' @param dir the directory to check
#' @details 
#'   from the root level, these files exist:
#'    mqpar.xml
#'    [[combined/]txt/]proteinGroups.txt
#' @return a list containing the info about MQ folder check


validMQFolder <- function(dir) {
  l <- list(valid = FALSE)
  i1 <- file.exists(file.path(dir, "mqpar.xml")) 
  i2 <- list.files(dir, pattern = "^proteinGroups.txt$", recursive = TRUE, full.names = TRUE)
  if (i1 & length(i2) == 1) {
    l$valid <- TRUE
    l$mqpar <- file.path(dir, "mqpar.xml")
    l$txt <- dirname(i2)
    l$filePath <- file.path(l$txt, "proteinGroups.txt")
  }
  l
}

#' Parse mqpar.xml file 
#' @description Getting the experimental informatione (TMT or label free) from mqpar.xml file.
#' @importFrom flatxml fxml_importXMLFlat
#' @param x the path to mqpar.xml file
#' @return a list of MQ paramters
#' 
getMQParams <- function(x) {
  
  mqpar <- fxml_importXMLFlat(x)
  
  # ftms <- which(mqpar$elem. %in% c("msmsParams", "Name") & mqpar$value. == "FTMS")
  # if (mqpar$level4[ftms+1] == "MatchTolerance")
  #   ms2Tol <- mqpar$value.[ftms+1] else
  #     ms2Tol <- NA
  method <- "LF"
  if (length(grep("TMT[0-9]*plex", mqpar$value.)) > 1)
    method <- "TMT"
  
  ifa <- which(mqpar$elem. == "fastaFilePath")
  if (length(ifa) == 0)
    ifa <- which(mqpar$elem. == "string" & mqpar$level2 == "fastaFiles")
  
  list(
    version = mqpar$value.[which(mqpar$elem. == "maxQuantVersion")],
    fasta = mqpar$value.[ifa],
    enzymes = na.omit(mqpar$value.[which(mqpar$level4 == "enzymes")]),
    varMod = na.omit(mqpar$value.[which(mqpar$level4 == "variableModifications")]),
    fixedMod = na.omit(mqpar$value.[which(mqpar$level4 == "fixedModifications")]),
    mainTol = mqpar$value.[which(mqpar$elem. == "mainSearchTol")],
    # ms2Tol = ms2Tol,
    MBR = mqpar$value.[which(mqpar$elem. == "matchBetweenRuns")],
    peptideFdr = mqpar$value.[which(mqpar$elem. == "peptideFdr")],
    proteinFdr = mqpar$value.[which(mqpar$elem. == "proteinFdr")],
    ibaq = mqpar$value.[which(mqpar$elem. == "ibaq")],
    lfqMode = mqpar$value.[which(mqpar$elem. == "lfqMode")],
    method = method
  )
}

#' @title Reading gene set .gmt file
#' @description Frequently the .gmt files are downloaed from MSigDB database
#' @param x the name/path of the gmt file to be read
#' @param id the id used in gene sets, if is not NA, it should be either "SYMBOL" 
#'   or "ENTREZ". Usually only used when reading the .gmt file downloaded from 
#'   MSigDB.
#' @param data.frame logical; whether to organize the data in \code{data.frame}
#'   format. Default is FALSE, a list will be returned.
#' @export
#' @examples 
#' file <- system.file("extdata", package = "omicsViewer")
#' file <- file.path(file, "geneset.gmt")
#' gs <- read_gmt(file)
#' @return a list or data frame of gene set. When data.frame = TRUE, the returned
#'   object is a \code{data.frame} with two columns: id and term. 

read_gmt <- function(x, id = NA, data.frame = FALSE) {
  
  x <- readLines(x)
  x <- strsplit(x, "\t")
  # names <- vapply(x, "[", 1, character(1))
  names <- vapply(x, function(xx) xx[1], character(1))
  x <- lapply(x, function(xx) {
    structure(xx[-(seq_len(2))], 
              name = xx[1],
              link = xx[2])
  })
  names(x) <- names
  
  if (!is.na(id)) {
    id <- match.arg(id, c("SYMBOL", "ENTREZ"))
  } else {
    if (any(grepl("[A-Z]", unlist(x[seq_len(10)]))))
      id <- "SYMBOL" else
        id <- "ENTREZ"
  }
  attr(x, "ID") <- id
  
  if (data.frame) {
    x <- data.frame(
      id = unlist(x),
      term = rep(names(x), vapply(x, length, integer(1))),
      stringsAsFactors = FALSE)
  }
  x
} 

# Adding extra columns to the pheno Data/colData or feature Data/rowData in \code{ExpressionSet} or \code{SummarizedExperiment}
#' @param object an object of \code{ExpressionSet-class}
#' @param newData a \code{data.frame} containing the data to be added
#' @param where where to add the extra columns, 
#'   should be one of "pData", "fData", "rowData" and "colData". 
#' @rdname extendMetaData
#' @importFrom Biobase pData fData
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom S4Vectors DataFrame
#' @export
#' @return an object of \code{ExpressionSet-class}
#' @note The attributes in the pheno data and feature data will be preserved
#' @examples
#' est <- Biobase::ExpressionSet(assayData=matrix(runif(1000), nrow=100, ncol=10))
#' Biobase::pData(est)
#' est <- extendMetaData(est, data.frame(letter = letters[1:10]), where = "pData")
#' Biobase::pData(est)
setGeneric("extendMetaData", function(object, newData, where) {
  standardGeneric("extendMetaData")
})

#' Add extra columns to the phenoData/colData or featureData/rowData in ExpressionSet/SummarizedExperiment
#' @rdname extendMetaData
#' @export
setMethod(
  "extendMetaData", 
  signature = c(object = "ExpressionSet", newData = "data.frame"), 
  function(object, newData, where = c("pData", "fData", "colData", "rowData")[1]) {
    
    .addPData <- function(x, y) {
      at <- names(attributes(x))
      at <- setdiff(at, c("names", "class", "row.names"))
      if (length(at) > 0) {
        oldAt <- lapply(at, function(an) attr(x, an))
        names(oldAt) <- at
      }
      x <- cbind(x, y)
      if (length(at) > 0) {
        for (i in at)
          attr(x, i) <- oldAt[[i]]
      }
      x
    }
    
    where <- match.arg(where, c("pData", "fData", "colData", "rowData"))
    if (where == "pData" || where == "colData") {
      d <- .addPData(pData(object), newData)
      Biobase::pData(object) <- d
    } else {
      d <- .addPData(fData(object), newData)
      Biobase::fData(object) <- d
    }
    object
  }
) 

#' Add extra columns to the phenoData/colData or featureData/rowData in ExpressionSet/SummarizedExperiment
#' @rdname extendMetaData
#' @export
setMethod(
  "extendMetaData", 
  signature = c(object = "SummarizedExperiment", newData = "data.frame"), 
  function(object, newData, where = c("pData", "fData", "colData", "rowData")[1]) {
    newData <-  DataFrame(newData, check.names = FALSE)
    extendMetaData(object, newData, where = where)
    })

#' Add extra columns to the phenoData/colData or featureData/rowData in ExpressionSet/SummarizedExperiment
#' @rdname extendMetaData
#' @export
setMethod(
  "extendMetaData", 
  signature = c(object = "SummarizedExperiment", newData = "DFrame"), 
  function(object, newData, where = c("pData", "fData", "colData", "rowData")[1]) {
    
    .addPData <- function(x, y) {
      at <- names(attributes(x))
      at <- setdiff(at, c("rownames", "nrows", "listData", "elementType", "elementMetadata", "metadata", "class"))
      if (length(at) > 0) {
        oldAt <- lapply(at, function(an) attr(x, an))
        names(oldAt) <- at
      }
      x <- cbind(x, y)
      if (length(at) > 0) {
        for (i in at)
          attr(x, i) <- oldAt[[i]]
      }
      x
    }    
    
    where <- match.arg(where, c("pData", "fData", "colData", "rowData"))
    if (where == "pData" || where == "colData") {
      d <- .addPData(colData(object), newData)
      SummarizedExperiment::colData(object) <- d
    } else {
      d <- .addPData(rowData(object), newData)
      SummarizedExperiment::rowData(object) <- d
    }
    object
  }
) 