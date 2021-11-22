#' Extract function annotation from uniprot .dat file
#' @param file the .dat or .dat.gz file
#' @param outputDir dir of output file
#' @param ... other parameters passed to readLines
#' @import stringr
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
    s <- lapply(1:(length(gname) - 1), function(i) {
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

#' row-wise normalization
#' @param x a matrix where rows are genes/proteins and columns are samples, log10 transoformed intensity
#' @param batch a factor or vecter has the same length as ncol(x) to indicate
#'   the batch assignment of samples
#' @param ref a logical vector has the same length as ncol(x) to indicated which
#'   columns are the common references among batches. If it is NULL (by default), 
#'   the mean of all channels will be used as batch reference.
#' @param useMean logical; whether to use means of batches, usually set to TRUE when no reference available
#' @import matrixStats
#' @export

rowshift <- function(x, batch, ref=NULL, useMean = FALSE) {
  
  if (is.data.frame(x))
    x <- apply(x, 2, as.numeric)
  if (is.null(ref))
    ref <- 1:ncol(x)
  
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

#' normalize soft quantile
#' @param x the input matrix, the columns of which were centered on quantiles, log10 transformed intensity
#' @param probs the quantiles used to controls/align. If probs is a length 1 numerical 
#'   vector, the quantiles will aligned. If probs' length is > 1, the quantiles are calculated
#'   and linear models are used to estimated slope and intersect between quantiles.
#' @param sharedProtein normalization based on the shared proteins between reference 
#'   experiment each individual experiment
#' @param ref the columns name or index to specify the reference experiment, only used 
#'   when sharedProtein = TRUE
#' @importFrom matrixStats colQuantiles

normalize.softQuantile <- function(x, probs = 0.5, sharedProtein = FALSE, ref = 1) {
  
  if (is.data.frame(x))
    x <- apply(x, 2, as.numeric)
  
  if (sharedProtein) {
    if (length(probs) == 1) {
      fac <- sapply(1:ncol(x), function(i) {
        ir <- !is.na(x[, ref]) & !is.na(x[, i])
        quantile(x[ir, ref], probs = probs) - quantile(x[ir, i], probs = probs)
      })
      x <- sweep(x, 2, fac, "+")
    } else {
      for (i in 1:ncol(x)) {
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
      for (i in 1:ncol(x)) {
        mod <- lm(grandquant ~ colquant[, i])
        x[, i] <- mod$coefficients[[2]] * x[, i] + mod$coefficients[[1]]
      }
    }
  }
  x
}

#' normalize total sum
#' @param x log10 transformed intensity 
#' @importFrom matrixStats colQuantiles
#' @importFrom stats median
normalize.totsum <- function(x) { 
  x <- 10^x
  xs <- colSums(x, na.rm = TRUE)
  x <- sweep(x, 2, xs, '/')
  x <- x*median(xs)
  log10(x)
}



#' Removing variance of reference samples
#' @param x expression matrix
#' @param ref the index of reference samples
#' @param ... if given, other parameters passed to normalize.softQuantile
#' @export
#' 
removeVarQC <- function(x, ref, ...) {
  ls <- list(...)
  if (length(ls) > 0)
    x <- normalize.softQuantile(x, ...)
  x0 <- x[, ref]
  decomp0 <- svd(x0)
  m0 <- decomp0$u %*% t(t(x) %*% decomp0$u)
  mm <- x - m0
  mm <- mm - rowMedians(mm)
  mm <- mm + rowMedians(x)
  mm[which(mm<0)] <- 0
  mm
}

#' normalize data row-wise or column-wise
#' @param x expression matrix where rows are features and columns are samples
#' @param reference index of reference samples
#' @param batch batch factor
#' @param colWise column wise normalization
#' @param rowWise row wise normalization method 

normalizeData <- function(
  x,
  colWise = c("None", "Median centering", "Median centering (shared ID)", "Total sum", "median centering + variance stablization")[1], 
  rowWise = c("None", "Reference", "Batch mean", "Batch reference")[1],
  reference = NULL, batch = NULL
) {
  d <- normalizeColWise( x, method = colWise )
  if (rowWise == "Reference") {
    if (is.null(reference))
      stop("Reference not given!")
    ina <- is.na(d)
    d[ina] <- min(d, na.rm = TRUE)
    d <- removeVarQC(d, reference)
    d[ina] <- NA
  } else if (grepl("Batch", rowWise)) {
    if (is.null(batch))
      stop("Batch not given!")
    if (is.null(reference) && rowWise == "Batch reference")
      stop("Reference not given!")
    d <- rowshift(d, batch = batch, ref=reference, useMean = rowWise == "Batch mean")
  } 
  d
}

normalizeColWise <- function(
  x, method  = c("Median centering", "Median centering (shared ID)", "Total sum", "median centering + variance stablization")[1]
) {
  
  nr <- which.max(colSums(!is.na(x)))
  
  normfun <- switch (
    method,
    "None" = function(x) x,
    "Median centering" = function(x) normalize.softQuantile(x, probs = 0.5, sharedProtein = FALSE, ref = nr), 
    "Median centering (shared ID)" = function(x) normalize.softQuantile(x, probs = 0.5, sharedProtein = TRUE, ref = nr), 
    "Total sum" = normalize.totsum, 
    "median centering + variance stablization" = function(x) 
      normalize.softQuantile(x, probs = seq(0.25, 0.75, by = 0.05), sharedProtein = FALSE, ref = nr)
  )
  
  normfun(x)
}

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
  ss <- names(which(unlist(sapply(v, nrow)) == nrow(v$annot)))
  for (i in ss)
    rownames(v[[i]]) <- rn
  v
}

#' Read protein groups output of maxquant output and split
#'   it to columns
#' @param file Maxquant proteinGroup.txt file path
#' @importFrom utils read.delim

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
  log = c(T, T, F, F, F, F, F, F, F, F), 
  stringsAsFactors = FALSE)
  
  vi <- sapply(df$val, function(x) length(grep(x, colnames(pg))) > 0)
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
    colnames(m) <- paste0("Var", 1:ncol(m))
    tab <- cbind(tab, m)
  }
  it <- which(sapply(tab, is.factor))
  if (length(it) > 0)
    tab[it] <- lapply(tab[it], as.character)
  tab
}


#' @return a logical vector where the TRUE means row to keep 
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
    sn <- sapply(unique(var), function(xx) {
      rowSums(!is.na(x[, var == xx, drop = FALSE]))
    })
    rsn <- rowMaxs(sn)
    f2 <- rsn < min.rep
  }
  print(table(f1))
  print(table(f2))
  !(f1 | f2)
}

#' MQ folder validator
#' Validate whether a folder is a MQ output folder
#' @param dir the directory to check
#' @details 
#'   from the root level, these files exist:
#'    mqpar.xml
#'    [[combined/]txt/]proteinGroups.txt

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
