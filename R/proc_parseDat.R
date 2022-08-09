#' Extract function annotation from uniprot .dat file
#' @param file the .dat or .dat.gz file
#' @param outputDir dir of output file
#' @param ... other parameters passed to readLines
#' @import stringr
#' @return a data.frame parse from .dat file
#' 
parseDatTerm <- function(file, outputDir = NULL, ...) {
  message('Reading dat file ...')
  d0 <- readLines(file, ...)
  d0 <- split(d0, cumsum(d0 == "//"))
  
  org <- trimws(sub("OS", "", grep("^OS", d0[[1]], value = TRUE)))
  org <- make.names(paste(org, collapse = ""))
  while(grepl("\\.\\.", org)) org <- gsub("\\.\\.", ".", org)
  
  fn <- basename(file)
  vn <- paste0("_", gsub("-", "", Sys.Date()), "_", org, ".annot", sep = "")
  
  if (is.null(outputDir))
    outputDir <- dirname(file)
  outputFile <- file.path(outputDir, sub("(.dat|.dat.gz)$", vn, fn))
  
  message('Processing ...')
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
    message('Writing table ...')
    write.table(dd, file = outputFile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  }
  invisible(dd)
}