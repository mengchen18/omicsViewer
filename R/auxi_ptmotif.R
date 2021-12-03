#' @importFrom ggseqlogo ggseqlogo
#' @note seq of length 0 will be removed

aaFreq <- function(x) {
  aa <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  x <- x[which(nchar(x) > 0)]
  if (length(unique(nchar(x))) > 1)
    stop("Different length in x!")
  s0 <- strsplit(x, "|")
  s0 <- do.call(rbind, s0)
  fg <- apply(s0, 2, function(x) table(x)[aa])
  fg[is.na(fg)] <- 0  
  s <- sweep(fg, 2, colSums(fg, na.rm = TRUE), "/")
  rownames(s) <- aa
  s
}


motifRF <- function(fg.seqs, bg.seqs, fg.pfm = NULL, bg.pfm=NULL) {
  aa <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  if (!is.null(fg.pfm))
    fg <- fg.pfm else
      fg <- aaFreq(fg.seqs)
  if (!is.null(bg.pfm))
    bk <- bg.pfm else
      bk <- aaFreq(bg.seqs)
  motif <- fg/bk
  motif[is.na(motif)] <- 0
  motif <- sweep(motif, 2, colSums(motif), "/")
  rownames(motif) <- aa
  motif
}
