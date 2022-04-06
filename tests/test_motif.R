library(unittest)

seqs <- c("PFERVGTATIDNLPT", "MCWPEYCHCLKPIAS", "HRPQTSVNMQDMDPC", "HRPQTSVNMQDMDPCUU")
aaFreq <- omicsViewer:::aaFreq
fq <- aaFreq(seqs[1:3])
ok(
  ut_cmp_equal(dim(fq), c(20, 15)),
  "aaFreq - dimension check"
  )
ok(
  ut_cmp_equal(colSums(fq), rep(1, 15)),
  "aaFreq - colSum check"
  )
ok(
  ut_cmp_error(aaFreq(seqs), "Different length in x!"),
  "aaFreq - length check"
  )


seqs <- c("PFERVGTATIDNLPT", "MCWPEYCHCLKPIAS", "HRPQTSVNMQDMDPC", "HRPQTSVNMQDMDPC", "HRPQTSVNMQDMDPC")
motifRF <- omicsViewer:::motifRF
sq <- motifRF(fg.seqs = seqs[3:5], bg.seqs = seqs, fg.pfm = NULL, bg.pfm=NULL)
ok(
  ut_cmp_equal(dim(sq), c(20, 15)),
  "motifRF - dim check"
  )
ok(
  ut_cmp_equal(unique(c(sq)), c(0, 1)),
  "motifRF - value check"
  )
ok(
  ut_cmp_equal(c(table(sq)), c(285, 15), check.attributes = FALSE),
  "motifRF - value count check"
  )