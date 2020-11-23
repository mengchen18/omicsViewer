
function(l) {
  
  # term list
  nn <- names(l)
  tt <- strsplit(nn, " ")
  nterm <- sapply(tt, function(x) {
    if (length(x) == 1)
      return(x)
    l <- lapply(2:length(x), function(n) {
      apply(cbind(1:(length(x)-n+1), n:length(x)), 1, function(j) {
        paste(x[j[1]:j[2]], collapse = " ")
      })
    })
    c(x, unlist(l))
  })
  
  ## binary matrix
  list2mat <- function(x) {
    ai <- unique(unlist(x))
    mat <- sapply(x, function(xx) fastmatch::`%fin%`(ai, xx))
    apply(mat, 2, as.integer)
  }
  bmat <- list2mat(x = l)
  
  
  ## cluster bmat
  bd <- as.dist(philentropy::distance(t(bmat), method = "jaccard"))
  hcl <- hclust(bd)
  cls <- cutree(hcl, h = 0.5)
  df <- data.frame(name = names(l), cluster = cls)
  
  options(stringsAsFactors = FALSE)
  lk <- lapply(unique(df$cluster), function(i) {
    i <- which(df$cluster == i)
    if (length(i) == 1)
      return(NULL)
    l1 <- table(unlist(nterm[i]))
    l2 <- table(unlist(nterm[-i]))
    dn <- l2[names(l1)]
    dn[is.na(dn)] <- 0
    dn <- dn+3
    q <- l1+3/dn
    data.frame(q, stringsAsFactors = FALSE)
  })
  names(lk) <- unique(df$cluster)
  
  lk1 <- sapply(names(lk), function(x) {

    lx <- lk[[x]]
    if (is.null(lx))
      return(NULL)
    mm <- max(lx$Freq)
    i <- which(lx$Freq > 2 & lx$Freq == mm)
    if (length(i) == 0)
      return(NULL)
    cd <- as.character(lx$Var1[i])
    cbind(x, cd[which.max(nchar(cd))], mm)
  })
  lk1 <- do.call(rbind, lk1)
  
  
  
}

