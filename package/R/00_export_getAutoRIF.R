#' Get genes associated with search terms and AutoRIF annotations
#' @param term a character vector of terms want to search
#' @param filter whether the result should be filtered. The least frequently 
#' mentioned genes (most like 1 or 2 times) will be removed.
#' @note https://amp.pharm.mssm.edu/geneshot/
#' @references Alexander Lachmann, Brian M Schilder, Megan L Wojciechowicz, 
#'   Denis Torre, Maxim V Kuleshov, Alexandra B Keenan, Avi Ma’ayan, Geneshot: 
#'   search engine for ranking genes from arbitrary text queries, Nucleic Acids 
#'   Research, Volume 47, Issue W1, 02 July 2019, Pages W571–W577, 
#'   https://doi.org/10.1093/nar/gkz393
#' @export
#' @importFrom jsonlite read_json
#' @examples 
#' a <- getAutoRIF("mtor signaling")

getAutoRIF <- function(term, filter = TRUE) {  
  term <- gsub(" ", "%20", term)
  term <- paste(term, collapse = ",")
  qry <- sprintf("http://amp.pharm.mssm.edu/geneshot/api/search/auto/%s", term)
  r <- jsonlite::read_json(qry)
  v <- sapply(r$gene_count, unlist)
  df <- data.frame(
    gene = colnames(v),
    n = v[1, ],
    perc = v[2, ],
    stringsAsFactors = FALSE
  )
  df$rank <- df$n * df$perc
  if (filter)
    df <- df[df$n > min(ceiling(nrow(df)/200), 3), ]
  attr(df, "term") <- r$search_term
  attr(df, "pubmedID_count") <- r$pubmedID_count
  df
}
