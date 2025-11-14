#' Get genes associated with search terms and AutoRIF annotations
#' @param term a character vector of terms want to search
#' @param rif either autorif or generif, see "https://maayanlab.cloud/geneshot/"
#' @param filter whether the result should be filtered. The least frequently 
#' mentioned genes (most like 1 or 2 times) will be removed.
#' @note https://amp.pharm.mssm.edu/geneshot/
#' @references Alexander Lachmann, Brian M Schilder, Megan L Wojciechowicz, 
#'   Denis Torre, Maxim V Kuleshov, Alexandra B Keenan, Avi Ma’ayan, Geneshot: 
#'   search engine for ranking genes from arbitrary text queries, Nucleic Acids 
#'   Research, Volume 47, Issue W1, 02 July 2019, Pages W571–W577, 
#'   https://doi.org/10.1093/nar/gkz393
#' @export
#' @references Alexander Lachmann, Brian M Schilder, Megan L Wojciechowicz, 
#'   Denis Torre, Maxim V Kuleshov, Alexandra B Keenan, Avi Ma’ayan, 
#'   Geneshot: search engine for ranking genes from arbitrary text queries, 
#'   Nucleic Acids Research, Volume 47, Issue W1, 02 July 2019, Pages W571–W577, 
#'   https://doi.org/10.1093/nar/gkz393
#' @importFrom httr POST content
#' @examples 
#' a <- getAutoRIF("mtor signaling")
#' @return a \code{data.frame} of 4 columns: gene, n, perc, rank. 

getAutoRIF <- function(term, rif = c("generif", "autorif")[1], filter = TRUE) {

  # Input validation
  if (is.null(term) || length(term) == 0 || all(nchar(trimws(term)) == 0)) {
    return(NULL)
  }

  term <- gsub(" ", "%20", term)
  term <- paste(term, collapse = ",")
  GENESHOT_URL <- 'https://maayanlab.cloud/geneshot/api/search'
  payload <- list("rif" = rif, "term" = term)

  # Use safe POST wrapper
  result <- safe_POST(
    url = GENESHOT_URL,
    body = payload,
    encode = "json",
    timeout = 30,
    api_name = "Geneshot"
  )

  # Check if request failed
  if (!result$success) {
    message("Geneshot API error: ", result$error)
    return(NULL)
  }

  # Extract response data
  r <- result$data

  # Validate response structure
  if (is.null(r) || !is.list(r) || is.null(r$gene_count)) {
    message("Geneshot returned unexpected data format.")
    return(NULL)
  }

  # Extract gene count data
  v <- tryCatch({
    vapply(r$gene_count, unlist, numeric(2))
  }, error = function(e) {
    message("Error parsing Geneshot results: ", conditionMessage(e))
    return(matrix(nrow = 0, ncol = 0))
  })

  if (length(v) == 0)
    return(NULL)

  df <- data.frame(
    gene = colnames(v),
    n = v[1, ],
    perc = v[2, ],
    stringsAsFactors = FALSE
  )
  df$rank <- df$n * df$perc

  if (filter)
    df <- df[df$n > min(ceiling(nrow(df)/AUTORIF_THRESHOLD_DIVISOR), AUTORIF_MIN_THRESHOLD), ]

  attr(df, "term") <- r$search_term
  attr(df, "pubmedID_count") <- r$pubmedID_count
  df
}