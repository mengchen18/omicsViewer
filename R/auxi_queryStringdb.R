#' @description Drawing network given network and gene set result
#' @param ntwk network result, often returned by "stringNetwork" function. 
#' @param gsa gene set result, often returned by "stringGSA" function
#' @param i row index of gsa, which should be highlighted in the network
#' @param label whether the point should be labelled in the network
#' @importFrom fastmatch fmatch
#' @importFrom networkD3 forceNetwork JS
stringD3Net <- function(ntwk, gsa, i, label = FALSE) {
  
  nd <- data.frame(name = unique(unlist(ntwk[, c('preferredName_A', 'preferredName_B')])), 
                   stringsAsFactors = FALSE)
  rownames(nd) <- nd$name
  nd$group <- 1
  i <- strsplit(gsa$preferredNames[i], ",")[[1]]
  nd$group[nd$name %in% i] <- 2
  links <- data.frame(
    source = fmatch(ntwk$preferredName_A, nd$name)-1,
    target = fmatch(ntwk$preferredName_B, nd$name)-1,
    value = (ntwk$score-0.4)^2*10
  )
  
  colorfunc <- networkD3::JS('colorfunc = function(i) { return i == 1 ? "#64A0C8" : "#E37222" };')
  lab <- ifelse(label, 1, 0)
  forceNetwork(Links = links, Nodes = nd, Source = "source", linkColour = "gray",
               Target = "target", Value = "value", NodeID = "name", charge = -5,
               Group = "group", colourScale = colorfunc, 
               opacity = 0.7, fontSize = 12, opacityNoHover = lab,
               zoom = TRUE, legend = TRUE)
}


#' @description Retrieve string network
#' @description get string network
#' @param genes the gene ids
#' @param taxid taxonomy ids
#' @param caller your identifier for string-db.org
#' @importFrom httr GET content

stringNetwork <- function(genes, taxid = 9606, caller = "omicsViewer") {
  
  string_api_url <- "https://string-db.org/api"
  output_format <- "tsv"
  method <- "network"
  if (inherits(genes, "list"))
    genes <- na.omit(unlist(genes))
  genes <- gsub("#|%", "", genes)
  
  request_url <- paste(string_api_url, output_format, method, sep = "/")
  
  params <- list(
    "identifiers" = paste(genes, collapse = "%0d"), 
    "species" = taxid, 
    "required_score" = 500,
    "caller_identity" = caller
  )
  params <- paste(
    vapply(names(params), function(x) paste(x, params[[x]], sep = "="), FUN.VALUE = character(1)),
    collapse = "&"
  )
  results <- GET( paste(request_url, params, sep = "?") )
  dd <- httr::content(results)
  dd <- unique(dd)
  if (!is.data.frame(dd)) {
    message('stringdb does not return valid results, return NULL!')
    return(dd)
  }
  as.data.frame(dd)
}


#' @description Mapping ids to string ids
#' @description the string ids can be used as background in the string enrichment analysis
#' @param genes a character vector of gene ids
#' @param taxid taxonomy ids
#' @param caller your identifier for string-db.org
#' @note https://string-db.org/cgi/help.pl?subpage=api%23mapping-identifiers
#' @importFrom httr GET content
#' @examples 
#' 1
#' # gg = c('P04637', 'P00533', 'P04626', "Q8IYB3", "O75494", "Q9Y696")
#' # getStringId(gg)

getStringId <- function(genes, taxid = 9606, caller = "omicsViewer") {
  string_api_url <- "https://string-db.org/api"
  output_format <- "tsv"
  method <- "get_string_ids"
  
  params <- list(
    "identifiers" = paste(genes, collapse = "%0d"),
    "species" = taxid, # species NCBI identifier 
    "limit" = 1, # only one (best) identifier per input protein
    "echo_query" = 1, # see your input identifiers in the output
    "caller_identity" = caller # your app name
  )
  params <- paste(
    vapply(names(params), function(x) paste(x, params[[x]], sep = "="), FUN.VALUE = character(1)),
    collapse = "&"
  )
  request_url <-  paste(string_api_url, output_format, method, sep = "/")
  ## Call STRING
  results <- httr::GET( paste(request_url, params, sep = "?") )
  httr::content(results)
}

#' @description Performing gene set analysis using string-db
#' @param genes a character vector of gene ids
#' @param taxid taxonomy ids (species NCBI identifier)
#' @param background the background genes
#' @param backgroundStringId logical value, whether the identifier given in background is 
#'   already stringid, only used when background != NULL. You can use \code{\link{getStringId}}
#'   to convet your id to string id.
#' @param caller your identifier for string-db.org
#' @note https://string-db.org/cgi/help.pl?subpage=api%23getting-functional-enrichment
#' @importFrom httr GET content
#' @examples
#' 1
#' # gg = c('P04637', 'P00533', 'P04626', "Q8IYB3", "O75494", "Q9Y696")
#' # u <- stringGSA(gg)

stringGSA <- function(genes, taxid = 9606, background = NULL, backgroundStringId = FALSE, caller = "omicsViewer") {
  
  string_api_url <- "https://string-db.org/api"
  output_format <- "tsv"
  method <- "enrichment"
  request_url <- paste(string_api_url, output_format, method, sep = "/")
  
  params <- list(
    "species" = taxid, 
    "caller_identity" = caller
  )
  
  if (!is.null(background)) {
    if (!backgroundStringId) {
      g <- getStringId(genes, taxid = taxid)
      if (inherits(g, "character"))
        return(g)
      sg <- g$stringId
      genes <- intersect(genes, g$queryItem)
    } else {
      sg <- background
      genes <- intersect(genes, background)
    }
    if (length(genes) == 0)
      stop("No genes exist in the current background!")
    params$background_string_identifiers <- paste(g$stringId, collapse = "%0d")
  }
  
  params$identifiers <- paste(genes, collapse = "%0d") 
  
  # call string
  response <- httr::GET(url = request_url, query = params) 
  dd <- httr::content(response) 
  if (!is.data.frame(dd)) {
    message('stringdb does not return valid results, return NULL!')
    return(dd)
  }
  as.data.frame(dd)
}


                      
                      
                      
                      