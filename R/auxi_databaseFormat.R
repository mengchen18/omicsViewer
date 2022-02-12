#' Save the xcmsViewer result object as sqlite database
#' @param obj an object of class ExpressionSet or SummarizedExperiment
#' @param db.file a character indicate file name of the database file 
#' @param overwrite logical. whether the database should be overwritten if exist already.
#' @return the directory where the database saved
#' @name saveOmicsViewerDb
#' @importFrom RSQLite SQLite dbConnect dbDisconnect dbGetQuery dbWriteTable dbListTables
#' @examples 
#' f <- system.file("extdata", "demo.RDS", package = "omicsViewer")
#' es <- readRDS(f)
#' # The following line will write a database file on your disk
#' # saveOmicsViewerDb(es, db.file = "./omicsViewerData.db")
#' @export

setGeneric("saveOmicsViewerDb", function(obj, db.file, overwrite = TRUE) {
  standardGeneric("saveOmicsViewerDb")
})

#' @rdname saveOmicsViewerDb
setMethod(
  "saveOmicsViewerDb", 
  signature(obj = "SummarizedExperiment", db.file = "character"),
  function(obj, db.file, overwrite = TRUE) {
    obj <- asEsetWithAttr(obj)
    saveOmicsViewerDb(obj = obj, db.file = db.file, overwrite = overwrite)
  })

#' @rdname saveOmicsViewerDb
setMethod(
  "saveOmicsViewerDb", 
  signature(obj = "ExpressionSet", db.file = "character"),
  function(obj, db.file, overwrite = TRUE) {
    
    expr <- exprs(obj)
    pd <- pData(obj)
    fd <- fData(obj)
    expr_impute <- exprsImpute(obj)
    
    con <- dbConnect(RSQLite::SQLite(), db.file)
    on.exit(dbDisconnect(con))
    
    dfrname <- function(x) {
      data.frame(x, rowname = rownames(x), check.names = FALSE)
    }
    
    dbWriteTable(conn = con, name = "exprs", value = dfrname(expr), overwrite = overwrite)
    dbWriteTable(conn = con, name = "feature", value = dfrname(fd), overwrite = overwrite)
    dbWriteTable(conn = con, name = "sample", value = dfrname(pd), overwrite = overwrite)
    if (!is.null(expr_impute))
      dbWriteTable(conn = con, name = "exprsimpute", value = dfrname(expr_impute), overwrite = overwrite)
    
    gs <- data.frame(featureId = character(0), gsId = character(0), weight = numeric(0))
    if (hasAttr(fd, "GS"))
      gs <- attr(fd, "GS")
    if (inherits(gs, c("dgCMatrix", "lgCMatrix")))
      gs <- csc2list(gs)
    dbWriteTable(conn = con, name = "GS", value = gs, overwrite = overwrite)
    
    defaultxy <- data.frame(axis = character(0), value = character(0))
    ax <- c("fx", "fy", "sx", "sy")
    ax <- ax[hasAttr(obj, attr.name = ax)]
    if (length(ax) > 0)
      defaultxy <- data.frame(
        axis = ax,
        value = vapply(ax, function(v) attr(obj, v), character(1))
      )
    dbWriteTable(conn = con, name = "axes", value = defaultxy, overwrite = overwrite)
    
    dend <- data.frame(name = character(0), object = character(0), dim = integer(0))
    if (hasAttr(obj, "dendrogram")) {
      dend <- attr(obj, "dendrogram")
      dend <- data.frame(
        name = names(dend),
        object = vapply(dend, function(x) hclust2str(as.hclust(x$hcl)), character(1)),
        dim = 1,
        row.names = NULL
      )
    }
    dbWriteTable(conn = con, name = "dendrogram", value = dend, overwrite = overwrite)
    dirname(db.file)
  }
)

#' Check whether an object has an attribute
#' @param x the object
#' @param attr.name a character vector containing the name of attributes to be checked
#' @return a logical value/vector has the same length as attr.name
#' 
hasAttr <- function(x, attr.name) {
  attr.name %in% names(attributes(x))
}

#' Convert hclust object to/from single character
#' @param x a character of length one or an hclust object
#' @note The $call element in hclust will not retained in the conversion. The conversion
#'   decrease the precision in $height element. 
#' @name hclust2str
#' @return a character stores the hclust object
#' @examples 
#' # not for end users
#' # m <- matrix(rnorm(50), 25)
#' # hc <- hclust(dist(m))
#' # plot(hc)
#' # te <- hclust2str(hc)
#' # hc2 <- str2hclust(te)
#' # plot(hc2)
#' 
hclust2str <- function(x) {
  cc <- c(
    merge = paste(paste(x$merge[, 1], x$merge[, 2], sep = "\t"), collapse = "\n"),
    height = paste(signif(x$height, digits = 5), collapse = ";"),
    order = paste(x$order, collapse = ";"),
    labels = paste(x$labels, collapse = "=;="),
    method = x$method,
    dist.method = x$dist.method
  )
  paste(cc, collapse = "__elementSplitter__")
}

#' @rdname hclust2str
#' @return a hclust object
str2hclust <- function(x) {
  x <- strsplit(x, split = "__elementSplitter__")[[1]]
  ll <- list(
    merge = apply(
      read.delim(
        file = textConnection(x[1]), stringsAsFactors = FALSE, header = FALSE
      ), 2, as.integer), 
    height = as.numeric(strsplit(x[2], split = ";")[[1]]),
    order = as.integer(strsplit(x[3], split = ";")[[1]]),
    labels = strsplit(x[4], split = "=;=")[[1]],
    call = list("convertedHclustObj", "str2hclust"),
    method = x[5],
    dist.method = x[6]
  )
  if (length(ll$labels) == 0)
    ll$labels <- NULL
  class(ll) <- "hclust"
  ll
}


#' Convert SummarizedExperiment to ExpressionSet retaining all attributes
#' @param x an object of class SummarizedExperiment
#' @return an object of class ExpressionSet
asEsetWithAttr <- function(x) {
  if (inherits(x, "SummarizedExperiment")) {
    eset <- as(x, "ExpressionSet")
    colnames(pData(eset)) <- colnames(colData(x))
    colnames(fData(eset)) <- colnames(rowData(x))
    
    DFattrs <- c("rownames", "nrows", "listData", "elementType", "elementMetadata", "metadata", "class")
    for (i in setdiff(names(attributes(colData(x))), DFattrs)) 
      attr(pData(eset), i) <- attr(colData(x), i)
    for (i in setdiff(names(attributes(rowData(x))), DFattrs))
      attr(fData(eset), i) <- attr(rowData(x), i)
    SEattrs <- c("assays", "colData", "NAMES", "elementMetadata", "metadata", "class")
    for (i in setdiff(names(attributes(x)), SEattrs))
      attr(eset, i) <- attr(x, i)
  } else if (inherits(x, "ExpressionSet")) {
    eset <- x
  } else
    stop("x should be either an SummarizedExperiment or ExpressionSet")
  eset
}

#' Read the object of SummarizedExperiment or ExpressetSet to be visualized using omicsViewer
#' @description This function accept a path to a sqlite database or RDS object. If an RDS file to be read, 
#'   The function is similar to \code{readRDS}. It reads the object to R working environment and perform extra two things. 
#'   1. If the loaded data an class of \code{SummarizedExperiment}, it will be converted to \code{ExpressionSet};
#'   2. If the gene set annotatio is in matrix format, the gene set annotation is converted to \code{data.frame} format.
#' @return an object of class \code{ExpressionSet} or \code{SummarizedExperiment} to be visualzied.
#' @param x the path of an object of \code{SummarizedExperiment} or \code{ExpressionSet}, passed to \link{readRDS}
#' @importFrom Biobase pData fData fData<- pData<-
#' @importFrom methods as
#' @export
#' @examples 
#' file <- system.file("extdata/demo.RDS", package = "omicsViewer")
#' obj <- readESVObj(file)
#' 
readESVObj <- function(x) {  
  if (grepl(".RDS$", x, ignore.case = TRUE)) {
    x <- asEsetWithAttr( readRDS(x) )
    x <- tallGS(x)
  } else if (grepl("(.db|.sqlite|.sqlite3)$", x, ignore.case = TRUE)) {
    x <- dbConnect(RSQLite::SQLite(), x)
  } else
    stop("readESVObj: Unkown input format!")
  x
}

getExprs <- function(x) {
  if (inherits(x, "SQLiteConnection")) {
    mat <- dbGetQuery(x, "SELECT * FROM exprs;")
    rn <- mat$rowname
    mat$rowname <- NULL
    mat <- apply(mat, 2, as.numeric)
    rownames(mat) <- rn
  } else if (inherits(x, "ExpressionSet"))
    mat <- exprs(x)
  mat
} 

getExprsImpute <- function(x) {
  if (inherits(x, "SQLiteConnection")) {
    if (!"exprsimpute" %in% dbListTables(x))
      return(NULL)
    mat <- dbGetQuery(x, "SELECT * FROM exprsimpute;")
    rn <- mat$rowname
    mat$rowname <- NULL
    mat <- apply(mat, 2, as.numeric)
    rownames(mat) <- rn
  } else if (inherits(x, "ExpressionSet"))
    mat <- exprsImpute(x)
  mat
} 

getPData <- function(x) {
  if (inherits(x, "SQLiteConnection")) {
    mat <- dbGetQuery(x, "SELECT * FROM sample;")
    rownames(mat) <- mat$rowname
    mat$rowname <- NULL
  } else if (inherits(x, "ExpressionSet")) {
    mat <- pData(x)
  }
  mat
}

getFData <- function(x) {
  if (inherits(x, "SQLiteConnection")) {
    mat <- dbGetQuery(x, "SELECT * FROM feature;")
    rownames(mat) <- mat$rowname
    mat$rowname <- NULL
    gs <- dbGetQuery(x, "SELECT * FROM GS;")
    if (nrow(gs) > 0) {
      gs$featureId <- as.factor(gs$featureId)
      gs$gsId <- as.factor(gs$gsId)
      attr(mat, "GS") <- gs
    }
  } else if (inherits(x, "ExpressionSet")) {
    mat <- fData(x)
  }
  mat
}

getAx <- function(x, what) {
  if (inherits(x, "SQLiteConnection")) {
    if (what == "dendrogram") {
      v <- getDend(x)
    } else {
      v <- dbGetQuery(x, "SELECT value FROM axes WHERE axis = :x;", params = list(x = what))[[1]]
      if (length(v) == 0)
        v <- NULL
    }
  } else if (inherits(x, "ExpressionSet")) {
    v <- attr(x, what)
  }
  v
}


getDend <- function(x) {
  dend <- dbGetQuery(x, "SELECT * FROM dendrogram;")
  if (nrow(dend) == 0)
    return(NULL)
  l <- lapply(dend$object, function(x) {
    o <- str2hclust(x)
    list(ord = o$order, hcl = as.dendrogram(o))
  })
  names(l) <- dend$name
  l
}
