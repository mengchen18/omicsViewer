#' Save the omicsViewer result object as sqlite database
#'
#' @description The object is normalized to the app-ready view before
#'   writing (see \code{\link{tallGS}}): wide \code{GS|*} gene-set columns
#'   are moved into the long \code{GS} table, so the \code{feature} table
#'   matches what the \code{.RDS} load path serves after
#'   \code{\link{readESVObj}}. Selected metadata attributes (user-defined
#'   quick views, \code{volcanoCategories}, \code{cormat}, \code{S6.6_drc},
#'   precomputed expression dendrograms, column type/level metadata) are
#'   preserved in the \code{attrs} table and restored by the getters
#'   (\code{getExprs}, \code{getPData}, \code{getFData}, \code{getAx}).
#'
#'   Tables written: \code{exprs}, \code{exprsimpute} (only when the object
#'   carries imputed values; a stale table from an earlier save is
#'   dropped), \code{feature}, \code{sample}, \code{GS}, \code{axes},
#'   \code{dendrogram}, \code{meta} (schema version marker; absent in
#'   pre-2.1.81 files) and \code{attrs}. All writes run in a single
#'   transaction, so a failed save never leaves a half-written database.
#'
#'   The \code{attrs.value} column stores base64-encoded serialized R
#'   objects - as with \code{\link{readRDS}}, only open databases you
#'   trust.
#'
#' @param obj an object of class ExpressionSet or SummarizedExperiment
#' @param db.file a character indicate file name of the database file
#' @param overwrite logical. whether the database should be overwritten if exist already.
#' @return the directory where the database saved
#' @name saveOmicsViewerDb
#' @importFrom RSQLite SQLite dbConnect dbDisconnect dbGetQuery dbWriteTable dbListTables dbBegin dbCommit dbRollback dbExecute SQLITE_RO
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

    # sqlite-revival WP-1: persist the app-ready view. tallGS() moves the
    # wide `GS|*` membership columns into the long gene-set annotation, so
    # the `feature` table matches what the .RDS load path serves after
    # readESVObj() (and the `GS` table is populated from either gene-set
    # representation). Idempotent for already-normalized objects.
    obj <- tallGS(obj)

    expr <- exprs(obj)
    pd <- pData(obj)
    fd <- fData(obj)
    expr_impute <- exprsImpute(obj)

    # sqlite-revival WP-5: validate column names BEFORE touching the file -
    # duplicate names cannot be represented in SQLite and 'rowname' collides
    # with the internal key column; unsupported classes fail deep inside
    # dbWriteTable otherwise.
    .db_guard_table(expr, "exprs")
    .db_guard_table(pd, "pData")
    .db_guard_table(fd, "fData")

    # sqlite-revival WP-5: SQLite cannot represent NaN (it degrades to
    # NULL); normalize NaN -> NA at save so the round trip is deterministic
    expr <- .db_nan_to_na(expr)
    pd <- .db_nan_to_na(pd)
    fd <- .db_nan_to_na(fd)
    if (!is.null(expr_impute))
      expr_impute <- .db_nan_to_na(expr_impute)

    con <- dbConnect(RSQLite::SQLite(), db.file)
    on.exit(dbDisconnect(con))

    if (!overwrite) {
      clash <- intersect(
        dbListTables(con),
        c("exprs", "exprsimpute", "feature", "sample", "GS", "axes",
          "dendrogram", "meta", "attrs"))
      if (length(clash) > 0)
        stop(
          "saveOmicsViewerDb: table(s) already exist (",
          paste(clash, collapse = ", "),
          ") - use overwrite = TRUE")
    }

    dfrname <- function(x) {
      data.frame(x, rowname = rownames(x), check.names = FALSE)
    }

    # sqlite-revival WP-1: one transaction around every write, so a failed
    # save never leaves a half-written database behind
    dbBegin(con)
    tryCatch({
      dbWriteTable(conn = con, name = "exprs", value = dfrname(expr), overwrite = TRUE)
      dbWriteTable(conn = con, name = "feature", value = dfrname(fd), overwrite = TRUE)
      dbWriteTable(conn = con, name = "sample", value = dfrname(pd), overwrite = TRUE)
      # the imputed matrix is written conditionally: a stale table from an
      # earlier save must not survive a re-save without imputed values
      if (!is.null(expr_impute)) {
        dbWriteTable(conn = con, name = "exprsimpute", value = dfrname(expr_impute), overwrite = TRUE)
      } else {
        dbExecute(con, "DROP TABLE IF EXISTS exprsimpute")
      }

      gs <- data.frame(featureId = character(0), gsId = character(0), weight = numeric(0))
      if (hasAttr(fd, "GS"))
        gs <- attr(fd, "GS")
      if (inherits(gs, c("dgCMatrix", "lgCMatrix")))
        gs <- csc2list(gs)
      dbWriteTable(conn = con, name = "GS", value = gs, overwrite = TRUE)

      defaultxy <- data.frame(axis = character(0), value = character(0))
      ax <- c("fx", "fy", "sx", "sy")
      ax <- ax[hasAttr(obj, attr.name = ax)]
      if (length(ax) > 0)
        defaultxy <- data.frame(
          axis = ax,
          value = vapply(ax, function(v) attr(obj, v), character(1))
        )
      dbWriteTable(conn = con, name = "axes", value = defaultxy, overwrite = TRUE)

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
      dbWriteTable(conn = con, name = "dendrogram", value = dend, overwrite = TRUE)

      # sqlite-revival WP-1: schema/version marker; absence marks a legacy
      # (pre-2.1.81) file for the reader
      dbWriteTable(conn = con, name = "meta", value = data.frame(
        key = c("schema_version", "pkg_version", "created",
                "n_feature", "n_sample"),
        value = as.character(c(
          1L,
          tryCatch(as.character(packageVersion("omicsViewer")),
                   error = function(e) "dev"),
          format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
          nrow(fd), ncol(expr))),
        stringsAsFactors = FALSE), overwrite = TRUE)

      # sqlite-revival WP-3: attribute payloads - quickViews/shortcut/
      # volcanoCategories (fData/pData), cormat/S6.6_drc (object),
      # rowDendrogram/colDendrogram (exprs), plus the column-type
      # metadata that lets the reader rebuild factors/logicals/Dates
      dbWriteTable(conn = con, name = "attrs", value = .db_attrs_save_frame(
        obj = obj, fd = fd, pd = pd, expr = expr), overwrite = TRUE)

      dbCommit(con)
    }, error = function(e) {
      tryCatch(dbRollback(con), error = function(e2) NULL)
      stop(e)
    })
    dirname(db.file)
  }
)

## -------------------------------------------------------------------
## sqlite-revival internal helpers (WP-1 .. WP-5)
## -------------------------------------------------------------------

# attribute whitelist per serialization scope; keys outside these lists are
# deliberately not persisted (they are either runtime-only or carried by
# their own table: fx/fy/sx/sy -> axes, dendrogram -> dendrogram, GS -> GS)
.db_attr_keys <- list(
  object = c("cormat", "S6.6_drc"),
  fdata  = c("quickViews", "shortcut", "volcanoCategories"),
  pdata  = c("quickViews", "shortcut"),
  exprs  = c("rowDendrogram", "colDendrogram")
)

# column classes dbWriteTable can represent without data loss; Date is
# stored as ISO character and rebuilt from the coltypes metadata
.db_allowed_col_classes <- c(
  "numeric", "double", "integer", "character", "factor", "ordered",
  "logical", "POSIXct", "Date")

.db_guard_table <- function(x, what) {
  cn <- colnames(x)
  if ("rowname" %in% cn)
    stop("saveOmicsViewerDb: ", what,
         " contains a column named 'rowname', which collides with the ",
         "internal database key column. Rename it before saving.")
  d <- anyDuplicated(cn)
  if (d > 0)
    stop("saveOmicsViewerDb: ", what, " contains duplicated column name(s): ",
         paste(unique(cn[duplicated(cn)]), collapse = ", "),
         ". SQLite requires unique column names; make the identifiers ",
         "unique before saving.")
  badcl <- cn[!vapply(x, function(col) class(col)[1] %in% .db_allowed_col_classes,
                      logical(1))]
  if (length(badcl) > 0)
    stop("saveOmicsViewerDb: ", what, " contains column(s) of unsupported ",
         "class (", paste(unique(vapply(x[, badcl, drop = FALSE],
                                        function(col) class(col)[1],
                                        character(1))), collapse = ", "),
         "): ", paste(badcl, collapse = ", "),
         ". Coerce to character/numeric/factor before saving.")
  invisible(TRUE)
}

.db_nan_to_na <- function(x) {
  if (is.matrix(x)) {
    if (any(idx <- is.nan(x)))
      x[idx] <- NA
    return(x)
  }
  for (i in which(vapply(x, is.numeric, logical(1)))) {
    v <- x[[i]]
    if (any(idx <- is.nan(v))) {
      v[idx] <- NA
      x[[i]] <- v
    }
  }
  x
}

.db_attr_serialize <- function(value) {
  jsonlite::base64_enc(serialize(value, connection = NULL, version = 2))
}

.db_attr_unserialize <- function(value) {
  unserialize(jsonlite::base64_dec(value))
}

.db_attrs_save_frame <- function(obj, fd, pd, expr) {
  rows <- list()
  add <- function(scope, key, value) {
    if (is.null(value))
      return(invisible(NULL))
    rows[[length(rows) + 1L]] <<- data.frame(
      scope = scope, key = key, format = "rds64",
      value = .db_attr_serialize(value), stringsAsFactors = FALSE)
  }
  for (k in .db_attr_keys$object) add("object", k, attr(obj, k))
  for (k in .db_attr_keys$fdata) add("fdata", k, attr(fd, k))
  add("fdata", "coltypes", .db_record_coltypes(fd))
  # gene-set factor levels: the .RDS path derives them from the wide
  # matrix dimnames (full rownames / all column names); record them so
  # the reader rebuilds bit-identical factors from the GS table
  gsdf <- attr(fd, "GS")
  if (!is.null(gsdf) && is.data.frame(gsdf) && nrow(gsdf) &&
      is.factor(gsdf$featureId) && is.factor(gsdf$gsId))
    add("fdata", "gsLevels", list(
      featureId = as.character(levels(gsdf$featureId)),
      gsId = as.character(levels(gsdf$gsId))))
  for (k in .db_attr_keys$pdata) add("pdata", k, attr(pd, k))
  add("pdata", "coltypes", .db_record_coltypes(pd))
  for (k in .db_attr_keys$exprs) add("exprs", k, attr(expr, k))
  if (!length(rows))
    return(data.frame(scope = character(0), key = character(0),
                      format = character(0), value = character(0),
                      stringsAsFactors = FALSE))
  do.call(rbind, rows)
}

# read all attribute payloads of one scope; tolerant of missing table
# (legacy files) and of unparsable payloads (skipped with a warning)
.db_read_attrs <- function(con, scope) {
  if (!"attrs" %in% dbListTables(con))
    return(list())
  tab <- tryCatch(
    dbGetQuery(con,
               "SELECT key, format, value FROM attrs WHERE scope = :s;",
               params = list(s = scope)),
    error = function(e) NULL)
  if (is.null(tab) || !nrow(tab))
    return(list())
  out <- list()
  for (i in seq_len(nrow(tab))) {
    v <- tryCatch(
      if (identical(tab$format[i], "rds64")) .db_attr_unserialize(tab$value[i]) else NULL,
      error = function(e) NULL)
    if (is.null(v)) {
      warning("getFData/getPData: skipping unparsable attribute '",
              tab$key[i], "' (scope '", scope, "')")
      next
    }
    out[[tab$key[i]]] <- v
  }
  out
}

.db_read_attr_one <- function(con, scope, key) {
  a <- .db_read_attrs(con, scope)
  if (is.null(a) || !length(a)) return(NULL)
  a[[key]]
}

# column-type metadata: classes SQLite cannot represent natively (factor
# levels, logical, Date) are recorded here and rebuilt on read
.db_record_coltypes <- function(df) {
  out <- lapply(df, function(col) {
    cl <- class(col)[1]
    if (identical(cl, "factor"))
      list(class = "factor", levels = as.character(levels(col)))
    else if (cl %in% c("logical", "integer", "Date", "ordered"))
      list(class = cl)
    else
      NULL
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) NULL else out
}

.db_apply_coltypes <- function(df, ct) {
  if (is.null(ct) || !length(ct))
    return(df)
  for (nm in intersect(names(ct), colnames(df))) {
    e <- ct[[nm]]
    if (identical(e$class, "factor")) {
      df[[nm]] <- factor(df[[nm]], levels = e$levels)
    } else if (identical(e$class, "logical")) {
      df[[nm]] <- as.logical(df[[nm]])
    } else if (identical(e$class, "integer")) {
      df[[nm]] <- as.integer(df[[nm]])
    } else if (identical(e$class, "Date")) {
      df[[nm]] <- as.Date(df[[nm]])
    }
  }
  df
}

# legacy (pre-2.1.81) files kept the wide `GS|*` membership columns in the
# feature table with an empty GS table - recover the tall annotation the
# same way tallGS() does on the .RDS path
.db_recover_wide_gs <- function(fd) {
  ir <- grep("^GS\\|", colnames(fd))
  if (!length(ir))
    return(NULL)
  igs <- fd[, ir, drop = FALSE]
  nm <- vapply(strsplit(colnames(igs), "|", fixed = TRUE),
               function(p) if (length(p) >= 3) p[3] else tail(p, 1),
               character(1))
  colnames(igs) <- make.names(nm)
  totall(igs)
}

# read a numeric matrix table (exprs/exprsimpute) without the historical
# apply(, 2, as.numeric) character round trip
.db_read_matrix <- function(con, table) {
  mat <- dbGetQuery(con, sprintf("SELECT * FROM %s;", table))
  rn <- mat$rowname
  mat$rowname <- NULL
  if (ncol(mat) == 0) {
    m <- matrix(numeric(0), nrow = length(rn), ncol = 0)
  } else if (all(vapply(mat, is.numeric, logical(1)))) {
    m <- as.matrix(mat)
  } else {
    m <- vapply(mat, as.numeric, numeric(nrow(mat)))
  }
  rownames(m) <- rn
  m
}

# close a dataset connection quietly; safe on non-connections and on
# already-closed connections (double disconnect after a dataset switch)
esv_close_connection <- function(x) {
  if (inherits(x, "SQLiteConnection"))
    tryCatch(dbDisconnect(x), error = function(e) NULL)
  invisible(NULL)
}

.db_legacy_notified <- new.env(parent = emptyenv())

.db_note_legacy <- function(path) {
  if (!is.null(.db_legacy_notified[[path]]))
    return(invisible(NULL))
  .db_legacy_notified[[path]] <- TRUE
  message("readESVObj: ", path, " carries no 'meta' table - legacy ",
          "(pre-2.1.81) omicsViewer database; reading with compatibility ",
          "recovery")
  invisible(NULL)
}

.db_connect_esv <- function(path) {
  if (!file.exists(path))
    stop("readESVObj: database file does not exist: ", path)
  # read-only: the app never writes to dataset connections, and a shared
  # read-only connection cannot corrupt the file or hold a write lock
  con <- dbConnect(RSQLite::SQLite(), path, flags = RSQLite::SQLITE_RO)
  ok <- tryCatch({
    tabs <- dbListTables(con)
    c("exprs", "feature", "sample") %in% tabs
  }, error = function(e) {
    tryCatch(dbDisconnect(con), error = function(e2) NULL)
    stop("readESVObj: not a valid SQLite database: ", path, " (",
         conditionMessage(e), ")")
  })
  if (!all(ok)) {
    missing <- c("exprs", "feature", "sample")[!ok]
    tryCatch(dbDisconnect(con), error = function(e) NULL)
    stop("readESVObj: database is missing required table(s): ",
         paste(missing, collapse = ", "), " - '", path,
         "' does not look like an omicsViewer database")
  }
  if (!"meta" %in% dbListTables(con))
    .db_note_legacy(path)
  con
}

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
#' @param x an object of class \code{SummarizedExperiment} or
#'   \code{ExpressionSet}, or any class with its own \code{asEsetWithAttr}
#'   method (dispatched automatically)
#' @return an object of class ExpressionSet
#' @name asEsetWithAttr
#' @export
setGeneric("asEsetWithAttr", function(x) standardGeneric("asEsetWithAttr"))

#' @rdname asEsetWithAttr
setMethod("asEsetWithAttr", "SummarizedExperiment", function(x) {
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
  eset
})

#' @rdname asEsetWithAttr
setMethod("asEsetWithAttr", "ExpressionSet", function(x) x)

#' @rdname asEsetWithAttr
setMethod("asEsetWithAttr", "ANY", function(x)
  stop("x should be either an SummarizedExperiment or ExpressionSet"))

#' Read the object of SummarizedExperiment or ExpressetSet to be visualized using omicsViewer
#' @description This function accept a path to a sqlite database or RDS object. If an RDS file to be read, 
#'   The function is similar to \code{readRDS}. It reads the object to R working environment and perform extra two things. 
#'   1. If the loaded data an class of \code{SummarizedExperiment}, it will be converted to \code{ExpressionSet};
#'   2. If the gene set annotatio is in matrix format, the gene set annotation is converted to \code{data.frame} format.
#'
#'   For SQLite files the function returns a live read-only
#'   \code{SQLiteConnection}; the app reads it through the getter methods
#'   (\code{getExprs}, \code{getPData}, \code{getFData}, \code{getAx}) and
#'   closes the connection when the dataset is switched or the session
#'   ends. Standalone users should close it themselves with
#'   \code{DBI::dbDisconnect()}. Databases missing the required tables are
#'   rejected with an explicit error; files without a \code{meta} table
#'   are legacy (pre-2.1.81) databases and are read with compatibility
#'   recovery (wide \code{GS|*} columns are converted back into the gene-set
#'   annotation).
#' @return an object of class \code{ExpressionSet} or an SQLite connection to be visualzied.
#' @param x the path of an object of \code{SummarizedExperiment} or \code{ExpressionSet}, passed to \link{readRDS}
#' @importFrom Biobase pData fData fData<- pData<-
#' @importFrom methods as
#' @export
#' @examples
#' file <- system.file("extdata/demo.RDS", package = "omicsViewer")
#' obj <- readESVObj(file)
#'
readESVObj <- function(x) {
  if (grepl("\\.rds$", x, ignore.case = TRUE)) {
    x <- asEsetWithAttr( readRDS(x) )
    x <- tallGS(x)
  } else if (grepl("\\.(db|sqlite|sqlite3)$", x, ignore.case = TRUE)) {
    x <- .db_connect_esv(x)
  } else
    stop("readESVObj: Unkown input format!")
  x
}

#' Extract the expression matrix of an omicsViewer dataset
#' @param x an SQLite connection, an \code{ExpressionSet}, or any class
#'   with its own \code{getExprs} method (dispatched automatically)
#' @return the (features x samples) expression matrix; for SQLite
#'   connections, matrix-level attributes (precomputed
#'   \code{rowDendrogram}/\code{colDendrogram}) are restored from the
#'   \code{attrs} table
#' @name getExprs
#' @export
setGeneric("getExprs", function(x) standardGeneric("getExprs"))

#' @rdname getExprs
setMethod("getExprs", "SQLiteConnection", function(x) {
  m <- .db_read_matrix(x, "exprs")
  a <- .db_read_attrs(x, "exprs")
  for (k in names(a))
    attr(m, k) <- a[[k]]
  m
})

#' @rdname getExprs
setMethod("getExprs", "ExpressionSet", function(x) exprs(x))

#' @rdname getExprs
setMethod("getExprs", "ANY", function(x)
  stop("getExprs: unsupported class '", class(x)[1], "' - expected a SQLite connection or an ExpressionSet (convert SummarizedExperiment via asEsetWithAttr)"))

#' Extract the imputed expression matrix of an omicsViewer dataset
#' @param x an SQLite connection, an \code{ExpressionSet}, or any class
#'   with its own \code{getExprsImpute} method (dispatched automatically)
#' @return the imputed expression matrix, or NULL if the dataset carries
#'   no imputed values
#' @name getExprsImpute
#' @export
setGeneric("getExprsImpute", function(x) standardGeneric("getExprsImpute"))

#' @rdname getExprsImpute
setMethod("getExprsImpute", "SQLiteConnection", function(x) {
  if (!"exprsimpute" %in% dbListTables(x))
    return(NULL)
  .db_read_matrix(x, "exprsimpute")
})

#' @rdname getExprsImpute
setMethod("getExprsImpute", "ExpressionSet", function(x) exprsImpute(x))

#' @rdname getExprsImpute
setMethod("getExprsImpute", "ANY", function(x)
  stop("getExprsImpute: unsupported class '", class(x)[1], "' - expected a SQLite connection or an ExpressionSet (convert SummarizedExperiment via asEsetWithAttr)"))

#' Extract the phenotype (sample) annotation of an omicsViewer dataset
#' @param x an SQLite connection, an \code{ExpressionSet}, or any class
#'   with its own \code{getPData} method (dispatched automatically)
#' @return the sample annotation data.frame; column types (factor levels,
#'   logicals) are rebuilt from the \code{attrs} metadata and
#'   user-defined quick-view shortcuts ride along as the
#'   \code{"quickViews"} attribute where available
#' @name getPData
#' @export
setGeneric("getPData", function(x) standardGeneric("getPData"))

#' @rdname getPData
setMethod("getPData", "SQLiteConnection", function(x) {
  mat <- dbGetQuery(x, "SELECT * FROM sample;")
  rownames(mat) <- mat$rowname
  mat$rowname <- NULL
  a <- .db_read_attrs(x, "pdata")
  mat <- .db_apply_coltypes(mat, a$coltypes)
  a$coltypes <- NULL
  for (k in names(a))
    attr(mat, k) <- a[[k]]
  mat
})

#' @rdname getPData
setMethod("getPData", "ExpressionSet", function(x) pData(x))

#' @rdname getPData
setMethod("getPData", "ANY", function(x)
  stop("getPData: unsupported class '", class(x)[1], "' - expected a SQLite connection or an ExpressionSet (convert SummarizedExperiment via asEsetWithAttr)"))

#' Extract the feature annotation of an omicsViewer dataset
#' @param x an SQLite connection, an \code{ExpressionSet}, or any class
#'   with its own \code{getFData} method (dispatched automatically)
#' @return the feature annotation data.frame; the gene-set annotation
#'   rides along as the \code{GS} attribute where available. Legacy
#'   (pre-2.1.81) databases are normalized on read: wide \code{GS|*}
#'   membership columns are converted back into the gene-set annotation,
#'   so both file generations serve the same app-ready view
#' @name getFData
#' @export
setGeneric("getFData", function(x) standardGeneric("getFData"))

#' @rdname getFData
setMethod("getFData", "SQLiteConnection", function(x) {
  mat <- dbGetQuery(x, "SELECT * FROM feature;")
  rownames(mat) <- mat$rowname
  mat$rowname <- NULL
  a <- .db_read_attrs(x, "fdata")
  gs <- dbGetQuery(x, "SELECT * FROM GS;")
  wide <- grep("^GS\\|", colnames(mat))
  # legacy recovery: wide GS|* columns + empty GS table
  if (nrow(gs) == 0 && length(wide) > 0)
    gs <- .db_recover_wide_gs(mat)
  if (length(wide) > 0)
    mat <- mat[, -wide, drop = FALSE]
  if (!is.null(gs) && nrow(gs) > 0 &&
      all(c("featureId", "gsId", "weight") %in% colnames(gs))) {
    gl <- a$gsLevels
    gs$featureId <- if (!is.null(gl))
      factor(gs$featureId, levels = gl$featureId) else as.factor(gs$featureId)
    gs$gsId <- if (!is.null(gl))
      factor(gs$gsId, levels = gl$gsId) else as.factor(gs$gsId)
    # normalize to the compact automatic row names the melt-derived GS
    # annotation carries on the .RDS path (identical() is strict here)
    rownames(gs) <- NULL
    attr(mat, "GS") <- gs
  }
  mat <- .db_apply_coltypes(mat, a$coltypes)
  a$coltypes <- NULL
  a$gsLevels <- NULL
  for (k in names(a))
    attr(mat, k) <- a[[k]]
  mat
})

#' @rdname getFData
setMethod("getFData", "ExpressionSet", function(x) fData(x))

#' @rdname getFData
setMethod("getFData", "ANY", function(x)
  stop("getFData: unsupported class '", class(x)[1], "' - expected a SQLite connection or an ExpressionSet (convert SummarizedExperiment via asEsetWithAttr)"))

#' Read a default axis of an omicsViewer dataset
#' @param x an SQLite connection, an \code{ExpressionSet}, or any class
#'   with its own \code{getAx} method (dispatched automatically)
#' @param what a character scalar naming the axis (\code{fx}, \code{fy},
#'   \code{sx}, \code{sy} or \code{dendrogram}); object-level attributes
#'   persisted through the SQLite path (e.g. \code{cormat},
#'   \code{S6.6_drc}) are served through the same seam
#' @return the axis value (usually a column name), or NULL when unset
#' @name getAx
#' @export
setGeneric("getAx", function(x, what) standardGeneric("getAx"))

#' @rdname getAx
setMethod("getAx", signature(x = "SQLiteConnection", what = "ANY"), function(x, what) {
  if (what == "dendrogram") {
    v <- getDend(x)
  } else {
    v <- dbGetQuery(x, "SELECT value FROM axes WHERE axis = :x;", params = list(x = what))[[1]]
    if (length(v) == 0) {
      # object-level attributes (cormat, S6.6_drc, ...) ride in the attrs
      # table; NULL when the dataset carries none
      v <- .db_read_attr_one(x, "object", what)
    }
  }
  v
})

#' @rdname getAx
setMethod("getAx", signature(x = "ExpressionSet", what = "ANY"), function(x, what)
  attr(x, what))

#' @rdname getAx
setMethod("getAx", signature(x = "ANY", what = "ANY"), function(x, what)
  stop("getAx: unsupported class '", class(x)[1], "' - expected a SQLite connection or an ExpressionSet (convert SummarizedExperiment via asEsetWithAttr)"))


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
