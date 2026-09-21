#' Internal state-format constants and helpers
#'
#' These helpers implement the versioned widget-state object used by omicsViewer
#' snapshots. The object is deliberately JSON-like: it contains widget settings
#' and stable identifiers, but no htmlwidgets, analysis results, browser-local
#' storage, session IDs or absolute temporary paths.
#'
#' @return Internal helpers return state objects, fingerprints, file names or
#'   (invisibly) side-effect results.
#'
#' @keywords internal
#' @name app_state_helpers
NULL

APP_STATE_FORMAT <- "omicsViewerState"
APP_STATE_SCHEMA_VERSION <- 1L
APP_STATE_POLICY <- "widget-only"
APP_STATE_GAPS <- list(
  stringdb = "STRING network widget/result state is not captured in this phase."
)

#' Create a versioned application state object
#'
#' @param dataset_id Character dataset identifier (usually the selected file name).
#' @param dataset The loaded omics object.
#' @param selection List with optional `features` and `samples` character vectors.
#' @param app List with app-level active-tab metadata.
#' @param data_space Data-space panel state.
#' @param result_space Result-space panel state.
#' @param label Optional user-facing snapshot label.
#' @param package_version Package version to store in the snapshot.
#' @param schema_version Integer state schema version.
#' @param created_at Snapshot creation timestamp.
#'
#' @keywords internal
#' @rdname app_state_helpers
new_app_state <- function(
    dataset_id = NA_character_,
    dataset = NULL,
    selection = list(features = character(), samples = character()),
    app = list(),
    data_space = list(),
    result_space = list(),
    label = NA_character_,
    package_version = as.character(utils::packageVersion("omicsViewer")),
    schema_version = APP_STATE_SCHEMA_VERSION,
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
) {
  ds <- dataset_state(dataset, id = dataset_id)
  selection <- normalize_selection(selection)

  structure(
    list(
      format = APP_STATE_FORMAT,
      schema_version = as.integer(schema_version),
      created_at = created_at,
      package_version = as.character(package_version),
      label = as.character(label),
      dataset = ds,
      selection = selection,
      app = as.list(app),
      panels = list(
        data_space = as.list(data_space),
        result_space = as.list(result_space)
      ),
      policy = APP_STATE_POLICY,
      gaps = APP_STATE_GAPS
    ),
    class = c("omicsViewerState", "list")
  )
}

#' Describe a dataset for snapshot compatibility checking
#'
#' @param dataset An `ExpressionSet`, `SummarizedExperiment`, or an object for
#'   which the package's expression/annotation getters work.
#' @param id Character dataset identifier.
#'
#' @keywords internal
#' @rdname app_state_helpers
dataset_state <- function(dataset, id = NA_character_) {
  if (is.null(dataset))
    return(list(id = as.character(id), class = NA_character_,
                dimensions = c(features = 0L, samples = 0L),
                fingerprint = NA_character_))

  feat <- tryCatch(
    rownames(getExprs(dataset)),
    error = function(e) rownames(dataset)
  )
  samp <- tryCatch(
    colnames(getExprs(dataset)),
    error = function(e) colnames(dataset)
  )
  fd <- tryCatch(colnames(getFData(dataset)), error = function(e) character())
  pd <- tryCatch(colnames(getPData(dataset)), error = function(e) character())
  gs <- tryCatch(
    {
      x <- unique(attr(getFData(dataset), "GS")$gsId)
      as.character(x[!is.na(x)])
    },
    error = function(e) character()
  )

  list(
    id = as.character(id),
    class = class(dataset)[1],
    dimensions = c(features = length(feat), samples = length(samp)),
    fingerprint = dataset_fingerprint(dataset, id = id)
  )
}

#' Compute a deterministic lightweight dataset fingerprint
#'
#' The fingerprint combines stable object structure and identifiers. It is not a
#' hash of expression values; expression revisions that do not alter IDs or
#' annotation columns are intentionally not treated as incompatible.
#'
#' @keywords internal
#' @rdname app_state_helpers
dataset_fingerprint <- function(dataset, id = NULL) {
  if (is.null(dataset))
    return(NA_character_)

  feat <- tryCatch(
    rownames(getExprs(dataset)),
    error = function(e) rownames(dataset)
  )
  samp <- tryCatch(
    colnames(getExprs(dataset)),
    error = function(e) colnames(dataset)
  )
  fd <- tryCatch(colnames(getFData(dataset)), error = function(e) character())
  pd <- tryCatch(colnames(getPData(dataset)), error = function(e) character())
  gs <- tryCatch(
    {
      x <- unique(attr(getFData(dataset), "GS")$gsId)
      as.character(x[!is.na(x)])
    },
    error = function(e) character()
  )

  parts <- c(
    as.character(id %.or_default% NA),
    class(dataset)[1],
    length(feat), length(samp),
    feat, samp, fd, pd, sort(gs)
  )
  .state_hash(parts)
}

`%.or_default%` <- function(x, y) if (is.null(x)) y else x

.state_hash <- function(x) {
  x <- paste0(x, collapse = "\u001f")
  # A deterministic base-R polynomial hash. It is not cryptographic, but it is
  # exact in double precision and sufficient for compatibility checks.
  rawv <- as.integer(charToRaw(enc2utf8(x)))
  hash <- 0
  for (b in rawv)
    hash <- (hash * 31 + b) %% 2147483647
  sprintf("%08x", hash)
}

#' Normalize semantic feature/sample selection
#'
#' @param selection List potentially containing `features` and `samples`.
#' @param allow_na Logical; keep NA for legacy states when true.
#'
#' @keywords internal
#' @rdname app_state_helpers
normalize_selection <- function(selection, allow_na = FALSE) {
  selection <- as.list(selection)
  clean <- function(x) {
    if (is.null(x))
      return(character())
    x <- as.character(x)
    keep <- !is.na(x) & nzchar(trimws(x))
    unique(x[keep])
  }
  list(
    features = clean(selection$features),
    samples = clean(selection$samples)
  )
}

#' Extract portable DataTable widget state
#'
#' Keep only the page position, page size, ordering, and column filters supplied
#' by DataTable's `input$<id>_state`. Callers add row selections separately using
#' stable row identifiers where available.
#'
#' @param state DataTable state input, usually `NULL` before the browser renders.
#'
#' @keywords internal
#' @rdname app_state_helpers
data_table_widget_state <- function(state) {
  if (!is.list(state) || length(state) == 0)
    return(NULL)

  columns <- NULL
  if (is.list(state$columns))
    columns <- lapply(state$columns, function(column) list(search = column$search))

  list(
    start = state$start,
    length = state$length,
    order = state$order,
    columns = columns
  )
}

#' Remove computed payloads from a module state
#'
#' Snapshot state is restricted to widget settings and semantic row identifiers.
#' Legacy snapshots can contain a few computed values; those fields are dropped so
#' downstream plots and analyses are recalculated after restoration.
#'
#' @param x Module state list.
#'
#' @keywords internal
#' @rdname app_state_helpers
.sanitize_widget_state <- function(x) {
  if (!is.list(x) || is.data.frame(x))
    return(x)

  dropped <- c(
    "external_cache", "htestV1", "htestV2", "rif",
    "rowDendrogram", "colDendrogram", "rowOrder", "colOrder"
  )
  nms <- names(x)
  if (is.null(nms))
    return(lapply(x, .sanitize_widget_state))
  x <- x[setdiff(seq_along(x), which(nms %in% dropped))]
  nms <- names(x)
  x[nms == "analyst_stringdb"] <- NULL
  x <- lapply(x, .sanitize_widget_state)
  names(x) <- names(x)
  x
}

#' Migrate a legacy flat snapshot into the versioned state object
#'
#' Legacy snapshots are plain lists whose names start with `eset_` or
#' `analyst_`, plus `active_feature`/`active_sample`. Migration only namespaces
#' those fields; child-module field names are left untouched.
#'
#' @param state A snapshot loaded from an `.ESS` file.
#'
#' @keywords internal
#' @rdname app_state_helpers
migrate_app_state <- function(state) {
  if (is.null(state))
    return(new_app_state())
  if (!is.list(state))
    stop("The selected file is not an omicsViewer snapshot.")

  if (identical(state$format, APP_STATE_FORMAT)) {
    if (!is.numeric(state$schema_version))
      stop("Snapshot has no schema version.")
    if (state$schema_version > APP_STATE_SCHEMA_VERSION)
      warning(sprintf(
        "Snapshot schema %s is newer than supported schema %s; some state may not be restored.",
        state$schema_version, APP_STATE_SCHEMA_VERSION
      ))
    if (state$schema_version < APP_STATE_SCHEMA_VERSION)
      state <- .migrate_v0_to_v1(state)
    state$panels <- lapply(state$panels, .sanitize_widget_state)
    state$policy <- APP_STATE_POLICY
    state$gaps <- APP_STATE_GAPS
    return(state)
  }

  # Legacy flat snapshot.
  nms <- names(state)
  data_space <- state[grepl("^eset_", nms)]
  result_space <- state[grepl("^analyst_", nms)]
  selection <- list(
    features = state$active_feature,
    samples = state$active_sample
  )
  new_app_state(
    dataset_id = NA_character_,
    dataset = NULL,
    selection = selection,
    app = list(
      data_active_tab = data_space$eset_active_tab,
      analysis_active_tab = result_space$analyst_active_tab
    ),
    data_space = .sanitize_widget_state(data_space),
    result_space = .sanitize_widget_state(result_space),
    label = state$label %.or_default% NA_character_,
    package_version = state$package_version %.or_default% NA_character_,
    schema_version = 1L,
    created_at = state$created_at %.or_default% NA_character_
  )
}

.migrate_v0_to_v1 <- function(state) {
  # Reserved for future migrations. Schema 1 is the first versioned format.
  state
}

#' Validate a state object against the active dataset
#'
#' Validation warns on incompatible metadata rather than blocking restoration.
#' Individual modules remain responsible for validating their own choices.
#'
#' @param state A migrated state object.
#' @param dataset Active dataset object.
#' @param dataset_id Active dataset identifier.
#'
#' @keywords internal
#' @rdname app_state_helpers
validate_app_state <- function(state, dataset = NULL, dataset_id = NA_character_) {
  state <- migrate_app_state(state)
  if (!identical(state$format, APP_STATE_FORMAT))
    warning("Unknown snapshot format; state restoration may be incomplete.")

  current <- dataset_state(dataset, id = dataset_id)
  old <- state$dataset
  if (!is.null(old) && !is.na(old$fingerprint) && !is.na(current$fingerprint) &&
      !identical(old$fingerprint, current$fingerprint)) {
    warning(sprintf(
      "Snapshot dataset fingerprint mismatch (snapshot: %s; current: %s). Missing features, samples and annotations will be ignored.",
      old$fingerprint, current$fingerprint
    ))
  }
  if (!is.null(old$id) && !is.na(old$id) && !is.na(dataset_id) &&
      !identical(as.character(old$id), as.character(dataset_id))) {
    warning(sprintf(
      "Snapshot was saved for dataset '%s' and is being restored into '%s'.",
      old$id, dataset_id
    ))
  }

  feat <- tryCatch(
    rownames(getExprs(dataset)),
    error = function(e) rownames(dataset)
  )
  samp <- tryCatch(
    colnames(getExprs(dataset)),
    error = function(e) colnames(dataset)
  )
  sel <- normalize_selection(state$selection)
  missing_f <- setdiff(sel$features, feat)
  missing_s <- setdiff(sel$samples, samp)
  if (length(missing_f))
    warning("Snapshot contains unknown features: ", paste(missing_f, collapse = ", "))
  if (length(missing_s))
    warning("Snapshot contains unknown samples: ", paste(missing_s, collapse = ", "))

  state
}

#' Sanitize a user-supplied snapshot name
#'
#' The returned name contains only portable filename characters and cannot
#' represent a path or parent-directory reference.
#'
#' @param name Character user input.
#' @param fallback Name used when input is empty.
#'
#' @keywords internal
#' @rdname app_state_helpers
sanitize_snapshot_name <- function(name, fallback = "snapshot") {
  name <- as.character(name)[1]
  if (is.na(name))
    name <- ""
  name <- gsub("[^A-Za-z0-9_.-]+", "_", name, perl = TRUE)
  name <- gsub("^[_.-]+", "", name)
  name <- gsub("[_.-]+$", "", name)
  if (!nzchar(name) || name == "." || name == "..")
    name <- fallback
  name
}

#' Generate a safe snapshot file name
#'
#' @param dataset_id Dataset identifier used in the legacy-compatible prefix.
#'
#' @keywords internal
#' @rdname app_state_helpers
snapshot_file_name <- function(name, dataset_id = "ESVObj.RDS", fallback = "snapshot") {
  name <- sanitize_snapshot_name(name, fallback = fallback)
  dataset_id <- sanitize_snapshot_name(dataset_id, fallback = "ESVObj.RDS")
  paste0("ESVSnapshot_", dataset_id, "_", name, ".ESS")
}

#' Build the top-level snapshot object from module states
#'
#' This pure helper is separated from the Shiny server so schema construction
#' can be tested without launching the full application.
#'
#' @param dataset Active dataset.
#' @param dataset_id Active dataset identifier.
#' @param data_status State returned by the data-space module.
#' @param result_status State returned by the result-space module.
#' @param selected_features Active semantic feature IDs.
#' @param selected_samples Active semantic sample IDs.
#' @param label Snapshot label/name.
#' @keywords internal
#' @rdname app_state_helpers
build_app_state <- function(
    dataset,
    dataset_id = NA_character_,
    data_status = list(),
    result_status = list(),
    selected_features = character(),
    selected_samples = character(),
    label = NA_character_,
    package_version = as.character(utils::packageVersion("omicsViewer"))
) {
  data_status <- as.list(data_status)
  result_status <- as.list(result_status)
  new_app_state(
    dataset_id = dataset_id,
    dataset = dataset,
    selection = list(features = selected_features, samples = selected_samples),
    app = list(
      data_active_tab = data_status$eset_active_tab,
      analysis_active_tab = result_status$analyst_active_tab
    ),
    data_space = .sanitize_widget_state(data_status),
    result_space = .sanitize_widget_state(result_status),
    label = label,
    package_version = package_version
  )
}

#' Atomically write a state object
#'
#' The object is first written to a temporary file in the destination directory
#' and then renamed. A failed write therefore never replaces a valid snapshot.
#'
#' @param state State object.
#' @param path Destination file path.
#'
#' @keywords internal
#' @rdname app_state_helpers
write_app_state <- function(state, path) {
  dir <- dirname(path)
  if (!dir.exists(dir))
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = ".ESS_", tmpdir = dir, fileext = ".tmp")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(state, tmp, compress = "xz")
  if (file.exists(path))
    stop("Snapshot file already exists: ", basename(path))
  if (!file.rename(tmp, path))
    stop("Could not finalize snapshot file: ", basename(path))
  invisible(path)
}

#' Restore a saved DataTable page length
#'
#' @param length Saved DataTable page length.
#' @param pageLength Module default when no valid saved value exists.
#'
#' @keywords internal
#' @rdname app_state_helpers
restore_table_page_length <- function(length, pageLength) {
  if (is.numeric(length) && length(length) == 1L && !is.na(length) && length > 0)
    as.integer(length)
  else
    pageLength
}
