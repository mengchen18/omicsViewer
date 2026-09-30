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
# Schema 2 (todo 2.6): the canonical widget store rides as a first-class
# section (values + explicit unset list); the semantic selection carries
# the full selection-bus records (per space: ids, clicked, origin, anchor,
# mirror); the panel sections keep only what the store cannot model (DT
# column order, heatmap zoom). Schema 1 snapshots are migrated on read.
APP_STATE_SCHEMA_VERSION <- 2L
APP_STATE_POLICY <- "widget-only"
APP_STATE_GAPS <- list(
  stringdb = "STRING network widget/result state is not captured in this phase."
)

#' Create a versioned application state object
#'
#' @param dataset_id Character dataset identifier (usually the selected file name).
#' @param dataset The loaded omics object.
#' @param selection List with optional `features` and `samples` character
#'   vectors and (schema 2) `records`, the selection-bus record per space.
#' @param app List with app-level active-tab metadata.
#' @param data_space Data-space panel state.
#' @param result_space Result-space panel state.
#' @param widget_store Optional canonical widget-store snapshot
#'   (`store_snapshot()` shape, plus an explicit `unset` list).
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
    widget_store = NULL,
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
      widget_store = if (is.list(widget_store)) widget_store else NULL,
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
# fingerprint cache for SQLite connections (M6): every fingerprint used to
# cost 2 full exprs-table reads + 2 full feature-table reads; saves and
# restores re-read the same unchanged file repeatedly. Cached by database
# path + modification time, so any file change re-computes.
.fingerprint_cache <- new.env(parent = emptyenv())

dataset_fingerprint <- function(dataset, id = NULL) {
  if (is.null(dataset))
    return(NA_character_)

  cache_key <- NULL
  if (inherits(dataset, "SQLiteConnection")) {
    info <- tryCatch(DBI::dbGetInfo(dataset), error = function(e) NULL)
    path <- if (is.null(info)) "" else (info$dbname %.or_default% "")
    mt <- tryCatch(file.mtime(path), error = function(e) NA_real_)[1]
    if (length(mt) == 1 && !is.na(mt)) {
      cache_key <- paste(path, format(mt), id %.or_default% "")
      hit <- tryCatch(.fingerprint_cache[[cache_key]], error = function(e) NULL)
      if (length(hit) == 1 && !is.na(hit))
        return(hit)
    }
  }

  # M6: one read per table (was: getExprs x2, getFData x2)
  ex <- tryCatch(getExprs(dataset), error = function(e) NULL)
  fd0 <- tryCatch(getFData(dataset), error = function(e) NULL)
  feat <- if (is.null(ex)) rownames(dataset) else rownames(ex)
  samp <- if (is.null(ex)) colnames(dataset) else colnames(ex)
  fd <- if (is.null(fd0)) character() else colnames(fd0)
  pd <- tryCatch(colnames(getPData(dataset)), error = function(e) character())
  gs <- if (is.null(fd0)) character() else {
    x <- unique(attr(fd0, "GS")$gsId)
    as.character(x[!is.na(x)])
  }

  parts <- c(
    as.character(id %.or_default% NA),
    class(dataset)[1],
    length(feat), length(samp),
    feat, samp, fd, pd, sort(gs)
  )
  fp <- .state_hash(parts)
  if (!is.null(cache_key))
    .fingerprint_cache[[cache_key]] <- fp
  fp
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
  out <- list(
    features = clean(selection$features),
    samples = clean(selection$samples)
  )
  # schema 2: the selection-bus records ride along untouched
  if (is.list(selection$records))
    out$records <- selection$records
  out
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
#' those fields; child-module field names are left untouched. A plain list
#' WITHOUT any of those names is rejected: a stray `.RDS` file that happens
#' to be a list must not "restore" as an empty state and silently clear the
#' selection (todo 2.7).
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
    if (state$schema_version < 1L)
      state <- .migrate_v0_to_v1(state)
    if (state$schema_version < 2L)
      state <- .migrate_v1_to_v2(state)
    state$panels <- lapply(state$panels, .sanitize_widget_state)
    state$policy <- APP_STATE_POLICY
    state$gaps <- APP_STATE_GAPS
    return(state)
  }

  # Legacy flat snapshot: require the naming convention before treating a
  # random list as one (todo 2.7)
  nms <- names(state)
  if (is.null(nms) ||
      !any(c(grepl("^eset_", nms), grepl("^analyst_", nms),
             nms %in% c("active_feature", "active_sample"))))
    stop("The selected file is not an omicsViewer snapshot.")
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
    schema_version = 2L,
    created_at = state$created_at %.or_default% NA_character_
  )
}

.migrate_v0_to_v1 <- function(state) {
  # Reserved for future migrations. Schema 1 is the first versioned format.
  state
}

#' Migrate a schema-1 snapshot to schema 2 (todo 2.6)
#'
#' Moves the store-owned scatter axis/mode/attribute state into the
#' `widget_store` section and derives the selection-bus records from the
#' legacy fields. The panel copies of the moved fields are dropped so a
#' restore reads them from exactly one place (the store).
.migrate_v1_to_v2 <- function(state) {
  if (!identical(state$schema_version, 1L)) {
    state$schema_version <- as.integer(state$schema_version)
  }
  values <- as.list(state$widget_store$values)
  if (is.null(values)) values <- list()

  fig_map <- c(
    eset_fdata_fig = "dataspace.feature_space",
    eset_pdata_fig = "dataspace.sample_space"
  )
  for (panel_key in names(fig_map)) {
    fig <- state$panels$data_space[[panel_key]]
    if (!is.list(fig)) next
    prefix <- fig_map[[panel_key]]
    if (identical(length(fig$xax), 3L))
      values <- .state_assign(values, paste0(prefix, ".", c("x_analysis", "x_subset", "x_variable")), as.list(fig$xax))
    if (identical(length(fig$yax), 3L))
      values <- .state_assign(values, paste0(prefix, ".", c("y_analysis", "y_subset", "y_variable")), as.list(fig$yax))
    if (!is.null(fig$axisMode) && fig$axisMode %in% c("quick", "custom"))
      values[[paste0(prefix, ".axis_mode")]] <- fig$axisMode
    a4_groups <- list(color = "selectColor", shape = "selectShape",
                      size = "selectSize", tooltip = "selectTooltip",
                      search = "searchOnCol")
    a4 <- fig$attr4
    if (is.list(a4)) {
      for (g in names(a4_groups)) {
        tr <- a4[[a4_groups[[g]]]]
        if (is.list(tr) && length(tr) == 3L &&
            !identical(tr$variable, "--select--") && nzchar(tr$variable %||% "")) {
          values <- .state_assign(values,
            paste0(prefix, ".attr4.", g, c("_analysis", "_subset", "_variable")),
            list(tr$analysis, tr$subset, tr$variable))
        }
      }
      if (!is.null(a4$xcut)) values[[paste0(prefix, ".attr4.xcut")]] <- a4$xcut
      if (!is.null(a4$ycut)) values[[paste0(prefix, ".attr4.ycut")]] <- a4$ycut
      if (!is.null(a4$acorner)) values[[paste0(prefix, ".attr4.scorner")]] <- a4$acorner
    }
    # drop the migrated axis/mode copies (the store is the single write
    # plane); attr4 stays - its searchValue field has no store key and the
    # panel restore writes the same values the store carries
    fig$xax <- NULL; fig$yax <- NULL; fig$axisMode <- NULL
    state$panels$data_space[[panel_key]] <- fig
  }

  if (!is.null(values) && length(values)) {
    if (is.null(state$widget_store) || !is.list(state$widget_store))
      state$widget_store <- list()
    state$widget_store$values <- values
    state$widget_store$unset <- as.character(
      state$widget_store$unset %||% names(values[vapply(values, is.null, logical(1))]))
  }

  # derive the selection-bus records (todo 2.2): the historic
  # selectByCorner display-authority flag is the best available origin
  # signal in a v1 snapshot; everything else restores as origin "restore"
  sel <- state$selection
  mk_rec <- function(ids, corner) list(
    ids = as.character(ids %||% character(0)),
    clicked = character(0),
    origin = if (isTRUE(corner)) "corner" else "restore",
    anchor = NULL, mirror = TRUE, epoch = 0L)
  if (is.null(sel$records) || !is.list(sel$records))
    sel$records <- list(
      feature = mk_rec(sel$features, state$panels$data_space$eset_fdata_fig$selectByCorner),
      sample = mk_rec(sel$samples, state$panels$data_space$eset_pdata_fig$selectByCorner)
    )
  state$selection <- sel
  state$schema_version <- 2L
  state
}

.state_assign <- function(values, keys, vals) {
  for (i in seq_along(keys))
    if (!is.null(vals[[i]])) values[[keys[[i]]]] <- vals[[i]]
  values
}

#' Validate a state object against the active dataset
#'
#' Validation warns on incompatible metadata rather than blocking
#' restoration; individual modules remain responsible for validating their
#' own choices. ALL warnings are collected and returned as the `warnings`
#' attribute (one character vector), so callers can surface every
#' adjustment instead of only the first (todo 2.7). The semantic selection
#' and the table-row mirrors are INTERSECTED with the current dataset ids
#' (reported, not dropped: the selection may legitimately name absent ids
#' after a dataset revision, but the tables would otherwise index out of
#' bounds on re-exported datasets).
#'
#' @param state A migrated state object.
#' @param dataset Active dataset object.
#' @param dataset_id Active dataset identifier.
#'
#' @keywords internal
#' @rdname app_state_helpers
validate_app_state <- function(state, dataset = NULL, dataset_id = NA_character_) {
  notes <- character()
  withCallingHandlers(
    state <- .validate_app_state_inner(state, dataset, dataset_id),
    warning = function(w) {
      notes <<- c(notes, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  attr(state, "warnings") <- notes
  state
}

.validate_app_state_inner <- function(state, dataset = NULL, dataset_id = NA_character_) {
  state <- migrate_app_state(state)
  if (!identical(state$format, APP_STATE_FORMAT))
    warning("Unknown snapshot format; state restoration may be incomplete.")

  current <- dataset_state(dataset, id = dataset_id)
  old <- state$dataset
  if (!is.null(old) && !is.na(old$fingerprint) && !is.na(current$fingerprint) &&
      !identical(old$fingerprint, current$fingerprint)) {
    warning(paste(
      "This snapshot was saved from a different version of the dataset",
      "(features, samples or annotation columns changed).",
      "Content that no longer exists will be ignored."))
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
    warning(sprintf(
      "%d of %d selected features are not in this dataset and will be ignored.",
      length(missing_f), length(sel$features)))
  if (length(missing_s))
    warning(sprintf(
      "%d of %d selected samples are not in this dataset and will be ignored.",
      length(missing_s), length(sel$samples)))
  # intersect the selection (and its bus records) with the current ids so
  # downstream consumers never see ids the dataset cannot resolve
  sel$features <- intersect(sel$features, feat)
  sel$samples <- intersect(sel$samples, samp)
  if (is.list(sel$records)) {
    if (is.list(sel$records$feature))
      sel$records$feature$ids <- intersect(as.character(sel$records$feature$ids %||% character(0)), feat)
    if (is.list(sel$records$sample))
      sel$records$sample$ids <- intersect(as.character(sel$records$sample$ids %||% character(0)), samp)
  }
  state$selection <- sel

  # intersect the table-row mirrors with the current ids (positional
  # mirrors would otherwise index out of bounds on re-exported/revised
  # datasets)
  intersect_ids <- function(ids, universe) {
    if (is.null(ids) || isTRUE(ids)) return(ids)
    keep <- as.character(ids) %in% universe
    if (!all(keep))
      warning(sprintf("%d mirrored ids are not in this dataset and were dropped.",
                      sum(!keep)))
    as.character(ids)[keep]
  }
  if (is.list(sel$records)) {
    if (is.list(sel$records$feature))
      sel$records$feature$mirror <- intersect_ids(sel$records$feature$mirror, feat)
    if (is.list(sel$records$sample))
      sel$records$sample$mirror <- intersect_ids(sel$records$sample$mirror, samp)
    state$selection$records <- sel$records
  }
  ds_panels <- state$panels$data_space
  if (is.list(ds_panels)) {
    if (!is.null(ds_panels$eset_fdata_tabrows) && !isTRUE(ds_panels$eset_fdata_tabrows))
      ds_panels$eset_fdata_tabrows <- intersect_ids(ds_panels$eset_fdata_tabrows, feat)
    if (!is.null(ds_panels$eset_pdata_tabrows) && !isTRUE(ds_panels$eset_pdata_tabrows))
      ds_panels$eset_pdata_tabrows <- intersect_ids(ds_panels$eset_pdata_tabrows, samp)
    state$panels$data_space <- ds_panels
  }
  state
}

#' Sanitize a user-supplied snapshot name
#'
#' The returned name contains only portable filename characters and cannot
#' represent a path or parent-directory reference. Unicode letters and
#' numbers are preserved (Perl `\p{L}`/`\p{N}` classes; todo 2.7); names
#' are capped at 80 characters so the generated file name stays comfortably
#' below the usual 255-byte file system limit even with the dataset-id
#' prefix (todo 1.5).
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
  name <- gsub("[^\\p{L}\\p{N}_.-]+", "_", name, perl = TRUE)
  name <- gsub("^[_.-]+", "", name, perl = TRUE)
  name <- gsub("[_.-]+$", "", name, perl = TRUE)
  if (!nzchar(name) || name == "." || name == "..")
    name <- fallback
  while (nchar(name, type = "bytes") > 80L && nzchar(name)) {
    name <- substr(name, 1L, max(1L, nchar(name, type = "chars") - 1L))
  }
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
#' @param widget_store Optional canonical widget-store snapshot (schema 2);
#'   augmented with the explicit `unset` key list when given.
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
    widget_store = NULL,
    label = NA_character_,
    package_version = as.character(utils::packageVersion("omicsViewer"))
) {
  data_status <- as.list(data_status)
  result_status <- as.list(result_status)
  if (is.list(widget_store) && is.list(widget_store$values) &&
      is.null(widget_store$unset)) {
    widget_store$unset <- names(
      widget_store$values[vapply(widget_store$values, is.null, logical(1))])
  }
  new_app_state(
    dataset_id = dataset_id,
    dataset = dataset,
    selection = list(features = selected_features,
                     samples = selected_samples,
                     records = data_status$eset_selection_records),
    app = list(
      data_active_tab = data_status$eset_active_tab,
      analysis_active_tab = result_status$analyst_active_tab
    ),
    data_space = .sanitize_widget_state(data_status),
    result_space = .sanitize_widget_state(result_status),
    widget_store = widget_store,
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
