#' Internal capability registry, discovery, and WP8 semantic tools
#'
#' One source of truth (AGENT_ACCURACY_PLAN.md sections 4 WP8/WP9 and 6.6)
#' for what the assistant can read and control. Capability records are
#' GENERATED: the widget half comes from the canonical widget-store
#' bindings (\code{\link{store_registry_view}}), the tool half from the
#' curated semantic-tool metadata below, so new widgets become
#' agent-discoverable by registering, and new semantic tools ship their
#' help text here from day one instead of being retrofitted.
#'
#' Every record carries: \code{id}, \code{panel}, \code{label}, \code{kind},
#' \code{description}/\code{help_text}, \code{writable},
#' \code{allowed_values} (bounded), and \code{operation} - the tool that
#' operates the capability (a semantic tool when one covers it, else
#' \code{set_widgets}).
#'
#' @keywords internal
#' @name agentCapabilityHelpers
NULL

############################################################################
### [1] semantic-tool metadata (the shared help-text structure)
############################################################################

# Tool records: id/panel/tier/label/description/help_text. tier is one of
# "read" (no side effects), "semantic" (curated write tool),
# "generic" (registry-driven widget tier), "discovery".
.agent_tool_metadata <- list(
  list(
    id = "get_omics_viewer_state", panel = "Global state", tier = "read",
    label = "Read application state",
    help_text = paste("Compact overview of the dataset, active tabs, selections,",
                      "quick views, and current scatter axes; optional sections",
                      "add the annotation catalog, quick-view records, panel",
                      "state, or the figure grammar in full detail.")
  ),
  list(
    id = "search_annotations", panel = "Annotations", tier = "read",
    label = "Search annotations",
    help_text = paste("Case-insensitive literal search over feature/sample IDs,",
                      "annotation column names, and bounded matching values.")
  ),
  list(
    id = "summarize_annotation", panel = "Annotations", tier = "read",
    label = "Summarize one annotation column",
    help_text = "Numeric quantiles or categorical counts for one exact annotation column."
  ),
  list(
    id = "set_omics_viewer_state", panel = "Global state", tier = "semantic",
    label = "Set tabs and selections",
    help_text = paste("Change the active data/analysis tab or the semantic",
                      "feature/sample selections after an explicit user request.")
  ),
  list(
    id = "set_scatter_view", panel = "Data space", tier = "semantic",
    label = "Set scatter view",
    help_text = paste("Point a feature-space or sample-space scatter at a quick",
                      "view or two exact annotation columns; never changes",
                      "selections or the quick/custom display mode.")
  ),
  list(
    id = "set_enrichment_parameters", panel = "Analysis space", tier = "semantic",
    label = "Set enrichment parameters",
    help_text = paste("Point the ORA or fGSEA panel at one collapse/ranking",
                      "annotation column and optionally select a pathway row",
                      "in the results table; opens the panel's tab. ORA tests",
                      "the currently selected features; fGSEA ranks all",
                      "features by a numeric column.")
  ),
  list(
    id = "set_table_view", panel = "Data space", tier = "semantic",
    label = "Set table view",
    help_text = paste("Configure the feature, sample, or expression table:",
                      "shown columns, multi-row selection, per-column filter",
                      "patterns, and the visible page; opens the table's tab.")
  ),
  list(
    id = "create_figure", panel = "Figures", tier = "semantic",
    label = "Create a figure",
    help_text = paste("Create a static ggplot2 figure from a template",
                      "(volcano/scatter/boxplot/histogram) or a full declarative",
                      "spec; the result carries the spec for later revision.")
  ),
  list(
    id = "update_figure", panel = "Figures", tier = "semantic",
    label = "Revise a figure",
    help_text = "Create a revised figure from the complete spec returned by the previous call."
  ),
  list(
    id = "list_widgets", panel = "Widgets", tier = "generic",
    label = "List widgets",
    help_text = paste("List every user-editable widget with canonical id, kind,",
                      "allowed values, and current value; optional prefix filter.")
  ),
  list(
    id = "get_widget", panel = "Widgets", tier = "generic",
    label = "Describe one widget",
    help_text = "Describe one widget by exact canonical id with its current value."
  ),
  list(
    id = "set_widgets", panel = "Widgets", tier = "generic",
    label = "Set widgets",
    help_text = paste("Apply a JSON-object patch of canonical widget ids to",
                      "values for controls the semantic tools do not cover.")
  ),
  list(
    id = "search_ui_capabilities", panel = "Capabilities", tier = "discovery",
    label = "Search capabilities",
    help_text = paste("Search by meaning across everything controllable:",
                      "semantic tools, figures, and every user-editable widget",
                      "with its panel, help text, and operating tool.")
  ),
  list(
    id = "get_ui_capability", panel = "Capabilities", tier = "discovery",
    label = "Describe one capability",
    help_text = "Describe one capability by exact id (a widget id or a semantic tool id)."
  )
)

#' Curated semantic-tool capability records
#'
#' The tool half of the shared metadata structure: every model-facing tool
#' registers its label, panel, and operating contract here (WP9's
#' one-source-of-truth principle; WP8 tools ship here from day one).
#'
#' @return A list of tool records.
#' @keywords internal
#' @rdname agentCapabilityHelpers
agent_tool_capabilities <- function() {
  lapply(.agent_tool_metadata, function(t) list(
    id = paste0("tool:", t$id),
    panel = t$panel,
    label = t$label,
    kind = paste0("tool_", t$tier),
    description = t$help_text,
    help_text = t$help_text,
    writable = t$tier %in% c("semantic", "generic"),
    operation = t$id,
    allowed_values = NULL,
    depends_on = NULL
  ))
}

############################################################################
### [2] widget-derived records (generated from bindings)
############################################################################

.agent_widget_panel_label <- function(id) {
  top <- strsplit(id, ".", fixed = TRUE)[[1]][1]
  switch(top,
    dataspace = "Data space",
    resultspace = "Analysis space",
    top)
}

# Which semantic tool covers a widget namespace? Everything else falls to
# the generic set_widgets tier.
.agent_widget_operation <- function(id) {
  if (startsWith(id, "dataspace.feature_space.") ||
      startsWith(id, "dataspace.sample_space."))
    return("set_scatter_view")
  if (startsWith(id, "resultspace.ora.") ||
      startsWith(id, "resultspace.fgsea."))
    return("set_enrichment_parameters")
  if (startsWith(id, "dataspace.tab_"))
    return("set_table_view")
  if (identical(id, "dataspace.active_tab") ||
      identical(id, "resultspace.analyst_tab"))
    return("set_omics_viewer_state")
  "set_widgets"
}

#' Generate every capability record
#'
#' Widget records come from the agent-visible store registry (labels, help
#' text, allowed values, dependencies); tool records from
#' \code{\link{agent_tool_capabilities}}. Together they answer "what can
#' you change in this app?" without any hand-maintained list.
#'
#' @param store Canonical widget store (or child view). NULL yields the
#'   tool records only.
#' @return A list of capability records.
#' @keywords internal
#' @rdname agentCapabilityHelpers
agent_capability_records <- function(store = NULL) {
  out <- list()
  if (!is.null(store)) {
    for (r in store_registry_view(store)) {
      help <- r$help %||% ""
      out[[length(out) + 1L]] <- list(
        id = r$id,
        panel = .agent_widget_panel_label(r$id),
        label = if (nzchar(r$label %||% "")) r$label else r$id,
        kind = r$kind,
        description = help,
        help_text = help,
        writable = TRUE,
        operation = .agent_widget_operation(r$id),
        allowed_values = r$allowed_values,
        depends_on = r$depends_on
      )
    }
  }
  c(out, agent_tool_capabilities())
}

#' Search capabilities by meaning
#'
#' Case-insensitive literal substring match over ids, labels, panels,
#' operations, help text, and descriptions. Bounded result set; the counts
#' keep the model informed about truncation.
#'
#' @param store Canonical widget store (or child view); may be NULL.
#' @param query Non-empty search query (1-128 characters).
#' @param max_results Maximum records returned (1-50).
#' @return \code{query}, \code{capability_count}, \code{match_count},
#'   \code{truncated}, and \code{capabilities} records.
#' @keywords internal
#' @rdname agentCapabilityHelpers
agent_capability_search <- function(store = NULL, query, max_results = 20L) {
  query <- .agent_trim_scalar(query)
  if (agent_sentinel_string(query) || !nzchar(query))
    stop("A non-empty capability search query is required.")
  if (nchar(query) > 128L)
    stop("Capability search query must be at most 128 characters.")
  max_results <- suppressWarnings(as.integer(max_results)[1])
  if (is.na(max_results) || max_results < 1L || max_results > 50L)
    stop("max_results must be an integer from 1 through 50.")

  records <- agent_capability_records(store)
  ql <- tolower(query)
  hit <- function(r) {
    hay <- tolower(c(r$id, r$label, r$panel, r$operation,
                     r$description %||% "", r$help_text %||% ""))
    any(grepl(ql, hay, fixed = TRUE))
  }
  hits <- Filter(hit, records)
  list(
    query = query,
    capability_count = length(records),
    match_count = length(hits),
    truncated = length(hits) > max_results,
    capabilities = if (length(hits))
      utils::head(hits, max_results) else list()
  )
}

#' Describe one capability by exact id
#'
#' Accepts a widget canonical id, a semantic tool name, or the
#' \code{"tool:"}-prefixed record id; unknown ids are rejected with
#' closest-match suggestions (WP2).
#'
#' @param store Canonical widget store (or child view); may be NULL.
#' @param id Exact capability id.
#' @return The matching capability record.
#' @keywords internal
#' @rdname agentCapabilityHelpers
agent_capability_get <- function(store = NULL, id) {
  id <- .agent_trim_scalar(id)
  if (!nzchar(id))
    stop("A capability id is required.")
  records <- agent_capability_records(store)
  keys <- vapply(records, function(r) sub("^tool:", "", r$id), character(1))
  want <- sub("^tool:", "", id)
  hit <- which(keys == want)
  if (!length(hit))
    stop("Unknown capability id: ", id, ".",
         .agent_suggest_text(want, keys),
         " Discover ids with search_ui_capabilities.")
  records[[hit[[1L]]]]
}

#' Compact capability counts for the state overview
#'
#' WP9 context hygiene: the \code{get_omics_viewer_state} overview carries
#' counts only (capability total, panels, semantic tool ids), never record
#' contents - discovery goes through \code{search_ui_capabilities}. Unlike
#' \code{\link{agent_capability_records}} this never invokes choices
#' providers, so it stays cheap.
#'
#' @param store Canonical widget store (or child view).
#' @return A bounded named list, or NULL when no store is given.
#' @keywords internal
#' @rdname agentCapabilityHelpers
agent_capability_summary <- function(store) {
  if (is.null(store))
    return(NULL)
  root <- if (is.null(store$parent)) store else store$parent
  ids <- names(root$bindings)
  ids <- ids[vapply(ids, function(k) isTRUE(root$bindings[[k]]$agent_writable),
                    logical(1))]
  panels <- sort(unique(vapply(ids, .agent_widget_panel_label, character(1))))
  semantic <- vapply(.agent_tool_metadata, function(t)
    if (identical(t$tier, "semantic")) t$id else "", character(1))
  list(
    capability_count = length(ids),
    panels = panels,
    semantic_tools = semantic[nzchar(semantic)]
  )
}

############################################################################
### [3] WP8 semantic-tool normalizers (thin validate + store patches)
############################################################################

# method -> analysis-space tab label (mirrors the analyst navbar choices)
.agent_enrichment_tabs <- c(ora = "ORA", fgsea = "fGSEA")

#' Validate an assistant-proposed enrichment update
#'
#' Thin semantic tier over the widget store: maps \code{method} to the
#' ORA/fGSEA panel, splits the collapse/ranking triple, and returns a
#' store patch (full canonical ids) plus the tab switch. Value-level
#' validation (cascade choices, pathway rows) happens transactionally in
#' \code{\link{store_apply}}, whose errors already carry suggestions.
#'
#' @param update Named list with \code{method} and optional
#'   \code{collapse} (full \code{Category|Subcategory|Variable} feature
#'   column) and \code{selected_pathway} (gene-set id from the results
#'   table).
#' @param feature_data Feature metadata (needs the \code{GS} attribute).
#' @return \code{list(method, tab, patch)}.
#' @keywords internal
#' @rdname agentCapabilityHelpers
agent_normalize_enrichment_update <- function(update, feature_data) {
  if (!is.list(update))
    stop("Enrichment update must be a named list.")
  allowed <- c("method", "collapse", "selected_pathway")
  unknown <- setdiff(names(update), allowed)
  if (length(unknown))
    stop("Unknown enrichment field(s): ", paste(unknown, collapse = ", "))

  method <- .agent_trim_scalar(update$method)
  if (!method %in% names(.agent_enrichment_tabs))
    stop("Unknown enrichment method: ", method, ".",
         .agent_suggest_text(method, names(.agent_enrichment_tabs)),
         " Use 'ora' or 'fgsea'.")
  tab <- .agent_enrichment_tabs[[method]]

  if (is.null(feature_data))
    stop("No feature metadata are loaded.")
  if (is.null(attr(feature_data, "GS")))
    stop("This dataset carries no gene-set annotations; ORA/fGSEA are unavailable.")

  prefix <- paste0("resultspace.", method)
  patch <- list()
  collapse <- .agent_nullable_scalar(update$collapse)
  if (nzchar(collapse)) {
    parts <- strsplit(collapse, "|", fixed = TRUE)[[1]]
    if (length(parts) != 3L || any(!nzchar(parts)))
      stop("collapse must be a full 'Category|Subcategory|Variable' feature",
           " annotation name, got: ", collapse, ".",
           " Use search_annotations to discover exact column names.")
    if (!collapse %in% colnames(feature_data))
      stop("Unknown feature annotation for the ", method, " input: ",
           collapse, ".",
           .agent_suggest_text(collapse, colnames(feature_data)),
           " Use search_annotations to confirm exact column names.")
    patch[[paste0(prefix, ".xax_analysis")]] <- parts[[1]]
    patch[[paste0(prefix, ".xax_subset")]] <- parts[[2]]
    patch[[paste0(prefix, ".xax_variable")]] <- parts[[3]]
  }
  pathway <- .agent_nullable_scalar(update$selected_pathway)
  if (nzchar(pathway))
    patch[[paste0(prefix, ".selected_row")]] <- pathway
  if (!length(patch))
    stop("Enrichment update contains no changes.")
  # make the effect visible: open the panel's tab (scatter-tool precedent)
  patch[["resultspace.analyst_tab"]] <- tab

  list(method = method, tab = tab, patch = patch)
}

# table -> (store key prefix, data-space tab label)
.agent_table_view_prefixes <- c(
  feature_table = "dataspace.tab_feature",
  sample_table = "dataspace.tab_pheno",
  expression_table = "dataspace.tab_expr"
)
.agent_table_view_tabs <- c(
  feature_table = "Feature table",
  sample_table = "Sample table",
  expression_table = "Expression"
)

# Sentinel-aware scalar-absent test: NULL, NA, and the literal sentinel
# strings all count as omitted (the AGENT_SENTINEL_STRINGS convention).
.agent_capability_absent <- function(x) {
  is.null(x) || (length(x) == 1L && is.na(x)) || agent_sentinel_string(x)
}

# Lenient logical coercion for boolean tool arguments (providers may send
# "true"/"false" strings).
.agent_capability_logical <- function(x) {
  if (is.logical(x) && length(x) == 1L && !is.na(x))
    return(x)
  txt <- tolower(trimws(suppressWarnings(as.character(x)[1])))
  if (is.na(txt)) return(FALSE)
  txt == "true" || txt == "1"
}

#' Normalize an optional column-filter mapping argument
#'
#' Accepts a named list (ellmer object) or a JSON-object string; sentinel
#' strings and empty inputs return NULL (treated as omitted - clearing is
#' the explicit \code{clear_filters} argument).
#'
#' @param value Named list or JSON object string.
#' @return Named character vector, or NULL when absent.
#' @keywords internal
#' @rdname agentCapabilityHelpers
.agent_normalize_column_filters <- function(value) {
  if (is.null(value))
    return(NULL)
  if (is.character(value) && length(value) == 1L) {
    txt <- trimws(value)
    if (!nzchar(txt) || agent_sentinel_string(txt))
      return(NULL)
    parsed <- tryCatch(
      jsonlite::fromJSON(txt, simplifyVector = FALSE),
      error = function(e)
        stop("column_filters is not valid JSON: ", conditionMessage(e))
    )
    value <- parsed
  }
  if (!is.list(value) && !is.character(value))
    stop("column_filters must be a JSON object of column name to search pattern.")
  if (is.list(value)) {
    if (!length(value))
      return(NULL)
    nms <- names(value)
    vals <- vapply(value, function(v) {
      v <- suppressWarnings(as.character(v)[1])
      if (is.na(v)) "" else v
    }, character(1), USE.NAMES = FALSE)
    if (is.null(nms) || any(is.na(nms)) || any(!nzchar(nms)))
      stop("column_filters must be a JSON object with named column keys.")
    value <- vals
    names(value) <- nms
  }
  if (!length(value))
    return(NULL)
  value
}

#' Validate an assistant-proposed table-view update
#'
#' Thin semantic tier over the widget store: resolves the table name to
#' the \code{dataspace.tab_*} namespace, normalizes every present field,
#' and returns the store patch (columns and multi-selection validate
#' transactionally against live choices in \code{\link{store_apply}}).
#'
#' @param update Named list with \code{table} plus optional
#'   \code{columns}, \code{multi_selection}, \code{column_filters},
#'   \code{page}, and \code{clear_filters}.
#' @return \code{list(table, tab, patch)}.
#' @keywords internal
#' @rdname agentCapabilityHelpers
agent_normalize_table_view_update <- function(update) {
  if (!is.list(update))
    stop("Table-view update must be a named list.")
  allowed <- c("table", "columns", "multi_selection", "column_filters",
               "page", "clear_filters")
  unknown <- setdiff(names(update), allowed)
  if (length(unknown))
    stop("Unknown table-view field(s): ", paste(unknown, collapse = ", "))

  table <- .agent_trim_scalar(update$table)
  if (!table %in% names(.agent_table_view_prefixes))
    stop("Unknown table: ", table, ".",
         .agent_suggest_text(table, names(.agent_table_view_prefixes)),
         " Available tables: ",
         paste(names(.agent_table_view_prefixes), collapse = ", "), ".")
  prefix <- .agent_table_view_prefixes[[table]]
  tab <- .agent_table_view_tabs[[table]]

  patch <- list()
  columns <- update$columns
  if (!is.null(columns)) {
    if (is.list(columns))
      columns <- unlist(columns, use.names = FALSE)
    if (is.character(columns) && length(columns) == 1L &&
        agent_sentinel_string(columns))
      columns <- NULL
    else if (!is.character(columns) || any(is.na(columns)) ||
             any(!nzchar(columns)))
      stop("columns must be an array of exact column names (at least one).")
    else if (!length(columns))
      stop("columns cannot be empty; at least one column must remain shown.")
    else
      patch[[paste0(prefix, ".columns")]] <- columns
  }
  if (!is.null(update$multi_selection) &&
      !.agent_capability_absent(update$multi_selection))
    patch[[paste0(prefix, ".multi_selection")]] <- update$multi_selection

  clear <- !.agent_capability_absent(update$clear_filters) &&
    isTRUE(.agent_capability_logical(update$clear_filters))
  filters <- .agent_normalize_column_filters(update$column_filters)
  if (clear)
    patch[[paste0(prefix, ".column_filters")]] <-
      setNames(character(0), character(0))
  else if (!is.null(filters))
    patch[[paste0(prefix, ".column_filters")]] <- filters

  page <- update$page
  if (!.agent_capability_absent(page)) {
    page_num <- suppressWarnings(as.integer(page)[1])
    if (is.na(page_num) || page_num < 1L)
      stop("page must be an integer of at least 1.")
    patch[[paste0(prefix, ".page")]] <- page_num
  }

  if (!length(patch))
    stop("Table-view update contains no changes.")
  # make the effect visible: open the table's tab (scatter-tool precedent)
  patch[["dataspace.active_tab"]] <- tab

  list(table = table, tab = tab, patch = patch)
}
