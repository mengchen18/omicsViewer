#' Internal AI-assistant state and validation helpers
#'
#' These helpers provide the compact, JSON-like state bridge between the Shiny
#' application and an optional ellmer-powered assistant. They deliberately
#' return bounded metadata and summaries rather than expression matrices or API
#' credentials.
#'
#' @return Provider configuration, compact state, annotation summaries, or a
#'   validated assistant update.
#'
#' @keywords internal
#' @name agentAssistantHelpers
NULL

.agent_trim_scalar <- function(x, fallback = "") {
  if (is.null(x) || length(x) == 0 || is.na(x))
    return(fallback)
  out <- trimws(as.character(x)[1])
  if (is.na(out) || !nzchar(out)) fallback else out
}

.agent_js_string <- function(x) {
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub("'", "\\\'", x, fixed = TRUE)
  x <- gsub("[\r\n]", "", x)
  paste0("'", x, "'")
}

.agent_shorten <- function(x, limit = 160L) {
  x <- as.character(x)
  too_long <- !is.na(x) & nchar(x, type = "chars", allowNA = TRUE) > limit
  if (any(too_long, na.rm = TRUE)) {
    x[too_long] <- paste0(
      substr(x[too_long], 1L, limit),
      " ... [truncated]"
    )
  }
  x
}

.agent_column_as_character <- function(x) {
  if (is.list(x)) {
    vapply(
      x,
      function(y) {
        if (is.null(y)) return("")
        paste(as.character(y), collapse = "; ")
      },
      character(1)
    )
  } else {
    as.character(x)
  }
}

#' Read the session request limit
#'
#' The limit applies to all provider requests in one Shiny session, including
#' intermediate requests generated while a model resolves tool calls.
#'
#' @return Integer request limit from 1 through 200. Default: 40.
#' @keywords internal
#' @rdname agentAssistantHelpers
agent_request_limit <- function() {
  value <- suppressWarnings(as.integer(Sys.getenv("OMICSVIEWER_LLM_MAX_REQUESTS")[1]))
  if (length(value) != 1L || is.na(value))
    value <- 40L
  max(1L, min(200L, value))
}

#' Read server-side assistant provider configuration
#'
#' Environment variables are intentionally limited to provider selection, model
#' selection, an optional OpenAI-compatible base URL, and credentials. The API
#' key is never placed in Shiny UI state, snapshots, tool results, or logs.
#'
#' @return A list with provider, model, base URL, API key, configuration
#'   status, and credential source. An empty API key means that the user must
#'   supply one in the current session.
#' @keywords internal
#' @rdname agentAssistantHelpers
agent_environment_config <- function() {
  provider <- tolower(.agent_trim_scalar(Sys.getenv("OMICSVIEWER_LLM_PROVIDER")))
  if (!provider %in% c("openai", "anthropic")) {
    provider <- if (nzchar(Sys.getenv("ANTHROPIC_API_KEY"))) "anthropic" else "openai"
  }

  api_key <- .agent_trim_scalar(Sys.getenv("OMICSVIEWER_LLM_API_KEY"))
  key_source <- if (nzchar(api_key)) "environment" else ""
  if (!nzchar(api_key)) {
    env_name <- if (provider == "anthropic") "ANTHROPIC_API_KEY" else "OPENAI_API_KEY"
    api_key <- .agent_trim_scalar(Sys.getenv(env_name))
    if (nzchar(api_key)) key_source <- "environment"
  }

  model <- .agent_trim_scalar(Sys.getenv("OMICSVIEWER_LLM_MODEL"))
  if (!nzchar(model)) {
    env_name <- if (provider == "anthropic") "OMICSVIEWER_ANTHROPIC_MODEL" else "OMICSVIEWER_OPENAI_MODEL"
    model <- .agent_trim_scalar(Sys.getenv(env_name))
  }

  base_url <- .agent_trim_scalar(Sys.getenv("OMICSVIEWER_LLM_BASE_URL"))
  if (!nzchar(base_url)) {
    env_name <- if (provider == "anthropic") "OMICSVIEWER_ANTHROPIC_BASE_URL" else "OMICSVIEWER_OPENAI_BASE_URL"
    base_url <- .agent_trim_scalar(Sys.getenv(env_name))
  }

  # Apply the same validation to administrator-provided settings as to settings
  # entered in the modal. In particular, do not silently send a session key to
  # a remote non-HTTPS endpoint.
  validated <- agent_validate_provider_config(
    provider = provider,
    model = model,
    api_key = api_key,
    base_url = base_url
  )
  list(
    provider = validated$provider,
    model = validated$model,
    base_url = validated$base_url,
    api_key = validated$api_key,
    configured = validated$configured,
    key_source = key_source
  )
}

#' Validate user-supplied assistant provider settings
#'
#' @param provider Single provider name: openai or anthropic.
#' @param model Optional model name. An empty value uses the provider default.
#' @param api_key Optional session API key.
#' @param base_url Optional HTTPS or local HTTP API endpoint.
#'
#' @return A normalized provider-configuration list without changing the
#'   credential source.
#' @keywords internal
#' @rdname agentAssistantHelpers
agent_validate_provider_config <- function(provider, model = "", api_key = "",
                                            base_url = "") {
  provider <- .agent_trim_scalar(provider)
  if (!provider %in% c("openai", "anthropic"))
    stop("LLM provider must be 'openai' or 'anthropic'.")

  model <- .agent_trim_scalar(model)
  api_key <- .agent_trim_scalar(api_key)
  base_url <- .agent_trim_scalar(base_url)

  if (nchar(model) > 128 || grepl("[\r\n]", model))
    stop("Model name is too long or contains control characters.")
  if (nchar(api_key) > 4096 || any(grepl("[\r\n[:space:]]", api_key)))
    stop("API key is invalid.")

  if (nzchar(base_url)) {
    if (nchar(base_url) > 2048 || grepl("[\r\n]", base_url))
      stop("API base URL is too long or contains control characters.")
    local_http <- grepl(
      "^http://(localhost|127\\.0\\.0\\.1)([/:]|$)",
      base_url
    )
    if (!startsWith(base_url, "https://") && !local_http)
      stop("API base URL must use HTTPS, or local HTTP for this machine only.")
  }

  list(
    provider = provider,
    model = model,
    api_key = api_key,
    base_url = base_url,
    configured = nzchar(api_key)
  )
}

#' Convert quick-view definitions to bounded records
#'
#' @param views A data.frame from \code{\link{prepare_quick_views}}, or NULL.
#' @return A list of JSON-like quick-view records.
#' @keywords internal
#' @rdname agentAssistantHelpers
agent_quick_view_records <- function(views) {
  if (is.null(views) || !is.data.frame(views) || nrow(views) == 0)
    return(list())
  lapply(seq_len(nrow(views)), function(i) as.list(views[i, , drop = FALSE]))
}

.agent_bound_value <- function(x, limit = 100L) {
  if (is.character(x)) {
    too_long <- !is.na(x) & nchar(x, type = "chars", allowNA = TRUE) > 160L
    if (any(too_long, na.rm = TRUE)) {
      x <- .agent_shorten(x)
    }
    if (length(x) > limit) {
      return(list(
        count = length(x),
        values = utils::head(x, limit),
        truncated = TRUE
      ))
    }
  }
  if (is.list(x)) {
    return(lapply(x, .agent_bound_value, limit = limit))
  }
  x
}

#' Build a compact model-facing application state
#'
#' @param state Versioned application state from \code{\link{build_app_state}}.
#' @param annotations Annotation catalog created by
#'   \code{\link{agent_annotation_catalog}}.
#' @param quick_views Named list containing feature and sample quick views.
#' @param available_tabs Named list of valid data-space and analysis-space tabs.
#'
#' @return A JSON-like list containing identifiers, active tabs, semantic
#'   selections, widget state, annotation metadata, and bounded quick-view
#'   definitions. Expression values, fingerprints, credentials, and timestamps
#'   are intentionally omitted.
#' @keywords internal
#' @rdname agentAssistantHelpers
agent_compact_state <- function(state, annotations = NULL, quick_views = NULL,
                                available_tabs = NULL, figure_grammar = NULL) {
  if (is.null(state))
    return(NULL)

  if (is.null(available_tabs))
    available_tabs <- list(data_space = character(), analysis_space = character())

  feature_views <- utils::head(agent_quick_view_records(quick_views$feature), 50L)
  sample_views <- utils::head(agent_quick_view_records(quick_views$sample), 50L)
  bound_selection <- function(ids) {
    list(
      count = length(ids),
      ids = utils::head(ids, 100L),
      truncated = length(ids) > 100L
    )
  }

  list(
    dataset = list(
      id = state$dataset$id,
      class = state$dataset$class,
      dimensions = state$dataset$dimensions
    ),
    active_tabs = list(
      data_space = state$app$data_active_tab,
      analysis_space = state$app$analysis_active_tab
    ),
    selection = list(
      features = bound_selection(state$selection$features),
      samples = bound_selection(state$selection$samples)
    ),
    available_tabs = available_tabs,
    annotations = annotations,
    quick_views = list(
      feature = feature_views,
      sample = sample_views
    ),
    panels = .agent_bound_value(state$panels),
    figure_grammar = figure_grammar,
    state_policy = state$policy
  )
}

#' Describe annotation structure without returning row-level data
#'
#' @param feature_data Feature metadata data.frame.
#' @param sample_data Sample metadata data.frame.
#' @return A named list containing dimensions, row ID labels, and column-level
#'   type/missingness/cardinality summaries.
#' @keywords internal
#' @rdname agentAssistantHelpers
agent_annotation_catalog <- function(feature_data, sample_data) {
  catalog_one <- function(x, row_label) {
    if (is.null(x)) {
      return(list(row_label = row_label, rows = 0L, columns = list()))
    }
    ids <- rownames(x)
    all_columns <- colnames(x)
    included_columns <- utils::head(all_columns, 200L)
    columns <- lapply(included_columns, function(nm) {
      values <- x[[nm]]
      missing <- if (is.list(values)) !lengths(values) else is.na(values)
      type <- if (is.numeric(values)) {
        "numeric"
      } else if (is.factor(values)) {
        "factor"
      } else if (is.logical(values)) {
        "logical"
      } else if (is.list(values)) {
        "list"
      } else {
        "character"
      }
      list(
        name = nm,
        type = type,
        missing = sum(missing),
        unique_count = length(unique(values[!missing]))
      )
    })
    names(columns) <- included_columns
    list(
      row_label = row_label,
      rows = nrow(x),
      row_id_count = if (is.null(ids)) 0L else length(ids),
      column_count = length(all_columns),
      columns_truncated = length(all_columns) > length(included_columns),
      columns = columns
    )
  }

  list(
    feature = catalog_one(feature_data, "featureId"),
    sample = catalog_one(sample_data, "sampleId")
  )
}

#' Search annotation IDs, column names, and bounded matching values
#'
#' @param space Either feature or sample.
#' @param query Case-insensitive query. Matching is literal, not regular-expression.
#' @param feature_data Feature metadata.
#' @param sample_data Sample metadata.
#' @param max_results Maximum number of row IDs to return.
#'
#' @return A JSON-like search result bounded by \code{max_results}. It never
#'   returns an entire annotation table or expression matrix.
#' @keywords internal
#' @rdname agentAssistantHelpers
agent_search_annotations <- function(space, query, feature_data, sample_data,
                                     max_results = 20L) {
  space <- match.arg(space, c("feature", "sample"))
  query <- .agent_trim_scalar(query)
  if (!nzchar(query))
    stop("Search query must not be empty.")
  if (nchar(query) > 128)
    stop("Search query must be at most 128 characters.")

  max_results <- suppressWarnings(as.integer(max_results)[1])
  if (is.na(max_results) || max_results < 1L || max_results > 50L)
    stop("max_results must be an integer from 1 through 50.")

  df <- if (space == "feature") feature_data else sample_data
  if (is.null(df))
    stop("No ", space, " metadata are loaded.")

  ids <- rownames(df)
  if (is.null(ids)) ids <- character()
  id_hits <- grepl(tolower(query), tolower(ids), fixed = TRUE)
  matching_ids <- ids[id_hits]

  column_names <- colnames(df)
  column_hit_flags <- grepl(tolower(query), tolower(column_names), fixed = TRUE)
  column_hits <- column_names[column_hit_flags]

  # Searching values is more expensive than searching column names. Keep the
  # scanned annotation payload bounded for very wide datasets.
  value_columns <- utils::head(column_names, 500L)
  value_hits <- list()
  value_match_column_count <- 0L
  for (nm in value_columns) {
    raw <- df[[nm]]
    values <- .agent_column_as_character(raw)
    valid <- !is.na(raw) & nzchar(values)
    if (is.list(raw)) valid <- lengths(raw) > 0 & nzchar(values)
    hit <- valid & grepl(tolower(query), tolower(values), fixed = TRUE)
    if (!any(hit)) next
    value_match_column_count <- value_match_column_count + 1L
    if (length(value_hits) >= 20L)
      next
    hit_ids <- ids[hit]
    hit_values <- unique(values[hit])
    n_examples <- min(5L, length(hit_ids))
    value_hits[[length(value_hits) + 1L]] <- list(
      column = nm,
      matching_rows = sum(hit),
      example_ids = if (n_examples) .agent_shorten(hit_ids[seq_len(n_examples)]) else character(),
      example_values = utils::head(.agent_shorten(hit_values), 5L)
    )
  }

  list(
    space = space,
    query = query,
    matching_id_count = length(matching_ids),
    matching_ids = utils::head(matching_ids, max_results),
    matching_column_count = sum(column_hit_flags),
    matching_columns = utils::head(column_hits, 50L),
    value_matches_truncated = value_match_column_count > length(value_hits),
    value_matches = value_hits
  )
}

#' Summarize one annotation column
#'
#' @param space Either feature or sample.
#' @param column Exact annotation column name.
#' @param feature_data Feature metadata.
#' @param sample_data Sample metadata.
#' @param max_values Maximum number of categorical values to return.
#'
#' @return A bounded numeric or categorical summary. Raw row-level values are
#'   not returned.
#' @keywords internal
#' @rdname agentAssistantHelpers
agent_summarize_annotation <- function(space, column, feature_data, sample_data,
                                       max_values = 12L) {
  space <- match.arg(space, c("feature", "sample"))
  column <- .agent_trim_scalar(column)
  max_values <- suppressWarnings(as.integer(max_values)[1])
  if (is.na(max_values) || max_values < 1L || max_values > 20L)
    stop("max_values must be an integer from 1 through 20.")

  df <- if (space == "feature") feature_data else sample_data
  if (is.null(df))
    stop("No ", space, " metadata are loaded.")
  if (!column %in% colnames(df))
    stop("Unknown ", space, " annotation column: ", column)

  values <- df[[column]]
  missing <- if (is.list(values)) !lengths(values) else is.na(values)
  observed <- values[!missing]

  out <- list(
    space = space,
    column = column,
    type = if (is.numeric(observed)) "numeric" else "categorical",
    total = length(values),
    missing = sum(missing),
    non_missing = length(observed)
  )

  if (!length(observed))
    return(out)

  if (is.numeric(observed)) {
    out$summary <- as.list(round(
      stats::quantile(observed, probs = c(0, .25, .5, .75, 1), na.rm = TRUE),
      digits = 6
    ))
    out$mean <- round(mean(observed), digits = 6)
    out$unique_count <- length(unique(observed))
  } else {
    observed <- .agent_shorten(.agent_column_as_character(observed), 120L)
    observed <- observed[nzchar(observed)]
    tab <- sort(table(observed), decreasing = TRUE)
    selected <- utils::head(tab, max_values)
    out$value_counts <- as.list(selected)
    out$other_count <- if (length(tab) > max_values) sum(tab[-seq_len(max_values)]) else 0L
    out$unique_count <- length(tab)
  }
  out
}

#' Validate an assistant-proposed application-state update
#'
#' @param update Named list with optional \code{data_space_tab},
#'   \code{analysis_space_tab}, \code{features}, and \code{samples}.
#' @param data_tabs Allowed data-space tab labels.
#' @param analysis_tabs Allowed analysis-space tab labels.
#' @param feature_ids Valid feature IDs.
#' @param sample_ids Valid sample IDs.
#' @param max_features Maximum number of selected features.
#' @param max_samples Maximum number of selected samples.
#'
#' @return A list containing only explicitly requested, validated changes. An
#'   empty character vector explicitly clears a selection.
#' @keywords internal
#' @rdname agentAssistantHelpers
agent_normalize_state_update <- function(update, data_tabs, analysis_tabs,
                                         feature_ids, sample_ids,
                                         max_features = 5000L,
                                         max_samples = 2000L) {
  if (!is.list(update))
    stop("Assistant state update must be a named list.")
  allowed <- c("data_space_tab", "analysis_space_tab", "features", "samples")
  unknown <- setdiff(names(update), allowed)
  if (length(unknown))
    stop("Unknown assistant state field(s): ", paste(unknown, collapse = ", "))

  normalize_ids <- function(x, ids, label, max_n) {
    if (is.null(x)) return(NULL)
    x <- as.character(x)
    if (any(is.na(x)) || any(!nzchar(x)))
      stop(label, " must contain non-empty IDs.")
    x <- x[!duplicated(x)]
    if (length(x) > max_n)
      stop("At most ", max_n, " ", label, "s can be selected by the assistant.")
    invalid <- setdiff(x, ids)
    if (length(invalid)) {
      example <- utils::head(invalid, 3L)
      stop("Unknown ", label, " ID(s): ", paste(example, collapse = ", "))
    }
    x
  }

  normalize_tab <- function(x, choices, label) {
    if (is.null(x)) return(NULL)
    x <- .agent_trim_scalar(x)
    if (!x %in% choices)
      stop("Unknown ", label, " tab: ", x)
    x
  }

  out <- list(
    data_space_tab = normalize_tab(update$data_space_tab, data_tabs, "data-space"),
    analysis_space_tab = normalize_tab(update$analysis_space_tab, analysis_tabs, "analysis-space"),
    features = normalize_ids(update$features, feature_ids, "feature", max_features),
    samples = normalize_ids(update$samples, sample_ids, "sample", max_samples)
  )
  nulls <- vapply(out, is.null, logical(1))
  if (all(nulls))
    stop("Assistant state update contains no changes.")
  out[!nulls]
}

#' Validate an assistant-proposed scatter view
#'
#' @param space Feature or sample scatter space.
#' @param quick_view_id Optional quick-view ID.
#' @param x_axis Optional exact annotation column for the X axis.
#' @param y_axis Optional exact annotation column for the Y axis.
#' @param quick_views Named list of feature/sample quick-view data.frames.
#' @param feature_columns Valid feature annotation columns.
#' @param sample_columns Valid sample annotation columns.
#'
#' @return A validated scatter-space and axis update.
#' @keywords internal
#' @rdname agentAssistantHelpers
agent_normalize_scatter_view <- function(space, quick_view_id = NULL,
                                         x_axis = NULL, y_axis = NULL,
                                         quick_views = NULL,
                                         feature_columns = character(),
                                         sample_columns = character()) {
  space <- match.arg(space, c("feature", "sample"))
  quick_view_id <- .agent_trim_scalar(quick_view_id)

  if (nzchar(quick_view_id)) {
    views <- quick_views[[space]]
    if (is.null(views) || !is.data.frame(views) || !nrow(views) ||
        !quick_view_id %in% views$id) {
      stop("Unknown ", space, " quick view: ", quick_view_id)
    }
    view <- views[views$id == quick_view_id, , drop = FALSE][1, ]
    return(list(space = space, mode = "quick", quick_view_id = view$id,
                x_axis = view$x, y_axis = view$y))
  }

  x_axis <- .agent_trim_scalar(x_axis)
  y_axis <- .agent_trim_scalar(y_axis)
  if (!nzchar(x_axis) || !nzchar(y_axis))
    stop("A scatter view requires either quick_view_id or both x_axis and y_axis.")

  columns <- if (space == "feature") feature_columns else sample_columns
  if (!x_axis %in% columns)
    stop("Unknown ", space, " X-axis annotation: ", x_axis)
  if (!y_axis %in% columns)
    stop("Unknown ", space, " Y-axis annotation: ", y_axis)

  axis_parts <- function(axis) {
    parts <- strsplit(axis, "|", fixed = TRUE)[[1]]
    if (length(parts) != 3L || any(!nzchar(parts)))
      return(NULL)
    parts
  }
  if (is.null(axis_parts(x_axis)) || is.null(axis_parts(y_axis)))
    stop("Custom axes must use the Category|Subcategory|Variable naming convention.")

  list(space = space, mode = "custom", quick_view_id = NULL,
       x_axis = x_axis, y_axis = y_axis)
}
