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

#' Suggest closest valid values for a rejected assistant input
#'
#' Used to make validation errors self-correcting: instead of a bare
#' rejection, the model receives the nearest valid candidates and can retry
#' without guessing. Scoring is deliberately cheap: exact case-insensitive
#' match first, then substring containment (either direction), then prefix,
#' and finally edit distance. Edit distance is skipped for very large
#' candidate sets (e.g. feature ID spaces on real datasets) where it would
#' dominate validation time.
#'
#' @param value Rejected value (single string).
#' @param candidates Character vector of valid values.
#' @param max Maximum number of suggestions to return.
#' @return A character vector of up to \code{max} suggestions (possibly empty).
#' @keywords internal
#' @rdname agentAssistantHelpers
.agent_suggest <- function(value, candidates, max = 3L) {
  value <- .agent_trim_scalar(value)
  if (!nzchar(value) || !length(candidates))
    return(character())
  candidates <- unique(as.character(candidates))
  value_lower <- tolower(value)
  cand_lower <- tolower(candidates)

  # never suggest the value itself on an exact (case-insensitive) hit; that
  # situation means the caller validated against a different set
  is_exact <- cand_lower == value_lower

  contains <- grepl(value_lower, cand_lower, fixed = TRUE) & !is_exact
  contained_by <- vapply(cand_lower, function(cl)
    grepl(cl, value_lower, fixed = TRUE), logical(1)) & !is_exact
  prefix <- startsWith(cand_lower, value_lower) & !is_exact

  pick <- function(keep, limit) {
    if (!any(keep) || limit <= 0L)
      return(list(idx = integer(), limit = limit))
    idx <- which(keep)
    list(idx = utils::head(idx, limit), limit = limit - length(utils::head(idx, limit)))
  }

  out <- character()
  remaining <- max
  for (keep in list(contains, prefix, contained_by)) {
    if (remaining <= 0L) break
    hit <- pick(keep, remaining)
    out <- c(out, candidates[hit$idx])
    remaining <- hit$limit
  }

  if (remaining > 0L && !any(is_exact) && length(candidates) <= 5000L) {
    # Compare edit distance against the candidate as a whole and against
    # each of its pipe-separated segments: annotation columns follow
    # Category|Subcategory|Variable, and typos land either in the variable
    # segment alone ("logg.fdrr" vs segment "log.fdr") or across the full
    # column name ("ttest|A_vs_B|logg.fdrr" vs "ttest|A_vs_B|log.fdr").
    segment_distance <- function(candidate) {
      segs <- strsplit(candidate, "|", fixed = TRUE)[[1]]
      if (length(segs) <= 1L)
        return(utils::adist(value_lower, tolower(candidate))[1, 1])
      min(utils::adist(value_lower, tolower(c(candidate, segs)))[1, ])
    }
    dist <- vapply(cand_lower, segment_distance, numeric(1))
    names(dist) <- NULL
    threshold <- max(2, floor(nchar(value_lower, type = "chars") / 3))
    near <- which(dist <= threshold & !is_exact)
    if (length(near)) {
      near <- near[order(dist[near])]
      out <- c(out, candidates[utils::head(near, remaining)])
    }
  }

  unique(utils::head(out, max))
}

.agent_suggest_text <- function(value, candidates, max = 3L,
                                 search_hint = NULL) {
  hits <- .agent_suggest(value, candidates, max = max)
  if (!length(hits))
    return("")
  paste0(
    " Closest matches: ", paste(hits, collapse = ", "), ".",
    if (!is.null(search_hint)) paste0(" Use ", search_hint, " to confirm exact values.") else ""
  )
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

#' Read the optional session cost/token budgets
#'
#' WP12 governance: administrators may cap the cumulative session spend in
#' USD and/or the cumulative token usage. Default unlimited-but-logged
#' (settled decision 11). Mirrors deputy's UsageLimits semantics without
#' adopting the package.
#'
#' @return List with optional \code{cost_usd} and \code{tokens} entries
#'   (NULL when unset or invalid).
#' @keywords internal
#' @rdname agentAssistantHelpers
agent_cost_limits <- function() {
  cost <- suppressWarnings(as.numeric(Sys.getenv("OMICSVIEWER_LLM_MAX_COST_USD")[1]))
  tokens <- suppressWarnings(as.numeric(Sys.getenv("OMICSVIEWER_LLM_MAX_TOKENS")[1]))
  list(
    cost_usd = if (length(cost) == 1L && !is.na(cost) && cost > 0) cost else NULL,
    tokens = if (length(tokens) == 1L && !is.na(tokens) && tokens > 0)
      as.integer(tokens) else NULL
  )
}

#' Extract one request's usage from a completed assistant turn
#'
#' @param turn An ellmer AssistantTurn (tokens = input/output/... counts;
#'   cost in USD when the provider reports it).
#' @return \code{list(tokens =, cost_usd =)} with zeros for unavailable
#'   fields.
#' @keywords internal
#' @rdname agentAssistantHelpers
agent_turn_usage <- function(turn) {
  tokens <- 0
  cost <- 0
  if (inherits(turn, "ellmer::AssistantTurn")) {
    toks <- suppressWarnings(as.numeric(turn@tokens))
    tokens <- sum(toks[!is.na(toks)])
    cst <- suppressWarnings(as.numeric(turn@cost)[1])
    cost <- if (!is.na(cst) && cst > 0) cst else 0
  }
  list(tokens = tokens, cost_usd = cost)
}

#' Check cumulative usage against the session budgets
#'
#' Pure helper evaluated before every provider request: the check must not
#' depend on reactive state.
#'
#' @param used_tokens Cumulative session tokens.
#' @param used_cost_usd Cumulative session cost in USD.
#' @param limits From \code{\link{agent_cost_limits}}.
#' @return An error message string when a budget is exhausted, else NULL.
#' @keywords internal
#' @rdname agentAssistantHelpers
agent_budget_violation <- function(used_tokens, used_cost_usd, limits) {
  if (!is.null(limits$tokens) && used_tokens > limits$tokens) {
    return(sprintf(
      paste("Assistant token budget reached for this session (%s of %s tokens).",
            "Start a new browser session or ask an administrator to adjust OMICSVIEWER_LLM_MAX_TOKENS."),
      format(used_tokens, big.mark = ","), format(limits$tokens, big.mark = ",")
    ))
  }
  if (!is.null(limits$cost_usd) && used_cost_usd > limits$cost_usd) {
    return(sprintf(
      paste("Assistant cost budget reached for this session ($%.4f of $%.2f).",
            "Start a new browser session or ask an administrator to adjust OMICSVIEWER_LLM_MAX_COST_USD."),
      used_cost_usd, limits$cost_usd
    ))
  }
  NULL
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

#' Normalize a requested state-section list
#'
#' Providers may serialize omitted array optionals as literal sentinel
#' strings (glm flash) or as JSON lists; both are normalized here so
#' \code{agent_compact_state} never sees an ambiguous request.
#'
#' @param sections Character vector (or list) of requested section names,
#'   or NULL for the overview only.
#' @return A deduplicated character vector of valid section names.
#' @keywords internal
#' @rdname agentAssistantHelpers
.agent_normalize_state_sections <- function(sections) {
  if (is.null(sections))
    return(character())
  if (agent_sentinel_string(sections))
    return(character())
  sections <- as.character(sections)
  sections <- trimws(sections)
  sections <- sections[!is.na(sections) & nzchar(sections)]
  sections <- sections[!duplicated(sections)]
  invalid <- setdiff(sections, AGENT_STATE_SECTIONS)
  if (length(invalid))
    stop("Unknown state section(s): ",
         paste(utils::head(invalid, 3L), collapse = ", "), ".",
         .agent_suggest_text(invalid[[1]], AGENT_STATE_SECTIONS),
         " Available sections: ", paste(AGENT_STATE_SECTIONS, collapse = ", "), ".")
  sections
}

#' Read the current data-space scatter views from the widget store
#'
#' Layer-0 overview anchors (plan section 3, WP1): the x/y axis triples and
#' display mode of both data-space scatters, read from the canonical widget
#' store, which is the single source of truth for axis state.
#'
#' @param store Canonical widget store (\code{\link{widget_store_new}}),
#'   typically the app's root store. NULL returns NULL.
#' @return A named list with \code{feature} and \code{sample} blocks (each
#'   with \code{x}/\code{y} triples and \code{axis_mode}), or NULL when the
#'   store carries no scatter-axis state.
#' @keywords internal
#' @rdname agentAssistantHelpers
agent_scatter_view_from_store <- function(store) {
  if (is.null(store))
    return(NULL)
  read_one <- function(prefix) {
    axes <- paste0(prefix, ".", rep(c("x", "y"), each = 3L), "_",
                   c("analysis", "subset", "variable"))
    vals <- store_read(store, unique(c(axes, paste0(prefix, ".axis_mode"))))
    triple <- function(axis) {
      parts <- vals[paste0(prefix, ".", axis, "_",
                          c("analysis", "subset", "variable"))]
      ok <- vapply(parts, function(p)
        is.character(p) && length(p) == 1L && nzchar(p), logical(1))
      if (!all(ok))
        return(NULL)
      list(analysis = parts[[1]], subset = parts[[2]], variable = parts[[3]],
           name = paste(unlist(parts, use.names = FALSE), collapse = "|"))
    }
    x <- triple("x")
    y <- triple("y")
    if (is.null(x) && is.null(y))
      return(NULL)
    list(x = x, y = y, axis_mode = vals[[paste0(prefix, ".axis_mode")]])
  }
  out <- list(
    feature = read_one("dataspace.feature_space"),
    sample = read_one("dataspace.sample_space")
  )
  if (is.null(out$feature) && is.null(out$sample))
    return(NULL)
  out
}

#' Build a compact model-facing application state
#'
#' Progressive disclosure (plan section 3, WP1): by default the payload is a
#' fixed-size overview - dataset, active and available tabs, selection
#' counts with at most 20 example IDs, quick-view id+label lists, the
#' widget-store scatter view, the \code{available_sections} menu, and the
#' state policy. Full-detail sections (annotation catalog, complete
#' quick-view records, bounded panel state, figure grammar) are returned by
#' the same call only when requested through \code{sections}.
#'
#' @param state Versioned application state from \code{\link{build_app_state}}.
#' @param annotations Annotation catalog created by
#'   \code{\link{agent_annotation_catalog}} (included only when the
#'   \code{annotations} section is requested).
#' @param quick_views Named list containing feature and sample quick views.
#' @param available_tabs Named list of valid data-space and analysis-space tabs.
#' @param figure_grammar Allowlisted figure grammar (included only when the
#'   \code{figure_grammar} section is requested).
#' @param sections Optional character vector requesting full-detail
#'   sections beyond the overview: \code{annotations}, \code{quick_views},
#'   \code{panels}, \code{figure_grammar}.
#' @param store Canonical widget store; when given, the overview gains a
#'   \code{scatter_view} block with the current x/y axis triples and
#'   axis mode of both data-space scatters, plus \code{capabilities}
#'   counts (WP9: the overview lists capability counts only, never record
#'   contents).
#'
#' @return A JSON-like list containing the overview (identifiers, active
#'   tabs, semantic selections, quick-view id/label lists, scatter view,
#'   available sections, state policy) plus any requested full-detail
#'   sections. Expression values, fingerprints, credentials, and timestamps
#'   are intentionally omitted.
#' @keywords internal
#' @rdname agentAssistantHelpers
agent_compact_state <- function(state, annotations = NULL, quick_views = NULL,
                                available_tabs = NULL, figure_grammar = NULL,
                                sections = NULL, store = NULL) {
  if (is.null(state))
    return(NULL)

  sections <- .agent_normalize_state_sections(sections)
  want <- function(section) section %in% sections

  if (is.null(available_tabs))
    available_tabs <- list(data_space = character(), analysis_space = character())

  feature_records <- utils::head(agent_quick_view_records(quick_views$feature), 50L)
  sample_records <- utils::head(agent_quick_view_records(quick_views$sample), 50L)
  bound_selection <- function(ids) {
    list(
      count = length(ids),
      ids = utils::head(ids, 20L),
      truncated = length(ids) > 20L
    )
  }
  overview_views <- function(records) {
    lapply(records, function(r) list(id = r$id, label = r$label))
  }

  out <- list(
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
    quick_views = if (want("quick_views")) {
      list(feature = feature_records, sample = sample_records)
    } else {
      list(
        feature = overview_views(feature_records),
        sample = overview_views(sample_records)
      )
    },
    available_sections = AGENT_STATE_SECTIONS,
    state_policy = state$policy
  )

  scatter_view <- agent_scatter_view_from_store(store)
  if (!is.null(scatter_view))
    out$scatter_view <- scatter_view

  # WP9: capability COUNTS only in the overview - record contents are
  # discovered through search_ui_capabilities (context hygiene).
  capabilities <- agent_capability_summary(store)
  if (!is.null(capabilities))
    out$capabilities <- capabilities

  if (want("annotations"))
    out$annotations <- annotations
  if (want("panels"))
    out$panels <- .agent_bound_value(state$panels)
  if (want("figure_grammar"))
    out$figure_grammar <- figure_grammar
  out
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

  out <- list(
    space = space,
    query = query,
    matching_id_count = length(matching_ids),
    matching_ids = utils::head(matching_ids, max_results),
    matching_column_count = sum(column_hit_flags),
    matching_columns = utils::head(column_hits, 50L),
    value_matches_truncated = value_match_column_count > length(value_hits),
    value_matches = value_hits
  )
  # A completely empty result is a dead end for the model; surface the
  # closest annotation column names so the next call can be a correction
  # instead of another guess.
  if (!out$matching_id_count && !out$matching_column_count && !length(value_hits))
    out$suggestions <- utils::head(
      .agent_suggest(query, column_names), 5L
    )
  out
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
    stop("Unknown ", space, " annotation column: ", column, ".",
         .agent_suggest_text(column, colnames(df)))

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
    # providers may serialize omitted array optionals as literal "null" /
    # "[]" strings (glm flash); treat them as absent
    if (is.character(x) && length(x) == 1L &&
        x %in% c("null", "NULL", "[]", "{}"))
      return(NULL)
    x <- as.character(x)
    if (any(is.na(x)) || any(!nzchar(x)))
      stop(label, " must contain non-empty IDs.")
    x <- x[!duplicated(x)]
    if (length(x) > max_n)
      stop("At most ", max_n, " ", label, "s can be selected by the assistant.")
    invalid <- setdiff(x, ids)
    if (length(invalid)) {
      example <- utils::head(invalid, 3L)
      hints <- vapply(example, function(v) {
        .agent_suggest_text(
          v, ids,
          search_hint = paste0('search_annotations(query="', v, '")')
        )
      }, character(1))
      stop(
        "Unknown ", label, " ID(s): ", paste(example, collapse = ", "), ".",
        paste(hints, collapse = "")
      )
    }
    x
  }

  normalize_tab <- function(x, choices, label) {
    if (is.null(x)) return(NULL)
    x <- .agent_nullable_scalar(x)
    if (!nzchar(x)) return(NULL)
    if (!x %in% choices)
      stop("Unknown ", label, " tab: ", x, ".",
           .agent_suggest_text(x, choices))
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
.agent_nullable_scalar <- function(x) {
  # Some providers serialize omitted optional string arguments as literal
  # sentinel strings instead of JSON null (see AGENT_SENTINEL_STRINGS);
  # normalize those artifacts to an empty string so downstream nzchar()
  # logic treats the argument as absent.
  x <- .agent_trim_scalar(x)
  if (agent_sentinel_string(x)) "" else x
}

agent_normalize_scatter_view <- function(space, quick_view_id = NULL,
                                         x_axis = NULL, y_axis = NULL,
                                         quick_views = NULL,
                                         feature_columns = character(),
                                         sample_columns = character()) {
  space <- match.arg(space, c("feature", "sample"))
  quick_view_id <- .agent_nullable_scalar(quick_view_id)

  if (nzchar(quick_view_id)) {
    views <- quick_views[[space]]
    if (is.null(views) || !is.data.frame(views) || !nrow(views) ||
        !quick_view_id %in% views$id) {
      known_ids <- if (is.null(views) || !is.data.frame(views) || !nrow(views))
        character() else views$id
      stop("Unknown ", space, " quick view: ", quick_view_id, ".",
           .agent_suggest_text(quick_view_id, known_ids),
           " Available quick views are listed under quick_views in get_omics_viewer_state.")
    }
    view <- views[views$id == quick_view_id, , drop = FALSE][1, ]
    return(list(space = space, mode = "quick", quick_view_id = view$id,
                x_axis = view$x, y_axis = view$y))
  }

  x_axis <- .agent_nullable_scalar(x_axis)
  y_axis <- .agent_nullable_scalar(y_axis)
  if (!nzchar(x_axis) || !nzchar(y_axis))
    stop("A scatter view requires either quick_view_id or both x_axis and y_axis.")

  columns <- if (space == "feature") feature_columns else sample_columns
  axis_parts <- function(axis) {
    parts <- strsplit(axis, "|", fixed = TRUE)[[1]]
    if (length(parts) != 3L || any(!nzchar(parts)))
      return(NULL)
    parts
  }
  convention_hint <-
    if (is.null(axis_parts(x_axis)) || is.null(axis_parts(y_axis)))
      " Axis names must use the Category|Subcategory|Variable naming convention."
    else ""
  for (axis_name in c("X", "Y")) {
    axis <- if (axis_name == "X") x_axis else y_axis
    if (!axis %in% columns)
      stop("Unknown ", space, " ", axis_name, "-axis annotation: ", axis, ".",
           .agent_suggest_text(axis, columns), convention_hint)
  }

  list(space = space, mode = "custom", quick_view_id = NULL,
       x_axis = x_axis, y_axis = y_axis)
}

############################################################################
### WP11: conversation-in-snapshot helpers
###
### The .ESS snapshot is the single, opt-in persistence path for the
### conversation (shinychat's own history stores are deliberately NOT
### enabled: file-based storage would persist transcripts outside the
### opt-in guardrail). Turns serialize with display-only payloads (base64
### figure previews) stripped and credential-like strings redacted; the
### figure registry rides along so update_figure revision keeps working
### after restore (specs re-validate against the CURRENT dataset at use -
### restored tool evidence never restores execution authority).
############################################################################

.agent_history_key_pattern <- paste(
  "sk-[A-Za-z0-9_-]{16,}",
  "sk-ant-[A-Za-z0-9_-]{16,}",
  "gsk_[A-Za-z0-9]{16,}",
  "xai-[A-Za-z0-9_-]{16,}",
  "Bearer[[:space:]]+[A-Za-z0-9._-]{16,}",
  sep = "|"
)

#' Redact credential-like strings from history text
#'
#' Belt-and-braces guardrail: API keys never legitimately enter chat turns
#' (configuration happens in a modal), but anything key-shaped is replaced
#' before transcript text is persisted or displayed.
#'
#' @param text Character vector.
#' @return Character vector with key-shaped substrings replaced by
#'   \code{"[redacted]"}.
#' @keywords internal
#' @rdname agentAssistantHelpers
agent_history_redact <- function(text) {
  gsub(.agent_history_key_pattern, "[redacted]", text)
}

#' Slim one turn for snapshot persistence
#'
#' Drops display-only payloads (embedded base64 figure previews inflate
#' every tool result by 100+ KB) from tool results; keeps values and errors
#' as inert context. Returns a NEW turn (the live session object is never
#' mutated).
#'
#' @param turn An ellmer Turn object.
#' @return A new Turn of the same class with slimmed contents.
#' @keywords internal
#' @rdname agentAssistantHelpers
.agent_slim_turn <- function(turn) {
  contents <- lapply(turn@contents, function(x) {
    if (inherits(x, "ellmer::ContentToolResult"))
      return(ellmer::ContentToolResult(value = x@value, error = x@error))
    x
  })
  if (inherits(turn, "ellmer::AssistantTurn")) {
    ellmer::AssistantTurn(
      contents = contents,
      tokens = turn@tokens, cost = turn@cost,
      duration = turn@duration, finish_reason = turn@finish_reason
    )
  } else if (inherits(turn, "ellmer::UserTurn")) {
    ellmer::UserTurn(contents = contents)
  } else {
    turn
  }
}

#' Turn list to bounded display records
#'
#' @param turns List of ellmer Turn objects.
#' @return List of \code{list(role, text)} records (text-only transcript;
#' assistant tool calls summarized by name in brackets).
#' @keywords internal
#' @rdname agentAssistantHelpers
agent_transcript_records <- function(turns) {
  lapply(turns, function(t) {
    texts <- vapply(t@contents, function(x)
      tryCatch(ellmer::contents_text(x) %||% "", error = function(e) ""),
      character(1))
    text <- paste(texts[nzchar(texts)], collapse = "\n")
    calls <- vapply(
      Filter(function(x) inherits(x, "ellmer::ContentToolRequest"), t@contents),
      function(x) x@name, character(1))
    if (length(calls))
      text <- paste(c(text, paste0("[called ", paste(calls, collapse = ", "), "]")),
                    collapse = "\n")
    list(role = if (identical(t@role, "user")) "user" else "assistant",
         text = agent_history_redact(text))
  })
}

#' Build the assistant snapshot payload
#'
#' Byte-capped (default 256 KB, measured by actual serialization size):
#' oldest turns drop first, at least the last two turns always survive.
#' Returns NULL when there is nothing to save.
#'
#' @param turns Client turns (\code{client$get_turns()}).
#' @param figures Figure registry (\code{figures()}); at most 20 entries.
#' @param max_bytes Serialized-size budget in bytes.
#' @return A versioned payload list, or NULL.
#' @keywords internal
#' @rdname agentAssistantHelpers
agent_history_payload <- function(turns, figures = list(),
                                  max_bytes = 256L * 1024L) {
  turns <- Filter(function(t) inherits(t, "ellmer::Turn"), turns)
  if (!length(turns))
    return(NULL)
  slim <- lapply(turns, .agent_slim_turn)
  records <- agent_transcript_records(slim)
  # object.size over-counts S7 objects (class metadata is charged to every
  # instance), so measure the real serialization footprint instead
  text_bytes <- c(0L, cumsum(nchar(vapply(records, function(r) r$text,
                                          character(1)), type = "bytes")))
  turn_bytes <- vapply(slim, function(t)
    length(serialize(t, connection = NULL)), numeric(1))
  cum_bytes <- text_bytes[-1] + cumsum(turn_bytes)
  keep <- length(slim)
  while (keep > 2L && cum_bytes[[keep]] > max_bytes)
    keep <- keep - 1L
  list(
    version = 1L,
    saved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    turns = if (keep == length(slim)) slim else slim[seq_len(keep)],
    transcript = if (keep == length(slim)) records else records[seq_len(keep)],
    figures = if (length(figures)) figures[seq_len(min(length(figures), 20L))] else list(),
    truncated = keep < length(slim),
    bytes = round(cum_bytes[[keep]])
  )
}

#' Validate a restored assistant payload
#'
#' Structure, size, and class checks at the restore boundary - snapshots
#' are untrusted files. Rebuilding is never attempted for malformed data.
#'
#' @param payload Candidate payload from a .ESS snapshot.
#' @param max_bytes Read-side budget (default 512 KB).
#' @return The validated payload (invisibly normalized), or an error.
#' @keywords internal
#' @rdname agentAssistantHelpers
agent_history_restore_payload <- function(payload, max_bytes = 512L * 1024L) {
  if (is.null(payload) || !is.list(payload) || is.null(payload$version))
    stop("Assistant payload is missing or malformed.")
  if (!identical(as.integer(payload$version), 1L))
    stop("Unsupported assistant payload version: ", payload$version, ".")
  turns <- payload$turns
  if (is.null(turns) || !is.list(turns) || !length(turns) ||
      !all(vapply(turns, function(t) inherits(t, "ellmer::Turn"), logical(1))))
    stop("Assistant payload turns are malformed.")
  transcript <- payload$transcript
  if (is.null(transcript) || !is.list(transcript) ||
      !all(vapply(transcript, function(r)
        is.list(r) && nzchar(r$role %||% "") && is.character(r$text),
        logical(1))))
    stop("Assistant payload transcript is malformed.")
  if (length(transcript) != length(turns))
    stop("Assistant payload transcript does not match its turns.")
  figures <- payload$figures %||% list()
  if (!is.list(figures) ||
      !all(vapply(figures, function(f) is.list(f) && nzchar(f$id %||% ""),
                  logical(1))))
    stop("Assistant payload figure registry is malformed.")
  size <- length(serialize(turns, connection = NULL)) +
    sum(nchar(vapply(transcript, function(r) r$text, character(1)),
                   type = "bytes"))
  if (size > max_bytes)
    stop(sprintf("Assistant payload is too large to restore (%.0f KB).", size / 1024))
  payload$figures <- figures[seq_len(min(length(figures), 20L))]
  payload
}
