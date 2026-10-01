#' Package-level tool registry for the AI assistant
#'
#' todo 4.4(b): every assistant tool definition lives here instead of
#' inline in \code{ai_assistant_module}'s moduleServer frame. The factory
#' takes the module's dependencies as EXPLICIT formals with the same names
#' the inline definitions closed over, so the moved code is verbatim and
#' closure scoping is deterministic by construction (ellmer review 7.3:
#' definitions used to resolve lazily in the moduleServer frame - budget
#' state like request_count/logger deliberately stays OUT and is wired by
#' the module's .ai_register_tools Governor hooks).
#'
#' @param state Compact-state builder (sections contract, WP1).
#' @param feature_data,sample_data,expression_data Dataset reactives.
#' @param selected_features,selected_samples Selection reactives.
#' @param apply_state,apply_scatter_view,apply_enrichment,apply_table_view
#'   The app-supplied apply callbacks (L0's one-write-plane functions).
#' @param store Canonical widget store (generic tier; NULL disables it).
#' @param figures Figure-registry reactiveVal (WP3/3.8).
#' @param figure_counter Figure-id counter reactiveVal.
#' @param figure_directory Session temp dir for rendered PNGs.
#' @param session_domain The module session (withReactiveDomain target).
#' @param logger Agent logger (figure eviction events).
#' @return A list with \code{register(client)} - registering every tool in
#'   the documented order, honouring the optional tiers - plus the
#'   individual ToolDefs for tests.
#' @keywords internal
#' @rdname agentToolRegistry
agent_tool_registry <- function(state, feature_data, sample_data,
                                expression_data, selected_features,
                                selected_samples, apply_state,
                                apply_scatter_view,
                                apply_enrichment = NULL,
                                apply_table_view = NULL, store = NULL,
                                figures, figure_counter, figure_directory,
                                session_domain, logger) {
  get_state_tool <- agent_tool(
    function(sections = NULL, `_intent`) {
      current <- state(sections = sections)
      if (is.null(current))
        stop("No dataset is currently available to the assistant.")
      .ai_tool_result(
        current,
        title = "Read current analysis state",
        label = paste(
          current$dataset$dimensions[["features"]], "features /",
          current$dataset$dimensions[["samples"]], "samples",
          if (length(sections))
            paste0("+ ", paste(sections, collapse = ", "))
          else "(overview)"
        ),
        preview = paste("data tab:", current$active_tabs$data_space,
                        "| analysis tab:", current$active_tabs$analysis_space)
      )
    },
    name = "get_omics_viewer_state",
    description = paste(
      "Read the current omicsViewer state. Call with no sections for a compact overview:",
      "dataset, active and available tabs, selection counts with up to 20 example IDs,",
      "quick-view id+label lists, and scatter_view (the current x/y axes and axis mode of",
      "both data-space scatters). Request sections only when needed: 'annotations' (full",
      "annotation column catalog), 'quick_views' (full quick-view records including axes),",
      "'panels' (bounded interface panel state), 'figure_grammar' (declarative figure",
      "grammar and templates). Requested sections are returned together with the overview."
    ),
    arguments = list(
      sections = ellmer::type_array(
        ellmer::type_enum(
          AGENT_STATE_SECTIONS,
          "State section to include in full detail beyond the overview."
        ),
        "Optional sections to include; omit or pass an empty array for the compact overview.",
        required = FALSE
      ),
      `_intent` = ellmer::type_string(
        "Short user-facing reason this state is needed."
      )
    ),
    annotations = ellmer::tool_annotations(
      title = "Reading application state",
      read_only_hint = TRUE,
      destructive_hint = FALSE,
      idempotent_hint = TRUE,
      open_world_hint = FALSE
    )
  )

  search_tool <- agent_tool(
    function(space, query, max_results = 20L, `_intent`) {
      result <- agent_search_annotations(
        space = space,
        query = query,
        feature_data = isolate(feature_data()),
        sample_data = isolate(sample_data()),
        max_results = max_results
      )
      .ai_tool_result(
        result,
        title = "Searched annotations",
        label = paste(result$space, ":", result$query),
        preview = paste(result$matching_id_count, "matching IDs")
      )
    },
    name = "search_annotations",
    description = paste(
      "Search feature/sample IDs, annotation column names, and bounded matching values.",
      "Returns at most 50 IDs. Use this before selecting IDs or discussing annotation names."
    ),
    arguments = list(
      space = ellmer::type_enum(c("feature", "sample"), "Annotation space to search."),
      query = ellmer::type_string("Case-insensitive literal query (1-128 characters)."),
      max_results = ellmer::type_integer("Maximum IDs to return, from 1 through 50.", required = FALSE),
      `_intent` = ellmer::type_string("Short user-facing reason for this search.")
    ),
    annotations = ellmer::tool_annotations(
      title = "Searching annotations",
      read_only_hint = TRUE,
      destructive_hint = FALSE,
      idempotent_hint = TRUE,
      open_world_hint = FALSE
    )
  )

  summary_tool <- agent_tool(
    function(space, column, max_values = 12L, `_intent`) {
      result <- agent_summarize_annotation(
        space = space,
        column = column,
        feature_data = isolate(feature_data()),
        sample_data = isolate(sample_data()),
        max_values = max_values
      )
      .ai_tool_result(
        result,
        title = "Summarized annotation",
        label = paste(result$space, ":", result$column),
        preview = paste(result$type, "|", result$non_missing, "observed")
      )
    },
    name = "summarize_annotation",
    description = paste(
      "Return bounded numeric quantiles or categorical counts for one exact annotation column.",
      "Does not return row-level data or expression values."
    ),
    arguments = list(
      space = ellmer::type_enum(c("feature", "sample"), "Annotation space."),
      column = ellmer::type_string("Exact annotation column name."),
      max_values = ellmer::type_integer("Maximum categorical values, from 1 through 20.", required = FALSE),
      `_intent` = ellmer::type_string("Short user-facing reason for this summary.")
    ),
    annotations = ellmer::tool_annotations(
      title = "Summarizing annotations",
      read_only_hint = TRUE,
      destructive_hint = FALSE,
      idempotent_hint = TRUE,
      open_world_hint = FALSE
    )
  )

  set_state_tool <- agent_tool(
    function(data_space_tab = NULL, analysis_space_tab = NULL,
             features = NULL, samples = NULL, `_intent`) {
      proposal <- list(
        data_space_tab = data_space_tab,
        analysis_space_tab = analysis_space_tab,
        features = features,
        samples = samples
      )
      result <- shiny::withReactiveDomain(
        session_domain,
        apply_state(proposal)
      )
      .ai_tool_result(
        result,
        title = "Updated analysis controls",
        label = paste(names(result), collapse = ", "),
        preview = paste0(
          length(result$features %||% character()), " features; ",
          length(result$samples %||% character()), " samples"
        )
      )
    },
    name = "set_omics_viewer_state",
    description = paste(
      "Change active data/analysis tabs or semantic feature/sample selections after an explicit user request.",
      "Omit fields that should remain unchanged; provide an empty array to clear a selection.",
      "Valid tab and ID values are reported by get_omics_viewer_state and search_annotations."
    ),
    arguments = list(
      data_space_tab = ellmer::type_string("Exact data-space tab label.", required = FALSE),
      analysis_space_tab = ellmer::type_string("Exact analysis-space tab label.", required = FALSE),
      features = ellmer::type_array(
        ellmer::type_string("Exact feature ID."),
        "Feature IDs to select; an empty array clears the feature selection.",
        required = FALSE
      ),
      samples = ellmer::type_array(
        ellmer::type_string("Exact sample ID."),
        "Sample IDs to select; an empty array clears the sample selection.",
        required = FALSE
      ),
      `_intent` = ellmer::type_string("Short user-facing reason for changing visible controls.")
    ),
    annotations = ellmer::tool_annotations(
      title = "Updating analysis controls",
      read_only_hint = FALSE,
      destructive_hint = FALSE,
      idempotent_hint = TRUE,
      open_world_hint = FALSE
    )
  )

  set_scatter_tool <- agent_tool(
    function(space, quick_view_id = NULL, x_axis = NULL, y_axis = NULL, `_intent`) {
      result <- shiny::withReactiveDomain(
        session_domain,
        apply_scatter_view(
          space = space,
          quick_view_id = quick_view_id,
          x_axis = x_axis,
          y_axis = y_axis
        )
      )
      .ai_tool_result(
        result,
        title = "Updated scatter view",
        label = paste(result$space, ":", result$x_axis, "vs", result$y_axis),
        preview = paste("mode:", result$mode)
      )
    },
    name = "set_scatter_view",
    description = paste(
      "Set a feature-space or sample-space scatter view after an explicit user request.",
      "Use quick_view_id from get_omics_viewer_state when possible; otherwise provide exact x_axis and y_axis names.",
      "quick_view_id is a shorthand for an axis pair: it changes the axes only and never switches the display mode.",
      "Does not change selections, tabs, or the quick/custom display mode, and runs no statistical tests."
    ),
    arguments = list(
      space = ellmer::type_enum(c("feature", "sample"), "Scatter space to update."),
      quick_view_id = ellmer::type_string("Exact quick-view ID.", required = FALSE),
      x_axis = ellmer::type_string("Exact X-axis annotation column.", required = FALSE),
      y_axis = ellmer::type_string("Exact Y-axis annotation column.", required = FALSE),
      `_intent` = ellmer::type_string("Short user-facing reason for changing the scatter view.")
    ),
    annotations = ellmer::tool_annotations(
      title = "Updating scatter view",
      read_only_hint = FALSE,
      destructive_hint = FALSE,
      idempotent_hint = TRUE,
      open_world_hint = FALSE
    )
  )

  # ----------------------------------------------------------------
  # S3 generic widget tier: registry-driven, thin wrappers over the
  # canonical widget store. isolate() is required: validation reads
  # choices providers, which read module reactives (triset() etc.).
  # ----------------------------------------------------------------
  list_widgets_tool <- agent_tool(
    function(section = NULL, `_intent`) {
      result <- shiny::isolate(shiny::withReactiveDomain(
        session_domain, agent_widget_list(store, section)))
      ids <- vapply(result$widgets, function(w) w$id, character(1))
      if (!length(ids)) ids <- "(none)"
      .ai_tool_result(
        result,
        title = "Listed controllable widgets",
        label = paste(result$widget_count, "widgets",
                     if (is.null(result$section)) "" else paste("in", result$section)),
        preview = paste(utils::head(ids, 5), collapse = ", ")
      )
    },
    name = "list_widgets",
    description = paste(
      "List the user-editable interface widgets you can control, with canonical ids, kinds,",
      "allowed values, and current values.",
      "Optional section filters by id prefix (e.g. 'dataspace' or 'dataspace.expr_heatmap').",
      "Use this before set_widgets and to answer questions about available interface controls."
    ),
    arguments = list(
      section = ellmer::type_string("Optional canonical id prefix to filter by.", required = FALSE),
      `_intent` = ellmer::type_string("Short user-facing reason for listing widgets.")
    ),
    annotations = ellmer::tool_annotations(
      title = "Listing interface widgets",
      read_only_hint = TRUE,
      destructive_hint = FALSE,
      idempotent_hint = TRUE,
      open_world_hint = FALSE
    )
  )

  get_widget_tool <- agent_tool(
    function(id, `_intent`) {
      result <- shiny::isolate(shiny::withReactiveDomain(
        session_domain, agent_widget_describe(store, id)))
      .ai_tool_result(
        result,
        title = "Described interface widget",
        label = result$id,
        preview = paste(result$kind, "|",
                        if (is.null(result$current_value)) "(unset)"
                        else as.character(result$current_value))
      )
    },
    name = "get_widget",
    description = paste(
      "Describe one user-editable widget by exact canonical id: kind, meaning, allowed",
      "values, dependencies, and its current value. Ids come from list_widgets."
    ),
    arguments = list(
      id = ellmer::type_string("Exact canonical widget id from list_widgets."),
      `_intent` = ellmer::type_string("Short user-facing reason for describing this widget.")
    ),
    annotations = ellmer::tool_annotations(
      title = "Describing interface widget",
      read_only_hint = TRUE,
      destructive_hint = FALSE,
      idempotent_hint = TRUE,
      open_world_hint = FALSE
    )
  )

  set_widgets_tool <- agent_tool(
    function(patch, `_intent`) {
      result <- shiny::isolate(shiny::withReactiveDomain(
        session_domain, agent_widget_apply(store, patch)))
      .ai_tool_result(
        result,
        title = "Updated interface widgets",
        label = paste(length(result$applied), "applied,",
                      length(result$rejected), "rejected"),
        preview = paste(
          if (length(result$applied))
            paste(result$applied, collapse = ", ") else "nothing applied",
          "; rejected:", if (length(result$rejected))
            paste(vapply(result$rejected, function(r) r$id, character(1)),
                  collapse = ", ") else "none"
        )
      )
    },
    name = "set_widgets",
    description = paste(
      "Set one or more user-editable interface widgets after an explicit user request,",
      "for controls not covered by set_scatter_view or set_omics_viewer_state.",
      "patch is a JSON object string mapping canonical widget ids to single values",
      "whose types match the widget kind (string/number/boolean; e.g.",
      "'{\"dataspace.expr_heatmap.heatmap_colors\": \"RdGy\"}').",
      "Only requested widgets change; invalid keys are reported per key with closest-match",
      "suggestions so you can correct and retry. Discover ids and allowed values with list_widgets."
    ),
    arguments = list(
      patch = ellmer::type_string(paste(
        "JSON object string of canonical widget id to value.",
        "Example: {\"dataspace.expr_heatmap.heatmap_colors\": \"RdGy\"}")),
      `_intent` = ellmer::type_string("Short user-facing reason for changing visible controls.")
    ),
    annotations = ellmer::tool_annotations(
      title = "Updating interface widgets",
      read_only_hint = FALSE,
      destructive_hint = FALSE,
      idempotent_hint = TRUE,
      open_world_hint = FALSE
    )
  )

  # ----------------------------------------------------------------
  # WP8 semantic tools: thin validate + store_apply wrappers over the
  # canonical widget store. Descriptions mirror the shared capability
  # metadata (auxi_agentCapabilities.R) - one source of truth.
  # ----------------------------------------------------------------
  set_enrichment_tool <- agent_tool(
    function(method, collapse = NULL, selected_pathway = NULL, `_intent`) {
      result <- shiny::withReactiveDomain(
        session_domain,
        apply_enrichment(list(
          method = method,
          collapse = collapse,
          selected_pathway = selected_pathway
        ))
      )
      .ai_tool_result(
        result,
        title = "Updated enrichment panel",
        label = paste0(result$method, " -> ", result$panel_tab),
        preview = paste(
          length(result$applied), "keys applied;",
          length(result$rejected), "rejected"
        )
      )
    },
    name = "set_enrichment_parameters",
    description = paste(
      "Configure the enrichment analysis panel after an explicit user request.",
      "method 'ora' runs over-representation analysis of the currently selected features;",
      "'fgsea' ranks all features and runs GSEA.",
      "collapse is the exact 'Category|Subcategory|Variable' feature annotation column",
      "features are collapsed on (ORA) or ranked by (fGSEA; must be numeric).",
      "selected_pathway optionally selects one gene-set row of the results table and must",
      "be a row of the CURRENT results (re-select in a follow-up call after changing the",
      "ranking). The panel's tab opens automatically; results recompute reactively.",
      "Valid columns come from search_annotations; valid pathway ids from the results tables"
    ),
    arguments = list(
      method = ellmer::type_enum(
        c("ora", "fgsea"),
        "Enrichment method: over-representation of the selection, or ranked GSEA."
      ),
      collapse = ellmer::type_string(
        "Exact Category|Subcategory|Variable feature annotation column.",
        required = FALSE
      ),
      selected_pathway = ellmer::type_string(
        "Exact gene-set id whose results row should be selected.",
        required = FALSE
      ),
      `_intent` = ellmer::type_string(
        "Short user-facing reason for changing the enrichment panel."
      )
    ),
    annotations = ellmer::tool_annotations(
      title = "Updating enrichment panel",
      read_only_hint = FALSE,
      destructive_hint = FALSE,
      idempotent_hint = TRUE,
      open_world_hint = FALSE
    )
  )

  set_table_view_tool <- agent_tool(
    function(table, columns = NULL, multi_selection = NULL,
             column_filters = NULL, page = NULL, clear_filters = NULL,
             `_intent`) {
      result <- shiny::withReactiveDomain(
        session_domain,
        apply_table_view(list(
          table = table,
          columns = columns,
          multi_selection = multi_selection,
          column_filters = column_filters,
          page = page,
          clear_filters = clear_filters
        ))
      )
      .ai_tool_result(
        result,
        title = "Updated table view",
        label = paste0(result$table, " -> ", result$panel_tab),
        preview = paste(
          length(result$applied), "keys applied;",
          length(result$rejected), "rejected"
        )
      )
    },
    name = "set_table_view",
    description = paste(
      "Configure one data table after an explicit user request: which columns are",
      "shown, whether multiple rows can be selected, per-column filter patterns, and",
      "the visible page. table is 'feature_table', 'sample_table', or 'expression_table';",
      "the table's tab opens automatically.",
      "column_filters is a JSON object of exact column name to substring pattern",
      "(e.g. '{\"group\": \"KO\"}'); pass clear_filters=true to clear all filters.",
      "Valid column names come from search_annotations or the columns widget",
      "(dataspace.<table>.columns) via get_widget."
    ),
    arguments = list(
      table = ellmer::type_enum(
        c("feature_table", "sample_table", "expression_table"),
        "Which data-space table to configure."
      ),
      columns = ellmer::type_array(
        ellmer::type_string("Exact column name."),
        "Exact set of columns to display; at least one must remain.",
        required = FALSE
      ),
      multi_selection = ellmer::type_boolean(
        "Allow selecting more than one row.",
        required = FALSE
      ),
      column_filters = ellmer::type_string(
        paste("JSON object of exact column name to substring pattern,",
              "e.g. {\"group\": \"KO\"}."),
        required = FALSE
      ),
      page = ellmer::type_integer("1-based page number to display.", required = FALSE),
      clear_filters = ellmer::type_boolean(
        "Clear every column filter (overrides column_filters).",
        required = FALSE
      ),
      `_intent` = ellmer::type_string(
        "Short user-facing reason for changing the table view."
      )
    ),
    annotations = ellmer::tool_annotations(
      title = "Updating table view",
      read_only_hint = FALSE,
      destructive_hint = FALSE,
      idempotent_hint = TRUE,
      open_world_hint = FALSE
    )
  )

  # ----------------------------------------------------------------
  # WP9 discovery tools: generated from the capability registry
  # (widget bindings + tool metadata) - one source of truth.
  # ----------------------------------------------------------------
  search_capabilities_tool <- agent_tool(
    function(query, max_results = 20L, `_intent`) {
      result <- shiny::isolate(shiny::withReactiveDomain(
        session_domain, agent_capability_search(store, query, max_results)))
      ids <- vapply(result$capabilities, function(r) r$id, character(1))
      if (!length(ids)) ids <- "(no matches)"
      .ai_tool_result(
        result,
        title = "Searched UI capabilities",
        label = paste0(result$query, " : ", result$match_count, " of ",
                       result$capability_count),
        preview = paste(utils::head(ids, 5), collapse = ", ")
      )
    },
    name = "search_ui_capabilities",
    description = paste(
      "Search everything you can read or control in this app by meaning:",
      "semantic tools (scatter, enrichment, tables, figures) and every user-editable",
      "widget with its panel, purpose, allowed values, and the tool that operates it.",
      "Case-insensitive substring match over ids, labels, help text, and panels.",
      "Use this before list_widgets when looking by purpose rather than exact id prefix;",
      "returns at most 50 records with match counts."
    ),
    arguments = list(
      query = ellmer::type_string("Case-insensitive search query (1-128 characters)."),
      max_results = ellmer::type_integer(
        "Maximum records to return, from 1 through 50.", required = FALSE),
      `_intent` = ellmer::type_string("Short user-facing reason for this search.")
    ),
    annotations = ellmer::tool_annotations(
      title = "Searching UI capabilities",
      read_only_hint = TRUE,
      destructive_hint = FALSE,
      idempotent_hint = TRUE,
      open_world_hint = FALSE
    )
  )

  get_capability_tool <- agent_tool(
    function(id, `_intent`) {
      result <- shiny::isolate(shiny::withReactiveDomain(
        session_domain, agent_capability_get(store, id)))
      .ai_tool_result(
        result,
        title = "Described UI capability",
        label = result$id,
        preview = paste(result$kind, "|", result$operation)
      )
    },
    name = "get_ui_capability",
    description = paste(
      "Describe one capability by exact id: a widget id (from list_widgets or",
      "search_ui_capabilities) or a semantic tool name. Returns the panel, purpose,",
      "allowed values, dependencies, and the tool that operates it."
    ),
    arguments = list(
      id = ellmer::type_string(
        "Exact capability id (widget id or semantic tool name)."),
      `_intent` = ellmer::type_string(
        "Short user-facing reason for describing this capability.")
    ),
    annotations = ellmer::tool_annotations(
      title = "Describing UI capability",
      read_only_hint = TRUE,
      destructive_hint = FALSE,
      idempotent_hint = TRUE,
      open_world_hint = FALSE
    )
  )

  # todo 4.2: the figure spec schema is generated from the single
  # grammar table (R/auxi_agentFigureGrammar.R) shared with the
  # validator and the prose grammar; filter nesting stops at depth 2
  # in the schema while the validator keeps accepting depth-3 trees
  # for old transcripts/snapshot replays.

  render_assistant_figure <- function(spec, parent_figure_id = NULL,
                                   template = NULL) {
    warning_messages <- character()
    rendered <- withCallingHandlers(
      {
        normalized <- agent_normalize_figure_spec(
          spec = spec,
          feature_data = isolate(feature_data()),
          sample_data = isolate(sample_data()),
          expression = isolate(expression_data()),
          selected_features = isolate(.agent_figure_selection(selected_features())),
          selected_samples = isolate(.agent_figure_selection(selected_samples()))
        )
        figure_data <- agent_build_figure_data(
          spec = normalized,
          feature_data = isolate(feature_data()),
          sample_data = isolate(sample_data()),
          expression = isolate(expression_data())
        )
        figure_plot <- agent_build_figure_plot(figure_data, normalized)
        current_count <- isolate(figure_counter()) + 1L
        figure_counter(current_count)
        figure_id <- paste0("fig_", current_count)
        files <- agent_render_figure(
          figure_plot,
          directory = figure_directory,
          figure_id = figure_id
        )
        list(
          normalized = normalized,
          data = figure_data,
          plot = figure_plot,
          files = files,
          figure_id = figure_id
        )
      },
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )

    spec <- rendered$normalized
    figure_id <- rendered$figure_id
    # todo 3.8: bounded registry - evict (LRU, superseded revisions
    # first) instead of failing once 20 figures exist; the evicted
    # PNGs are transient (agent_render_figure unlinks its files), so
    # only the registry entry and its revision base are dropped.
    registry <- isolate(figures())
    if (!is.null(parent_figure_id) && !is.null(registry[[parent_figure_id]])) {
      touched <- registry[[parent_figure_id]]
      touched$last_used_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
      registry[[parent_figure_id]] <- touched
    }
    room <- agent_figure_registry_evict(
      registry, exclude = parent_figure_id %||% character())
    if (length(room$evicted)) {
      agent_logger_event(
        logger, "figure_evicted",
        list(
          evicted = room$evicted,
          reason = "registry_full",
          registry_size_before = length(registry)
        )
      )
      warning_messages <- c(
        warning_messages,
        paste0(
          "Figure registry full: evicted ",
          paste(room$evicted, collapse = ", "),
          " (least recently used; superseded revisions first). ",
          "Their specs are no longer revisable; all other figures are unchanged."
        )
      )
    }
    registry <- room$registry
    registry[[figure_id]] <- list(
      id = figure_id,
      parent_id = parent_figure_id,
      template = template,
      created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      last_used_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      row_count = nrow(rendered$data),
      geoms = vapply(spec$layers, function(x) x$geom, character(1)),
      # WP3/WP7: the normalized spec is the canonical revision base;
      # patch-mode update_figure merges onto exactly this copy.
      spec = spec
    )
    figures(registry)

    metadata <- list(
      figure_id = figure_id,
      parent_figure_id = parent_figure_id,
      template = template,
      data_source = spec$data_source,
      row_count = nrow(rendered$data),
      # resolved counts (what was actually plotted), independent of
      # whether the id array is echoed (todo 3.1 elides defaults)
      feature_count = if (!is.null(spec$features)) length(spec$features) else
        if (identical(spec$data_source, "expression")) 0L else
          nrow(isolate(feature_data())),
      sample_count = if (!is.null(spec$samples)) length(spec$samples) else
        if (identical(spec$data_source, "expression")) 0L else
          nrow(isolate(sample_data())),
      layers = lapply(spec$layers, function(x)
        list(geom = x$geom, mappings = x$mappings, filter = x$filter)),
      facet_by = spec$facet_by,
      x_transform = spec$x_transform,
      y_transform = spec$y_transform,
      theme = spec[["theme"]],
      palette = spec$palette,
      labels = spec$labels,
      preview_dimensions = rendered$files$preview_dimensions,
      full_dimensions = rendered$files$full_dimensions,
      full_png_bytes = rendered$files$full_bytes,
      download_filename = paste0("omicsviewer-", figure_id, "-2400x1800.png"),
      # WP3 round-trip (todo 3.1: bounded): the compact echo carries
      # the documented (flat-aesthetic) input shape but never a
      # dataset-sized id array - small explicit subsets ride along,
      # defaults/large sets are elided in favour of the counts above.
      # Revise through update_figure(figure_id, changes) or read the
      # compact spec back with get_figure(figure_id).
      spec = agent_figure_spec_echo(
        spec,
        all_features = rownames(isolate(feature_data())),
        all_samples = rownames(isolate(sample_data()))
      ),
      warnings = utils::head(unique(warning_messages), 10L)
    )

    list(
      metadata = metadata,
      display = shinychat::tool_result_display(
        title = paste("Created AI figure", figure_id),
        label = paste(spec$data_source, "|", nrow(rendered$data), "rows"),
        value_preview = paste(
          length(spec$layers), "layer(s);",
          rendered$files$full_bytes, "byte full PNG"
        ),
        html = agent_figure_html(rendered$files, spec, figure_id),
        show_request = FALSE,
        open = TRUE,
        full_screen = TRUE,
        open_style = "framed"
      )
    )
  }

  # WP6: template arguments expand server-side into a full validated
  # spec (agent_figure_template_spec) and then flow through the exact
  # same render/round-trip path as hand-written specs, so template
  # figures stay revisable through update_figure.
  expand_figure_template <- function(template, x, y, color, label_top_n,
                                    title, space, features, samples) {
    agent_figure_template_spec(
      template = template,
      x = x,
      y = y,
      color = color,
      label_top_n = label_top_n,
      title = title,
      space = space,
      features = features,
      samples = samples,
      feature_data = isolate(feature_data()),
      sample_data = isolate(sample_data()),
      expression = isolate(expression_data()),
      selected_features = isolate(.agent_figure_selection(selected_features())),
      selected_samples = isolate(.agent_figure_selection(selected_samples()))
    )
  }

  create_figure_tool <- agent_tool(
    function(spec = NULL, template = NULL, x = NULL, y = NULL, color = NULL,
             label_top_n = NULL, title = NULL, space = NULL, features = NULL,
             samples = NULL, `_intent`) {
      template_name <- NULL
      spec_absent <- is.null(spec) || length(spec) == 0L
      template_present <- !.agent_param_absent(template)
      template_ignored <- FALSE
      if (template_present) {
        if (!spec_absent) {
          # WP6b (2026-09-24 benchmark evidence): models echo the WP3
          # round-trip spec AND fill template shorthand in one call.
          # A hard "not both" error was never recovered from (16
          # identical retries in one task). The spec is authoritative —
          # it is the full previous figure — so honor it and surface a
          # warning instead of rejecting the call.
          template_ignored <- TRUE
        } else {
          expanded <- shiny::withReactiveDomain(
            session_domain,
            expand_figure_template(template, x, y, color, label_top_n,
                                  title, space, features, samples)
          )
          spec <- expanded$spec
          template_name <- expanded$template
        }
      } else if (spec_absent) {
        stop("create_figure requires either a template (with x and optional y, color, label_top_n, features, samples, title, space) or a full spec.")
      }
      result <- shiny::withReactiveDomain(
        session_domain,
        render_assistant_figure(spec, template = template_name)
      )
      if (template_ignored)
        result$metadata$warnings <- c(
          result$metadata$warnings,
          "Template arguments ignored: a full spec was also provided and takes precedence. To use template shorthand, resend WITHOUT the spec."
        )
      ellmer::ContentToolResult(
        value = result$metadata,
        extra = list(display = result$display)
      )
    },
    name = "create_figure",
    description = paste(
      "Create a static ggplot2 figure.",
      "Preferred for common plots: pass template ('volcano', 'scatter', 'boxplot', 'histogram') with a few exact column names (x; y; color; label_top_n; features; samples; title; space) and the server expands it to a validated figure.",
      "Templates: volcano = fold change vs log-scale significance (x, y required; label_top_n labels the most significant features); boxplot = numeric column by group, or without y the expression of the selected features (or an explicit features subset) by a sample column; scatter and histogram = two or one annotation columns.",
      "Never send a template AND a spec together: when both are present the spec is used and the template arguments are ignored.",
      "The chat displays a small preview and a high-resolution PNG download.",
      "The result carries figure_id and a compact spec (large id sets are elided - see feature_count/sample_count). Revise the figure with update_figure(figure_id, changes={partial spec}); read the current compact spec with get_figure(figure_id).",
      "For advanced multi-layer figures pass the full declarative spec instead; for expression data use feature__ and sample__ prefixed metadata columns described by the figure grammar.",
      "Highlighting and labeling: any layer accepts a structured filter (column/operator/value; all/any combinators) and constant params.color/params.fill hex colors, so e.g. label only rows where a column exceeds a threshold in a custom color; text/label layers take params.order_by to rank rows before the max_labels cap; scale sets explicit per-category colors and theme_options tunes legend/rotation/grid.",
      "Use exact columns returned by get_omics_viewer_state/search_annotations and never invent R code."
    ),
    arguments = list(
      template = ellmer::type_enum(
        .agent_figure_templates,
        "Template name; expands server-side to a full figure specification.",
        required = FALSE
      ),
      x = ellmer::type_string(
        "Exact x column (fold change for volcano; grouping column for boxplot).",
        required = FALSE
      ),
      y = ellmer::type_string(
        "Exact y column (log-scale significance for volcano; omit for expression-mode boxplot).",
        required = FALSE
      ),
      color = ellmer::type_string(
        "Exact column mapped to point color or box fill.",
        required = FALSE
      ),
      label_top_n = ellmer::type_integer(
        "Volcano/scatter only: label the top n features (0-50).",
        required = FALSE
      ),
      title = ellmer::type_string("Figure title (at most 200 characters).", required = FALSE),
      space = ellmer::type_string(
        "'feature' or 'sample'; disambiguates column names present in both spaces.",
        required = FALSE
      ),
      features = ellmer::type_array(
        ellmer::type_string("Exact feature ID."),
        paste("Optional: restrict the plotted features (e.g. the first N selected genes);",
              "defaults to the current selection for expression figures and all features otherwise."),
        required = FALSE
      ),
      samples = ellmer::type_array(
        ellmer::type_string("Exact sample ID."),
        "Optional: restrict the plotted samples of expression figures.",
        required = FALSE
      ),
      spec = agent_figure_spec_schema(required = FALSE),
      `_intent` = ellmer::type_string("Short user-facing description of the requested figure.")
    ),
    annotations = ellmer::tool_annotations(
      title = "Creating analysis figure",
      read_only_hint = TRUE,
      destructive_hint = FALSE,
      idempotent_hint = FALSE,
      open_world_hint = FALSE
    )
  )

  get_figure_tool <- agent_tool(
    function(figure_id, `_intent`) {
      current <- isolate(figures())
      fid <- .agent_figure_scalar(figure_id, arg = "figure_id")
      entry <- current[[fid]]
      if (is.null(entry)) {
        stop("Unknown figure ID: ", fid, ".", .agent_suggest_text(
          fid, names(current)),
          " Current figures: ",
          if (length(current)) paste(names(current), collapse = ", ")
          else "(none)", ".")
      }
      touched <- entry
      touched$last_used_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
      current[[fid]] <- touched
      figures(current)
      .ai_tool_result(
        list(
          figure_id = entry$id,
          parent_figure_id = entry$parent_id,
          template = entry$template,
          created_at = entry$created_at,
          row_count = entry$row_count,
          geoms = entry$geoms,
          spec = agent_figure_spec_echo(
            entry$spec,
            all_features = rownames(isolate(feature_data())),
            all_samples = rownames(isolate(sample_data()))
          ),
          revision_note = paste(
            "Revise with update_figure(figure_id, changes): a partial",
            "spec of only the fields to change; unmentioned fields keep",
            "their current values. Large id sets are elided from the",
            "echo - pass an explicit features/samples array only when",
            "changing it."
          )
        ),
        title = "Read figure",
        label = fid,
        preview = paste(entry$geoms, collapse = "+")
      )
    },
    name = "get_figure",
    description = paste(
      "Read one figure's current compact specification by exact figure_id:",
      "data source, layers with mappings/filters, labels, theme, and id",
      "counts. Large feature/sample id sets are elided from the echo",
      "(see feature_count/sample_count in create_figure results).",
      "Use it before update_figure when you need the current spec, and",
      "revise through changes rather than reconstructing the full spec."
    ),
    arguments = list(
      figure_id = ellmer::type_string("Exact figure ID returned by create_figure."),
      `_intent` = ellmer::type_string("Short user-facing reason for reading this figure.")
    ),
    annotations = ellmer::tool_annotations(
      title = "Reading figure",
      read_only_hint = TRUE,
      destructive_hint = FALSE,
      idempotent_hint = TRUE,
      open_world_hint = FALSE
    )
  )

  update_figure_tool <- agent_tool(
    function(figure_id, spec = NULL, changes = NULL, `_intent`) {
      current <- isolate(figures())
      fid <- .agent_figure_scalar(figure_id, arg = "figure_id")
      entry <- current[[fid]]
      if (is.null(entry)) {
        stop("Unknown figure ID: ", fid, ".", .agent_suggest_text(
          fid, names(current)),
          " Current figures: ",
          if (length(current)) paste(names(current), collapse = ", ")
          else "(none)", ".")
      }
      spec_absent <- .agent_param_absent(spec)
      changes_absent <- .agent_param_absent(changes)
      if (spec_absent && changes_absent)
        stop(paste(
          "update_figure requires either changes (a partial spec of only",
          "the fields to change; preferred) or spec (a complete",
          "specification)."
        ))
      full_spec_ignored <- FALSE
      if (!changes_absent && !spec_absent) {
        # WP6b posture: never hard-reject an over-eager double fill;
        # the complete spec is authoritative and the warning says so.
        full_spec_ignored <- TRUE
        merged <- spec
      } else if (!spec_absent) {
        merged <- spec
      } else {
        merged <- agent_figure_spec_patch(entry$spec, changes)
      }
      result <- shiny::withReactiveDomain(
        session_domain,
        render_assistant_figure(merged, parent_figure_id = fid)
      )
      if (full_spec_ignored)
        result$metadata$warnings <- c(
          result$metadata$warnings,
          "changes ignored: a complete spec was also provided and takes precedence. To patch, resend WITHOUT the spec."
        )
      ellmer::ContentToolResult(
        value = result$metadata,
        extra = list(display = result$display)
      )
    },
    name = "update_figure",
    description = paste(
      "Revise an existing figure into a new revision.",
      "Send changes - a partial spec of only the fields to change",
      "(e.g. {\"labels\":{\"title\":\"New\"},\"theme\":\"classic\"}) - merged onto the",
      "stored spec; unmentioned fields keep their current values.",
      "Scalars and arrays (layers, features, samples) replace wholesale;",
      "labels/theme_options/scale merge per key (explicit null removes a key);",
      "null/empty features or samples clears back to the default set.",
      "Call get_figure(figure_id) first when you need the current spec.",
      "The new result keeps figure_id as its parent. Use create_figure for a new, unrelated figure.",
      "For expression data, use feature__ and sample__ prefixed metadata columns described by the figure grammar."
    ),
    arguments = list(
      figure_id = ellmer::type_string("Exact figure ID returned by create_figure."),
      changes = agent_figure_spec_schema(required = FALSE, layers_required = FALSE),
      # todo 4.2: the FULL spec schema lives only on create_figure's
      # advanced path; this legacy argument stays dispatch-compatible
      # (formals must match the schema exactly) as an opaque object so
      # pre-narrowing transcripts and echoed complete specs still
      # revise — the documented path is `changes`.
      spec = ellmer::type_object(
        paste("Complete specification (legacy; prefer changes).",
              "Same fields as create_figure's spec argument."),
        .required = FALSE
      ),
      `_intent` = ellmer::type_string("Short user-facing description of the requested change.")
    ),
    annotations = ellmer::tool_annotations(
      title = "Updating analysis figure",
      read_only_hint = TRUE,
      destructive_hint = FALSE,
      idempotent_hint = FALSE,
      open_world_hint = FALSE
    )
  )

  list(
    register = function(client) {
      client$register_tool(get_state_tool)
      client$register_tool(search_tool)
      client$register_tool(summary_tool)
      client$register_tool(set_state_tool)
      client$register_tool(set_scatter_tool)
      if (!is.null(store)) {
        client$register_tool(list_widgets_tool)
        client$register_tool(get_widget_tool)
        client$register_tool(set_widgets_tool)
        if (!is.null(apply_enrichment))
          client$register_tool(set_enrichment_tool)
        if (!is.null(apply_table_view))
          client$register_tool(set_table_view_tool)
        client$register_tool(search_capabilities_tool)
        client$register_tool(get_capability_tool)
      }
      client$register_tool(create_figure_tool)
      client$register_tool(get_figure_tool)
      client$register_tool(update_figure_tool)
      client
    },
    tools = list(
      get_state = get_state_tool, search = search_tool,
      summary = summary_tool, set_state = set_state_tool,
      set_scatter = set_scatter_tool, list_widgets = list_widgets_tool,
      get_widget = get_widget_tool, set_widgets = set_widgets_tool,
      set_enrichment = set_enrichment_tool,
      set_table_view = set_table_view_tool,
      search_capabilities = search_capabilities_tool,
      get_capability = get_capability_tool,
      create_figure = create_figure_tool, get_figure = get_figure_tool,
      update_figure = update_figure_tool
    )
  )
}
