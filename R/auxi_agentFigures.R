#' Internal constrained AI figure helpers
#'
#' These helpers implement a declarative, allowlisted ggplot2 subset for the
#' optional AI assistant. The model composes JSON-like figure specifications;
#' omicsViewer validates every field and performs all rendering itself. The
#' assistant never supplies R code, function names, file paths, raw HTML, or
#' arbitrary expression strings.
#'
#' @return Figure grammar/catalog objects, normalized specifications, bounded
#'   plotting data, ggplot objects, rendered files, or server-generated HTML.
#'
#' @importFrom base64enc dataURI
#' @keywords internal
#' @name agentFigureHelpers
NULL

.agent_figure_sources <- c("feature_annotation", "sample_annotation", "expression")
.agent_figure_reserved <- c("__feature_id__", "__sample_id__", "__expression__")
.agent_figure_geoms <- c(
  "point", "line", "path", "bar", "boxplot", "violin", "histogram", "density",
  "text", "label", "smooth", "errorbar", "ribbon", "hline", "vline"
)
.agent_figure_aesthetics <- c(
  "x", "y", "color", "fill", "group", "size", "alpha", "shape", "linetype",
  "label", "ymin", "ymax"
)
.agent_figure_transforms <- c("identity", "log2", "log10", "reverse", "sqrt")
.agent_figure_themes <- c("minimal", "classic", "light", "grey", "bw")
.agent_figure_palettes <- c("default", "colorblind", "sequential", "diverging", "grey")
# WP6 figure templates (plan decision 4: volcano/boxplot/histogram/scatter
# ship now; density/barplot stay log-gated; pca is deferred)
.agent_figure_templates <- c("volcano", "scatter", "boxplot", "histogram")

#' Describe the allowlisted model-facing figure grammar
#'
#' @return A JSON-like description of supported data sources, geoms, aesthetic
#'   mappings, transformations, themes, palettes, templates (WP6), and hard
#'   limits.
#' @keywords internal
#' @rdname agentFigureHelpers
agent_figure_grammar <- function() {
  list(
    data_sources = list(
      feature_annotation = "One row per feature; feature metadata columns plus __feature_id__.",
      sample_annotation = "One row per sample; sample metadata columns plus __sample_id__.",
      expression = paste(
        "Long-format expression values for selected features/samples;",
        "feature metadata is prefixed with feature__, sample metadata with sample__,",
        "and reserved values are __feature_id__, __sample_id__, and __expression__."
      )
    ),
    reserved_columns = .agent_figure_reserved,
    expression_metadata_prefixes = c(feature = "feature__", sample = "sample__"),
    geoms = .agent_figure_geoms,
    aesthetics = .agent_figure_aesthetics,
    transforms = .agent_figure_transforms,
    themes = .agent_figure_themes,
    palettes = .agent_figure_palettes,
    templates = agent_figure_templates(),
    limits = list(
      max_layers = 12L,
      max_annotation_rows = 20000L,
      max_expression_features = 50L,
      max_expression_samples = 200L,
      max_text_labels = 50L,
      max_full_png_mb = 10
    )
  )
}

.agent_figure_scalar <- function(x, fallback = NULL, max_chars = 200L) {
  if (is.null(x) || length(x) == 0 || is.na(x))
    return(fallback)
  out <- trimws(as.character(x)[1])
  # providers may serialize omitted optional strings as literal sentinels
  # (see AGENT_SENTINEL_STRINGS); treat them like absent values
  if (is.na(out) || agent_sentinel_string(out)) return(fallback)
  if (nzchar(out) && nchar(out, type = "chars", allowNA = TRUE) > max_chars)
    out <- paste0(substr(out, 1L, max_chars), " ...")
  out
}

.agent_figure_numeric_param <- function(value, name, min, max, default) {
  # ellmer converts JSON null to NA and some providers (glm flash) echo
  # omitted optionals as empty objects or literal sentinels; treat NA,
  # length-0, empty-list, and sentinel-string values like an omitted value.
  if (.agent_param_absent(value)) return(default)
  value <- suppressWarnings(as.numeric(value)[1])
  if (is.na(value) || value < min || value > max)
    stop("Figure parameter ", name, " must be between ", min, " and ", max, ".")
  value
}

.agent_figure_integer_param <- function(value, name, min, max, default) {
  # ellmer converts JSON null to NA and some providers (glm flash) echo
  # omitted optionals as empty objects or literal sentinels; treat NA,
  # length-0, empty-list, and sentinel-string values like an omitted value.
  if (.agent_param_absent(value)) return(default)
  value <- suppressWarnings(as.integer(value)[1])
  if (is.na(value) || value < min || value > max)
    stop("Figure parameter ", name, " must be an integer between ", min, " and ", max, ".")
  value
}

.agent_param_absent <- function(value) {
  is.null(value) || length(value) == 0L || is.na(value[1]) ||
    (is.character(value) && agent_sentinel_string(value[1]))
}

.agent_figure_ids_param <- function(value, name) {
  # WP6b: optional ID-array template argument (features/samples). Normalizes
  # the provider-artifact shapes seen in tool traffic — ellmer tibbles,
  # record lists, NA entries, literal sentinel strings ("null", "[]") — to
  # a trimmed character vector, or NULL when effectively omitted.
  if (.agent_param_absent(value)) return(NULL)
  if (inherits(value, "data.frame"))
    value <- unlist(lapply(as.list(value), as.character), use.names = FALSE)
  else if (is.list(value))
    value <- unlist(lapply(value, function(v) as.character(v)[1]), use.names = FALSE)
  else
    value <- as.character(value)
  value <- trimws(value[!is.na(value)])
  value <- value[nzchar(value) & !(value %in% AGENT_SENTINEL_STRINGS)]
  if (!length(value)) return(NULL)
  if (anyDuplicated(value))
    stop("Figure ", name, " must be unique IDs; duplicate entries found.")
  value
}

.agent_figure_choice <- function(value, choices, name, fallback = NULL) {
  if (is.null(value) || length(value) == 0 || is.na(value) ||
      agent_sentinel_string(value))
    return(fallback)
  value <- .agent_figure_scalar(value, max_chars = 100L)
  if (is.null(value) || !value %in% choices)
    stop("Unsupported figure ", name, ": ", value, ".",
         .agent_suggest_text(value, choices))
  value
}

.agent_figure_selection <- function(x) {
  if (is.null(x) || length(x) == 0) return(character())
  if (agent_sentinel_string(x)) return(character())
  x <- as.character(x)
  unique(x[!is.na(x) & nzchar(x)])
}

.agent_figure_ids <- function(x, valid_ids, label, max_n, allow_default = TRUE) {
  if (is.null(x) || agent_sentinel_string(x)) {
    if (!allow_default)
      stop("Figure requires at least one ", label, " ID.")
    return(NULL)
  }
  x <- as.character(x)
  if (any(is.na(x)) || any(!nzchar(x)))
    stop("Figure ", label, " IDs must be non-empty strings.")
  x <- x[!duplicated(x)]
  if (length(x) > max_n)
    stop("Figure can plot at most ", max_n, " ", label, "s.")
  invalid <- setdiff(x, valid_ids)
  if (length(invalid)) {
    example <- utils::head(invalid, 3L)
    stop("Unknown figure ", label, " ID(s): ", paste(example, collapse = ", "))
  }
  x
}

#' Convert a normalized figure spec to its re-submittable echo shape
#'
#' Tool results must carry the spec in the exact shape the tool schema
#' documents (aesthetics as direct layer fields): providers and ellmer's
#' schema-driven argument conversion drop properties that are not in the
#' declared schema, so a model echoing the result verbatim would otherwise
#' lose the axis mappings. The session registry keeps the normalized shape
#' (the canonical base for a future patch-mode update_figure).
#'
#' @param spec Normalized specification from
#'   \code{\link{agent_normalize_figure_spec}}.
#' @return A JSON-like spec in the documented input shape; normalizing it
#'   again reproduces \code{spec} exactly.
#' @keywords internal
#' @rdname agentFigureHelpers
agent_figure_spec_echo <- function(spec) {
  layers <- lapply(spec$layers, function(layer) {
    out <- c(list(geom = layer$geom), layer$mappings)
    if (!is.null(layer$params)) out$params <- layer$params
    out
  })
  out <- spec
  out$layers <- layers
  out
}

#' Normalize and validate a model-proposed figure specification
#'
#' @param spec Model-proposed declarative figure specification.
#' @param feature_data Feature metadata.
#' @param sample_data Sample metadata.
#' @param expression Expression matrix.
#' @param selected_features Current semantic feature selection.
#' @param selected_samples Current semantic sample selection.
#'
#' @return A normalized specification with only allowlisted keys and values.
#' @keywords internal
#' @rdname agentFigureHelpers
agent_normalize_figure_spec <- function(spec, feature_data, sample_data, expression,
                                       selected_features = character(),
                                       selected_samples = character()) {
  if (!is.list(spec))
    stop("Figure specification must be an object.")

  reserved_metadata <- intersect(
    c(colnames(feature_data), colnames(sample_data)),
    .agent_figure_reserved
  )
  if (length(reserved_metadata))
    stop("Metadata column names are reserved for AI figures: ",
         paste(utils::head(reserved_metadata, 3L), collapse = ", "))

  allowed <- c(
    "data_source", "features", "samples", "layers", "facet_by", "facet_ncol",
    "x_transform", "y_transform", "theme", "palette", "labels"
  )
  unknown <- setdiff(names(spec), allowed)
  if (length(unknown))
    stop("Unknown figure specification field(s): ", paste(unknown, collapse = ", "))

  data_source <- .agent_figure_choice(
    spec$data_source, .agent_figure_sources, "data source", "feature_annotation"
  )

  feature_ids <- rownames(feature_data)
  sample_ids <- rownames(sample_data)
  if (is.null(feature_ids)) feature_ids <- character()
  if (is.null(sample_ids)) sample_ids <- character()

  features <- .agent_figure_ids(
    spec$features, feature_ids, "feature",
    if (data_source == "expression") 50L else 20000L,
    allow_default = TRUE
  )
  # WP3 round-trip: expression figures normalize to the FULL sample set by
  # default, so a spec echoed back from a tool result may explicitly carry
  # every sample. The ad-hoc cap exists to bound hand-written subsets; a
  # verbatim full set must survive re-submission on larger datasets.
  samples_are_full_set <- !is.null(spec$samples) &&
    !agent_sentinel_string(spec$samples) &&
    identical(as.character(spec$samples), sample_ids)
  samples <- .agent_figure_ids(
    spec$samples, sample_ids, "sample",
    if (data_source == "expression" && !samples_are_full_set) 200L else 20000L,
    allow_default = TRUE
  )
  if (data_source == "expression") {
    if (is.null(features) || !length(features)) {
      if (is.null(selected_features) || !length(selected_features))
        stop("Expression figures require selected feature IDs or an explicit feature array.")
      features <- .agent_figure_ids(
        selected_features, feature_ids, "feature", 50L, allow_default = TRUE
      )
    }
    if (is.null(samples) || !length(samples))
      samples <- sample_ids
  }

  # ellmer converts type_array(of type_object) arguments into tibbles, so a
  # figure spec arriving through the chat tool path has layers as an
  # N-row/15-column data.frame (where length() counts columns!) and params
  # as df-columns. Coerce both to plain row-lists BEFORE counting, and keep
  # plain-list specs (unit tests, programmatic callers) working unchanged.
  if (is.data.frame(spec$layers)) {
    spec$layers <- lapply(seq_len(nrow(spec$layers)), function(i)
      as.list(spec$layers[i, , drop = FALSE]))
  }
  spec$layers <- lapply(spec$layers, function(layer) {
    if (is.data.frame(layer$params))
      layer$params <- as.list(layer$params[1, , drop = FALSE])
    layer
  })

  if (!is.list(spec$layers) || length(spec$layers) < 1L)
    stop("Figure requires between 1 and 12 layers.")
  if (length(spec$layers) > 12L)
    stop("Figure supports at most 12 layers. Merge same-geom layers using the color/fill aesthetic instead of one layer per group, or split the request into two figures.")

  layers <- lapply(spec$layers, function(layer) {
    if (!is.list(layer))
      stop("Every figure layer must be an object.")
    # WP3 round-trip: normalized specs (returned in tool results and stored
    # in the session figure registry) carry aesthetics under a `mappings`
    # named list. Expand it into the documented flat shape so a normalized
    # spec re-submitted through update_figure validates identically;
    # explicitly given flat aesthetics win over expanded mappings.
    if (!is.null(layer$mappings)) {
      if (!is.list(layer$mappings))
        stop("Figure layer mappings must be an object.")
      unknown_mappings <- setdiff(names(layer$mappings), .agent_figure_aesthetics)
      if (length(unknown_mappings))
        stop("Unknown figure layer mapping(s): ",
             paste(unknown_mappings, collapse = ", "))
      for (nm in names(layer$mappings)) {
        if (is.null(layer[[nm]]))
          layer[[nm]] <- layer$mappings[[nm]]
      }
      layer$mappings <- NULL
    }
    allowed_layer <- c("geom", .agent_figure_aesthetics, "params")
    unknown_layer <- setdiff(names(layer), allowed_layer)
    if (length(unknown_layer))
      stop("Unknown figure layer field(s): ", paste(unknown_layer, collapse = ", "))
    geom <- .agent_figure_choice(layer$geom, .agent_figure_geoms, "geom")
    if (is.null(geom))
      stop("Every figure layer requires an allowlisted geom.")

    mappings <- lapply(.agent_figure_aesthetics, function(aesthetic) {
      .agent_figure_scalar(layer[[aesthetic]], max_chars = 500L)
    })
    names(mappings) <- .agent_figure_aesthetics
    mappings <- mappings[!vapply(mappings, is.null, logical(1))]

    required <- switch(
      geom,
      point = , line = , path = , smooth = c("x", "y"),
      boxplot = , violin = c("x", "y"),
      bar = c("x"),
      histogram = , density = c("x"),
      text = , label = c("x", "y", "label"),
      errorbar = , ribbon = c("x", "ymin", "ymax"),
      hline = character(),
      vline = character()
    )
    missing <- setdiff(required, names(mappings))
    if (length(missing))
      stop("Figure layer ", geom, " requires mapping(s): ", paste(missing, collapse = ", "))

    if (!is.list(layer$params))
      params <- list()
    else
      params <- layer$params
    allowed_params <- c(
      "alpha", "size", "linewidth", "bins", "method", "se", "position",
      "xintercept", "yintercept", "max_labels"
    )
    unknown_params <- setdiff(names(params), allowed_params)
    if (length(unknown_params))
      stop("Unknown figure layer parameter(s): ", paste(unknown_params, collapse = ", "))

    params$alpha <- .agent_figure_numeric_param(params$alpha, "alpha", 0, 1, 0.85)
    params$size <- .agent_figure_numeric_param(params$size, "size", 0.05, 12, 1.8)
    params$linewidth <- .agent_figure_numeric_param(params$linewidth, "linewidth", 0.05, 6, 0.8)
    params$bins <- .agent_figure_integer_param(params$bins, "bins", 5L, 100L, 30L)
    params$method <- .agent_figure_choice(params$method, c("auto", "lm", "loess"), "method", "auto")
    if (is.null(params$se) || length(params$se) == 0L || is.na(params$se[1])) params$se <- TRUE
    if (!isTRUE(params$se %in% c(TRUE, FALSE)))
      stop("Figure layer parameter se must be true or false.")
    params$position <- .agent_figure_choice(
      params$position, c("stack", "dodge", "fill", "jitter"), "position", "stack"
    )
    params$xintercept <- .agent_figure_numeric_param(params$xintercept, "xintercept", -1e9, 1e9, NULL)
    params$yintercept <- .agent_figure_numeric_param(params$yintercept, "yintercept", -1e9, 1e9, NULL)
    params$max_labels <- .agent_figure_integer_param(params$max_labels, "max_labels", 0L, 50L, 20L)
    if (geom %in% c("hline", "vline")) {
      if (geom == "hline" && is.null(params$yintercept))
        stop("Figure layer hline requires numeric yintercept.")
      if (geom == "vline" && is.null(params$xintercept))
        stop("Figure layer vline requires numeric xintercept.")
    }

    list(geom = geom, mappings = mappings, params = params)
  })

  facet_by <- .agent_figure_scalar(spec$facet_by, max_chars = 500L)
  facet_ncol <- .agent_figure_integer_param(spec$facet_ncol, "facet_ncol", 1L, 6L, NULL)
  x_transform <- .agent_figure_choice(
    spec$x_transform, .agent_figure_transforms, "x transform", "identity"
  )
  y_transform <- .agent_figure_choice(
    spec$y_transform, .agent_figure_transforms, "y transform", "identity"
  )
  theme <- .agent_figure_choice(spec$theme, .agent_figure_themes, "theme", "minimal")
  palette <- .agent_figure_choice(spec$palette, .agent_figure_palettes, "palette", "default")

  labels <- if (is.list(spec$labels)) spec$labels else list()
  unknown_labels <- setdiff(names(labels), c("title", "subtitle", "x", "y", "caption"))
  if (length(unknown_labels))
    stop("Unknown figure label field(s): ", paste(unknown_labels, collapse = ", "))
  labels <- lapply(labels, function(x) .agent_figure_scalar(x, max_chars = 200L))
  labels <- labels[!vapply(labels, is.null, logical(1))]

  list(
    data_source = data_source,
    features = features,
    samples = samples,
    layers = layers,
    facet_by = facet_by,
    facet_ncol = facet_ncol,
    x_transform = x_transform,
    y_transform = y_transform,
    theme = theme,
    palette = palette,
    labels = labels
  )
}

#' Describe the allowlisted figure templates
#'
#' Templates (WP6) are the concise path through \code{create_figure}: a few
#' well-named columns expand server-side into a full validated specification.
#' The generic grammar remains the advanced path for multi-layer figures.
#'
#' @return A JSON-like description of every supported template with its named
#'   arguments and their semantics.
#' @keywords internal
#' @rdname agentFigureHelpers
agent_figure_templates <- function() {
  list(
    volcano = list(
      description = paste(
        "Volcano plot of feature-level differential-analysis results:",
        "fold change on x, log-scale significance on y (higher = more",
        "significant, e.g. a log.pvalue or log.fdr column), a zero",
        "fold-change reference line, and optional labels on the most",
        "significant features."
      ),
      arguments = list(
        x = "Required numeric feature column: fold change (e.g. a ttest mean.diff column).",
        y = "Required numeric feature column: log-scale significance where higher = more significant (e.g. a log.fdr or log.pvalue column).",
        color = "Optional feature column mapped to point color.",
        label_top_n = "Optional integer 0-50: label the n features ranked by y, highest first. Default 0 (no labels).",
        title = "Optional title; defaults to 'Volcano: <x> vs <y>'."
      )
    ),
    scatter = list(
      description = paste(
        "Scatter plot of one annotation column against another in either",
        "the feature or the sample space, with optional point coloring",
        "and optional labels."
      ),
      arguments = list(
        x = "Required column for the x axis.",
        y = "Required column for the y axis.",
        color = "Optional column mapped to point color (same space as x/y).",
        label_top_n = "Optional integer 0-50: label the first n rows in data order. Default 0 (no labels).",
        title = "Optional title; defaults to '<x> vs <y>'.",
        space = "Optional 'feature' or 'sample'; required to disambiguate when a column name exists in both spaces."
      )
    ),
    boxplot = list(
      description = paste(
        "Boxplot of a numeric column grouped by a categorical column.",
        "Without y it plots the expression distribution of the currently",
        "selected features grouped by a sample annotation column; pass an",
        "explicit features array to plot a subset (e.g. the first N",
        "selected genes)."
      ),
      arguments = list(
        x = "Required grouping column: a categorical annotation column, or (when y is omitted) a sample annotation column.",
        y = "Optional numeric column; omit it to plot the expression of the selected features (__expression__) by a sample grouping.",
        color = "Optional column mapped to fill (defaults to the grouping column x).",
        features = "Optional array of exact feature IDs; restricts expression mode to these features (at most 50; the current selection is the default).",
        samples = "Optional array of exact sample IDs; restricts expression mode to these samples (at most 200).",
        title = "Optional title; defaults to '<y> by <x>' (or 'Expression by <x>')."
      )
    ),
    histogram = list(
      description = "Histogram of one numeric annotation column.",
      arguments = list(
        x = "Required numeric column.",
        color = "Optional column mapped to fill for grouped histograms.",
        title = "Optional title; defaults to 'Distribution of <x>'."
      )
    ),
    shared_arguments = list(
      features = "Accepted by every template: optional array of exact feature IDs restricting the plotted rows (validated against the dataset; expression figures allow at most 50).",
      samples = "Accepted by every template: optional array of exact sample IDs restricting expression figures (at most 200)."
    )
  )
}

.agent_template_space <- function(space) {
  value <- .agent_figure_scalar(space, max_chars = 100L)
  if (is.null(value)) return(NULL)
  if (!value %in% c("feature", "sample"))
    stop("Figure template space must be 'feature' or 'sample', not: ", value, ".")
  value
}

# Resolve a template column argument against both annotation spaces.
# Returns list(column, space) or NULL for an absent optional column; errors
# carry closest-match suggestions (WP2 control surface) and name the space
# conflict when a column exists in both spaces and none was pinned.
.agent_template_column <- function(value, role, feature_data, sample_data,
                                   space = NULL, required = TRUE) {
  column <- .agent_figure_scalar(value, max_chars = 500L)
  if (is.null(column)) {
    if (required)
      stop("Figure template requires a ", role, " column.")
    return(NULL)
  }
  feature_columns <- colnames(feature_data)
  sample_columns <- colnames(sample_data)
  in_feature <- column %in% feature_columns
  in_sample <- column %in% sample_columns
  if (!in_feature && !in_sample)
    stop("Unknown figure template ", role, " column: ", column, ".",
         .agent_suggest_text(column, c(feature_columns, sample_columns),
                             search_hint = "search_annotations"))
  space <- .agent_template_space(space)
  if (!is.null(space)) {
    columns <- if (identical(space, "feature")) feature_columns else sample_columns
    if (!column %in% columns)
      stop("Unknown ", space, "-space figure template ", role, " column: ", column, ".",
           .agent_suggest_text(column, columns, search_hint = "search_annotations"))
    return(list(column = column, space = space))
  }
  if (in_feature && in_sample)
    stop("Figure template ", role, " column exists in both the feature and the sample annotations: ",
         column, ". Pass space='feature' or space='sample' to disambiguate.")
  list(column = column, space = if (in_feature) "feature" else "sample")
}

.agent_template_volcano <- function(x, y, color, label_top_n, title,
                                    feature_data, sample_data,
                                    features = NULL) {
  x_col <- .agent_template_column(
    x, "x (fold change)", feature_data, sample_data, space = "feature")
  y_col <- .agent_template_column(
    y, "y (log-scale significance)", feature_data, sample_data, space = "feature")
  if (!is.numeric(feature_data[[x_col$column]]))
    stop("Volcano x column must be numeric (fold change): ", x_col$column, ".")
  if (!is.numeric(feature_data[[y_col$column]]))
    stop("Volcano y column must be numeric (log-scale significance): ", y_col$column, ".")
  point <- list(geom = "point", x = x_col$column, y = y_col$column)
  color_col <- .agent_template_column(
    color, "color", feature_data, sample_data, space = "feature", required = FALSE)
  if (!is.null(color_col)) point$color <- color_col$column
  spec <- list(
    data_source = "feature_annotation",
    layers = list(point, list(geom = "vline", params = list(xintercept = 0))),
    labels = list(
      title = title %||% paste("Volcano:", x_col$column, "vs", y_col$column),
      x = x_col$column,
      y = y_col$column
    )
  )
  if (label_top_n > 0L) {
    # reorder the plotted rows so the capped label layer marks exactly the
    # most significant features: higher y = more significant (the app's
    # log.pvalue/log.fdr convention; see module_meta_scatter volcano detection).
    # An explicit features subset is ranked within itself.
    significance <- feature_data[[y_col$column]]
    ids <- if (!is.null(features)) features else rownames(feature_data)
    spec$features <- ids[order(
      significance[match(ids, rownames(feature_data))],
      decreasing = TRUE, na.last = TRUE
    )]
    spec$layers <- c(spec$layers, list(list(
      geom = "label",
      x = x_col$column,
      y = y_col$column,
      label = "__feature_id__",
      params = list(max_labels = label_top_n, size = 3)
    )))
  }
  spec
}

.agent_template_scatter <- function(x, y, color, label_top_n, title,
                                    feature_data, sample_data, space) {
  x_col <- .agent_template_column(x, "x", feature_data, sample_data, space)
  y_col <- .agent_template_column(y, "y", feature_data, sample_data, space)
  if (!identical(x_col$space, y_col$space))
    stop("Scatter x and y columns must come from the same annotation space: x is ",
         x_col$space, "-space, y is ", y_col$space, "-space.")
  color_col <- .agent_template_column(
    color, "color", feature_data, sample_data, x_col$space, required = FALSE)
  id_column <- if (identical(x_col$space, "feature")) "__feature_id__" else "__sample_id__"
  point <- list(geom = "point", x = x_col$column, y = y_col$column)
  if (!is.null(color_col)) point$color <- color_col$column
  spec <- list(
    data_source = paste0(x_col$space, "_annotation"),
    layers = list(point),
    labels = list(
      title = title %||% paste(x_col$column, "vs", y_col$column),
      x = x_col$column,
      y = y_col$column
    )
  )
  if (label_top_n > 0L) {
    spec$layers <- c(spec$layers, list(list(
      geom = "label",
      x = x_col$column,
      y = y_col$column,
      label = id_column,
      params = list(max_labels = label_top_n, size = 3)
    )))
  }
  spec
}

.agent_template_boxplot <- function(x, y, color, title,
                                    feature_data, sample_data, space,
                                    selected_features, features = NULL) {
  x_col <- .agent_template_column(x, "x (grouping)", feature_data, sample_data, space)
  y_value <- .agent_figure_scalar(y, max_chars = 500L)
  if (is.null(y_value)) {
    # expression mode: distribution of the selected features' expression
    # grouped by a sample annotation column
    if (!identical(x_col$space, "sample"))
      stop("Boxplot without a y column plots the expression of selected features grouped by a sample annotation; x column '",
           x_col$column, "' is feature-space. Pass a numeric y column for a feature-space boxplot.")
    if (!length(selected_features) && is.null(features))
      stop("Boxplot expression mode requires plotted features; select features in the app first, or pass explicit feature IDs through the features argument.")
    color_col <- .agent_template_column(
      color, "color", feature_data, sample_data, "sample", required = FALSE)
    grouped <- paste0("sample__", x_col$column)
    fill <- if (!is.null(color_col)) paste0("sample__", color_col$column) else grouped
    return(list(
      data_source = "expression",
      layers = list(list(
        geom = "boxplot", x = grouped, y = "__expression__", fill = fill,
        params = list(alpha = 0.65)
      )),
      labels = list(
        title = title %||% paste("Expression by", x_col$column),
        x = x_col$column,
        y = "Expression"
      )
    ))
  }
  y_col <- .agent_template_column(y, "y", feature_data, sample_data, space)
  if (!identical(x_col$space, y_col$space))
    stop("Boxplot x and y columns must come from the same annotation space: x is ",
         x_col$space, "-space, y is ", y_col$space, "-space.")
  frame <- if (identical(x_col$space, "feature")) feature_data else sample_data
  if (!is.numeric(frame[[y_col$column]]))
    stop("Boxplot y column must be numeric: ", y_col$column, ".")
  color_col <- .agent_template_column(
    color, "color", feature_data, sample_data, x_col$space, required = FALSE)
  list(
    data_source = paste0(x_col$space, "_annotation"),
    layers = list(list(
      geom = "boxplot",
      x = x_col$column,
      y = y_col$column,
      fill = if (!is.null(color_col)) color_col$column else x_col$column,
      params = list(alpha = 0.65)
    )),
    labels = list(
      title = title %||% paste(y_col$column, "by", x_col$column),
      x = x_col$column,
      y = y_col$column
    )
  )
}

.agent_template_histogram <- function(x, color, title,
                                      feature_data, sample_data, space) {
  x_col <- .agent_template_column(x, "x", feature_data, sample_data, space)
  frame <- if (identical(x_col$space, "feature")) feature_data else sample_data
  if (!is.numeric(frame[[x_col$column]]))
    stop("Histogram column must be numeric: ", x_col$column, ".")
  color_col <- .agent_template_column(
    color, "color", feature_data, sample_data, x_col$space, required = FALSE)
  layer <- list(geom = "histogram", x = x_col$column)
  if (!is.null(color_col)) layer$fill <- color_col$column
  list(
    data_source = paste0(x_col$space, "_annotation"),
    layers = list(layer),
    labels = list(
      title = title %||% paste("Distribution of", x_col$column),
      x = x_col$column
    )
  )
}

#' Expand a figure template into a full validated specification
#'
#' WP6: templates are the concise \code{create_figure} path. A few well-named
#' columns (\code{x}, \code{y}, \code{color}, \code{label_top_n}, \code{title},
#' \code{space}) expand server-side into a complete specification that is
#' normalized and validated by \code{\link{agent_normalize_figure_spec}}, so
#' template figures carry the full WP3 round-trip spec and remain revisable
#' through \code{update_figure}. The generic grammar stays the advanced path.
#'
#' @param template Template name: one of volcano, scatter, boxplot, histogram.
#' @param x Required x column (grouping column for boxplot).
#' @param y Optional y column (fold-change significance for volcano; omitted y
#'   selects expression mode for boxplot).
#' @param color Optional column mapped to point color or fill.
#' @param label_top_n Optional integer 0-50: label top features (volcano,
#'   ranked by y) or the first rows (scatter); unsupported elsewhere.
#' @param title Optional plot title.
#' @param space Optional "feature" or "sample" disambiguation when a column
#'   name exists in both annotation spaces.
#' @param features Optional character vector of exact feature IDs restricting
#'   the plotted rows (e.g. a subset of the current selection); the current
#'   selection remains the default.
#' @param samples Optional character vector of exact sample IDs restricting
#'   the plotted rows of expression figures.
#' @param feature_data Feature metadata.
#' @param sample_data Sample metadata.
#' @param expression Expression matrix (passed through to normalization).
#' @param selected_features Current semantic feature selection (expression-mode
#'   boxplot source).
#' @param selected_samples Current semantic sample selection.
#'
#' @return A list with \code{template} (the validated template name) and
#'   \code{spec} (the normalized specification).
#' @keywords internal
#' @rdname agentFigureHelpers
agent_figure_template_spec <- function(template = NULL, x = NULL, y = NULL,
                                       color = NULL, label_top_n = NULL,
                                       title = NULL, space = NULL,
                                       features = NULL, samples = NULL,
                                       feature_data, sample_data,
                                       expression = NULL,
                                       selected_features = character(),
                                       selected_samples = character()) {
  template_value <- .agent_figure_scalar(template, max_chars = 100L)
  if (is.null(template_value))
    stop("Figure template is required. Available templates: ",
         paste(.agent_figure_templates, collapse = ", "), ".")
  if (!template_value %in% .agent_figure_templates)
    stop("Unsupported figure template: ", template_value,
         ". Available templates: ", paste(.agent_figure_templates, collapse = ", "), ".",
         .agent_suggest_text(template_value, .agent_figure_templates))
  label_top_n <- .agent_figure_integer_param(label_top_n, "label_top_n", 0L, 50L, 0L)
  if (label_top_n > 0L && template_value %in% c("boxplot", "histogram"))
    stop("label_top_n is supported by the volcano and scatter templates only.")
  title_value <- .agent_figure_scalar(title)
  # WP6b: normalize the optional ID subsets up front so the volcano
  # template can rank within the subset.
  features_value <- .agent_figure_ids_param(features, "features")
  samples_value <- .agent_figure_ids_param(samples, "samples")

  spec <- switch(
    template_value,
    volcano = .agent_template_volcano(
      x, y, color, label_top_n, title_value, feature_data, sample_data,
      features_value),
    scatter = .agent_template_scatter(
      x, y, color, label_top_n, title_value, feature_data, sample_data,
      .agent_template_space(space)),
    boxplot = .agent_template_boxplot(
      x, y, color, title_value, feature_data, sample_data,
      .agent_template_space(space), selected_features, features_value),
    histogram = .agent_template_histogram(
      x, color, title_value, feature_data, sample_data,
      .agent_template_space(space))
  )

  # WP6b: attach the optional explicit ID subsets (e.g. "the first 20
  # selected genes" — previously inexpressible on the template path, which
  # forced models into the template+spec collision). The volcano template
  # already set spec$features when ranking for label_top_n; IDs are
  # validated downstream by agent_normalize_figure_spec against the live
  # rownames, and each subset only lands on data sources that plot it.
  if (!is.null(features_value) &&
      spec$data_source %in% c("feature_annotation", "expression") &&
      is.null(spec$features))
    spec$features <- features_value
  if (!is.null(samples_value) &&
      spec$data_source %in% c("sample_annotation", "expression"))
    spec$samples <- samples_value

  list(
    template = template_value,
    spec = agent_normalize_figure_spec(
      spec,
      feature_data = feature_data,
      sample_data = sample_data,
      expression = expression,
      selected_features = selected_features,
      selected_samples = selected_samples
    )
  )
}

#' Build bounded plotting data for a validated figure specification
#'
#' @param spec Normalized specification from
#'   \code{\link{agent_normalize_figure_spec}}.
#'
#' @return A data.frame with reserved ID/expression columns and original
#'   annotation column names.
#' @keywords internal
#' @rdname agentFigureHelpers
agent_build_figure_data <- function(spec, feature_data, sample_data, expression) {
  feature_ids <- rownames(feature_data)
  sample_ids <- rownames(sample_data)

  if (identical(spec$data_source, "feature_annotation")) {
    ids <- if (length(spec$features)) spec$features else feature_ids
    if (length(ids) > 20000L)
      stop("Feature annotation figures can plot at most 20000 rows.")
    out <- data.frame(`__feature_id__` = ids, check.names = FALSE, stringsAsFactors = FALSE)
    if (length(ids)) {
      idx <- match(ids, feature_ids)
      out <- cbind(out, feature_data[idx, , drop = FALSE])
    }
    if (anyDuplicated(colnames(out)))
      stop("Figure data contains duplicate annotation column names.")
    rownames(out) <- NULL
    return(out)
  }

  if (identical(spec$data_source, "sample_annotation")) {
    ids <- if (length(spec$samples)) spec$samples else sample_ids
    if (length(ids) > 20000L)
      stop("Sample annotation figures can plot at most 20000 rows.")
    out <- data.frame(`__sample_id__` = ids, check.names = FALSE, stringsAsFactors = FALSE)
    if (length(ids)) {
      idx <- match(ids, sample_ids)
      out <- cbind(out, sample_data[idx, , drop = FALSE])
    }
    if (anyDuplicated(colnames(out)))
      stop("Figure data contains duplicate annotation column names.")
    rownames(out) <- NULL
    return(out)
  }

  if (is.null(spec$features) || !length(spec$features))
    stop("Expression figures require selected feature IDs.")
  if (is.null(spec$samples) || !length(spec$samples))
    stop("Expression figures require at least one sample.")

  mat <- expression[spec$features, spec$samples, drop = FALSE]
  long <- data.frame(
    `__feature_id__` = rep(spec$features, times = ncol(mat)),
    `__sample_id__` = rep(spec$samples, each = nrow(mat)),
    `__expression__` = as.numeric(mat),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  fi <- match(long[["__feature_id__"]], feature_ids)
  si <- match(long[["__sample_id__"]], sample_ids)
  feature_metadata <- feature_data[fi, , drop = FALSE]
  sample_metadata <- sample_data[si, , drop = FALSE]
  colnames(feature_metadata) <- paste0("feature__", colnames(feature_metadata))
  colnames(sample_metadata) <- paste0("sample__", colnames(sample_metadata))
  out <- cbind(long, feature_metadata, sample_metadata)
  if (anyDuplicated(colnames(out)))
    stop("Expression figure metadata contain duplicate names after namespacing.")
  out
}

#' Build a ggplot from a validated figure specification
#'
#' @param data Plotting data from \code{\link{agent_build_figure_data}}.
#' @param spec Normalized figure specification.
#'
#' @return A ggplot object composed exclusively from allowlisted geoms,
#'   mappings, scales, themes, and palettes.
#' @keywords internal
#' @rdname agentFigureHelpers
agent_build_figure_plot <- function(data, spec) {
  columns <- colnames(data)
  all_mappings <- unique(unlist(lapply(spec$layers, function(x) as.character(x$mappings))))
  all_mappings <- unname(all_mappings[!is.na(all_mappings) & nzchar(all_mappings)])
  invalid <- setdiff(all_mappings, columns)
  if (length(invalid)) {
    hints <- vapply(utils::head(invalid, 3L), function(v)
      .agent_suggest_text(v, columns), character(1))
    stop("Figure mappings refer to unavailable column(s): ",
         paste(utils::head(invalid, 3L), collapse = ", "), ".",
         paste(hints, collapse = ""))
  }
  if (!is.null(spec$facet_by) && !spec$facet_by %in% columns)
    stop("Figure facet column is unavailable: ", spec$facet_by, ".",
         .agent_suggest_text(spec$facet_by, columns))

  mapping_symbol <- function(name) as.name(name)
  make_mapping <- function(mappings) {
    args <- lapply(mappings, mapping_symbol)
    do.call(ggplot2::aes, args)
  }

  plot <- ggplot2::ggplot(data)
  for (layer in spec$layers) {
    geom <- layer$geom
    mappings <- layer$mappings
    params <- layer$params
    mapping <- make_mapping(mappings)

    if (geom %in% c("text", "label") && params$max_labels > 0L) {
      layer_data <- utils::head(data, params$max_labels)
    } else {
      layer_data <- NULL
    }

    new_layer <- switch(
      geom,
      point = if (identical(params$position, "jitter")) {
        ggplot2::geom_jitter(mapping = mapping, alpha = params$alpha, size = params$size)
      } else {
        ggplot2::geom_point(mapping = mapping, alpha = params$alpha, size = params$size)
      },
      line = ggplot2::geom_line(mapping = mapping, alpha = params$alpha, linewidth = params$linewidth),
      path = ggplot2::geom_path(mapping = mapping, alpha = params$alpha, linewidth = params$linewidth),
      bar = {
        args <- list(mapping = mapping, alpha = params$alpha)
        if (!is.null(mappings$y)) args$stat <- "identity"
        if (params$position %in% c("stack", "dodge", "fill"))
          args$position <- params$position
        do.call(ggplot2::geom_bar, args)
      },
      boxplot = ggplot2::geom_boxplot(mapping = mapping, alpha = params$alpha),
      violin = ggplot2::geom_violin(mapping = mapping, alpha = params$alpha),
      histogram = ggplot2::geom_histogram(
        mapping = mapping, alpha = params$alpha, bins = params$bins,
        position = if (params$position %in% c("stack", "dodge", "fill")) params$position else "stack"
      ),
      density = ggplot2::geom_density(mapping = mapping, alpha = params$alpha, linewidth = params$linewidth),
      text = do.call(
        ggplot2::geom_text,
        list(data = layer_data, mapping = mapping, size = params$size, alpha = params$alpha)
      ),
      label = do.call(
        ggplot2::geom_label,
        list(data = layer_data, mapping = mapping, size = params$size, alpha = params$alpha)
      ),
      smooth = ggplot2::geom_smooth(
        mapping = mapping, method = params$method, se = params$se,
        alpha = params$alpha, linewidth = params$linewidth
      ),
      errorbar = ggplot2::geom_errorbar(
        mapping = mapping, alpha = params$alpha, linewidth = params$linewidth
      ),
      ribbon = ggplot2::geom_ribbon(mapping = mapping, alpha = params$alpha),
      hline = ggplot2::geom_hline(yintercept = params$yintercept, alpha = params$alpha, linewidth = params$linewidth),
      vline = ggplot2::geom_vline(xintercept = params$xintercept, alpha = params$alpha, linewidth = params$linewidth),
      stop("Unsupported figure geom")
    )
    plot <- plot + new_layer
  }

  add_scale <- function(aesthetic, transform) {
    if (identical(transform, "identity"))
      return(NULL)
    fields <- unique(unlist(lapply(spec$layers, function(layer) layer$mappings[[aesthetic]])))
    fields <- fields[!is.na(fields) & nzchar(fields)]
    if (!length(fields))
      return(NULL)
    field <- fields[1]
    values <- data[[field]]
    if (!is.numeric(values))
      stop("Log, square-root, and reverse transforms require a numeric ", aesthetic, " axis.")
    if (identical(transform, "reverse")) {
      if (aesthetic == "x") return(ggplot2::scale_x_reverse())
      return(ggplot2::scale_y_reverse())
    }
    trans <- switch(transform, log2 = "log2", log10 = "log10", sqrt = "sqrt")
    if (identical(transform, "log2") && any(values <= 0, na.rm = TRUE))
      stop("log2 axis transformation requires positive values.")
    if (identical(transform, "log10") && any(values <= 0, na.rm = TRUE))
      stop("log10 axis transformation requires positive values.")
    if (identical(transform, "sqrt") && any(values < 0, na.rm = TRUE))
      stop("sqrt axis transformation requires non-negative values.")
    if (aesthetic == "x")
      return(ggplot2::scale_x_continuous(trans = trans))
    ggplot2::scale_y_continuous(trans = trans)
  }
  plot <- plot + add_scale("x", spec$x_transform)
  plot <- plot + add_scale("y", spec$y_transform)

  color_field <- unique(unlist(lapply(spec$layers, function(x) {
    c(x$mappings$color, x$mappings$fill)
  })))
  color_field <- color_field[!is.na(color_field) & nzchar(color_field)]
  if (length(color_field) && !identical(spec$palette, "default")) {
    primary <- color_field[1]
    categorical <- !is.numeric(data[[primary]])
    if (identical(spec$palette, "grey")) {
      if (any(vapply(spec$layers, function(x) identical(x$mappings$color, primary), logical(1))))
        plot <- plot + ggplot2::scale_color_grey()
      if (any(vapply(spec$layers, function(x) identical(x$mappings$fill, primary), logical(1))))
        plot <- plot + ggplot2::scale_fill_grey()
    } else if (categorical) {
      brewer <- switch(
        spec$palette,
        colorblind = "Set2",
        sequential = "Blues",
        diverging = "RdBu"
      )
      if (any(vapply(spec$layers, function(x) identical(x$mappings$color, primary), logical(1))))
        plot <- plot + ggplot2::scale_color_brewer(palette = brewer)
      if (any(vapply(spec$layers, function(x) identical(x$mappings$fill, primary), logical(1))))
        plot <- plot + ggplot2::scale_fill_brewer(palette = brewer)
    } else {
      colors <- switch(
        spec$palette,
        colorblind = c("#2166ac", "#B2182B"),
        sequential = c("#F7FBFF", "#08519C"),
        diverging = c("#B2182B", "#2166ac")
      )
      if (any(vapply(spec$layers, function(x) identical(x$mappings$color, primary), logical(1))))
        plot <- plot + ggplot2::scale_color_gradient(low = colors[1], high = colors[2])
      if (any(vapply(spec$layers, function(x) identical(x$mappings$fill, primary), logical(1))))
        plot <- plot + ggplot2::scale_fill_gradient(low = colors[1], high = colors[2])
    }
  }

  if (!is.null(spec$facet_by)) {
    args <- list(facets = as.name(spec$facet_by))
    if (!is.null(spec$facet_ncol)) args$ncol <- spec$facet_ncol
    plot <- plot + do.call(ggplot2::facet_wrap, args)
  }

  label_args <- spec$labels[names(spec$labels) %in% c("title", "subtitle", "x", "y", "caption")]
  if (length(label_args))
    plot <- plot + do.call(ggplot2::labs, label_args)

  plot + switch(
    spec$theme,
    minimal = ggplot2::theme_minimal(base_size = 11),
    classic = ggplot2::theme_classic(base_size = 11),
    light = ggplot2::theme_light(base_size = 11),
    grey = ggplot2::theme_grey(base_size = 11),
    bw = ggplot2::theme_bw(base_size = 11)
  )
}

#' Render low-resolution and full-resolution PNG figures
#'
#' @param plot A ggplot object.
#' @param directory Session-specific output directory.
#' @param figure_id Sanitized figure identifier.
#'
#' @return Data URIs, dimensions, and encoded file sizes for both previews and
#'   high-resolution downloads. Temporary files are removed after encoding.
#' @keywords internal
#' @rdname agentFigureHelpers
agent_render_figure <- function(plot, directory, figure_id) {
  if (!dir.exists(directory))
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  preview <- file.path(directory, paste0(figure_id, "_preview.png"))
  full <- file.path(directory, paste0(figure_id, "_full.png"))

  on.exit(unlink(c(preview, full)), add = TRUE)
  ggplot2::ggsave(
    preview, plot = plot, width = 5.2, height = 3.9, units = "in",
    dpi = 110, bg = "white", limitsize = FALSE
  )
  ggplot2::ggsave(
    full, plot = plot, width = 8, height = 6, units = "in",
    dpi = 300, bg = "white", limitsize = FALSE
  )

  full_bytes <- unname(file.size(full))
  preview_bytes <- unname(file.size(preview))
  if (full_bytes > 10 * 1024^2)
    stop("High-resolution figure exceeds the 10 MB AI-figure limit; reduce layers, labels, or plotted rows.")
  if (preview_bytes > 2 * 1024^2)
    stop("Figure preview exceeds the 2 MB AI-figure limit; reduce layers, labels, or plotted rows.")

  preview_uri <- base64enc::dataURI(file = preview, mime = "image/png")
  full_uri <- base64enc::dataURI(file = full, mime = "image/png")
  list(
    preview_uri = preview_uri,
    preview_bytes = preview_bytes,
    full_uri = full_uri,
    full_bytes = full_bytes,
    preview_dimensions = c(width = 572L, height = 429L),
    full_dimensions = c(width = 2400L, height = 1800L)
  )
}

#' Create trusted HTML for a rendered assistant figure
#'
#' @param figure Figure result from \code{\link{agent_render_figure}}.
#' @param spec Normalized figure specification.
#' @param figure_id Figure identifier.
#'
#' @return htmltools tags containing a low-resolution preview and a
#'   high-resolution PNG download button. Model text is escaped by htmltools.
#' @keywords internal
#' @rdname agentFigureHelpers
agent_figure_html <- function(figure, spec, figure_id) {
  title <- if (nzchar(spec$labels$title %||% "")) spec$labels$title else paste("AI figure", figure_id)
  alt <- paste("AI-generated", spec$data_source, "figure:", title)
  filename <- paste0("omicsviewer-", figure_id, "-2400x1800.png")

  htmltools::div(
    class = "omicsviewer-ai-figure",
    htmltools::img(
      src = figure$preview_uri,
      alt = alt,
      style = "width:100%; height:auto; border:1px solid #d8dce0; border-radius:4px; background:#fff;"
    ),
    htmltools::div(
      style = "display:flex; align-items:center; justify-content:space-between; gap:8px; margin-top:6px;",
      htmltools::span(
        style = "font-size:12px; color:#56606a;",
        sprintf("%d x %d preview | %0.1f MB full PNG",
                figure$preview_dimensions[["width"]], figure$preview_dimensions[["height"]],
                figure$full_bytes / 1024^2)
      ),
      htmltools::a(
        href = figure$full_uri,
        download = filename,
        class = "btn btn-primary btn-sm",
        `aria-label` = paste("Download high-resolution figure", title),
        htmltools::HTML("&#8681; High-res PNG")
      )
    )
  )
}
