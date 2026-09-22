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

#' Describe the allowlisted model-facing figure grammar
#'
#' @return A JSON-like description of supported data sources, geoms, aesthetic
#'   mappings, transformations, themes, palettes, and hard limits.
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
  # omitted optionals as empty objects; treat NA, length-0, and empty
  # list values like an omitted value.
  if (is.null(value) || length(value) == 0L || is.na(value[1])) return(default)
  value <- suppressWarnings(as.numeric(value)[1])
  if (is.na(value) || value < min || value > max)
    stop("Figure parameter ", name, " must be between ", min, " and ", max, ".")
  value
}

.agent_figure_integer_param <- function(value, name, min, max, default) {
  # ellmer converts JSON null to NA and some providers (glm flash) echo
  # omitted optionals as empty objects; treat NA, length-0, and empty
  # list values like an omitted value.
  if (is.null(value) || length(value) == 0L || is.na(value[1])) return(default)
  value <- suppressWarnings(as.integer(value)[1])
  if (is.na(value) || value < min || value > max)
    stop("Figure parameter ", name, " must be an integer between ", min, " and ", max, ".")
  value
}

.agent_figure_choice <- function(value, choices, name, fallback = NULL) {
  if (is.null(value) || length(value) == 0 || is.na(value) ||
      agent_sentinel_string(value))
    return(fallback)
  value <- .agent_figure_scalar(value, max_chars = 100L)
  if (is.null(value) || !value %in% choices)
    stop("Unsupported figure ", name, ": ", value)
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
