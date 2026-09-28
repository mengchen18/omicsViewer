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
# Widened grammar (2026-09-26): structured row predicates (layer `filter`),
# constant color/fill params, explicit scale overrides, and bounded theme
# tweaks. Still fully declarative — a fixed operator table interpreted
# server-side, never parsed or evaluated as code.
.agent_figure_where_ops <- c(
  ">", ">=", "<", "<=", "abs>", "abs>=", "==", "!=", "between",
  "in", "not_in", "is_na", "not_na", "starts_with", "ends_with"
)
.agent_figure_where_max_depth <- 3L
.agent_figure_where_max_leaves <- 8L
.agent_figure_where_max_values <- 50L
.agent_figure_scale_max_values <- 26L
.agent_figure_hex_pattern <- "^#[0-9a-fA-F]{3,8}$"
.agent_figure_legend_positions <- c("top", "bottom", "left", "right", "none")

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
    layer_filters = list(
      description = paste(
        "Optional row filter on any layer: only matching rows are drawn by",
        "that layer (applied before max_labels). Combine a highlight overlay,",
        "a labeled subset, or per-group coloring by stacking filtered layers",
        "with constant params.color/params.fill."
      ),
      form = paste(
        "One leaf condition {column, op, ...} or one combinator",
        "{all: [leaf, ...]} / {any: [leaf, ...]}; at most",
        .agent_figure_where_max_depth, "nesting levels and",
        .agent_figure_where_max_leaves, "leaves per filter."
      ),
      ops = list(
        `>` = "numeric column greater than value.",
        `>=` = "numeric column greater than or equal to value.",
        `<` = "numeric column less than value.",
        `<=` = "numeric column less than or equal to value.",
        `abs>` = "absolute numeric column value greater than value.",
        `abs>=` = "absolute numeric column value greater than or equal to value.",
        `==` = "equal to value (numeric columns) or text (text columns).",
        `!=` = "not equal to value (numeric columns) or text (text columns).",
        between = "numeric column between min and max (inclusive).",
        `in` = "column value one of values (send numbers as strings on numeric columns).",
        not_in = "column value none of values.",
        is_na = "column value is missing.",
        not_na = "column value is present.",
        starts_with = "text column value starts with text.",
        ends_with = "text column value ends with text."
      ),
      leaf_fields = list(
        column = paste(
          "Exact data column (feature__/sample__ prefixed in expression figures),",
          "or 'x'/'y' for this layer's own axis mapping."
        ),
        value = "Number, for > >= < <= abs> abs>= == !=.",
        text = "String, for == != starts_with ends_with on text columns.",
        min = "Number, between lower bound.",
        max = "Number, between upper bound.",
        values = paste("1-", .agent_figure_where_max_values,
                       " strings, for in/not_in.", sep = ""),
        not = "true negates the leaf condition."
      ),
      semantics = paste(
        "Rows with NA values never match a comparison; use is_na/not_na",
        "to match them explicitly."
      )
    ),
    constant_layer_params = list(
      color = "Hex color (e.g. #b2182b) for every row of the layer; mutually exclusive with mapping the color aesthetic.",
      fill = "Hex color for every row of the layer; mutually exclusive with mapping the fill aesthetic."
    ),
    scale_overrides = list(
      description = paste(
        "Explicit per-channel scale control; beats the palette preset for",
        "that channel. Discrete mappings take values (category-to-hex pairs);",
        "continuous mappings take limits (two numbers) or midpoint",
        "(diverging gradient around the midpoint)."
      ),
      channels = c("color", "fill")
    ),
    theme_options = list(
      description = "Bounded tweaks applied on top of the chosen theme.",
      fields = list(
        base_size = "Integer 8-24; base font size.",
        legend_position = "'top', 'bottom', 'left', 'right', or 'none'.",
        rotate_x_labels = "Number 0-90; degrees to rotate x-axis labels.",
        show_grid = "true/false; hide panel grid lines when false."
      )
    ),
    limits = list(
      max_layers = 12L,
      max_annotation_rows = 20000L,
      max_expression_features = 50L,
      max_expression_samples = 200L,
      max_text_labels = 50L,
      max_filter_depth = .agent_figure_where_max_depth,
      max_filter_leaves = .agent_figure_where_max_leaves,
      max_in_values = .agent_figure_where_max_values,
      max_scale_values = .agent_figure_scale_max_values,
      max_full_png_mb = 10
    )
  )
}

.agent_figure_scalar <- function(x, fallback = NULL, max_chars = 200L) {
  if (is.null(x) || length(x) == 0 || is.na(x))
    return(fallback)
  out <- trimws(as.character(x)[1])
  if (is.na(out)) return(fallback)
  if (nzchar(out) && nchar(out, type = "chars", allowNA = TRUE) > max_chars)
    out <- paste0(substr(out, 1L, max_chars), " ...")
  out
}

.agent_figure_numeric_param <- function(value, name, min, max, default) {
  # Absent (NULL / empty) means omitted and takes the default; every
  # other out-of-range or non-numeric value is a self-correctable error.
  if (.agent_param_absent(value)) return(default)
  value <- suppressWarnings(as.numeric(value)[1])
  if (is.na(value) || value < min || value > max)
    stop("Figure parameter ", name, " must be between ", min, " and ", max, ".")
  value
}

.agent_figure_integer_param <- function(value, name, min, max, default) {
  if (.agent_param_absent(value)) return(default)
  value <- suppressWarnings(as.integer(value)[1])
  if (is.na(value) || value < min || value > max)
    stop("Figure parameter ", name, " must be an integer between ", min, " and ", max, ".")
  value
}

.agent_param_absent <- function(value) {
  # Absent in the canonical tool-argument document: NULL, empty, or a
  # list (object or array) whose every element is itself absent -- the
  # model may echo an all-null object explicitly. isTRUE() guards
  # is.na()'s logical(0) on empty-list elements so the scalar branches
  # can never error on exotic shapes.
  if (is.null(value) || length(value) == 0L)
    return(TRUE)
  if (is.list(value))
    return(all(vapply(value, function(v) .agent_param_absent(v), logical(1))))
  isTRUE(is.na(value[1]))
}

.agent_figure_ids_param <- function(value, name) {
  # WP6b: optional ID-array template argument (features/samples).
  # Canonical transport delivers a character vector or NULL; trim and
  # drop empty entries, or return NULL when effectively omitted.
  if (.agent_param_absent(value)) return(NULL)
  if (is.list(value))
    value <- unlist(lapply(value, function(v) as.character(v)[1]), use.names = FALSE)
  else
    value <- as.character(value)
  value <- trimws(value[!is.na(value)])
  value <- value[nzchar(value)]
  if (!length(value)) return(NULL)
  if (anyDuplicated(value))
    stop("Figure ", name, " must be unique IDs; duplicate entries found.")
  value
}

.agent_figure_choice <- function(value, choices, name, fallback = NULL) {
  if (is.null(value) || length(value) == 0 || is.na(value))
    return(fallback)
  value <- .agent_figure_scalar(value, max_chars = 100L)
  if (is.null(value) || !value %in% choices)
    stop("Unsupported figure ", name, ": ", value, ".",
         .agent_suggest_text(value, choices))
  value
}

.agent_figure_selection <- function(x) {
  if (is.null(x) || length(x) == 0) return(character())
  x <- as.character(x)
  unique(x[!is.na(x) & nzchar(x)])
}

.agent_figure_ids <- function(x, valid_ids, label, max_n, allow_default = TRUE) {
  if (is.null(x)) {
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

# ---- structured row predicates (widened grammar) -------------------------
#
# A filter is a small JSON tree: leaves carry a column, an operator from a
# closed table, and typed value fields; combinators are all/any. Nothing is
# ever parsed or evaluated as code — agent_where_normalize canonicalizes the
# tree at spec-normalization time (validating operators, value types, caps,
# and x/y aesthetic references) and agent_where_eval interprets the canonical
# tree against the built plotting data. Field names are deliberately typed
# (value = number, text = string, min/max = numbers, values = string array)
# because ellmer NA-converts type-mismatched scalars inside arrays of
# objects — a polymorphic value field would silently arrive as NA.

.agent_where_number <- function(value, name, min = -1e12, max = 1e12) {
  if (.agent_param_absent(value)) stop("Figure filter ", name, " is required for this operator.")
  value <- suppressWarnings(as.numeric(value)[1])
  if (is.na(value) || value < min || value > max)
    stop("Figure filter ", name, " must be a number between ", min, " and ", max, ".")
  value
}

.agent_where_text <- function(value, name) {
  if (.agent_param_absent(value)) stop("Figure filter ", name, " is required for this operator.")
  out <- .agent_figure_scalar(value, max_chars = 200L)
  if (is.null(out)) stop("Figure filter ", name, " is required for this operator.")
  out
}

.agent_where_values <- function(values, name) {
  if (.agent_param_absent(values)) stop("Figure filter ", name, " is required for this operator.")
  values <- .agent_flatten_strings(values)
  values <- trimws(values[!is.na(values)])
  values <- values[nzchar(values)]
  if (!length(values))
    stop("Figure filter ", name, " needs at least one value ",
         "(send numbers as strings on numeric columns).")
  if (length(values) > .agent_figure_where_max_values)
    stop("Figure filter ", name, " accepts at most ",
         .agent_figure_where_max_values, " values.")
  values
}

# Flatten a values argument into a plain character vector: atomic
# vectors, lists of scalars, or nested lists.
.agent_flatten_strings <- function(values) {
  if (!is.list(values))
    return(as.character(values))
  unlist(lapply(values, function(v) {
    if (is.null(v) || length(v) == 0L) character(0)
    else if (is.list(v)) .agent_flatten_strings(v)
    else as.character(v)
  }), use.names = FALSE)
}

.agent_where_flag <- function(value, name) {
  if (.agent_param_absent(value)) return(FALSE)
  if (is.logical(value)) return(isTRUE(value[1]))
  value <- tolower(.agent_figure_scalar(value, max_chars = 10L) %||% "")
  if (!value %in% c("true", "false"))
    stop("Figure filter ", name, " must be true or false.")
  identical(value, "true")
}

# Canonicalize one filter object. mappings is the layer's own aesthetic
# mapping table, used to resolve the 'x'/'y' column shorthands.
agent_where_normalize <- function(where, mappings = NULL) {
  if (.agent_param_absent(where)) return(NULL)
  if (!is.list(where))
    stop("Figure layer filter must be an object.")

  all_value <- where$all
  any_value <- where$any
  has_all <- !.agent_param_absent(all_value)
  has_any <- !.agent_param_absent(any_value)
  if (has_all || has_any) {
    if (has_all && has_any)
      stop("Figure layer filter accepts either 'all' or 'any', not both.")
    children <- if (has_all) all_value else any_value
    if (!is.list(children) || !length(children))
      stop("Figure layer filter combinator must be a non-empty array of conditions.")
    out <- lapply(children, agent_where_normalize, mappings = mappings)
    out <- out[!vapply(out, is.null, logical(1))]
    # an all/any whose every child is absent is itself absent (the model
    # may echo an all-null combinator array explicitly)
    if (!length(out)) return(NULL)
    return(setNames(list(out), if (has_all) "all" else "any"))
  }

  allowed <- c("column", "op", "value", "text", "min", "max", "values", "not")
  # all/any were already deemed absent by the guard above; other keys count
  # only when their value is not absent
  present_names <- names(where)[vapply(names(where), function(n) {
    !identical(n, "all") && !identical(n, "any") &&
      !.agent_param_absent(where[[n]])
  }, logical(1))]
  unknown <- setdiff(present_names, allowed)
  if (length(unknown))
    stop("Unknown figure filter field(s): ", paste(unknown, collapse = ", "))

  column <- .agent_figure_scalar(where$column, max_chars = 500L)
  if (is.null(column))
    stop("Figure filter requires a column.")
  if (column %in% c("x", "y")) {
    mapped <- if (!is.null(mappings)) mappings[[column]] else NULL
    mapped <- .agent_figure_scalar(mapped, max_chars = 500L)
    if (is.null(mapped))
      stop("Figure filter column '", column, "' refers to this layer's ", column,
           " aesthetic, which the layer does not map.")
    column <- mapped
  }

  op <- .agent_figure_choice(where$op, .agent_figure_where_ops, "filter operator")
  if (is.null(op))
    stop("Figure filter requires an operator from: ",
         paste(.agent_figure_where_ops, collapse = ", "), ".")

  # presence must be tested with the absence helper, not names(): the
  # model may echo an all-null object with every field explicitly null
  given <- Filter(function(n) !.agent_param_absent(where[[n]]),
                  c("value", "text", "min", "max", "values"))
  needed <- switch(
    op,
    `>` = , `>=` = , `<` = , `<=` = , `abs>` = , `abs>=` = "value",
    between = c("min", "max"),
    `==` = , `!=` = c("value", "text"),
    `in` = , `not_in` = "values",
    `is_na` = , `not_na` = character(),
    starts_with = , ends_with = "text"
  )
  if (!length(needed)) {
    if (length(given))
      stop("Figure filter operator ", op, " takes no comparison value.")
  } else if (identical(needed, c("value", "text"))) {
    if (length(given) != 1L || !given %in% needed)
      stop("Figure filter operator ", op, " takes exactly one of value (number) or text (string).")
  } else {
    missing <- setdiff(needed, given)
    if (length(missing))
      stop("Figure filter operator ", op, " requires field(s): ",
           paste(missing, collapse = ", "), ".")
    extra <- setdiff(given, needed)
    if (length(extra))
      stop("Figure filter operator ", op, " does not accept field(s): ",
           paste(extra, collapse = ", "), ".")
  }

  leaf <- list(column = column, op = op)
  if (identical(op, "between")) {
    leaf$min <- .agent_where_number(where$min, "min")
    leaf$max <- .agent_where_number(where$max, "max")
    if (leaf$min > leaf$max)
      stop("Figure filter between requires min <= max.")
  } else if (op %in% c("==", "!=")) {
    if ("value" %in% given) leaf$value <- .agent_where_number(where$value, "value")
    else leaf$text <- .agent_where_text(where$text, "text")
  } else if (op %in% c(">", ">=", "<", "<=", "abs>", "abs>=")) {
    leaf$value <- .agent_where_number(where$value, "value")
  } else if (op %in% c("in", "not_in")) {
    leaf$values <- .agent_where_values(where$values, "values")
  } else if (op %in% c("starts_with", "ends_with")) {
    leaf$text <- .agent_where_text(where$text, "text")
  }
  if (.agent_where_flag(where$not, "not")) leaf$not <- TRUE
  leaf
}

.agent_where_depth <- function(node) {
  if (!is.null(node$all) || !is.null(node$any))
    return(1L + max(vapply(.agent_where_children(node), .agent_where_depth, integer(1))))
  1L
}

.agent_where_children <- function(node) {
  if (!is.null(node$all)) return(node$all)
  node$any
}

.agent_where_leaves <- function(node) {
  if (!is.null(node$all) || !is.null(node$any))
    return(sum(vapply(.agent_where_children(node), .agent_where_leaves, integer(1))))
  1L
}

.agent_where_validate_caps <- function(node) {
  if (.agent_where_depth(node) > .agent_figure_where_max_depth)
    stop("Figure layer filter nests at most ", .agent_figure_where_max_depth, " levels deep.")
  if (.agent_where_leaves(node) > .agent_figure_where_max_leaves)
    stop("Figure layer filter accepts at most ", .agent_figure_where_max_leaves, " conditions.")
  invisible(node)
}

# Interpret a canonical filter against the built plotting data. Returns a
# settled logical vector (no NAs): rows with NA values never match a
# comparison — is_na/not_na are the explicit way to match them.
agent_where_eval <- function(where, data) {
  if (is.null(where)) return(rep(TRUE, nrow(data)))
  hit <- .agent_where_eval_node(where, data)
  hit[is.na(hit)] <- FALSE
  hit
}

.agent_where_eval_node <- function(node, data) {
  if (!is.null(node$all) || !is.null(node$any)) {
    hits <- lapply(.agent_where_children(node), .agent_where_eval_node, data)
    if (any(vapply(hits, function(h) length(h) != nrow(data), logical(1))))
      stop("Figure filter produced an invalid result.")
    if (!is.null(node$all)) return(Reduce(`&`, hits))
    return(Reduce(`|`, hits))
  }

  column <- .agent_figure_scalar(node$column, max_chars = 500L)
  if (is.null(column) || !column %in% colnames(data))
    stop("Figure filter column is unavailable: ", column, ".",
         .agent_suggest_text(column, colnames(data)))
  vals <- data[[column]]
  op <- node$op

  require_numeric <- function() {
    if (!is.numeric(vals))
      stop("Figure filter operator ", op, " requires a numeric column, but '",
           column, "' is not numeric.")
  }
  numeric_text <- function(value, text) {
    out <- suppressWarnings(as.numeric(if (is.null(value)) text else value)[1])
    if (is.na(out))
      stop("Figure filter compares numeric column '", column,
           "' with a non-numeric value; pass it as a number.")
    out
  }

  hit <- switch(
    op,
    `>` = , `>=` = , `<` = , `<=` = , `abs>` = , `abs>=` = {
      require_numeric()
      lhs <- if (op %in% c("abs>", "abs>=")) abs(vals) else vals
      switch(op,
             `>` = lhs > node$value,
             `>=` = lhs >= node$value,
             `<` = lhs < node$value,
             `<=` = lhs <= node$value,
             `abs>` = lhs > node$value,
             `abs>=` = lhs >= node$value)
    },
    between = {
      require_numeric()
      vals >= node$min & vals <= node$max
    },
    `==` = , `!=` = {
      if (is.numeric(vals)) {
        rhs <- numeric_text(node$value, node$text)
        if (identical(op, "==")) vals == rhs else vals != rhs
      } else {
        rhs <- as.character(if (is.null(node$text)) node$value else node$text)
        if (identical(op, "==")) as.character(vals) == rhs else as.character(vals) != rhs
      }
    },
    `in` = , `not_in` = {
      rhs <- if (is.numeric(vals)) {
        nums <- suppressWarnings(as.numeric(node$values))
        if (anyNA(nums))
          stop("Figure filter in/not_in on numeric column '", column,
               "' needs numeric values (send them as strings, e.g. \"1.5\").")
        nums
      } else {
        node$values
      }
      hit <- vals %in% rhs
      if (identical(op, "not_in")) hit <- !hit & !is.na(vals)
      hit
    },
    is_na = is.na(vals),
    not_na = !is.na(vals),
    starts_with = , ends_with = {
      if (!is.character(vals) && !is.factor(vals))
        stop("Figure filter operator ", op, " requires a text column, but '",
             column, "' is not text.")
      fun <- if (identical(op, "starts_with")) startsWith else endsWith
      fun(as.character(vals), node$text)
    },
    stop("Unsupported figure filter operator: ", op, ".")
  )

  if (isTRUE(node$not)) hit <- !hit
  hit
}

# ---- constant colors, scale overrides, theme tweaks --------------------

.agent_figure_color_scalar <- function(value, what) {
  out <- .agent_figure_scalar(value, max_chars = 9L)
  if (is.null(out)) return(NULL)
  if (!grepl(.agent_figure_hex_pattern, out))
    stop("Figure ", what, " must be a hex color like '#b2182b', not: ", out, ".")
  tolower(out)
}

.agent_figure_clean_strings <- function(value, name, max_n, empty_hint = NULL) {
  if (.agent_param_absent(value)) return(NULL)
  value <- .agent_flatten_strings(value)
  value <- trimws(value[!is.na(value)])
  value <- value[nzchar(value)]
  if (!length(value))
    stop("Figure ", name, " needs at least one entry",
         if (!is.null(empty_hint)) paste0(" ", empty_hint) else "", ".")
  if (length(value) > max_n)
    stop("Figure ", name, " accepts at most ", max_n, " entries.")
  value
}

# Discrete per-category colors: schema shape is an array of
# {category, color} pairs (ellmer type_object cannot declare dynamic-key
# maps); named lists/vectors from programmatic callers are accepted too.
# Returns a named character vector (category -> lowercase hex).
.agent_figure_scale_values <- function(values) {
  if (.agent_param_absent(values)) return(NULL)
  if (is.list(values)) {
    scalar <- vapply(values, function(v) is.atomic(v), logical(1))
    if (all(scalar)) {
      category <- names(values)
      color <- vapply(values, function(v) as.character(v)[1], character(1))
    } else {
      category <- vapply(values, function(v) as.character(v$category)[1], character(1))
      color <- vapply(values, function(v) as.character(v$color)[1], character(1))
    }
  } else {
    category <- names(values)
    color <- as.character(values)
  }
  keep <- !is.na(category) & nzchar(trimws(category)) & !is.na(color)
  category <- trimws(category[keep])
  color <- color[keep]
  if (!length(category))
    stop("Figure scale values need at least one category/color pair.")
  if (anyDuplicated(category))
    stop("Figure scale values contain duplicate categories: ",
         paste(utils::head(category[duplicated(category)], 3L), collapse = ", "), ".")
  if (length(category) > .agent_figure_scale_max_values)
    stop("Figure scale values accept at most ", .agent_figure_scale_max_values,
         " categories; merge rare levels or use a mapped palette instead.")
  color <- vapply(seq_along(color), function(i)
    .agent_figure_color_scalar(color[i], paste0("scale color for '", category[i], "'")),
    character(1))
  setNames(color, category)
}

.agent_figure_scale_normalize <- function(scale) {
  if (.agent_param_absent(scale)) return(NULL)
  if (!is.list(scale))
    stop("Figure scale override must be an object.")
  unknown <- setdiff(names(scale), c("color", "fill"))
  if (length(unknown))
    stop("Unknown figure scale field(s): ", paste(unknown, collapse = ", "))
  out <- list()
  for (channel in c("color", "fill")) {
    channel_spec <- scale[[channel]]
    if (.agent_param_absent(channel_spec)) next
    if (!is.list(channel_spec))
      stop("Figure scale ", channel, " override must be an object.")
    unknown_channel <- setdiff(names(channel_spec), c("values", "limits", "midpoint"))
    if (length(unknown_channel))
      stop("Unknown figure scale ", channel, " field(s): ",
           paste(unknown_channel, collapse = ", "))
    channel_out <- list(
      values = .agent_figure_scale_values(channel_spec$values),
      limits = .agent_figure_clean_strings(
        channel_spec$limits, paste0("scale ", channel, " limits"),
        .agent_figure_scale_max_values)
    )
    if (!.agent_param_absent(channel_spec$midpoint)) {
      midpoint <- suppressWarnings(as.numeric(channel_spec$midpoint)[1])
      if (is.na(midpoint))
        stop("Figure scale ", channel, " midpoint must be a number.")
      channel_out$midpoint <- midpoint
    }
    if (is.null(channel_out$values) && is.null(channel_out$limits) &&
        is.null(channel_out$midpoint))
      stop("Figure scale ", channel,
           " override requires values, limits, or midpoint.")
    if (!is.null(channel_out$values) && !is.null(channel_out$midpoint))
      stop("Figure scale ", channel,
           " values apply to discrete mappings while midpoint applies to continuous ones; pick one.")
    channel_out <- channel_out[!vapply(channel_out, is.null, logical(1))]
    out[[channel]] <- channel_out
  }
  if (!length(out))
    stop("Figure scale override requires a color or fill channel.")
  out
}

.agent_figure_theme_options_normalize <- function(options) {
  if (.agent_param_absent(options)) return(NULL)
  if (!is.list(options))
    stop("Figure theme options must be an object.")
  unknown <- setdiff(names(options),
                     c("base_size", "legend_position", "rotate_x_labels", "show_grid"))
  if (length(unknown))
    stop("Unknown figure theme option(s): ", paste(unknown, collapse = ", "))
  out <- list()
  base_size <- .agent_figure_integer_param(
    options$base_size, "theme base_size", 8L, 24L, NULL)
  if (!is.null(base_size)) out$base_size <- base_size
  legend_position <- .agent_figure_choice(
    options$legend_position, .agent_figure_legend_positions,
    "theme legend_position")
  if (!is.null(legend_position)) out$legend_position <- legend_position
  if (!.agent_param_absent(options$rotate_x_labels)) {
    rot <- options$rotate_x_labels
    if (is.logical(rot)) rot <- if (isTRUE(rot[1])) 45 else 0
    rot <- suppressWarnings(as.numeric(rot)[1])
    if (is.na(rot) || rot < 0 || rot > 90)
      stop("Figure theme option rotate_x_labels must be between 0 and 90 degrees.")
    out$rotate_x_labels <- rot
  }
  if (!.agent_param_absent(options$show_grid)) {
    flag <- options$show_grid
    parsed <- tolower(as.character(flag)[1])
    if (is.logical(flag)) parsed <- if (isTRUE(flag[1])) "true" else "false"
    if (!parsed %in% c("true", "false"))
      stop("Figure theme option show_grid must be true or false.")
    out$show_grid <- identical(parsed, "true")
  }
  if (!length(out)) return(NULL)
  out
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
    if (!is.null(layer$filter)) out$filter <- layer$filter
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
    "x_transform", "y_transform", "theme", "palette", "labels",
    "scale", "theme_options"
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

  # Layers arrive as a list of layer objects in the canonical
  # tool-argument document (agent_args_sanitize at the tool seam).
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
    allowed_layer <- c("geom", .agent_figure_aesthetics, "params", "filter")
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
      "xintercept", "yintercept", "max_labels", "color", "fill"
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

    # widened grammar: constant per-layer colors (hex-validated) and a
    # structured row filter. A constant color is mutually exclusive with
    # mapping the same aesthetic — otherwise the constant silently wins
    # and the model cannot tell why its mapped legend disappeared.
    params$color <- .agent_figure_color_scalar(params$color, "layer params.color")
    params$fill <- .agent_figure_color_scalar(params$fill, "layer params.fill")
    if (!is.null(params$color) && !is.null(mappings$color))
      stop("Figure layer maps color to '", mappings$color,
           "' and also sets constant params.color; remove one of the two.")
    if (!is.null(params$fill) && !is.null(mappings$fill))
      stop("Figure layer maps fill to '", mappings$fill,
           "' and also sets constant params.fill; remove one of the two.")

    filter <- agent_where_normalize(layer$filter, mappings = mappings)
    if (!is.null(filter)) .agent_where_validate_caps(filter)

    list(geom = geom, mappings = mappings, params = params, filter = filter)
  })

  facet_by <- .agent_figure_scalar(spec$facet_by, max_chars = 500L)
  facet_ncol <- .agent_figure_integer_param(spec$facet_ncol, "facet_ncol", 1L, 6L, NULL)
  x_transform <- .agent_figure_choice(
    spec$x_transform, .agent_figure_transforms, "x transform", "identity"
  )
  y_transform <- .agent_figure_choice(
    spec$y_transform, .agent_figure_transforms, "y transform", "identity"
  )
  # exact access: list `$` partial matching would resolve the absent
  # `theme` field to `theme_options`
  theme <- .agent_figure_choice(
    spec[["theme"]], .agent_figure_themes, "theme", "minimal")
  palette <- .agent_figure_choice(spec[["palette"]], .agent_figure_palettes, "palette", "default")

  labels <- if (is.list(spec$labels)) spec$labels else list()
  unknown_labels <- setdiff(names(labels), c("title", "subtitle", "x", "y", "caption"))
  if (length(unknown_labels))
    stop("Unknown figure label field(s): ", paste(unknown_labels, collapse = ", "))
  labels <- lapply(labels, function(x) .agent_figure_scalar(x, max_chars = 200L))
  labels <- labels[!vapply(labels, is.null, logical(1))]

  scale <- .agent_figure_scale_normalize(spec[["scale"]])
  theme_options <- .agent_figure_theme_options_normalize(spec[["theme_options"]])

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
    labels = labels,
    scale = scale,
    theme_options = theme_options
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

    # widened grammar: the structured row filter composites before the
    # max_labels cap (filter first, then cap the surviving rows)
    layer_data <- NULL
    if (!is.null(layer$filter)) {
      hit <- agent_where_eval(layer$filter, data)
      if (!any(hit))
        warning("Figure layer ", geom, " filter matches no rows; the layer is empty.")
      layer_data <- data[hit, , drop = FALSE]
    }
    if (geom %in% c("text", "label") && params$max_labels > 0L) {
      layer_data <- utils::head(
        if (is.null(layer_data)) data else layer_data, params$max_labels
      )
    }

    # constant per-layer colors (hex-validated at normalization time)
    constant <- list()
    if (!is.null(params$color)) constant$color <- params$color
    if (!is.null(params$fill)) constant$fill <- params$fill
    # filtered layers draw their subset explicitly; unfiltered ones inherit
    # the plot data (passing data = NULL would be equivalent)
    inherit <- if (is.null(layer_data)) list() else list(data = layer_data)

    new_layer <- switch(
      geom,
      point = if (identical(params$position, "jitter")) {
        do.call(ggplot2::geom_jitter, c(inherit, list(
          mapping = mapping, alpha = params$alpha, size = params$size
        ), constant))
      } else {
        do.call(ggplot2::geom_point, c(inherit, list(
          mapping = mapping, alpha = params$alpha, size = params$size
        ), constant))
      },
      line = do.call(ggplot2::geom_line, c(inherit, list(
        mapping = mapping, alpha = params$alpha, linewidth = params$linewidth
      ), constant)),
      path = do.call(ggplot2::geom_path, c(inherit, list(
        mapping = mapping, alpha = params$alpha, linewidth = params$linewidth
      ), constant)),
      bar = {
        args <- c(inherit, list(mapping = mapping, alpha = params$alpha), constant)
        if (!is.null(mappings$y)) args$stat <- "identity"
        if (params$position %in% c("stack", "dodge", "fill"))
          args$position <- params$position
        do.call(ggplot2::geom_bar, args)
      },
      boxplot = do.call(ggplot2::geom_boxplot, c(inherit, list(
        mapping = mapping, alpha = params$alpha
      ), constant)),
      violin = do.call(ggplot2::geom_violin, c(inherit, list(
        mapping = mapping, alpha = params$alpha
      ), constant)),
      histogram = do.call(ggplot2::geom_histogram, c(inherit, list(
        mapping = mapping, alpha = params$alpha, bins = params$bins,
        position = if (params$position %in% c("stack", "dodge", "fill")) params$position else "stack"
      ), constant)),
      density = do.call(ggplot2::geom_density, c(inherit, list(
        mapping = mapping, alpha = params$alpha, linewidth = params$linewidth
      ), constant)),
      text = do.call(
        ggplot2::geom_text,
        c(list(data = layer_data, mapping = mapping, size = params$size, alpha = params$alpha), constant)
      ),
      label = do.call(
        ggplot2::geom_label,
        c(list(data = layer_data, mapping = mapping, size = params$size, alpha = params$alpha), constant)
      ),
      smooth = do.call(ggplot2::geom_smooth, c(inherit, list(
        mapping = mapping, method = params$method, se = params$se,
        alpha = params$alpha, linewidth = params$linewidth
      ), constant)),
      errorbar = do.call(ggplot2::geom_errorbar, c(inherit, list(
        mapping = mapping, alpha = params$alpha, linewidth = params$linewidth
      ), constant)),
      ribbon = do.call(ggplot2::geom_ribbon, c(inherit, list(
        mapping = mapping, alpha = params$alpha
      ), constant)),
      hline = do.call(ggplot2::geom_hline, c(list(
        yintercept = params$yintercept, alpha = params$alpha, linewidth = params$linewidth
      ), constant)),
      vline = do.call(ggplot2::geom_vline, c(list(
        xintercept = params$xintercept, alpha = params$alpha, linewidth = params$linewidth
      ), constant)),
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

  mapped_fields <- function(channel) {
    fields <- unique(unlist(lapply(spec$layers, function(x) x$mappings[[channel]])))
    fields[!is.na(fields) & nzchar(fields)]
  }
  color_field <- unique(unlist(lapply(spec$layers, function(x) {
    c(x$mappings$color, x$mappings$fill)
  })))
  color_field <- color_field[!is.na(color_field) & nzchar(color_field)]
  gradient_colors <- switch(
    spec$palette,
    colorblind = c("#2166ac", "#B2182B"),
    sequential = c("#F7FBFF", "#08519C"),
    diverging = c("#B2182B", "#2166ac"),
    c("#132b43", "#56b1f7")
  )

  for (channel in c("color", "fill")) {
    fields <- mapped_fields(channel)
    override <- if (is.null(spec$scale)) NULL else spec$scale[[channel]]
    if (!is.null(override)) {
      # widened grammar: an explicit scale override beats the palette
      # preset for this channel
      if (!length(fields))
        stop("Figure scale ", channel,
             " override requires a layer that maps ", channel, ".")
      plot <- plot + .agent_scale_override(
        override, channel, data[[fields[1]]], fields[1], gradient_colors
      )
      next
    }
    if (!length(fields) || identical(spec$palette, "default") ||
        !identical(fields[1], color_field[1]))
      next
    primary <- fields[1]
    categorical <- !is.numeric(data[[primary]])
    if (identical(spec$palette, "grey")) {
      if (any(vapply(spec$layers, function(x) identical(x$mappings[[channel]], primary), logical(1))))
        plot <- plot + if (identical(channel, "color"))
          ggplot2::scale_color_grey() else ggplot2::scale_fill_grey()
    } else if (categorical) {
      brewer <- switch(
        spec$palette,
        colorblind = "Set2",
        sequential = "Blues",
        diverging = "RdBu"
      )
      if (any(vapply(spec$layers, function(x) identical(x$mappings[[channel]], primary), logical(1))))
        plot <- plot + if (identical(channel, "color"))
          ggplot2::scale_color_brewer(palette = brewer) else
          ggplot2::scale_fill_brewer(palette = brewer)
    } else {
      if (any(vapply(spec$layers, function(x) identical(x$mappings[[channel]], primary), logical(1))))
        plot <- plot + if (identical(channel, "color"))
          ggplot2::scale_color_gradient(low = gradient_colors[1], high = gradient_colors[2]) else
          ggplot2::scale_fill_gradient(low = gradient_colors[1], high = gradient_colors[2])
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

  theme_opts <- spec$theme_options
  base_size <- if (is.null(theme_opts) || is.null(theme_opts$base_size)) 11 else theme_opts$base_size
  plot <- plot + switch(
    spec$theme,
    minimal = ggplot2::theme_minimal(base_size = base_size),
    classic = ggplot2::theme_classic(base_size = base_size),
    light = ggplot2::theme_light(base_size = base_size),
    grey = ggplot2::theme_grey(base_size = base_size),
    bw = ggplot2::theme_bw(base_size = base_size)
  )
  tweaks <- list()
  if (!is.null(theme_opts)) {
    if (!is.null(theme_opts$legend_position))
      tweaks$legend.position <- theme_opts$legend_position
    if (identical(theme_opts$show_grid, FALSE)) {
      tweaks$panel.grid.major <- ggplot2::element_blank()
      tweaks$panel.grid.minor <- ggplot2::element_blank()
    }
    if (!is.null(theme_opts$rotate_x_labels) && theme_opts$rotate_x_labels > 0)
      tweaks$axis.text.x <- ggplot2::element_text(
        angle = theme_opts$rotate_x_labels, hjust = 1
      )
  }
  if (length(tweaks))
    plot <- plot + do.call(ggplot2::theme, tweaks)
  plot
}

# Apply one explicit scale override channel. Discrete mappings take
# per-category hex values (unmatched names warn; missing levels render
# grey) and optional legend limits; numeric mappings take a limits pair
# or a diverging midpoint. warnings() raised here are collected by the
# tool-level handler and surfaced in the figure result.
.agent_scale_override <- function(override, channel, values, field, gradient_colors) {
  if (is.numeric(values)) {
    if (!is.null(override$values))
      stop("Figure scale ", channel, " values require a discrete mapping; '",
           field, "' is numeric. Use limits or midpoint instead.")
    if (!is.null(override$midpoint)) {
      gradient2 <- if (identical(channel, "color"))
        ggplot2::scale_color_gradient2 else ggplot2::scale_fill_gradient2
      return(gradient2(
        low = gradient_colors[2], mid = "#f7f7f7", high = gradient_colors[1],
        midpoint = override$midpoint
      ))
    }
    limits <- suppressWarnings(as.numeric(override$limits))
    if (is.null(override$limits) || length(limits) != 2L || anyNA(limits))
      stop("Figure scale ", channel, " override on numeric '", field,
           "' requires limits as two numbers or a midpoint.")
    gradient <- if (identical(channel, "color"))
      ggplot2::scale_color_gradient else ggplot2::scale_fill_gradient
    return(gradient(low = gradient_colors[1], high = gradient_colors[2], limits = limits))
  }
  if (!is.null(override$midpoint))
    stop("Figure scale ", channel, " midpoint requires a numeric mapping; '",
         field, "' is discrete.")
  levels_ <- if (is.factor(values)) levels(values) else
    unique(as.character(values[!is.na(values)]))
  if (!is.null(override$limits)) {
    bad <- setdiff(override$limits, levels_)
    if (length(bad))
      stop("Unknown figure scale ", channel, " limits: ",
           paste(utils::head(bad, 3L), collapse = ", "), ".",
           .agent_suggest_text(bad[1], levels_))
  }
  if (is.null(override$values)) {
    # limits-only discrete override: reorder/filter the legend, keep the
    # default colors
    discrete <- if (identical(channel, "color"))
      ggplot2::scale_color_discrete else ggplot2::scale_fill_discrete
    return(discrete(limits = override$limits))
  }
  unmatched <- setdiff(names(override$values), levels_)
  if (length(unmatched))
    warning("Figure scale ", channel, " values name levels not present in '",
            field, "': ", paste(utils::head(unmatched, 3L), collapse = ", "), ".")
  missing_levels <- setdiff(levels_, names(override$values))
  if (length(missing_levels))
    warning("Figure scale ", channel, " has no explicit color for: ",
            paste(utils::head(missing_levels, 3L), collapse = ", "),
            "; those levels render grey.")
  manual <- if (identical(channel, "color"))
    ggplot2::scale_color_manual else ggplot2::scale_fill_manual
  args <- list(values = override$values)
  if (!is.null(override$limits)) args$limits <- override$limits
  do.call(manual, args)
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
