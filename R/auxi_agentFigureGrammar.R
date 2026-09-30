#' Single-source AI figure grammar table and schema builders
#'
#' todo 4.2: one grammar table (geom -> mappings/params/filters) is the
#' single source for the three artifacts that used to duplicate every
#' enum, range and default by hand:
#' \itemize{
#'   \item the ellmer tool-argument schema (\code{agent_figure_spec_schema}),
#'   \item the spec validator (\code{agent_normalize_figure_spec} reads the
#'     same tables through the helpers below),
#'   \item the model-facing prose (\code{agent_figure_grammar}).
#' }
#' The schema narrows filter nesting to
#' \code{.agent_figure_where_schema_depth} (2) while the validator keeps
#' accepting depth-3 trees (idempotent-normalizer precedent) so old
#' transcripts and snapshot replays still parse; the prose advertises the
#' schema contract.
#'
#' @return Grammar tables and ellmer type objects for the figure tools.
#' @keywords internal
#' @name agentFigureGrammar
NULL

# Model-facing filter contract. The SCHEMA and the prose grammar stop at
# depth 2 (one combinator around leaves); the VALIDATOR keeps accepting
# depth-3 trees (below) so pre-narrowing transcripts and golden replays
# still parse — acceptance, not advertisement.
.agent_figure_where_schema_depth <- 2L

#' Compact prose view of the layer-parameter table
#'
#' One string per parameter for the model-facing grammar section
#' (\code{agent_figure_grammar()$layer_params}); generated from
#' \code{\link{.agent_figure_param_specs}} so the prose can never drift
#' from the schema/validator.
#'
#' @return Named character vector.
#' @keywords internal
#' @rdname agentFigureGrammar
.agent_figure_param_descriptions <- function() {
  specs <- .agent_figure_param_specs()
  vapply(names(specs), function(nm) {
    spec <- specs[[nm]]
    default <- if (is.null(spec$default)) ""
      else paste0("; default ", spec$default)
    switch(spec$kind,
      number = paste0("Number ", spec$min, "-", spec$max, default, "."),
      integer = paste0("Integer ", spec$min, "-", spec$max, default, "."),
      enum = paste0("One of ", paste(spec$values, collapse = "/"), default, "."),
      boolean = paste0("True or false", default, "."),
      hex = "Constant hex color like '#b2182b'; mutually exclusive with mapping the same aesthetic.",
      order_by = "{column, decreasing}: rank rows before the max_labels cap; 'x'/'y' refer to the layer's own axis mapping."
    )
  }, character(1))
}

#' The layer-parameter grammar table
#'
#' One row per allowlisted \code{params} field of a figure layer. Order is
#' load-bearing: the schema emits properties in this order and the
#' validator normalizes in this order, so a spec with several invalid
#' parameters reports the same first error from either artifact.
#' \code{kind} dispatches both consumers:
#' \describe{
#'   \item{number / integer}{\code{min}, \code{max}, \code{default};
#'     absent means the default, out-of-range errors.}
#'   \item{enum}{\code{values}, \code{default}.}
#'   \item{boolean}{\code{default}.}
#'   \item{hex}{constant colour string, validated against
#'     \code{.agent_figure_hex_pattern}; default NULL.}
#'   \item{order_by}{structured render-time ordering (see
#'     \code{.agent_figure_order_by_normalize}).}
#' }
#' @return Named list of parameter specifications.
#' @keywords internal
#' @rdname agentFigureGrammar
.agent_figure_param_specs <- function() {
  list(
    alpha = list(
      kind = "number", min = 0, max = 1, default = 0.85,
      schema = "Transparency from 0 through 1."),
    size = list(
      kind = "number", min = 0.05, max = 12, default = 1.8,
      schema = "Point/text size from 0.05 through 12."),
    linewidth = list(
      kind = "number", min = 0.05, max = 6, default = 0.8,
      schema = "Line width from 0.05 through 6."),
    bins = list(
      kind = "integer", min = 5L, max = 100L, default = 30L,
      schema = "Histogram bins from 5 through 100."),
    method = list(
      kind = "enum", values = c("auto", "lm", "loess"), default = "auto",
      schema = "Smooth method."),
    se = list(
      kind = "boolean", default = TRUE,
      schema = "Show confidence interval for smooth."),
    position = list(
      kind = "enum", values = c("stack", "dodge", "fill", "jitter"),
      default = "stack",
      schema = "Position adjustment."),
    xintercept = list(
      kind = "number", min = -1e9, max = 1e9, default = NULL,
      schema = "Numeric vertical-line intercept."),
    yintercept = list(
      kind = "number", min = -1e9, max = 1e9, default = NULL,
      schema = "Numeric horizontal-line intercept."),
    max_labels = list(
      kind = "integer", min = 0L, max = 50L, default = 20L,
      schema = "Maximum text/label rows from 0 through 50."),
    order_by = list(
      kind = "order_by",
      schema = "Render-time row ordering for text/label layers, applied before max_labels."),
    color = list(
      kind = "hex", default = NULL,
      schema = "Constant hex color for every row of this layer, e.g. '#b2182b'."),
    fill = list(
      kind = "hex", default = NULL,
      schema = "Constant hex fill color for every row of this layer.")
  )
}

#' Required aesthetic mappings per geom
#'
#' Replaces the hand-maintained \code{switch(geom, ...)} in the validator;
#' the schema documents the same requirement through each aesthetic's
#' description and the prose grammar exposes it as
#' \code{geom_requirements}.
#'
#' @return Named list: geom -> character vector of required aesthetics.
#' @keywords internal
#' @rdname agentFigureGrammar
.agent_figure_geom_required <- function() {
  list(
    point = c("x", "y"),
    line = c("x", "y"),
    path = c("x", "y"),
    bar = "x",
    boxplot = c("x", "y"),
    violin = c("x", "y"),
    histogram = "x",
    density = "x",
    text = c("x", "y", "label"),
    label = c("x", "y", "label"),
    smooth = c("x", "y"),
    errorbar = c("x", "ymin", "ymax"),
    ribbon = c("x", "ymin", "ymax"),
    hline = character(),
    vline = character()
  )
}

#' The theme-option grammar table
#'
#' Fields of the bounded \code{theme_options} object; the schema builder
#' and \code{.agent_figure_theme_options_normalize} share these ranges.
#'
#' @return Named list of option specifications.
#' @keywords internal
#' @rdname agentFigureGrammar
.agent_figure_theme_option_specs <- function() {
  list(
    base_size = list(
      kind = "integer", min = 8L, max = 24L, default = NULL,
      schema = "Base font size from 8 through 24."),
    legend_position = list(
      kind = "enum", values = .agent_figure_legend_positions, default = NULL,
      schema = "Legend placement."),
    rotate_x_labels = list(
      kind = "number", min = 0, max = 90, default = NULL,
      schema = "Degrees to rotate x-axis labels, 0 through 90."),
    show_grid = list(
      kind = "boolean", default = NULL,
      schema = "Show panel grid lines; false hides them.")
  )
}

# ---- ellmer schema builders (todo 4.2: moved from module_aiAssistant) ----
# Pure functions of the grammar tables above — package level so the module
# only wires them into tool registrations (and 4.4's ToolRegistry split
# inherits them as-is). Field types are deliberately non-polymorphic
# (value = number, text = string, ...): ellmer NA-converts type-mismatched
# scalars inside arrays of objects, so a union field would arrive as NA.

#' Schema for one row-filter node
#'
#' Recursion stops at \code{.agent_figure_where_schema_depth} (2): one
#' combinator around leaf conditions. The validator still accepts one
#' more level for pre-narrowing replays.
#'
#' @param depth Current recursion depth (internal).
#' @return ellmer type object.
#' @keywords internal
#' @rdname agentFigureGrammar
agent_figure_where_schema <- function(depth = 1L) {
  args <- list(
    column = ellmer::type_string(
      "Exact data column, or 'x'/'y' for this layer's own axis mapping.",
      required = FALSE
    ),
    op = ellmer::type_enum(
      .agent_figure_where_ops, "Comparison operator.", required = FALSE
    ),
    value = ellmer::type_number(
      "Numeric comparison value, for > >= < <= abs> abs>= == !=.",
      required = FALSE
    ),
    text = ellmer::type_string(
      "Text comparison value, for == != starts_with ends_with.",
      required = FALSE
    ),
    min = ellmer::type_number("between lower bound.", required = FALSE),
    max = ellmer::type_number("between upper bound.", required = FALSE),
    values = ellmer::type_array(
      ellmer::type_string(
        "One value; numbers as strings on numeric columns."),
      paste0("1-", .agent_figure_where_max_values,
             " values, for in/not_in."),
      required = FALSE
    ),
    not = ellmer::type_boolean("Negate this condition.", required = FALSE)
  )
  if (depth < .agent_figure_where_schema_depth) {
    args$all <- ellmer::type_array(
      agent_figure_where_schema(depth + 1L),
      "ALL conditions must match.",
      required = FALSE
    )
    args$any <- ellmer::type_array(
      agent_figure_where_schema(depth + 1L),
      "ANY condition matches.",
      required = FALSE
    )
  }
  do.call(ellmer::type_object, c(list(
    "Row filter: one leaf condition or one all/any combinator.",
    .required = FALSE
  ), args))
}

.agent_figure_scale_channel_schema <- function(channel) {
  ellmer::type_object(
    paste0("Explicit ", channel, " scale override."),
    values = ellmer::type_array(
      ellmer::type_object(
        "One category-to-color pair.",
        category = ellmer::type_string("Exact category level."),
        color = ellmer::type_string("Hex color like '#b2182b'.")
      ),
      paste0("Discrete per-category colors (at most ",
             .agent_figure_scale_max_values, ")."),
      required = FALSE
    ),
    limits = ellmer::type_array(
      ellmer::type_string("Level name, or a number as a string."),
      "Discrete legend order, or exactly two numbers bounding a continuous range.",
      required = FALSE
    ),
    midpoint = ellmer::type_number(
      "Continuous diverging-gradient midpoint.",
      required = FALSE
    ),
    .required = FALSE
  )
}

#' Schema for the spec-level scale override
#' @return ellmer type object.
#' @keywords internal
#' @rdname agentFigureGrammar
agent_figure_scale_schema <- function() {
  ellmer::type_object(
    "Explicit scale overrides; beat the palette preset per channel.",
    color = .agent_figure_scale_channel_schema("color"),
    fill = .agent_figure_scale_channel_schema("fill"),
    .required = FALSE
  )
}

#' Schema for the full declarative figure specification
#'
#' Generated from the grammar tables; with patch-mode
#' \code{update_figure(figure_id, changes)} this full schema is advertised
#' only on \code{create_figure}'s advanced path (todo 3.1/4.2). The
#' \code{layers} element and \code{order_by} sub-object are shared with the
#' patch schema through \code{layers_required}.
#'
#' @param required Whether the whole object is a required argument.
#' @param layers_required Whether \code{layers} must be present (full-spec
#'   mode; the patch mode omits it).
#' @return ellmer type object.
#' @keywords internal
#' @rdname agentFigureGrammar
agent_figure_spec_schema <- function(required = TRUE, layers_required = TRUE) {
  param_args <- lapply(.agent_figure_param_specs(), function(spec) {
    switch(spec$kind,
      number = ellmer::type_number(spec$schema, required = FALSE),
      integer = ellmer::type_integer(spec$schema, required = FALSE),
      enum = ellmer::type_enum(spec$values, spec$schema, required = FALSE),
      boolean = ellmer::type_boolean(spec$schema, required = FALSE),
      hex = ellmer::type_string(spec$schema, required = FALSE),
      order_by = ellmer::type_object(
        spec$schema,
        column = ellmer::type_string(
          "Exact data column, or 'x'/'y' for this layer's own axis mapping.",
          required = FALSE
        ),
        decreasing = ellmer::type_boolean(
          "Sort descending (default true).", required = FALSE
        ),
        .required = FALSE
      ),
      stop("Unknown figure param kind: ", spec$kind)
    )
  })
  aesthetics <- lapply(.agent_figure_aesthetics, function(aesthetic) {
    ellmer::type_string(
      paste0("Exact data column for ", aesthetic, "."), required = FALSE
    )
  })
  names(aesthetics) <- .agent_figure_aesthetics
  theme_option_args <- lapply(.agent_figure_theme_option_specs(), function(spec) {
    switch(spec$kind,
      integer = ellmer::type_integer(spec$schema, required = FALSE),
      enum = ellmer::type_enum(spec$values, spec$schema, required = FALSE),
      number = ellmer::type_number(spec$schema, required = FALSE),
      boolean = ellmer::type_boolean(spec$schema, required = FALSE),
      stop("Unknown figure theme option kind: ", spec$kind)
    )
  })
  ellmer::type_object(
    "Declarative allowlisted ggplot2 figure specification. Fields map to validated omicsViewer rendering code, never arbitrary R.",
    data_source = ellmer::type_enum(
      .agent_figure_sources,
      "Plot data source.",
      required = FALSE
    ),
    features = ellmer::type_array(
      ellmer::type_string("Exact feature ID."),
      "Feature IDs to plot; expression figures default to the current semantic selection.",
      required = FALSE
    ),
    samples = ellmer::type_array(
      ellmer::type_string("Exact sample ID."),
      "Sample IDs to plot; expression figures default to all samples.",
      required = FALSE
    ),
    layers = ellmer::type_array(
      ellmer::type_object(
        "One allowlisted ggplot2 layer.",
        geom = ellmer::type_enum(.agent_figure_geoms, "Allowlisted geom."),
        x = ellmer::type_string("Exact data column for x.", required = FALSE),
        y = ellmer::type_string("Exact data column for y.", required = FALSE),
        color = ellmer::type_string("Exact data column for color.", required = FALSE),
        fill = ellmer::type_string("Exact data column for fill.", required = FALSE),
        group = ellmer::type_string("Exact data column for group.", required = FALSE),
        size = ellmer::type_string("Exact data column for size.", required = FALSE),
        alpha = ellmer::type_string("Exact data column for alpha.", required = FALSE),
        shape = ellmer::type_string("Exact data column for shape.", required = FALSE),
        linetype = ellmer::type_string("Exact data column for line type.", required = FALSE),
        label = ellmer::type_string("Exact data column for text labels.", required = FALSE),
        ymin = ellmer::type_string("Exact data column for minimum y.", required = FALSE),
        ymax = ellmer::type_string("Exact data column for maximum y.", required = FALSE),
        params = do.call(ellmer::type_object, c(
          list("Validated numeric/display parameters."), param_args,
          list(.required = FALSE)
        )),
        filter = agent_figure_where_schema()
      ),
      "One to twelve validated figure layers.",
      required = layers_required
    ),
    facet_by = ellmer::type_string("Exact facet column.", required = FALSE),
    facet_ncol = ellmer::type_integer("Facet columns from 1 through 6.", required = FALSE),
    x_transform = ellmer::type_enum(.agent_figure_transforms, "Allowlisted x-axis transform.", required = FALSE),
    y_transform = ellmer::type_enum(.agent_figure_transforms, "Allowlisted y-axis transform.", required = FALSE),
    theme = ellmer::type_enum(.agent_figure_themes, "Allowlisted ggplot2 theme.", required = FALSE),
    palette = ellmer::type_enum(.agent_figure_palettes, "Allowlisted color palette.", required = FALSE),
    scale = agent_figure_scale_schema(),
    theme_options = do.call(ellmer::type_object, c(
      list("Bounded tweaks applied on top of the chosen theme."),
      theme_option_args,
      list(.required = FALSE)
    )),
    labels = ellmer::type_object(
      "Escaped plot labels.",
      title = ellmer::type_string("Title (at most 200 characters).", required = FALSE),
      subtitle = ellmer::type_string("Subtitle (at most 200 characters).", required = FALSE),
      x = ellmer::type_string("X-axis label (at most 200 characters).", required = FALSE),
      y = ellmer::type_string("Y-axis label (at most 200 characters).", required = FALSE),
      caption = ellmer::type_string("Caption (at most 200 characters).", required = FALSE),
      .required = FALSE
    ),
    .required = required
  )
}
