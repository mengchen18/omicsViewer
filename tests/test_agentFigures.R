library(omicsViewer)
library(unittest, quietly = TRUE)

agent_figure_grammar <- omicsViewer:::agent_figure_grammar
agent_normalize_figure_spec <- omicsViewer:::agent_normalize_figure_spec
agent_figure_spec_echo <- omicsViewer:::agent_figure_spec_echo
agent_build_figure_data <- omicsViewer:::agent_build_figure_data
agent_build_figure_plot <- omicsViewer:::agent_build_figure_plot
agent_render_figure <- omicsViewer:::agent_render_figure
agent_figure_html <- omicsViewer:::agent_figure_html

fd <- data.frame(
  score = 1:10,
  category = rep(c("kinase", "phosphatase"), 5),
  row.names = paste0("F", 1:10),
  check.names = FALSE
)
pd <- data.frame(
  group = rep(c("WT", "KO"), 5),
  row.names = paste0("S", 1:10),
  check.names = FALSE
)
mat <- matrix(rnorm(100), nrow = 10, dimnames = list(rownames(fd), rownames(pd)))

spec <- list(
  data_source = "expression",
  features = paste0("F", 1:3),
  layers = list(
    list(
      geom = "boxplot",
      x = "sample__group",
      y = "__expression__",
      fill = "sample__group",
      params = list(alpha = 0.65)
    )
  ),
  labels = list(title = "Selected feature expression", x = "Group", y = "Expression"),
  theme = "minimal",
  palette = "colorblind"
)

grammar <- agent_figure_grammar()
ok(
  ut_cmp_identical(grammar$limits$max_expression_features, 50L),
  "figure grammar exposes bounded expression features"
)
ok(
  ut_cmp_identical("__expression__" %in% grammar$reserved_columns, TRUE),
  "figure grammar reserves the expression value column"
)

normalized <- agent_normalize_figure_spec(
  spec, fd, pd, mat, character(), character()
)
ok(
  ut_cmp_identical(normalized$data_source, "expression"),
  "expression figure source is normalized"
)

# provider sentinel sweep (glm flash serializes omitted optionals as
# literal "null"/"{}"/"[]" strings): every sentinel falls back to the
# omitted-value default instead of erroring or rendering literally
sentinel_spec <- list(
  data_source = "expression",
  features = "null",
  samples = "[]",
  layers = list(
    list(geom = "point", x = "feature__score", y = "__expression__",
         color = "null")
  ),
  theme = "null",
  palette = "{}",
  facet_by = "null",
  x_transform = "null",
  y_transform = "[]",
  labels = list(title = "null", caption = "{}", x = "score")
)
sentinel_normalized <- agent_normalize_figure_spec(
  sentinel_spec, fd, pd, mat, paste0("F", 1:3), character()
)
ok(
  ut_cmp_identical(sentinel_normalized$theme, "minimal") &&
    ut_cmp_identical(sentinel_normalized$palette, "default") &&
    ut_cmp_identical(sentinel_normalized$x_transform, "identity") &&
    ut_cmp_identical(sentinel_normalized$y_transform, "identity") &&
    is.null(sentinel_normalized$facet_by),
  "sentinel figure choices fall back to defaults"
)
ok(
  is.null(sentinel_normalized$labels$title) &&
    is.null(sentinel_normalized$labels$caption) &&
    ut_cmp_identical(sentinel_normalized$labels$x, "score"),
  "sentinel figure labels are dropped, real labels kept"
)
ok(
  is.null(sentinel_normalized$layers[[1]]$mappings$color) &&
    ut_cmp_identical(sentinel_normalized$layers[[1]]$mappings$x, "feature__score"),
  "sentinel aesthetic mappings are dropped"
)

# glm flash also echoes omitted optionals as EMPTY OBJECTS on the spec
# round-trip path (observed live 2026-09-24: facet_ncol = {} rejected
# twice with "must be an integer between 1 and 6"): length-0 values must
# behave exactly like omitted values in every param helper
empty_object_spec <- list(
  data_source = "feature_annotation",
  layers = list(list(
    geom = "point", x = "score", y = "score",
    params = list(alpha = list(), size = integer(0), se = list())
  )),
  facet_ncol = list(),
  facet_by = character(0)
)
empty_object_normalized <- agent_normalize_figure_spec(
  empty_object_spec, fd, pd, mat, character(), character()
)
ok(
  is.null(empty_object_normalized$facet_ncol) &&
    is.null(empty_object_normalized$facet_by) &&
    ut_cmp_identical(empty_object_normalized$layers[[1]]$params$alpha, 0.85) &&
    ut_cmp_identical(empty_object_normalized$layers[[1]]$params$size, 1.8) &&
    ut_cmp_identical(empty_object_normalized$layers[[1]]$params$se, TRUE),
  "empty-object sentinels fall back to param defaults (live regression)"
)
ok(
  ut_cmp_identical(sentinel_normalized$features, paste0("F", 1:3)) &&
    ut_cmp_identical(sentinel_normalized$samples, rownames(pd)),
  "sentinel feature/sample arrays resolve to selection/all samples"
)
ok(
  ut_cmp_identical(normalized$features, paste0("F", 1:3)),
  "explicit figure feature IDs are retained"
)
ok(
  ut_cmp_identical(normalized$layers[[1]]$params$alpha, 0.65),
  "validated layer parameters are retained"
)

figure_data <- agent_build_figure_data(normalized, fd, pd, mat)
ok(ut_cmp_identical(nrow(figure_data), 30L), "expression figure data has one row per feature/sample")
ok(
  ut_cmp_identical(
    colnames(figure_data),
    c(
      "__feature_id__", "__sample_id__", "__expression__",
      "feature__score", "feature__category", "sample__group"
    )
  ),
  "expression figure data joins feature and sample annotations"
)

figure_plot <- agent_build_figure_plot(figure_data, normalized)
ok(inherits(figure_plot, "ggplot"), "validated specification builds a ggplot")
rendered <- agent_render_figure(
  figure_plot,
  directory = file.path(tempdir(), paste0("agent-figures-", Sys.getpid())),
  figure_id = "fig_unit"
)
ok(
  ut_cmp_identical(startsWith(rendered$preview_uri, "data:image/png;base64,"), TRUE),
  "figure preview is encoded as PNG"
)
ok(
  ut_cmp_identical(startsWith(rendered$full_uri, "data:image/png;base64,"), TRUE),
  "high-resolution figure is encoded as PNG"
)
ok(
  ut_cmp_identical(rendered$full_bytes > rendered$preview_bytes, TRUE),
  "full-resolution figure is larger than its preview"
)
html <- as.character(agent_figure_html(rendered, normalized, "fig_unit"))
ok(
  ut_cmp_identical(grepl("<img", html, fixed = TRUE), TRUE),
  "figure card embeds an image directly"
)
ok(
  ut_cmp_identical(grepl("High-res PNG", html, fixed = TRUE), TRUE),
  "figure card includes a high-resolution download control"
)

ok(
  ut_cmp_error(
    agent_normalize_figure_spec(
      list(layers = list(list(geom = "eval", x = "score", y = "score"))),
      fd, pd, mat, character(), character()
    ),
    "Unsupported figure geom: eval"
  ),
  "arbitrary R geoms are rejected"
)
unknown_param_error <- tryCatch(
  agent_normalize_figure_spec(
    list(
      data_source = "feature_annotation",
      layers = list(list(geom = "point", x = "score", y = "score", params = list(path = "/tmp/x")))
    ),
    fd, pd, mat, character(), character()
  ),
  error = function(e) conditionMessage(e)
)
ok(
  ut_cmp_identical(
    grepl("Unknown figure layer parameter", unknown_param_error, fixed = TRUE),
    TRUE
  ),
  "arbitrary layer parameters are rejected"
)
transform_spec <- agent_normalize_figure_spec(
  list(
    data_source = "feature_annotation",
    layers = list(list(geom = "point", x = "category", y = "score")),
    x_transform = "log10"
  ),
  fd, pd, mat, character(), character()
)
transform_error <- tryCatch(
  {
    transform_data <- agent_build_figure_data(transform_spec, fd, pd, mat)
    agent_build_figure_plot(transform_data, transform_spec)
    NULL
  },
  error = function(e) conditionMessage(e)
)
ok(
  ut_cmp_identical(
    grepl("require a numeric x axis", transform_error, fixed = TRUE),
    TRUE
  ),
  "invalid numeric transformations are rejected"
)

# layer-count cap: 12 allowed, 13 rejected with guidance
mk_layers <- function(n)
  lapply(seq_len(n), function(i) list(geom = "point", x = "score", y = "score"))

ok(
  ut_cmp_identical(
    length(agent_normalize_figure_spec(
      list(layers = mk_layers(12)), fd, pd, mat, character(), character()
    )$layers),
    12L
  ),
  "figures accept up to 12 layers"
)
ok(
  ut_cmp_error(
    agent_normalize_figure_spec(
      list(layers = mk_layers(13)), fd, pd, mat, character(), character()
    ),
    "Figure supports at most 12 layers"
  ),
  "thirteen layers are rejected with merge guidance"
)

# ellmer converts type_array(type_object) tool arguments into tibbles: the
# chat-path spec has layers as an N-row/15-column data.frame (length() counts
# columns) and params as df-columns, with JSON nulls turned into NA. This
# shape previously failed the layer-count check (15 > cap) even though the
# model sent four layers.
if (requireNamespace("ellmer", quietly = TRUE)) {
  layer_type <- ellmer::type_object(
    geom = ellmer::type_string("geom"),
    x = ellmer::type_string("x", required = FALSE),
    y = ellmer::type_string("y", required = FALSE),
    color = ellmer::type_string("color", required = FALSE),
    fill = ellmer::type_string("fill", required = FALSE),
    params = ellmer::type_object(
      alpha = ellmer::type_number(required = FALSE),
      yintercept = ellmer::type_number(required = FALSE),
      .required = FALSE
    ),
    .required = FALSE
  )
  spec_type <- ellmer::type_object(
    layers = ellmer::type_array(layer_type, "layers", required = TRUE)
  )
  raw <- jsonlite::fromJSON(jsonlite::toJSON(list(
    layers = list(
      list(geom = "point", x = "score", y = "score",
           params = list(alpha = NULL)),
      list(geom = "hline", params = list(yintercept = 1.5))
    ),
    auto_unbox = FALSE
  ), auto_unbox = TRUE), simplifyVector = FALSE)
  converted <- ellmer:::convert_from_type(raw, spec_type)
  ok(
    ut_cmp_identical(is.data.frame(converted$layers), TRUE),
    "ellmer array-of-object conversion yields a data.frame (shape replicated)"
  )
  normalized_tibble <- agent_normalize_figure_spec(
    converted, fd, pd, mat, character(), character()
  )
  ok(
    ut_cmp_identical(length(normalized_tibble$layers), 2L),
    "tibble-shaped layers are coerced to row-lists and counted correctly"
  )
  ok(
    ut_cmp_identical(normalized_tibble$layers[[2]]$params$yintercept, 1.5),
    "numeric params pass through the tibble coercion"
  )
  ok(
    ut_cmp_identical(is.null(normalized_tibble$layers[[1]]$params$alpha), FALSE),
    "NA alpha from JSON null falls back to its default"
  )

  # WP3 live-path replication: the echo spec is serialized to JSON (tool
  # result), re-parsed and schema-converted by ellmer (tool arguments),
  # then normalized - and must reproduce the original normalized spec.
  echo_raw <- jsonlite::fromJSON(
    jsonlite::toJSON(agent_figure_spec_echo(normalized), auto_unbox = TRUE),
    simplifyVector = FALSE
  )
  echo_converted <- ellmer:::convert_from_type(
    list(layers = echo_raw$layers), spec_type
  )
  echo_renormalized <- agent_normalize_figure_spec(
    echo_converted, fd, pd, mat, character(), character()
  )
  ok(
    ut_cmp_identical(echo_renormalized$layers, normalized$layers),
    "echo spec survives the JSON -> ellmer -> normalize round-trip"
  )
} else {
  ok(TRUE, "ellmer conversion-shape test skipped: ellmer unavailable")
}

# ---- WP3: normalized specs are re-submittable (round-trip) -------------
renormalized <- agent_normalize_figure_spec(
  normalized, fd, pd, mat, character(), character()
)
ok(
  ut_cmp_identical(renormalized, normalized),
  "normalizing a normalized spec is the identity (WP3 round-trip)"
)

# the echo shape (what tool results carry) uses flat layer aesthetics and
# re-normalizes to the identical normalized spec
echo <- agent_figure_spec_echo(normalized)
ok(
  ut_cmp_identical("mappings" %in% names(echo$layers[[1]]), FALSE) &&
    ut_cmp_identical(echo$layers[[1]]$y, "__expression__"),
  "echo shape carries flat aesthetics, not mappings keys"
)
ok(
  ut_cmp_identical(
    agent_normalize_figure_spec(echo, fd, pd, mat, character(), character()),
    normalized
  ),
  "echo-shaped specs re-normalize to the identical normalized spec"
)

# the mappings-keyed layer shape (what tool results and the registry
# carry) is accepted and expanded into the documented flat aesthetics
mappings_shape <- agent_normalize_figure_spec(
  list(
    data_source = "feature_annotation",
    layers = list(list(
      geom = "point",
      mappings = list(x = "score", y = "score", color = "category")
    ))
  ),
  fd, pd, mat, character(), character()
)
ok(
  ut_cmp_identical(
    mappings_shape$layers[[1]]$mappings,
    list(x = "score", y = "score", color = "category")
  ),
  "mappings-keyed layers normalize to the same mappings"
)
ok(
  ut_cmp_identical(
    agent_normalize_figure_spec(
      list(layers = list(list(
        geom = "point",
        x = "score", y = "score",
        mappings = list(x = "category")
      ))),
      fd, pd, mat, character(), character()
    )$layers[[1]]$mappings$x,
    "score"
  ),
  "flat aesthetics win over expanded mappings"
)
ok(
  ut_cmp_error(
    agent_normalize_figure_spec(
      list(layers = list(list(geom = "point", mappings = list(zz = "a")))),
      fd, pd, mat, character(), character()
    ),
    "Unknown figure layer mapping"
  ),
  "unknown mappings keys are rejected"
)

# expression figures default to the FULL sample set; datasets larger than
# the ad-hoc 200-sample cap must still round-trip the echoed spec
pd_big <- data.frame(
  group = rep(c("WT", "KO"), 125),
  row.names = paste0("S", 1:250),
  check.names = FALSE
)
mat_big <- matrix(
  rnorm(10 * 250), nrow = 10,
  dimnames = list(rownames(fd), rownames(pd_big))
)
big <- agent_normalize_figure_spec(
  list(
    data_source = "expression",
    features = "F1",
    layers = list(list(geom = "boxplot", x = "sample__group", y = "__expression__"))
  ),
  fd, pd_big, mat_big, character(), character()
)
ok(
  ut_cmp_identical(length(big$samples), 250L),
  "expression default fills the full sample set beyond the 200 cap"
)
ok(
  ut_cmp_identical(
    agent_normalize_figure_spec(big, fd, pd_big, mat_big, character(), character()),
    big
  ) &&
    ut_cmp_identical(
      agent_normalize_figure_spec(
        agent_figure_spec_echo(big), fd, pd_big, mat_big,
        character(), character()
      ),
      big
    ),
  "full-set sample specs re-normalize identically (WP3 round-trip)"
)
ok(
  ut_cmp_error(
    agent_normalize_figure_spec(
      list(
        data_source = "expression",
        features = "F1",
        samples = paste0("S", 1:201),
        layers = list(list(geom = "boxplot", x = "sample__group", y = "__expression__"))
      ),
      fd, pd_big, mat_big, character(), character()
    ),
    "Figure can plot at most 200 samples"
  ),
  "hand-written sample subsets keep the 200 cap"
)

# ---- WP6: figure templates -----------------------------------------

agent_figure_templates <- omicsViewer:::agent_figure_templates
agent_figure_template_spec <- omicsViewer:::agent_figure_template_spec

ok(
  ut_cmp_identical(
    names(agent_figure_grammar()$templates),
    c("volcano", "scatter", "boxplot", "histogram")
  ) &&
    all(vapply(agent_figure_templates(), function(t)
      is.character(t$description) && is.list(t$arguments), logical(1))),
  "figure grammar advertises the four WP6 templates with argument help"
)

fd2 <- data.frame(
  logFC = c(-2, -1, 0, 1, 2, 3, -0.5, 0.5, 1.5, -1.5),
  logFdr = c(5, 0.1, 0, 4, 6, 0.2, 0, 3, 0, 1),
  category = rep(c("kinase", "phosphatase"), 5),
  shared = 1:10,
  row.names = paste0("F", 1:10),
  check.names = FALSE
)
pd2 <- data.frame(
  group = rep(c("WT", "KO"), 5),
  score2 = (1:10) * 1.5,
  shared = 1:10,
  row.names = paste0("S", 1:10),
  check.names = FALSE
)
pd2$shared <- as.character(pd2$shared)
mat2 <- matrix(rnorm(100), nrow = 10,
               dimnames = list(rownames(fd2), rownames(pd2)))

# volcano: full validated expansion with defaults
volcano <- agent_figure_template_spec(
  template = "volcano", x = "logFC", y = "logFdr",
  feature_data = fd2, sample_data = pd2, expression = mat2
)
ok(
  ut_cmp_identical(volcano$template, "volcano") &&
    ut_cmp_identical(volcano$spec$data_source, "feature_annotation") &&
    ut_cmp_identical(volcano$spec$layers[[1]]$geom, "point") &&
    ut_cmp_identical(volcano$spec$layers[[1]]$mappings$x, "logFC") &&
    ut_cmp_identical(volcano$spec$layers[[2]]$geom, "vline") &&
    ut_cmp_identical(volcano$spec$layers[[2]]$params$xintercept, 0),
  "volcano template expands to point + zero fold-change reference line"
)
ok(
  grepl("logFC", volcano$spec$labels$title, fixed = TRUE) &&
    ut_cmp_identical(volcano$spec$labels$x, "logFC") &&
    is.null(volcano$spec$features),
  "volcano template derives default labels and plots all features"
)

# volcano label_top_n: rows reordered by significance (y, highest first)
# so the capped label layer marks the most significant features
volcano_labeled <- agent_figure_template_spec(
  template = "volcano", x = "logFC", y = "logFdr", color = "category",
  label_top_n = 3L, feature_data = fd2, sample_data = pd2, expression = mat2
)
ok(
  ut_cmp_identical(
    volcano_labeled$spec$features,
    c("F5", "F1", "F4", "F8", "F10", "F6", "F2", "F3", "F7", "F9")
  ),
  "volcano label_top_n ranks features by significance, highest first"
)
ok(
  ut_cmp_identical(length(volcano_labeled$spec$features), 10L) &&
    ut_cmp_identical(volcano_labeled$spec$layers[[3]]$geom, "label") &&
    ut_cmp_identical(volcano_labeled$spec$layers[[3]]$mappings$label, "__feature_id__") &&
    ut_cmp_identical(volcano_labeled$spec$layers[[3]]$params$max_labels, 3L) &&
    ut_cmp_identical(volcano_labeled$spec$layers[[1]]$mappings$color, "category"),
  "volcano label layer is capped at the requested count and carries color"
)

ok(
  ut_cmp_error(
    agent_figure_template_spec(
      template = "volcano", x = "logFC", y = "category",
      feature_data = fd2, sample_data = pd2
    ),
    "must be numeric"
  ),
  "volcano rejects non-numeric significance columns"
)
ok(
  ut_cmp_error(
    agent_figure_template_spec(
      template = "histogram", x = "logFdr", label_top_n = 2L,
      feature_data = fd2, sample_data = pd2
    ),
    "label_top_n is supported by the volcano and scatter templates"
  ),
  "label_top_n on boxplot/histogram is rejected honestly"
)

# scatter: space resolution, ambiguity handling
scatter_f <- agent_figure_template_spec(
  template = "scatter", x = "score2", y = "score2",
  feature_data = fd2, sample_data = pd2
)
ok(
  ut_cmp_identical(scatter_f$spec$data_source, "sample_annotation") &&
    ut_cmp_identical(scatter_f$spec$layers[[1]]$mappings$y, "score2"),
  "scatter resolves the annotation space from the columns"
)
ok(
  ut_cmp_error(
    agent_figure_template_spec(
      template = "scatter", x = "shared", y = "score2",
      feature_data = fd2, sample_data = pd2
    ),
    "exists in both the feature and the sample"
  ),
  "ambiguous scatter column without space is rejected with guidance"
)
scatter_disambiguated <- agent_figure_template_spec(
  template = "scatter", x = "shared", y = "shared", space = "feature",
  feature_data = fd2, sample_data = pd2
)
ok(
  ut_cmp_identical(scatter_disambiguated$spec$data_source, "feature_annotation"),
  "explicit space disambiguates shared column names"
)
ok(
  ut_cmp_error(
    agent_figure_template_spec(
      template = "scatter", x = "score", y = "score",
      feature_data = fd2, sample_data = pd2
    ),
    "Closest matches"
  ),
  "unknown template columns surface closest-match suggestions"
)

# boxplot: annotation mode (y given) and expression mode (y omitted)
box_annotation <- agent_figure_template_spec(
  template = "boxplot", x = "category", y = "logFC",
  feature_data = fd2, sample_data = pd2
)
ok(
  ut_cmp_identical(box_annotation$spec$data_source, "feature_annotation") &&
    ut_cmp_identical(box_annotation$spec$layers[[1]]$geom, "boxplot") &&
    ut_cmp_identical(box_annotation$spec$layers[[1]]$mappings$fill, "category"),
  "annotation-mode boxplot groups a numeric column with fill defaulting to x"
)
box_expression <- agent_figure_template_spec(
  template = "boxplot", x = "group",
  feature_data = fd2, sample_data = pd2, expression = mat2,
  selected_features = c("F1", "F2", "F3")
)
ok(
  ut_cmp_identical(box_expression$spec$data_source, "expression") &&
    ut_cmp_identical(box_expression$spec$layers[[1]]$mappings$x, "sample__group") &&
    ut_cmp_identical(box_expression$spec$layers[[1]]$mappings$y, "__expression__") &&
    ut_cmp_identical(box_expression$spec$features, c("F1", "F2", "F3")),
  "expression-mode boxplot namespaces the sample grouping and uses the selection"
)
ok(
  ut_cmp_error(
    agent_figure_template_spec(
      template = "boxplot", x = "group",
      feature_data = fd2, sample_data = pd2, expression = mat2,
      selected_features = character()
    ),
    "requires selected features"
  ),
  "expression-mode boxplot needs a feature selection"
)
ok(
  ut_cmp_error(
    agent_figure_template_spec(
      template = "boxplot", x = "logFC",
      feature_data = fd2, sample_data = pd2
    ),
    "sample annotation"
  ),
  "expression-mode boxplot rejects feature-space grouping columns"
)

# histogram
histogram <- agent_figure_template_spec(
  template = "histogram", x = "logFdr", color = "category",
  feature_data = fd2, sample_data = pd2
)
ok(
  ut_cmp_identical(histogram$spec$data_source, "feature_annotation") &&
    ut_cmp_identical(histogram$spec$layers[[1]]$geom, "histogram") &&
    ut_cmp_identical(histogram$spec$layers[[1]]$mappings$fill, "category") &&
    grepl("logFdr", histogram$spec$labels$title, fixed = TRUE),
  "histogram template maps the optional color argument to fill"
)

# template validation + sentinel hardening (glm flash omitted-optional class)
ok(
  ut_cmp_error(
    agent_figure_template_spec(
      template = "heatmap", x = "logFC", y = "logFdr",
      feature_data = fd2, sample_data = pd2
    ),
    "Available templates"
  ),
  "unknown templates are rejected with the available list"
)
ok(
  ut_cmp_error(
    agent_figure_template_spec(feature_data = fd2, sample_data = pd2),
    "Figure template is required"
  ),
  "missing template is rejected"
)
sentinel_template <- agent_figure_template_spec(
  template = "volcano", x = "logFC", y = "logFdr",
  color = "null", label_top_n = "null", title = "{}", space = "[]",
  feature_data = fd2, sample_data = pd2
)
ok(
  is.null(sentinel_template$spec$layers[[1]]$mappings$color) &&
    length(sentinel_template$spec$layers) == 2L &&
    grepl("logFC", sentinel_template$spec$labels$title, fixed = TRUE),
  "sentinel template optionals behave exactly like omitted values"
)
empty_object_template <- agent_figure_template_spec(
  template = "volcano", x = "logFC", y = "logFdr",
  label_top_n = list(),
  feature_data = fd2, sample_data = pd2
)
ok(
  length(empty_object_template$spec$layers) == 2L,
  "empty-object label_top_n behaves exactly like an omitted value"
)
ok(
  ut_cmp_identical(
    agent_normalize_figure_spec(
      volcano$spec, fd2, pd2, mat2, character(), character()
    ),
    volcano$spec
  ),
  "template expansion result is normalized and re-normalizes identically"
)

# expanded template specs build and render through the generic path
volcano_data <- agent_build_figure_data(
  volcano_labeled$spec, fd2, pd2, mat2
)
ok(
  ut_cmp_identical(volcano_data[["__feature_id__"]][1], "F5") &&
    ut_cmp_identical(nrow(volcano_data), 10L),
  "labeled volcano data keeps all rows with the most significant first"
)
ok(
  inherits(agent_build_figure_plot(volcano_data, volcano_labeled$spec), "ggplot") &&
    inherits(
      agent_build_figure_plot(
        agent_build_figure_data(histogram$spec, fd2, pd2, mat2), histogram$spec
      ),
      "ggplot"
    ),
  "expanded template specs build ggplots through the generic path"
)
