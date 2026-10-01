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
# literal "null"/"{}"/"[]" strings): the WP15 seam
# (agent_args_sanitize at tool entry) turns every sentinel into an
# absent value BEFORE validation, so the omitted-value defaults apply
# instead of erroring or rendering literally
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
  omicsViewer:::agent_args_sanitize(sentinel_spec), fd, pd, mat,
  paste0("F", 1:3), character()
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
# twice with "must be an integer between 1 and 6"): the seam collapses
# empty objects/arrays to NULL, so length-0 values behave exactly like
# omitted values in every param helper
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
  omicsViewer:::agent_args_sanitize(empty_object_spec), fd, pd, mat,
  character(), character()
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

# WP15: tool arguments travel on the canonical JSON-document transport
# (agent_tool/agent_args_sanitize at the seam; ellmer convert = FALSE), so
# the chat-path spec is the parsed JSON document itself -- named lists,
# atomic vectors, NULL for absent -- and never an ellmer tibble. Replicate
# the wire path: serialize, re-parse as the provider payload, sanitize at
# the seam, normalize.
wire_spec <- jsonlite::fromJSON(jsonlite::toJSON(list(
  layers = list(
    list(geom = "point", x = "score", y = "score",
         params = list(alpha = NULL)),
    list(geom = "hline", params = list(yintercept = 1.5))
  )
), auto_unbox = TRUE), simplifyVector = FALSE)
sanitized_spec <- omicsViewer:::agent_args_sanitize(wire_spec)
ok(
  ut_cmp_identical(is.data.frame(sanitized_spec$layers), FALSE) &&
    ut_cmp_identical(length(sanitized_spec$layers), 2L),
  "the seam delivers layers as a plain list, not an ellmer tibble"
)
normalized_wire <- agent_normalize_figure_spec(
  sanitized_spec, fd, pd, mat, character(), character()
)
ok(
  ut_cmp_identical(length(normalized_wire$layers), 2L),
  "canonical layers are counted correctly"
)
ok(
  ut_cmp_identical(normalized_wire$layers[[2]]$params$yintercept, 1.5),
  "numeric params pass through the canonical transport"
)
ok(
  ut_cmp_identical(is.null(normalized_wire$layers[[1]]$params$alpha), FALSE),
  "explicit JSON null params fall back to their default"
)

# 2026-09-28 post-mortem: a decorative hline without yintercept is the
# conventional no-change reference line; deriving the default (with a
# surfaced warning) beats rejecting the figure - a hard stop cost a full
# provider round trip just to drop one decorative layer.
no_intercept <- agent_normalize_figure_spec(list(
  layers = list(
    list(geom = "hline", params = list(color = "#555555")),
    list(geom = "vline", params = list())
  )
), fd, pd, mat, character(), character())
ok(
  ut_cmp_identical(no_intercept$layers[[1]]$params$yintercept, 0) &&
    ut_cmp_identical(no_intercept$layers[[2]]$params$xintercept, 0),
  "hline/vline without intercept default to the conventional 0 reference"
)
warned <- character()
withCallingHandlers(
  normalized_warn <- agent_normalize_figure_spec(list(
    layers = list(list(geom = "hline", params = list()))
  ), fd, pd, mat, character(), character()),
  warning = function(w) {
    warned <<- c(warned, conditionMessage(w)); invokeRestart("muffleWarning")
  }
)
ok(
  length(warned) == 1L && grepl("yintercept", warned),
  "the derived default surfaces a warning the model can correct on revision"
)

# WP3 live-path replication: the echo spec is serialized to JSON (tool
# result), re-parsed as the provider payload (tool arguments), sanitized
# at the seam, then normalized - and must reproduce the original.
echo_raw <- jsonlite::fromJSON(
  jsonlite::toJSON(agent_figure_spec_echo(normalized), auto_unbox = TRUE),
  simplifyVector = FALSE
)
echo_renormalized <- agent_normalize_figure_spec(
  omicsViewer:::agent_args_sanitize(echo_raw), fd, pd, mat,
  character(), character()
)
ok(
  ut_cmp_identical(echo_renormalized$layers, normalized$layers),
  "echo spec survives the JSON -> seam -> normalize round-trip"
)

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
# todo 3.1: the echo of a default (full-set) sample resolution is ELIDED
# (bounded independent of dataset size); re-normalization re-derives it
echo_default <- agent_figure_spec_echo(
  big, all_features = rownames(fd), all_samples = rownames(pd_big))
ok(
  is.null(echo_default$samples) &&
    nchar(jsonlite::toJSON(echo_default, auto_unbox = TRUE), type = "bytes") < 4000L,
  "full-set id arrays are elided from the echo (bounded tool results)"
)
# a small explicit subset rides along verbatim and still round-trips
small <- agent_normalize_figure_spec(
  list(data_source = "expression", features = "F1",
       samples = c("S3", "S1"),
       layers = list(list(geom = "boxplot", x = "sample__group", y = "__expression__"))),
  fd, pd_big, mat_big, character(), character())
echo_small <- agent_figure_spec_echo(
  small, all_features = rownames(fd), all_samples = rownames(pd_big))
ok(
  ut_cmp_identical(echo_small$samples, c("S3", "S1")) &&
    ut_cmp_identical(
      agent_normalize_figure_spec(echo_small, fd, pd_big, mat_big,
                                  character(), character()),
      small
    ),
  "small explicit id subsets ride the echo verbatim (WP3 round-trip)"
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
    c("volcano", "scatter", "boxplot", "histogram", "shared_arguments")
  ) &&
    all(vapply(agent_figure_templates()[c("volcano", "scatter", "boxplot", "histogram")], function(t)
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

# volcano label_top_n: rows ranked at RENDER time (y, highest first)
# so the capped label layer marks the most significant features while
# spec$features stays NULL (todo 3.1: the echo never grows with the set)
volcano_labeled <- agent_figure_template_spec(
  template = "volcano", x = "logFC", y = "logFdr", color = "category",
  label_top_n = 3L, feature_data = fd2, sample_data = pd2, expression = mat2
)
ok(
  is.null(volcano_labeled$spec$features),
  "volcano label_top_n no longer resolves the reordered feature vector into the spec"
)
volcano_data <- omicsViewer:::agent_build_figure_data(
  volcano_labeled$spec, fd2, pd2, mat2)
volcano_plot <- agent_build_figure_plot(volcano_data, volcano_labeled$spec)
volcano_labels <- volcano_plot$layers[[3]]$data[["__feature_id__"]]
ok(
  ut_cmp_identical(volcano_labels,
                   c("F5", "F1", "F4")) &&
    ut_cmp_identical(length(volcano_labels), 3L),
  "volcano label layer marks the top-n by significance at render time"
)
ok(
  ut_cmp_identical(volcano_labeled$spec$layers[[3]]$geom, "label") &&
    ut_cmp_identical(volcano_labeled$spec$layers[[3]]$mappings$label, "__feature_id__") &&
    ut_cmp_identical(volcano_labeled$spec$layers[[3]]$params$max_labels, 3L) &&
    ut_cmp_identical(volcano_labeled$spec$layers[[3]]$params$order_by,
                     list(column = "logFdr", decreasing = TRUE)) &&
    ut_cmp_identical(volcano_labeled$spec$layers[[1]]$mappings$color, "category"),
  "volcano label layer carries color and a render-time order_by"
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
    "requires plotted features"
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
sentinel_template <- do.call(
  agent_figure_template_spec,
  c(
    omicsViewer:::agent_args_sanitize(list(
      template = "volcano", x = "logFC", y = "logFdr",
      color = "null", label_top_n = "null", title = "{}", space = "[]"
    )),
    list(feature_data = fd2, sample_data = pd2)
  )
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
  ut_cmp_identical(volcano_data[["__feature_id__"]][1], "F1") &&
    ut_cmp_identical(nrow(volcano_data), 10L),
  "labeled volcano data keeps the natural row order (ranking is render-time)"
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

# ---- WP6b: template features/samples subsets + id-param normalization ----

# explicit features subset on a feature-space template
volcano_subset <- agent_figure_template_spec(
  template = "volcano", x = "logFC", y = "logFdr", features = c("F3", "F1", "F7"),
  feature_data = fd2, sample_data = pd2, expression = mat2
)
ok(
  ut_cmp_identical(sort(volcano_subset$spec$features), c("F1", "F3", "F7")) &&
    ut_cmp_identical(volcano_subset$spec$data_source, "feature_annotation"),
  "WP6b: volcano accepts an explicit features subset"
)

# features subset is ranked within itself when label_top_n is requested
# (todo 3.1: the subset rides the spec unchanged; ranking is render-time)
volcano_subset_labeled <- agent_figure_template_spec(
  template = "volcano", x = "logFC", y = "logFdr", features = c("F3", "F1", "F7"),
  label_top_n = 2L, feature_data = fd2, sample_data = pd2, expression = mat2
)
subset_data <- omicsViewer:::agent_build_figure_data(
  volcano_subset_labeled$spec, fd2, pd2, mat2)
subset_plot <- agent_build_figure_plot(subset_data, volcano_subset_labeled$spec)
subset_labels <- subset_plot$layers[[3]]$data[["__feature_id__"]]
ok(
  ut_cmp_identical(volcano_subset_labeled$spec$features, c("F3", "F1", "F7")) &&
    ut_cmp_identical(subset_labels, c("F1", "F3")) &&
    ut_cmp_identical(volcano_subset_labeled$spec$layers[[3]]$params$max_labels, 2L),
  "WP6b: label_top_n ranks the explicit subset at render time (y: F1=5 leads; F3/F7 tie at 0)"
)

# unknown IDs are rejected (validated downstream by the normalizer)
ok(
  ut_cmp_error(
    agent_figure_template_spec(
      template = "volcano", x = "logFC", y = "logFdr", features = c("F1", "nope"),
      feature_data = fd2, sample_data = pd2
    ),
    "Unknown figure feature"
  ),
  "WP6b: template features are validated against the dataset rownames"
)

# boxplot expression mode works from an explicit subset with NO selection
boxplot_subset <- agent_figure_template_spec(
  template = "boxplot", x = "group", features = c("F2", "F3"),
  feature_data = fd2, sample_data = pd2, expression = mat2,
  selected_features = character()
)
ok(
  ut_cmp_identical(boxplot_subset$spec$data_source, "expression") &&
    ut_cmp_identical(sort(boxplot_subset$spec$features), c("F2", "F3")) &&
    ut_cmp_identical(sort(boxplot_subset$spec$samples), sort(rownames(pd2))),
  "WP6b: expression boxplot from explicit features without a selection"
)

# explicit samples subset rides along on expression figures
boxplot_samples <- agent_figure_template_spec(
  template = "boxplot", x = "group", features = c("F2", "F3"),
  samples = c("S1", "S3", "S5", "S7", "S9"),
  feature_data = fd2, sample_data = pd2, expression = mat2
)
ok(
  ut_cmp_identical(sort(boxplot_samples$spec$samples),
                   c("S1", "S3", "S5", "S7", "S9")),
  "WP6b: expression boxplot honors an explicit samples subset"
)

# sentinel shapes are treated as omitted across array providers: the
# seam collapses sentinel strings, empty arrays, and all-null arrays to
# NULL before the id-array helper runs
agent_figure_ids_param <- omicsViewer:::.agent_figure_ids_param
sanitize <- omicsViewer:::agent_args_sanitize
ok(
  is.null(agent_figure_ids_param(NULL, "features")) &&
    is.null(agent_figure_ids_param(character(), "features")) &&
    is.null(agent_figure_ids_param(NA_character_, "features")) &&
    is.null(agent_figure_ids_param(sanitize(list()), "features")) &&
    ut_cmp_identical(
      agent_figure_ids_param(sanitize(list("null", "[]")), "features"),
      c("null", "[]")
    ) &&
    ut_cmp_identical(
      agent_figure_ids_param(sanitize(list("F1", "F2")), "features"), c("F1", "F2")
    ) &&
    ut_cmp_identical(
      agent_figure_ids_param(c(" F1 ", "", "F2"), "features"), c("F1", "F2")
    ),
  "WP6b: id-array params keep sentinel-equal elements as data (3.3)"
)

# subsets only attach to data sources that actually plot them
scatter_sample_space <- agent_figure_template_spec(
  template = "scatter", x = "score2", y = "group", features = c("F1"),
  feature_data = fd2, sample_data = pd2
)
ok(
  ut_cmp_identical(scatter_sample_space$spec$data_source, "sample_annotation") &&
    is.null(scatter_sample_space$spec$features),
  "WP6b: feature ids do not attach to sample-space figures"
)

# ---- widened grammar: structured layer filters ---------------------------
# Layer `filter` predicates: canonicalization, x/y shorthand, caps, sentinel
# tolerance, and the closed-operator eval interpreter.
agent_where_normalize <- omicsViewer:::agent_where_normalize
agent_where_eval <- omicsViewer:::agent_where_eval

wnum <- agent_where_normalize(list(column = "score", op = ">", value = 5))
ok(
  ut_cmp_identical(wnum, list(column = "score", op = ">", value = 5)),
  "where leaf numeric op canonicalizes"
)
# numeric strings (provider artifacts) coerce like numbers
ok(
  ut_cmp_identical(
    agent_where_normalize(list(column = "score", op = "<=", value = "2.5")),
    list(column = "score", op = "<=", value = 2.5)
  ),
  "where leaf coerces numeric-string values"
)
# 'x'/'y' resolve to the layer's own mapping at normalize time
wxy <- agent_where_normalize(
  list(column = "x", op = "between", min = 1, max = 3),
  mappings = list(x = "score", y = "category")
)
ok(
  ut_cmp_identical(wxy$column, "score"),
  "where x/y shorthand resolves to the layer mapping"
)
ut_fails <- function(expr, pattern) {
  tryCatch({ force(expr); FALSE },
           error = function(e) grepl(pattern, conditionMessage(e), fixed = TRUE))
}
ok(
  ut_fails(
    agent_where_normalize(list(column = "y", op = ">", value = 1), mappings = list(x = "score")),
    "does not map"
  ),
  "where x/y shorthand without that mapping errors clearly"
)
ok(
  ut_fails(
    agent_where_normalize(list(column = "score", op = "between", min = 5, max = 1)),
    "min <= max"
  ) &&
    ut_fails(
      agent_where_normalize(list(column = "score", op = ">", value = 1, text = "a")),
      "does not accept"
    ) &&
    ut_fails(
      agent_where_normalize(list(column = "score", op = "==")),
      "exactly one of value (number) or text (string)"
    ) &&
    ut_fails(
      agent_where_normalize(list(column = "score", op = "is_na", value = 1)),
      "takes no comparison value"
    ),
  "where operator field contracts are enforced"
)
ok(
  ut_fails(
    agent_where_normalize(list(column = "score", op = "much_greater", value = 1)),
    "Unsupported figure filter operator"
  ),
  "unknown where operator is rejected with a suggestion"
)
wall <- agent_where_normalize(list(all = list(
  list(column = "score", op = "abs>=", value = 2),
  list(any = list(
    list(column = "category", op = "==", text = "kinase"),
    list(column = "__feature_id__", op = "starts_with", text = "F1")
  ))
)))
ok(
  ut_cmp_identical(names(wall), "all") &&
    ut_cmp_identical(wall$all[[2]]$any[[1]]$text, "kinase"),
  "where combinators canonicalize recursively"
)
ok(
  ut_fails(
    agent_where_normalize(list(all = list(list(column = "score", op = ">", value = 1)),
                               any = list(list(column = "score", op = "<", value = 9)))),
    "either 'all' or 'any'"
  ) &&
    is.null(agent_where_normalize(list(all = list(list())))) &&
    ut_fails(
      agent_where_normalize(list(all = 5)),
      "non-empty array"
    ),
  "where combinator misuse errors clearly"
)
# caps: depth <= 3, leaves <= 8
deep <- list(column = "score", op = ">", value = 1)
for (i in 1:5) deep <- list(all = list(deep))
ok(
  ut_fails(
    omicsViewer:::.agent_where_validate_caps(deep),
    "nests at most 3 levels"
  ),
  "where depth cap is enforced"
)
many <- list(all = rep(list(list(column = "score", op = ">", value = 1)), 9))
ok(
  ut_fails(
    omicsViewer:::.agent_where_validate_caps(many),
    "at most 8 conditions"
  ),
  "where leaf cap is enforced"
)
# sentinels behave like omitted values everywhere: the seam turns
# sentinel strings and empty/all-null structures into NULL before the
# filter canonicalizer runs
ok(
  is.null(agent_where_normalize(sanitize("null"))) &&
    is.null(agent_where_normalize(sanitize(list()))) &&
    is.null(agent_where_normalize(NA)),
  "where sentinels are treated as omitted"
)
ok(
  ut_cmp_identical(
    agent_where_normalize(list(column = "category", op = "in",
                               values = sanitize(list("kinase", "null", "", NA))))[["values"]],
    c("kinase", "null")
  ),
  "where in-values drop NA/empty entries; sentinel-equal strings are data (3.3)"
)
# idempotency: canonical filters re-normalize identically
ok(
  ut_cmp_identical(agent_where_normalize(wnum), wnum) &&
    ut_cmp_identical(agent_where_normalize(wall), wall),
  "where normalization is idempotent"
)

# ---- where eval against built plotting data ------------------------------
edata <- agent_build_figure_data(
  agent_normalize_figure_spec(
    list(data_source = "feature_annotation",
         layers = list(list(geom = "point", x = "score", y = "score"))),
    fd, pd, mat, character(), character()
  ),
  fd, pd, mat
)
ok(
  ut_cmp_identical(
    as.integer(agent_where_eval(list(column = "score", op = ">", value = 7), edata)),
    c(rep(0L, 7), rep(1L, 3L))
  ),
  "where eval applies numeric comparisons rowwise"
)
ok(
  ut_cmp_identical(
    as.integer(agent_where_eval(
      list(column = "score", op = "between", min = 2, max = 4), edata)),
    c(0L, 1L, 1L, 1L, rep(0L, 6L))
  ) &&
    ut_cmp_identical(
      as.integer(agent_where_eval(
        list(column = "category", op = "==", text = "kinase"), edata)),
      rep(c(1L, 0L), 5L)
    ),
  "where eval handles between and text equality"
)
ok(
  ut_cmp_identical(
    as.integer(agent_where_eval(
      list(column = "score", op = "in", values = c("1", "10")), edata)),
    c(1L, rep(0L, 8L), 1L)
  ),
  "where eval coerces number-strings on numeric columns"
)
ok(
  ut_cmp_identical(
    as.integer(agent_where_eval(
      list(column = "category", op = "not_in", values = "kinase"), edata)),
    rep(c(0L, 1L), 5L)
  ),
  "where eval applies not_in on text columns"
)
ok(
  ut_cmp_identical(
    as.integer(agent_where_eval(
      list(column = "score", op = "<", value = 5, not = TRUE), edata)),
    c(rep(0L, 4L), rep(1L, 6L))
  ),
  "where eval honors the not flag"
)
ok(
  ut_cmp_identical(
    agent_where_eval(list(column = "score", op = ">", value = 0), edata),
    rep(TRUE, 10L)
  ) &&
    isTRUE(all(
      agent_where_eval(
        list(all = list(
          list(column = "score", op = ">", value = 100),
          list(column = "score", op = "<", value = 0)
        )),
        edata
      ) == FALSE
    )),
  "where eval returns settled logicals for impossible matches"
)
# NA semantics: comparisons never match NA; is_na/not_na are explicit
fd_na <- fd
fd_na$score[c(2, 5)] <- NA
edata_na <- agent_build_figure_data(
  agent_normalize_figure_spec(
    list(data_source = "feature_annotation",
         layers = list(list(geom = "point", x = "score", y = "score"))),
    fd_na, pd, mat, character(), character()
  ),
  fd_na, pd, mat
)
ok(
  ut_cmp_identical(
    as.integer(agent_where_eval(list(column = "score", op = ">", value = -1), edata_na)),
    c(1L, 0L, 1L, 1L, 0L, rep(1L, 5L))
  ) &&
    ut_cmp_identical(
      as.integer(agent_where_eval(list(column = "score", op = "is_na"), edata_na)),
      c(0L, 1L, 0L, 0L, 1L, rep(0L, 5L))
    ) &&
    ut_cmp_identical(
      as.integer(agent_where_eval(
        list(column = "score", op = ">", value = -1, not = TRUE), edata_na)),
      rep(0L, 10L)
    ),
  "where NA rows never match comparisons, is_na and not are explicit"
)
ok(
  ut_fails(
    agent_where_eval(list(column = "scor", op = ">", value = 1), edata),
    "Closest matches"
  ) &&
    ut_fails(
      agent_where_eval(list(column = "category", op = ">", value = 1), edata),
      "requires a numeric column"
    ) &&
    ut_fails(
      agent_where_eval(list(column = "score", op = "starts_with", text = "F"), edata),
      "requires a text column"
    ),
  "where eval errors carry suggestions and type diagnostics"
)

# ---- widened spec: filters, constant colors, scale, theme options --------
wide_spec <- list(
  data_source = "feature_annotation",
  layers = list(
    list(geom = "point", x = "score", y = "score", fill = "category",
         params = list(color = "#B0B0B0")),
    list(geom = "point", x = "score", y = "score",
         filter = list(column = "x", op = ">", value = 7),
         params = list(color = "#B2182B", size = 2.6)),
    list(geom = "label", x = "score", y = "score", label = "__feature_id__",
         filter = list(column = "score", op = ">", value = 7),
         params = list(color = "#b2182b", size = 3, max_labels = 2))
  ),
  scale = list(fill = list(values = list(
    list(category = "kinase", color = "#2166ac"),
    list(category = "phosphatase", color = "#B2182B")
  ))),
  theme_options = list(legend_position = "bottom", rotate_x_labels = 30,
                       base_size = 12L, show_grid = FALSE)
)
wide_norm <- agent_normalize_figure_spec(
  wide_spec, fd, pd, mat, character(), character()
)
ok(
  ut_cmp_identical(wide_norm$layers[[1]]$params$color, "#b0b0b0") &&
    ut_cmp_identical(wide_norm$layers[[2]]$params$color, "#b2182b") &&
    ut_cmp_identical(wide_norm$layers[[2]]$filter$column, "score") &&
    is.null(wide_norm$layers[[1]]$filter),
  "constant colors lowercase and x-shorthand resolves inside normalize"
)
ok(
  ut_cmp_identical(
    wide_norm$scale$fill$values,
    c(kinase = "#2166ac", phosphatase = "#b2182b")
  ) &&
    ut_cmp_identical(wide_norm$theme_options$legend_position, "bottom") &&
    ut_cmp_identical(wide_norm$theme_options$rotate_x_labels, 30) &&
    ut_cmp_identical(wide_norm$theme_options$show_grid, FALSE),
  "scale values pairs and theme options canonicalize"
)
ok(
  ut_fails(
    agent_normalize_figure_spec(
      list(data_source = "feature_annotation",
           layers = list(list(geom = "point", x = "score", y = "score",
                              color = "score", params = list(color = "#ff0000")))),
      fd, pd, mat, character(), character()),
    "remove one of the two"
  ),
  "constant color conflicts with a mapped aesthetic loudly"
)
ok(
  ut_fails(
    agent_normalize_figure_spec(
      list(data_source = "feature_annotation",
           layers = list(list(geom = "point", x = "score", y = "score",
                              params = list(fill = "not-a-color")))),
      fd, pd, mat, character(), character()),
    "must be a hex color"
  ),
  "invalid hex colors are rejected with guidance"
)
ok(
  ut_fails(
    agent_normalize_figure_spec(
      list(data_source = "feature_annotation",
           layers = list(list(geom = "point", x = "score", y = "score")),
           scale = list(color = list(values = list(
             list(category = "kinase", color = "#000000"),
             list(category = "kinase", color = "#ffffff"))))),
      fd, pd, mat, character(), character()),
    "duplicate categories"
  ),
  "duplicate scale categories are rejected"
)
# sentinel sweep on the new spec fields (through the WP15 seam)
sentinel_wide <- agent_normalize_figure_spec(
  omicsViewer:::agent_args_sanitize(list(
    data_source = "feature_annotation",
    layers = list(list(geom = "point", x = "score", y = "score",
                       filter = "null")),
    scale = "{}", theme_options = "null"
  )),
  fd, pd, mat, character(), character()
)
ok(
  is.null(sentinel_wide$layers[[1]]$filter) &&
    is.null(sentinel_wide$scale) &&
    is.null(sentinel_wide$theme_options),
  "sentinel filter/scale/theme_options behave like omitted values"
)
# round-trip: the echo shape re-normalizes identically (WP3 contract)
wide_echo <- agent_figure_spec_echo(wide_norm)
wide_again <- agent_normalize_figure_spec(
  wide_echo, fd, pd, mat, character(), character()
)
ok(
  ut_cmp_identical(wide_again$layers[[2]]$filter, wide_norm$layers[[2]]$filter) &&
    ut_cmp_identical(wide_again$scale, wide_norm$scale) &&
    ut_cmp_identical(wide_again$theme_options, wide_norm$theme_options) &&
    ut_cmp_identical(wide_again$layers[[1]]$params, wide_norm$layers[[1]]$params),
  "widened spec survives the echo round-trip"
)

# ---- widened plots render ------------------------------------------------
wide_plot <- agent_build_figure_plot(
  agent_build_figure_data(wide_norm, fd, pd, mat), wide_norm
)
ok(
  inherits(wide_plot, "ggplot") && length(wide_plot$layers) == 3L,
  "filtered multi-layer highlight spec builds a ggplot"
)
# label layer composite: filter (score > 7 -> F8..F10) then max_labels 2
label_layer <- wide_plot$layers[[3]]
label_rows <- label_layer$data
ok(
  nrow(label_rows) == 2L &&
    all(label_rows[["__feature_id__"]] %in% c("F8", "F9", "F10")),
  "filter composites before the max_labels cap"
)
# filtered NON-text layers draw their subset explicitly (regression: the
# highlight point layer initially dropped the filter silently)
highlight_layer <- wide_plot$layers[[2]]
ok(
  nrow(highlight_layer$data) == 3L &&
    identical(highlight_layer$data[["__feature_id__"]], c("F8", "F9", "F10")),
  "filtered point layers draw only the matching rows"
)
ok(
  !is.data.frame(wide_plot$layers[[1]]$data),
  "unfiltered layers keep inheriting the plot data"
)
# scale fill override waits for a layer that maps fill
fill_spec <- agent_normalize_figure_spec(
  list(data_source = "feature_annotation",
       layers = list(list(geom = "point", x = "score", y = "score",
                          fill = "category")),
       scale = list(fill = list(values = list(
         list(category = "kinase", color = "#2166ac"),
         list(category = "phosphatase", color = "#B2182B")
       )))),
  fd, pd, mat, character(), character()
)
fill_plot <- agent_build_figure_plot(
  agent_build_figure_data(fill_spec, fd, pd, mat), fill_spec
)
manual_hit <- vapply(fill_plot$scales$scales, function(s) {
  is.function(s$palette) &&
    identical(sort(unique(s$palette(2))), sort(c("#2166ac", "#b2182b")))
}, logical(1))
ok(
  any(manual_hit),
  "discrete scale override installs a manual scale"
)
ok(
  ut_fails(
    agent_build_figure_plot(
      agent_build_figure_data(
        agent_normalize_figure_spec(
          list(data_source = "feature_annotation",
               layers = list(list(geom = "point", x = "score", y = "score",
                                  fill = "category")),
               scale = list(fill = list(limits = "nope"))),
          fd, pd, mat, character(), character()),
        fd, pd, mat),
      agent_normalize_figure_spec(
        list(data_source = "feature_annotation",
             layers = list(list(geom = "point", x = "score", y = "score",
                                fill = "category")),
             scale = list(fill = list(limits = "nope"))),
        fd, pd, mat, character(), character())
    ),
    "Unknown figure scale fill limits"
  ),
  "unknown discrete scale limits error at build time"
)
# continuous overrides: limits pair and diverging midpoint
cont_limits <- agent_normalize_figure_spec(
  list(data_source = "feature_annotation",
       layers = list(list(geom = "point", x = "score", y = "score", color = "score")),
       scale = list(color = list(limits = c("2", "8")))),
  fd, pd, mat, character(), character()
)
cont_plot <- agent_build_figure_plot(
  agent_build_figure_data(cont_limits, fd, pd, mat), cont_limits
)
cont_classes <- vapply(cont_plot$scales$scales, function(s) class(s)[[1]], character(1))
ok(
  any(grepl("ScaleContinuous", cont_classes, fixed = TRUE)),
  "continuous limits override installs a gradient scale"
)
mid_spec <- agent_normalize_figure_spec(
  list(data_source = "feature_annotation",
       layers = list(list(geom = "point", x = "score", y = "score", color = "score")),
       scale = list(color = list(midpoint = 5))),
  fd, pd, mat, character(), character()
)
mid_plot <- agent_build_figure_plot(
  agent_build_figure_data(mid_spec, fd, pd, mat), mid_spec
)
mid_classes <- vapply(mid_plot$scales$scales, function(s) class(s)[[1]], character(1))
ok(
  any(grepl("ScaleContinuous", mid_classes, fixed = TRUE)),
  "midpoint override installs a diverging gradient"
)
# unmatched scale levels warn (collected by the tool handler), empty filter
# warns but still renders
warn_msgs <- character()
withCallingHandlers(
  {
    warn_spec <- agent_normalize_figure_spec(
      list(data_source = "feature_annotation",
           layers = list(
             list(geom = "point", x = "score", y = "score", fill = "category"),
             list(geom = "point", x = "score", y = "score",
                  filter = list(column = "score", op = ">", value = 99))
           ),
           scale = list(fill = list(values = list(
             list(category = "kinase", color = "#2166ac"),
             list(category = "ghost-level", color = "#000000")
           )))),
      fd, pd, mat, character(), character()
    )
    agent_build_figure_plot(
      agent_build_figure_data(warn_spec, fd, pd, mat), warn_spec
    )
  },
  warning = function(w) {
    warn_msgs <<- c(warn_msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
ok(
  any(grepl("levels not present", warn_msgs, fixed = TRUE)) &&
    any(grepl("no explicit color for", warn_msgs, fixed = TRUE)) &&
    any(grepl("filter matches no rows", warn_msgs, fixed = TRUE)),
  "unmatched scale levels and empty filters surface as warnings"
)
# theme options land on the plot
themed <- agent_build_figure_plot(
  agent_build_figure_data(wide_norm, fd, pd, mat), wide_norm
)
ok(
  identical(themed$theme$legend.position, "bottom") &&
    !is.null(themed$theme$axis.text.x$angle) &&
    identical(themed$theme$axis.text.x$angle, 30) &&
    inherits(themed$theme$panel.grid.major, "element_blank"),
  "theme options apply legend, rotation, and grid tweaks"
)
# grammar advertises the new surface
grammar_wide <- agent_figure_grammar()
ok(
  ut_cmp_identical(grammar_wide$limits$max_filter_leaves, 8L) &&
    identical(sort(names(grammar_wide$layer_filters$ops)),
              sort(omicsViewer:::.agent_figure_where_ops)) &&
    identical(grammar_wide$scale_overrides$channels, c("color", "fill")) &&
    !is.null(grammar_wide$theme_options$fields$legend_position),
  "figure grammar documents filters, scale overrides, and theme options"
)

## ---- Stage 2 (todo 3.1/3.8): bounded echo, patch mode, registry LRU ----
agent_figure_spec_patch <- omicsViewer:::agent_figure_spec_patch
agent_figure_registry_evict <- omicsViewer:::agent_figure_registry_evict

# todo 3.1 size budget: volcano + labels + echo on a synthetic 20k-feature
# dataset stays kilobytes (the pre-Stage-2 spec resolved all 20k ids into
# features AND echoed them in every tool result)
fd20k <- data.frame(
  logFC = rnorm(20000), logFdr = abs(rnorm(20000)),
  row.names = paste0("G", 1:20000), check.names = FALSE)
v20k <- agent_figure_template_spec(
  template = "volcano", x = "logFC", y = "logFdr", label_top_n = 10L,
  feature_data = fd20k, sample_data = pd)
echo20k <- agent_figure_spec_echo(
  v20k$spec, all_features = rownames(fd20k), all_samples = rownames(pd))
ok(
  is.null(v20k$spec$features) && is.null(echo20k$features) &&
    nchar(jsonlite::toJSON(echo20k, auto_unbox = TRUE), type = "bytes") < 4000L,
  "volcano label_top_n on a 20k-feature dataset yields a <4 KB echo"
)

# patch semantics (todo 3.1): unmentioned fields keep the base value
patch_base <- agent_normalize_figure_spec(
  list(data_source = "feature_annotation",
       layers = list(list(geom = "point", x = "score", y = "score")),
       theme = "minimal", palette = "colorblind",
       labels = list(title = "Base", x = "score")),
  fd, pd, mat, character(), character())
patched <- agent_figure_spec_patch(patch_base, list(
  theme = "classic",
  labels = list(title = "Revised", y = NULL),
  layers = list(list(geom = "point", x = "score", y = "category",
                     params = list(alpha = 0.4)))
))
ok(
  ut_cmp_identical(patched$theme, "classic") &&
    ut_cmp_identical(patched$palette, "colorblind") &&
    ut_cmp_identical(patched$labels$title, "Revised") &&
    ut_cmp_identical(patched$labels$x, "score") &&
    is.null(patched$labels$y) &&
    ut_cmp_identical(patched$layers[[1]]$params$alpha, 0.4),
  "patch replaces mentioned fields, keeps unmentioned, nulls clear keys"
)
patched_norm <- agent_normalize_figure_spec(
  patched, fd, pd, mat, character(), character())
ok(
  ut_cmp_identical(patched_norm$theme, "classic") &&
    ut_cmp_identical(patched_norm$layers[[1]]$mappings$y, "category"),
  "patched specs re-normalize against the live dataset"
)
ok(
  ut_cmp_error(
    agent_figure_spec_patch(patch_base, list(bogus = 1)),
    "Unknown figure changes field"
  ),
  "unknown patch fields are rejected for model self-correction"
)
# clearing features back to the default set
patch_clear <- agent_figure_spec_patch(
  agent_normalize_figure_spec(
    list(data_source = "feature_annotation",
         features = c("F1", "F2"),
         layers = list(list(geom = "point", x = "score", y = "score"))),
    fd, pd, mat, character(), character()),
  list(features = NULL))
ok(
  is.null(patch_clear$features),
  "explicit null on an array field clears it (re-normalization re-defaults)"
)

# registry LRU eviction (todo 3.8)
mk_registry <- function(n, parent_of = NULL, used = NULL) {
  out <- list()
  for (i in seq_len(n)) {
    out[[paste0("fig_", i)]] <- list(
      id = paste0("fig_", i),
      parent_id = if (i %in% parent_of) "fig_1" else NULL,
      created_at = sprintf("2026-01-%02dT00:00:00Z", i),
      last_used_at = used[[i]]
    )
  }
  out
}
heads_only <- mk_registry(20)
ok(
  ut_cmp_identical(
    agent_figure_registry_evict(heads_only, exclude = "fig_20")$evicted, "fig_1"
  ),
  "registry at capacity evicts the oldest lineage head"
)
ok(
  ut_cmp_identical(
    agent_figure_registry_evict(mk_registry(19))$evicted, character()),
  "registry below capacity evicts nothing"
)
superseded <- mk_registry(20, parent_of = 2:20)
ok(
  ut_cmp_identical(
    agent_figure_registry_evict(superseded, exclude = "fig_20")$evicted, "fig_1"
  ),
  "superseded revisions evict before lineage heads (fig_1 has children)"
)
stale_touch <- mk_registry(20)
stale_touch[["fig_15"]]$last_used_at <- "2025-06-01T00:00:00Z"
ok(
  ut_cmp_identical(
    agent_figure_registry_evict(stale_touch)$evicted, "fig_15"
  ),
  "LRU uses last_used_at, not creation order"
)
ok(
  ut_cmp_identical(
    agent_figure_registry_evict(
      mk_registry(2, parent_of = 2), capacity = 2L, exclude = "fig_1")$evicted,
    "fig_2"
  ),
  "the figure under revision is never evicted"
)
ok(
  ut_cmp_identical(
    length(agent_figure_registry_evict(heads_only)$registry), 19L),
  "eviction leaves room for exactly one new entry"
)
# grammar advertises the revision contract + label ordering
grammar_s2 <- agent_figure_grammar()
ok(
  !is.null(grammar_s2$revision$description) &&
    !is.null(grammar_s2$label_layers$example$params$order_by),
  "figure grammar documents patch-mode revision and label ordering"
)

## ---- Stage 3 (todo 4.2): single-source grammar + narrowed schema --------
agent_figure_spec_schema <- omicsViewer:::agent_figure_spec_schema
agent_figure_where_schema <- omicsViewer:::agent_figure_where_schema
agent_where_normalize <- omicsViewer:::agent_where_normalize

# The model-facing contract narrows filter nesting to depth 2 (one
# combinator around leaves); the validator keeps accepting depth 3 so
# pre-narrowing transcripts and golden replays still parse.
where_t1 <- agent_figure_where_schema()
ok(
  ut_cmp_identical(
    "all" %in% names(where_t1@properties), TRUE) &&
    ut_cmp_identical(
      "all" %in% names(
        where_t1@properties$all@items@properties), FALSE),
  "filter schema offers one combinator level only (depth 2)"
)
deep3 <- list(column = "score", op = ">", value = 1)
for (i in 1:2) deep3 <- list(all = list(deep3))
norm_deep3 <- agent_where_normalize(deep3)
ok(
  ut_cmp_identical(
    isTRUE(tryCatch(
      {omicsViewer:::.agent_where_validate_caps(norm_deep3); TRUE},
      error = function(e) FALSE)), TRUE),
  "validator still accepts depth-3 filters (replay compatibility)"
)
deep4 <- list(column = "score", op = ">", value = 1)
for (i in 1:3) deep4 <- list(all = list(deep4))
ok(
  ut_fails(
    omicsViewer:::.agent_where_validate_caps(
      agent_where_normalize(deep4)),
    "nests at most 3 levels"),
  "validator still rejects filters beyond depth 3"
)

# The spec schema derives from the grammar tables (single source):
# layer params in table order, geom enum == the geoms constant.
patch_type <- agent_figure_spec_schema(required = FALSE, layers_required = FALSE)
ok(
  ut_cmp_identical(
    names(patch_type@properties$layers@items@properties$params@properties),
    names(omicsViewer:::.agent_figure_param_specs())),
  "layer params schema comes from the grammar table (table order)"
)
ok(
  ut_cmp_identical(
    patch_type@properties$layers@items@properties$geom@values,
    omicsViewer:::.agent_figure_geoms),
  "geom enum comes from the geoms constant"
)
ok(
  ut_cmp_identical(
    names(patch_type@properties$theme_options@properties),
    names(omicsViewer:::.agent_figure_theme_option_specs())),
  "theme option schema comes from the grammar table"
)

# Validator defaults come from the same table.
plain <- agent_normalize_figure_spec(
  list(layers = list(list(geom = "point", x = "score", y = "score"))),
  fd, pd, mat, character(), character())
specs_table <- omicsViewer:::.agent_figure_param_specs()
ok(
  ut_cmp_identical(plain$layers[[1]]$params$size, specs_table$size$default) &&
    ut_cmp_identical(plain$layers[[1]]$params$alpha, specs_table$alpha$default) &&
    ut_cmp_identical(plain$layers[[1]]$params$max_labels,
                     specs_table$max_labels$default),
  "validator parameter defaults come from the grammar table"
)

# Geom requirements are table-driven for both validator and prose.
reqs <- omicsViewer:::.agent_figure_geom_required()
text_missing_label <- tryCatch(
  agent_normalize_figure_spec(
    list(layers = list(list(geom = "text", x = "score", y = "score"))),
    fd, pd, mat, character(), character()),
  error = function(e) conditionMessage(e))
ok(
  ut_cmp_identical(
    grepl("requires mapping", text_missing_label), TRUE) &&
    ut_cmp_identical(
      identical(reqs$text, c("x", "y", "label")) &&
        identical(reqs$bar, "x") && identical(reqs$hline, character()),
      TRUE),
  "geom requirements come from the grammar table"
)

# The prose grammar advertises the schema contract, not the looser
# validator cap, and its params/geoms sections derive from the tables.
grammar_42 <- agent_figure_grammar()
ok(
  ut_cmp_identical(grammar_42$limits$max_filter_depth, 2L) &&
    grepl("2 nesting levels", grammar_42$layer_filters$form, fixed = TRUE),
  "prose grammar advertises the depth-2 filter contract"
)
ok(
  ut_cmp_identical(names(grammar_42$layer_params),
                   names(specs_table)) &&
    ut_cmp_identical(names(grammar_42$geom_requirements),
                     omicsViewer:::.agent_figure_geoms),
  "prose grammar derives params and geom requirements from the tables"
)
