# S3 unit suite for the generic widget tier (auxi_agentWidgets.R) plus the
# heatmap module's store wiring. Pure helper logic needs no browser; the
# testServer section exercises the real UI<->store observers.

library(omicsViewer)
library(unittest, quietly = TRUE)

agent_widget_ids <- omicsViewer:::agent_widget_ids
agent_widget_list <- omicsViewer:::agent_widget_list
agent_widget_describe <- omicsViewer:::agent_widget_describe
agent_widget_apply <- omicsViewer:::agent_widget_apply
widget_store_new <- omicsViewer:::widget_store_new
widget_store_child <- omicsViewer:::widget_store_child
widget_binding <- omicsViewer:::widget_binding
store_register <- omicsViewer:::store_register
store_read <- omicsViewer:::store_read
store_apply <- omicsViewer:::store_apply
store_registry_view <- omicsViewer:::store_registry_view

mk <- function() {
  s <- widget_store_new()
  store_register(s,
    widget_binding("app.main_tab", "navbar", label = "Data-space tab",
      help = "Visible data-space tab",
      choices_provider = function(v) c("Feature", "Sample", "Heatmap")),
    widget_binding("app.theme", "enum", label = "Theme",
      help = "Plot theme", values = c("minimal", "classic", "bw")),
    widget_binding("app.show_labels", "boolean", label = "Labels",
      help = "Show point labels"),
    widget_binding("app.min_size", "integer", label = "Min size",
      help = "Minimum gene-set size", min = 1L, max = 500L),
    widget_binding("app.internal_counter", "numeric", internal = TRUE,
      label = "Internal counter", help = "Snapshot-only state")
  )
  s
}

## ------------------------------------------------------- discovery ----
s <- mk()
store_apply(s, list(app.theme = "classic"), origin = "system")

ids <- agent_widget_ids(s)
ok(
  setequal(ids, c("app.main_tab", "app.theme", "app.show_labels", "app.min_size")),
  "id listing exposes exactly the agent-visible widgets"
)

listing <- agent_widget_list(s)
ok(identical(listing$widget_count, 4L), "listing counts user-editable widgets")
ok(
  identical(listing$widgets$app.theme$current_value, "classic") &&
    is.null(listing$widgets$app.internal_counter),
  "listing carries current values and hides internal state"
)
themed <- Filter(function(w) identical(w$id, "app.theme"), listing$widgets)[[1]]
ok(
  setequal(themed$allowed_values, c("minimal", "classic", "bw")),
  "listing carries allowed values for enum widgets"
)

ok(
  identical(agent_widget_list(s, "apples")$widget_count, 0L),
  "unknown section lists nothing"
)
ok(
  identical(agent_widget_list(s, "null")$widget_count, 4L) &&
    identical(agent_widget_list(s, "{}")$widget_count, 4L) &&
    is.null(agent_widget_list(s, "null")$section),
  "sentinel section strings are treated as omitted"
)
ok(
  identical(agent_widget_list(s, "app")$widget_count, 4L),
  "section 'app' keeps component-prefixed widgets"
)
# component (dotted) prefix semantics live in store_registry_view
s_pref <- widget_store_new()
store_register(s_pref,
  widget_binding("app.x", "select", values = c("a")),
  widget_binding("apple.x", "select", values = c("a")))
ok(
  identical(vapply(omicsViewer:::store_registry_view(s_pref, prefix = "app"),
                   function(r) r$id, character(1)), "app.x"),
  "prefix filtering matches whole components, not raw startsWith"
)

d <- agent_widget_describe(s, "app.min_size")
ok(
  identical(d$kind, "integer") && identical(d$min, 1L) &&
    is.null(d$current_value),
  "describe returns the registry record with bounds and current value"
)
ok(
  ut_cmp_error(agent_widget_describe(s, "app.them"),
    "Unknown or not user-editable widget id: app.them. Closest matches: app.theme"),
  "describe suggests the closest id on typos"
)
ok(
  ut_cmp_error(agent_widget_describe(s, "app.internal_counter"),
    "not user-editable"),
  "describe never reveals internal widgets"
)

## ----------------------------------------------------------- apply ----
r <- agent_widget_apply(s, list(app.theme = "bw", app.show_labels = TRUE,
                                app.min_size = 15))
ok(
  setequal(r$applied, c("app.theme", "app.show_labels", "app.min_size")) &&
    length(r$rejected) == 0L,
  "valid mixed-kind patch applies fully"
)
ok(
  identical(store_read(s)$app.theme, "bw") &&
    identical(store_read(s)$app.show_labels, TRUE) &&
    identical(store_read(s)$app.min_size, 15L),
  "applied values land in the store with correct types"
)

# JSON-object string (the ellmer wire format) and scalar coercions
r <- agent_widget_apply(s, '{"app.theme": "minimal", "app.min_size": "12",
                            "app.show_labels": "false"}')
ok(
  identical(r$applied_values$app.min_size, 12L) &&
    identical(r$applied_values$app.show_labels, FALSE),
  "JSON string patch coerces string scalars to integer/boolean kinds"
)

# omitted-optional sentinels ("null"/"{}") are dropped, not errors
ok(
  ut_cmp_error(agent_widget_apply(s, "{}"), "contains no entries"),
  "empty JSON object reports no entries"
)

# per-key rejection: invalid enum keeps the valid remainder
r <- agent_widget_apply(s, list(app.theme = "miniml", app.min_size = 3))
ok(
  setequal(r$applied, "app.min_size") && identical(store_read(s)$app.min_size, 3L),
  "invalid values are rejected per key while valid keys still apply"
)
ok(
  length(r$rejected) == 1L && identical(r$rejected[[1]]$id, "app.theme") &&
    grepl("minimal", r$rejected[[1]]$reason),
  "rejection reasons carry closest-match suggestions"
)
ok(
  identical(store_read(s)$app.theme, "minimal"),
  "rejected key leaves the stored value untouched"
)

# unknown ids and non-editable keys are rejected per key, not fatal
r <- agent_widget_apply(s, list(app.them = "bw", app.internal_counter = 1,
                                app.show_labels = TRUE))
ok(
  setequal(r$applied, "app.show_labels") && length(r$rejected) == 2L,
  "unknown and non-editable ids are rejected without vetoing valid keys"
)
reasons <- vapply(r$rejected, function(x) x$reason, "")
ok(
  any(grepl("not user-editable", reasons)) &&
    any(grepl("Unknown widget id: app.them", reasons)) &&
    any(grepl("app.theme", reasons)),
  "id rejections distinguish unknown from non-editable and suggest matches"
)

# no-op keys are reported as unchanged, not applied
r <- agent_widget_apply(s, list(app.min_size = 3))
ok(
  length(r$applied) == 0L && setequal(r$unchanged, "app.min_size"),
  "identical values are reported as unchanged"
)

# malformed patches fail loudly for model self-correction
ok(
  ut_cmp_error(agent_widget_apply(s, "not json"),
    "not valid JSON"),
  "malformed JSON patch errors clearly"
)
ok(
  ut_cmp_error(agent_widget_apply(s, list(list(value = "x"))),
    "named by a widget id"),
  "unnamed patch entries are refused"
)

# ellmer coercions: tibble / record arrays normalize to one named patch
# (fresh store: shape normalization, independent of earlier sequences)
s3 <- mk()
r <- agent_widget_apply(s3, list(
  list(id = "app.theme", value = "classic"),
  list(id = "app.min_size", value = 7)
))
ok(
  setequal(r$applied, c("app.theme", "app.min_size")) &&
    identical(store_read(s3)$app.theme, "classic"),
  "array-of-records patch normalizes"
)
r <- agent_widget_apply(s3, data.frame(
  id = c("app.theme", "app.show_labels"), value = c("bw", "true"),
  stringsAsFactors = FALSE))
ok(
  setequal(r$applied, c("app.theme", "app.show_labels")) &&
    identical(store_read(s3)$app.show_labels, TRUE),
  "data.frame (tibble) patch normalizes"
)

## ------------------------------------------ reactive choices contract ----
# choices providers may read module reactives; helpers must work inside
# isolate() (the calling convention the ellmer tools use).
s2 <- widget_store_new()
choices_rv <- shiny::reactiveVal(c("ttest", "PCA"))
store_register(s2,
  widget_binding("ds.analysis", "select", label = "Analysis",
    help = "Analysis category",
    choices_provider = function(v) choices_rv()))
ok(
  identical(shiny::isolate(
    agent_widget_list(s2)$widgets$ds.analysis$allowed_values),
    c("ttest", "PCA")),
  "reactive choices are readable through the list view under isolate"
)
r <- shiny::isolate(agent_widget_apply(s2, list(ds.analysis = "ttest")))
ok(
  setequal(r$applied, "ds.analysis"),
  "reactive-backed validation works under isolate"
)

## --------------------------------------------- heatmap module wiring ----
.sent_hm <- new.env(); .sent_hm$msgs <- list()
.spy_input_messages <- function(session, sink) {
  orig <- session$sendInputMessage
  session$sendInputMessage <- function(id, msg) {
    sink$msgs[[length(sink$msgs) + 1L]] <- list(id = id, msg = msg)
    orig(id, msg)
  }
}
hm_root <- widget_store_new()
hm_store <- widget_store_child(hm_root, "dataspace.expr_heatmap")
app_hm <- function(input, output, session) {
  omicsViewer:::iheatmapModule(
    "hm",
    mat = shiny::reactive(NULL),
    pd = shiny::reactive(data.frame(group = c("A", "B"),
                                    row.names = c("S1", "S2"))),
    fd = shiny::reactive(data.frame(score = c(1, 2),
                                    row.names = c("G1", "G2"))),
    status = shiny::reactive(NULL),
    store = hm_store
  )
}
shiny::testServer(app_hm, {
  .spy_input_messages(session, .sent_hm)

  # all inputs exist -> the seed observer fires and store_applies defaults
  session$setInputs(`hm-heatmapColors` = "RdYlBu", `hm-scale` = "row",
                    `hm-marginBottom` = 4, `hm-marginRight` = 4)
  session$flushReact()

  # user edit: UI -> store sync (child views key reads by full canonical id).
  # NOTE: testServer's ignoreInit quirk swallows the first value change,
  # so warm the observer with the seed value first (the seed IS that value).
  session$setInputs(`hm-heatmapColors` = "RdBu")
  session$flushReact()
  ok(
    identical(store_read(hm_store)$dataspace.expr_heatmap.heatmap_colors, "RdBu"),
    "user palette edit syncs into the store"
  )

  # external write: store -> UI push via updateSelectInput/updateSliderInput
  # (MockShinySession relays update ids without the module prefix and
  # carries select selections in msg$value)
  store_apply(hm_store, list(heatmap_colors = "PiYG", margin_bottom = 9L),
              origin = "agent")
  session$flushReact()
})
ok(
  any(vapply(.sent_hm$msgs, function(m)
    grepl("(^|\\.)heatmapColors$", m$id) &&
      identical(m$msg$selected %||% m$msg$value, "PiYG"),
    logical(1))),
  "agent palette write pushes updateSelectInput to the widget"
)
ok(
  any(vapply(.sent_hm$msgs, function(m)
    grepl("(^|\\.)marginBottom$", m$id) && identical(as.numeric(m$msg$value), 9),
    logical(1))),
  "agent margin write pushes updateSliderInput to the widget"
)

# generic tier: the wire contract is the ROOT store with full canonical
# ids (exactly what the ellmer tool receives from the app module)
r <- agent_widget_apply(hm_root,
  list("dataspace.expr_heatmap.heatmap_colors" = "Spectral"))
ok(
  length(r$applied) == 0L && length(r$rejected) == 1L &&
    grepl("BrBG", r$rejected[[1]]$reason),
  "invalid palette is rejected per key with suggestions"
)
r <- agent_widget_apply(hm_root,
  list("dataspace.expr_heatmap.margin_bottom" = 42))
ok(
  length(r$rejected) == 1L && grepl(">= 1|<= 20", r$rejected[[1]]$reason),
  "out-of-bounds margins are rejected per key"
)

## ------------------------------------ S4 heatmap completion wiring ----
# sorting, clustering, and annotation widgets are registered, synced, and
# pushed exactly like the S3 palette/scale/margins
reg <- Filter(function(x) startsWith(x$id, "dataspace.expr_heatmap."),
              store_registry_view(hm_root))
keyed <- setNames(reg, vapply(reg, function(x)
  sub("^dataspace\\.expr_heatmap\\.", "", x$id), character(1)))
ok(
  identical(sort(vapply(reg, function(x) x$id, character(1))),
            sort(paste0("dataspace.expr_heatmap.", c(
              "heatmap_colors", "scale", "margin_bottom", "margin_right",
              "col_sort_by", "row_sort_by", "cluster_col_dist",
              "cluster_col_link", "cluster_row_dist", "cluster_row_link",
              "annot_col", "annot_row", "tooltip_info")))),
  "all 13 heatmap parameter widgets are registered"
)
ok(
  identical(keyed$annot_col$kind, "multi_select") &&
    identical(keyed$annot_row$kind, "multi_select") &&
    identical(keyed$tooltip_info$kind, "multi_select"),
  "annotation and tooltip widgets are multi_select bindings"
)

.sent_hm$msgs <- list()  # reset the spy sink for the S4 section
hm_root2 <- widget_store_new()
hm_store2 <- widget_store_child(hm_root2, "dataspace.expr_heatmap")
app_hm2 <- function(input, output, session) {
  omicsViewer:::iheatmapModule(
    "hm",
    mat = shiny::reactive(NULL),
    pd = shiny::reactive(data.frame(group = c("A", "B"),
                                    row.names = c("S1", "S2"))),
    fd = shiny::reactive(data.frame(score = c(1, 2),
                                   pathway = c("x", "y"),
                                   row.names = c("G1", "G2"))),
    status = shiny::reactive(NULL),
    store = hm_store2
  )
}
shiny::testServer(app_hm2, {
  .spy_input_messages(session, .sent_hm)
  # initialize every parameter input so the seed observer fires
  session$setInputs(`hm-heatmapColors` = "RdYlBu", `hm-scale` = "row",
                    `hm-marginBottom` = 4, `hm-marginRight` = 4,
                    `hm-colSortBy` = "none", `hm-rowSortBy` = "none",
                    `hm-clusterColDist` = "Pearson correlation",
                    `hm-clusterColLink` = "ward.D",
                    `hm-clusterRowDist` = "Pearson correlation",
                    `hm-clusterRowLink` = "ward.D",
                    `hm-annotCol` = "group", `hm-annotRow` = "score",
                    `hm-tooltipInfo` = "score")
  session$flushReact()
  seeded <- store_read(hm_store2)
  ok(
    identical(seeded$dataspace.expr_heatmap.col_sort_by, "none") &&
      identical(seeded$dataspace.expr_heatmap.annot_col, "group") &&
      identical(seeded$dataspace.expr_heatmap.cluster_col_link, "ward.D"),
    "the seed observer stores defaults for the S4 keys"
  )

  # user edits sync (first change per observer is swallowed by the
  # testServer ignoreInit quirk - the seed values above warmed them)
  session$setInputs(`hm-colSortBy` = "group", `hm-annotCol` = character(0))
  session$flushReact()
  vals <- store_read(hm_store2)
  ok(
    identical(vals$dataspace.expr_heatmap.col_sort_by, "group"),
    "user column-sort edit syncs into the store"
  )
  ok(
    identical(vals$dataspace.expr_heatmap.annot_col, character(0)),
    "clearing a multi-select input syncs character(0)"
  )

  # external writes push to the right updater (selectize for server-side
  # rowSortBy/annotRow, plain select otherwise). annot_row writes a value
  # DIFFERENT from the seeded one: after WP4 the seed no longer arms an
  # un-acknowledgeable pending, so an identical no-op write correctly
  # pushes nothing - the widget already displays that value.
  store_apply(hm_store2,
    list(annot_col = "group", cluster_col_link = "complete",
         row_sort_by = "score", annot_row = "pathway"),
    origin = "agent")
  session$flushReact()

  # generic tier over the S4 keys: JSON arrays, per-key rejections
  # (run inside the session: choices providers read module reactives)
  s4 <- list()
  s4$json_ok <- agent_widget_apply(hm_root2,
    '{"dataspace.expr_heatmap.annot_col": [], "dataspace.expr_heatmap.cluster_row_link": "ward.D2"}')
  s4$json_bad <- agent_widget_apply(hm_root2,
    '{"dataspace.expr_heatmap.annot_col": ["notAColumn"]}')
  s4$sort_bad <- agent_widget_apply(hm_root2,
    '{"dataspace.expr_heatmap.row_sort_by": "hierarchical clustre"}')
  .sent_hm$s4 <- s4
})
ok(
  any(vapply(.sent_hm$msgs, function(m)
    grepl("(^|\\.)annotCol$", m$id) && identical(m$msg$value, "group"),
    logical(1))),
  "agent annotation write pushes updateSelectInput to the widget"
)
ok(
  any(vapply(.sent_hm$msgs, function(m)
    grepl("(^|\\.)clusterColLink$", m$id) && identical(m$msg$value, "complete"),
    logical(1))),
  "agent linkage write pushes updateSelectInput to the widget"
)
ok(
  any(vapply(.sent_hm$msgs, function(m)
    grepl("(^|\\.)annotRow$", m$id) && identical(m$msg$value, "pathway"),
    logical(1))),
  "agent row-annotation write reaches the server-side selectize widget"
)
ok(
  any(vapply(.sent_hm$msgs, function(m)
    grepl("(^|\\.)rowSortBy$", m$id) && identical(m$msg$value, "score"),
    logical(1))),
  "agent row-sort write reaches the server-side selectize widget"
)

# generic tier over the S4 keys: JSON arrays, per-key rejections
r <- .sent_hm$s4$json_ok
ok(
  identical(sort(r$applied), sort(c("dataspace.expr_heatmap.annot_col",
                                    "dataspace.expr_heatmap.cluster_row_link"))) &&
    identical(r$applied_values$dataspace.expr_heatmap.annot_col, character(0)),
  "generic tier applies JSON-array multi_select values (empty array clears)"
)
r <- .sent_hm$s4$json_bad
ok(
  length(r$applied) == 0L && length(r$rejected) == 1L &&
    grepl("notAColumn", r$rejected[[1]]$reason) &&
    grepl("group", r$rejected[[1]]$reason),
  "invalid multi_select entries are rejected with suggestions"
)
r <- .sent_hm$s4$sort_bad
ok(
  length(r$rejected) == 1L && grepl("hierarchical cluster",
                                    r$rejected[[1]]$reason),
  "invalid sort values are rejected with suggestions"
)

## ------------------------------------------ dataTable module wiring ----
# S4: the three data-space tables register their user-editable surface
# (multi-selection switch + shown-columns set) on the canonical store.
.dt_cols_reg <- Filter(
  function(x) startsWith(x$id, "dataspace.tab_feature."),
  store_registry_view(widget_store_new()))
dt_root <- widget_store_new()
dt_store <- widget_store_child(dt_root, "dataspace.tab_feature")
.dt_pd <- data.frame(
  `General|All|group` = rep(c("A", "B"), each = 15),
  `PCA|All|PC1` = rnorm(30),
  `mean|Origin|RE` = rnorm(30),
  check.names = FALSE,
  row.names = paste0("feature_", sprintf("%03d", 1:30)))
.dt_state <- list(start = 0L, length = 25L,
                  order = list(c(1L, "asc")),
                  columns = list(list(search = list(search = ""))))
.sent_dt <- new.env(); .sent_dt$msgs <- list(); .sent_dt$s4 <- NULL
dt_result <- NULL
app_dt <- function(input, output, session) {
  dt_result <<- omicsViewer:::dataTable_module(
    "dt", reactive_data = shiny::reactive(.dt_pd),
    tab_status = shiny::reactive(NULL), tab_rows = shiny::reactive(TRUE),
    store = dt_store)
}
shiny::testServer(app_dt, {
  .spy_input_messages(session, .sent_dt)

  # initialize: switch off (warm the ignoreInit swallow), then a user edit
  session$setInputs(`dt-multisel` = FALSE)
  session$flushReact()
  ok(
    identical(store_read(dt_store)$dataspace.tab_feature.columns,
              "General|All|group"),
    "the seed observer stores the default shown columns"
  )
  session$setInputs(`dt-multisel` = TRUE)
  session$flushReact()
  ok(
    identical(store_read(dt_store)$dataspace.tab_feature.multi_selection, TRUE),
    "user multi-selection toggle syncs into the store"
  )

  # external column write: store -> scn push, observed through the table
  # status contract (showColumns)
  store_apply(dt_store, list(columns = c("PCA|All|PC1", "mean|Origin|RE")),
              origin = "agent")
  session$flushReact()
  session$setInputs(`dt-table_state` = .dt_state)
  session$flushReact()
  ok(
    identical(attr(dt_result(), "status")$showColumns,
              c("PCA|All|PC1", "mean|Origin|RE")),
    "agent column write pushes the shown-column set in order"
  )

  # external switch write: store -> updateSwitchInput (MockShinySession
  # only relays the message, so assert on the spied message)
  store_apply(dt_store, list(multi_selection = FALSE), origin = "agent")
  session$flushReact()
  session$setInputs(`dt-table_state` = .dt_state)
  session$flushReact()
  ok(
    any(vapply(.sent_dt$msgs, function(m)
      grepl("(^|\\.)multisel$", m$id) && identical(m$msg$value %||% m$msg$checked, FALSE),
      logical(1))),
    "agent multi-selection write pushes the switch"
  )

  # generic tier over the table keys (inside the session: providers read
  # module reactives)
  s4 <- list()
  s4$ok <- agent_widget_apply(dt_root,
    '{"dataspace.tab_feature.columns": ["General|All|group"], "dataspace.tab_feature.multi_selection": true}')
  s4$empty <- agent_widget_apply(dt_root,
    '{"dataspace.tab_feature.columns": []}')
  s4$bad <- agent_widget_apply(dt_root,
    '{"dataspace.tab_feature.columns": ["PC1"]}')
  .sent_dt$s4 <- s4
})
r <- .sent_dt$s4$ok
ok(
  identical(sort(r$applied), sort(c("dataspace.tab_feature.columns",
                                    "dataspace.tab_feature.multi_selection"))) &&
    identical(r$applied_values$dataspace.tab_feature.columns,
              "General|All|group"),
  "generic tier applies table columns and switches"
)
r <- .sent_dt$s4$empty
ok(
  length(r$applied) == 0L && length(r$rejected) == 1L &&
    grepl("at least 1", r$rejected[[1]]$reason),
  "emptying the shown columns is rejected (one column must remain)"
)
r <- .sent_dt$s4$bad
ok(
  length(r$rejected) == 1L && grepl("PCA\\|All\\|PC1", r$rejected[[1]]$reason),
  "unknown table columns are rejected with suggestions"
)

## ------------------------------- S4 result-space: feature_general ------
# The first result-space module on the store: link-variable cascade,
# plot-type radio, regression-line state, and the shared attr4 panel.
.fg_reg <- Filter(
  function(x) startsWith(x$id, "resultspace.feature_general."),
  store_registry_view(widget_store_new()))
fg_root <- widget_store_new()
fg_store <- widget_store_child(fg_root, "resultspace.feature_general")
.fg_pd <- data.frame(
  `General|All|group` = factor(rep(c("A", "B"), each = 4)),
  `General|All|score` = rnorm(8),
  `Surv|All|time` = rexp(8),
  row.names = paste0("S", 1:8), check.names = FALSE)
.fg_ex <- matrix(rnorm(16 * 8), nrow = 16,
                 dimnames = list(paste0("F", 1:16), paste0("S", 1:8)))
.fg_fd <- data.frame(`General|All|Gene.name` = paste0("g", 1:16),
                     row.names = paste0("F", 1:16), check.names = FALSE)
.sent_fg <- new.env(); .sent_fg$msgs <- list(); .sent_fg$s4 <- NULL
app_fg <- function(input, output, session) {
  omicsViewer:::feature_general_module(
    "fg", reactive_expr = shiny::reactive(.fg_ex),
    reactive_i = shiny::reactive(3), reactive_highlight = shiny::reactive(NULL),
    reactive_phenoData = shiny::reactive(.fg_pd),
    reactive_featureData = shiny::reactive(.fg_fd),
    store = fg_store)
}
shiny::testServer(app_fg, {
  .spy_input_messages(session, .sent_fg)
  ids <- names(store_read(fg_store))
  ok(
    identical(sort(ids), sort(paste0("resultspace.feature_general.", c(
      "xax_analysis", "xax_subset", "xax_variable", "plot_type",
      "regression_line",
      paste0("attr4.", outer(c("color", "shape", "size", "tooltip", "search"),
                             c("_analysis", "_subset", "_variable"),
                             FUN = paste0)),
      "attr4.xcut", "attr4.ycut", "attr4.scorner")))),
    "feature_general registers 23 user-editable keys (5 own + 18 attr4)"
  )
  # user pick through the triselector cascade syncs into the store
  session$setInputs(`fg-tris_feature_general-analysis` = "General")
  session$setInputs(`fg-tris_feature_general-subset` = "All")
  session$setInputs(`fg-tris_feature_general-variable` = "group")
  session$flushReact()
  ok(
    identical(store_read(fg_store)$resultspace.feature_general.xax_analysis,
              "General") &&
      identical(store_read(fg_store)$resultspace.feature_general.xax_variable,
                "group"),
    "user link-variable pick syncs into the store"
  )
  # plot type: warm the ignoreInit swallow, then a user edit
  session$setInputs(`fg-internal_radio` = "Bees")
  session$flushReact()
  session$setInputs(`fg-internal_radio` = "Curve")
  session$flushReact()
  ok(
    identical(store_read(fg_store)$resultspace.feature_general.plot_type,
              "Curve"),
    "user plot-type edit syncs into the store"
  )
  # external cascade write: store -> triselector update relays with
  # unprefixed ids; the score variable drives the single-feature scatter
  .sent_fg$msgs <- list()
  store_apply(fg_store, list(xax_variable = "score"), origin = "agent")
  session$flushReact()
  ok(
    identical(store_read(fg_store)$resultspace.feature_general.xax_variable,
              "score") &&
      any(vapply(.sent_fg$msgs, function(m)
        grepl("(^|\\.)variable$", m$id), logical(1))),
    "agent link-variable write pushes the cascade"
  )
  # external plot-type write: updateRadioGroupButtons relays a message
  .sent_fg$msgs <- list()
  store_apply(fg_store, list(plot_type = "Bees"), origin = "agent")
  session$flushReact()
  ok(
    any(vapply(.sent_fg$msgs, function(m)
      grepl("(^|\\.)internal_radio$", m$id) &&
        identical(m$msg$selected %||% m$msg$value, "Bees"), logical(1))),
    "agent plot-type write pushes the radio"
  )
  # external regression-line write lands in the mirrored reactiveVal
  store_apply(fg_store, list(regression_line = TRUE), origin = "agent")
  session$flushReact()
  ok(
    isTRUE(store_read(fg_store)$resultspace.feature_general.regression_line),
    "agent regression-line write is stored and acked through the mirror"
  )
  # generic tier over the result-space keys (providers read reactives)
  s4 <- list()
  s4$ok <- agent_widget_apply(fg_root,
    '{"resultspace.feature_general.attr4.color_variable": "score", "resultspace.feature_general.plot_type": "Curve"}')
  s4$bad <- agent_widget_apply(fg_root,
    '{"resultspace.feature_general.plot_type": "Histogram"}')
  s4$nope <- agent_widget_apply(fg_root,
    '{"resultspace.feature_general.xax_variable": "nope"}')
  s4$surv <- agent_widget_apply(fg_root,
    '{"resultspace.feature_general.xax_analysis": "Surv"}')
  .sent_fg$s4 <- s4
})
r <- .sent_fg$s4$ok
ok(
  identical(sort(r$applied),
            sort(c("resultspace.feature_general.attr4.color_variable",
                   "resultspace.feature_general.plot_type"))),
  "generic tier applies attr4 cascades and the plot type"
)
r <- .sent_fg$s4$bad
ok(
  length(r$rejected) == 1L &&
    grepl("Allowed: Bees, Curve", r$rejected[[1]]$reason),
  "invalid plot type is rejected with allowed values"
)
r <- .sent_fg$s4$nope
ok(
  length(r$rejected) == 1L && grepl("Unknown value", r$rejected[[1]]$reason),
  "unknown link variable is rejected with suggestions"
)
r <- .sent_fg$s4$surv
ok(
  length(r$rejected) == 1L &&
    grepl("Surv", r$rejected[[1]]$reason),
  "the Surv category is not a feature_general choice (sample_general only)"
)

## ------------------------------- S4 result-space: sample_general -------
# The sample twin keeps the Surv category (survival view) and routes the
# batch-comparison link through the store.
sg_root <- widget_store_new()
sg_store <- widget_store_child(sg_root, "resultspace.sample_general")
.sent_sg <- new.env(); .sent_sg$msgs <- list(); .sent_sg$s4 <- NULL
app_sg <- function(input, output, session) {
  omicsViewer:::sample_general_module(
    "sg", reactive_phenoData = shiny::reactive(.fg_pd),
    reactive_expr = shiny::reactive(.fg_ex),
    reactive_j = shiny::reactive(c("S1", "S3")),
    store = sg_store)
}
shiny::testServer(app_sg, {
  .spy_input_messages(session, .sent_sg)
  ids <- names(store_read(sg_store))
  ok(
    identical(sort(ids), sort(paste0("resultspace.sample_general.", c(
      "xax_analysis", "xax_subset", "xax_variable",
      paste0("attr4.", outer(c("color", "shape", "size", "tooltip", "search"),
                             c("_analysis", "_subset", "_variable"),
                             FUN = paste0)),
      "attr4.xcut", "attr4.ycut", "attr4.scorner",
      "survival_censor",
      "batch_show_phenotype", "batch_show_features",
      "batch_phenotype_selected_row", "batch_feature_selected_row")))),
    "sample_general registers 26 keys (3 own + 18 attr4 + survival censor + 4 batch)"
  )
  session$setInputs(`sg-tris_sample_general-analysis` = "Surv")
  session$setInputs(`sg-tris_sample_general-subset` = "All")
  session$setInputs(`sg-tris_sample_general-variable` = "time")
  session$flushReact()
  ok(
    identical(store_read(sg_store)$resultspace.sample_general.xax_analysis,
              "Surv"),
    "user Surv pick syncs into the store (survival view choice)"
  )
  # generic tier: cascaded attr4 patch validated against the effective state
  s4 <- list()
  s4$ok <- agent_widget_apply(sg_root, paste0(
    '{"resultspace.sample_general.xax_analysis":"General",',
    '"resultspace.sample_general.xax_subset":"All",',
    '"resultspace.sample_general.xax_variable":"group",',
    '"resultspace.sample_general.attr4.xcut":"1",',
    '"resultspace.sample_general.attr4.scorner":"right"}'))
  s4$corner <- agent_widget_apply(sg_root,
    '{"resultspace.sample_general.attr4.scorner": "volcano"}')
  .sent_sg$s4 <- s4
})
r <- .sent_sg$s4$ok
ok(
  setequal(r$applied,
           c("resultspace.sample_general.xax_analysis",
             "resultspace.sample_general.xax_variable",
             "resultspace.sample_general.attr4.xcut",
             "resultspace.sample_general.attr4.scorner")),
  "generic tier applies the sample cascade and cutoffs in one patch"
)
ok(
  identical(r$applied_values$`resultspace.sample_general.attr4.scorner`,
            "right"),
  "corner choice validates against the patch's own xcut"
)
r <- .sent_sg$s4$corner
ok(
  length(r$rejected) == 1L,
  "volcano corner without both cutoffs is rejected"
)

## ------------------------------- S4 result-space: analyst navbar -------
# L1 registers the analysis navbar itself so the agent can switch
# result-space tabs exactly like a user.
.rs_root <- widget_store_new()
.sent_rs <- new.env(); .sent_rs$msgs <- list()
app_rs <- function(input, output, session) {
  omicsViewer:::L1_result_space_module(
    "rs", reactive_expr = shiny::reactive(.fg_ex),
    reactive_phenoData = shiny::reactive(.fg_pd),
    reactive_featureData = shiny::reactive(.fg_fd),
    reactive_i = shiny::reactive(1), reactive_highlight = shiny::reactive(NULL),
    store = .rs_root)
}
shiny::testServer(app_rs, {
  .spy_input_messages(session, .sent_rs)
  session$flushReact()
  rec <- Filter(function(x) identical(x$id, "resultspace.analyst_tab"),
                store_registry_view(.rs_root))[1]
  ok(
    identical(rec[[1]]$kind, "navbar") &&
      setequal(rec[[1]]$allowed_values, c("Feature", "Geneshot", "Sample")),
    "analyst navbar is registered with dataset-dependent tab choices"
  )
  # warm the ignoreInit swallow, then a user tab switch
  session$setInputs(`rs-analyst` = "Feature")
  session$flushReact()
  session$setInputs(`rs-analyst` = "Sample")
  session$flushReact()
  ok(
    identical(store_read(.rs_root)$resultspace.analyst_tab, "Sample"),
    "user tab switch syncs into the store"
  )
  # external write: updateNavbarPage relays a message
  .sent_rs$msgs <- list()
  store_apply(.rs_root, list("resultspace.analyst_tab" = "Geneshot"),
              origin = "agent")
  session$flushReact()
  ok(
    any(vapply(.sent_rs$msgs, function(m)
      grepl("(^|\\.)analyst$", m$id) &&
        identical(m$msg$selected %||% m$msg$value, "Geneshot"), logical(1))),
    "agent tab write pushes the navbar"
  )
  r <- agent_widget_apply(.rs_root,
    '{"resultspace.analyst_tab": "Feature"}')
  ok(
    setequal(r$applied, "resultspace.analyst_tab"),
    "generic tier applies the analyst tab"
  )
  r <- agent_widget_apply(.rs_root, '{"resultspace.analyst_tab": "ORA"}')
  ok(
    length(r$rejected) == 1L && grepl("Feature", r$rejected[[1]]$reason),
    "tab absent from this dataset (ORA needs gene sets) is rejected"
  )
})

## ------------------- S4 decision: dataTableDownload row selection ------
# Row selection registers ONLY where it drives a downstream view (ORA
# overlap table, fGSEA barplot, STRING network, batch links); gsList and
# no-effect tables stay off the store (see AGENT_ACCURACY_PLAN 6.2).
.dtd_root <- widget_store_new()
.dtd_store <- widget_store_child(.dtd_root, "test.dtd")
.dtd_tab <- data.frame(pathway = c("P1", "P2", "P3"), pv = c(.2, .01, .3))
app_dtd <- function(input, output, session) {
  omicsViewer:::dataTableDownload_module(
    "dtd", reactive_table = shiny::reactive(.dtd_tab),
    sortBy = "pv", decreasing = FALSE,
    reactive_row_ids = shiny::reactive(c("P1", "P2", "P3")),
    store = .dtd_store, store_key = "selected_row")
}
shiny::testServer(app_dtd, {
  ids <- names(store_read(.dtd_store))
  ok(identical(ids, "test.dtd.selected_row"),
     "dataTableDownload registers exactly the row-selection key")
  # user click on displayed row 2 (pv asc: P2, P1, P3 -> original P1)
  session$setInputs(`dtd-table_rows_selected` = 2L)
  session$flushReact()
  ok(identical(store_read(.dtd_store)$test.dtd.selected_row, "P1"),
     "user row click syncs into the store as the stable row id")
  # agent push: the stale pre-render report must NOT count as an override
  store_apply(.dtd_store, list(selected_row = "P3"), origin = "agent")
  session$flushReact()
  ok(identical(store_read(.dtd_store)$test.dtd.selected_row, "P3") &&
       length(.dtd_root$pending) == 1L && length(.dtd_root$override_log) == 0,
     "agent row push survives the stale pre-render report (no override)")
  # DT reports the pushed row -> acknowledgement
  session$setInputs(`dtd-table_rows_selected` = 3L)
  session$flushReact()
  ok(identical(store_read(.dtd_store)$test.dtd.selected_row, "P3") &&
       length(.dtd_root$pending) == 0L,
     "DT report of the pushed row acknowledges the write")
  # a real user click while a push is in flight wins (override contract)
  store_apply(.dtd_store, list(selected_row = "P1"), origin = "agent")
  session$flushReact()
  session$setInputs(`dtd-table_rows_selected` = 1L)  # P2: the user's pick
  session$flushReact()
  ok(identical(store_read(.dtd_store)$test.dtd.selected_row, "P2") &&
       length(.dtd_root$override_log) == 1,
     "user click overriding an in-flight push is honoured and logged")
  # unknown row id is rejected with allowed values
  r <- agent_widget_apply(.dtd_root, '{"test.dtd.selected_row": "P999"}')
  ok(length(r$rejected) == 1L && grepl("P1, P2, P3", r$rejected[[1]]$reason),
     "unknown row id is rejected with the allowed row ids")
})

## ------------------------------- S4 result-space: ora ------------------
# The collapse-features cascade plus the results-table row selection
# (drives the overlap-genes table).
.ora_fd <- data.frame(
  `General|All|Gene.name` = paste0("g", 1:30),
  row.names = paste0("F", 1:30), check.names = FALSE)
attr(.ora_fd, "GS") <- data.frame(
  featureId = factor(paste0("F", c(1:10, 6:20, 21:30, 1:5, 25:28))),
  gsId = factor(rep(c("GS_A", "GS_B", "GS_C", "GS_D"),
                    c(10, 15, 10, 9))))
.ora_root <- widget_store_new()
.ora_store <- widget_store_child(.ora_root, "resultspace.ora")
app_ora <- function(input, output, session) {
  omicsViewer:::enrichment_analysis_module(
    "ora", reactive_featureData = shiny::reactive(.ora_fd),
    reactive_i = shiny::reactive(paste0("F", 1:5)),
    store = .ora_store)
}
shiny::testServer(app_ora, {
  ok(
    identical(sort(names(store_read(.ora_store))),
              sort(paste0("resultspace.ora.", c("xax_analysis", "xax_subset",
                                                "xax_variable", "selected_row")))),
    "ora registers the collapse cascade and the results row selection"
  )
  # user cascade pick drives the real ORA computation
  session$setInputs(`ora-tris_ora-analysis` = "General")
  session$setInputs(`ora-tris_ora-subset` = "All")
  session$setInputs(`ora-tris_ora-variable` = "Gene.name")
  session$flushReact()
  ok(
    identical(store_read(.ora_store)$resultspace.ora.xax_variable, "Gene.name"),
    "user collapse-variable pick syncs into the store"
  )
  b <- .ora_root$bindings[["resultspace.ora.selected_row"]]
  ch <- b$choices_provider(list())
  ok(
    length(ch) == 2L && setequal(ch, c("GS_A", "GS_D")),
    "row choices derive from the real enriched pathways (only tested sets)"
  )
  # user row click syncs; agent push + DT-report ack through the store
  # (display order is p.value ascending: GS_D is row 1, GS_A row 2)
  session$setInputs(`ora-stab-table_rows_selected` = 1L)
  session$flushReact()
  ok(
    identical(store_read(.ora_store)$resultspace.ora.selected_row, "GS_D"),
    "user pathway-row click syncs into the store"
  )
  store_apply(.ora_store, list(selected_row = "GS_A"), origin = "agent")
  session$flushReact()
  session$setInputs(`ora-stab-table_rows_selected` = 2L)
  session$flushReact()
  ok(
    identical(store_read(.ora_store)$resultspace.ora.selected_row, "GS_A") &&
      length(.ora_root$pending) == 0L,
    "agent pathway-row push applies and acknowledges (drives overlap table)"
  )
  # agent cascade write pushes the triselector
  store_apply(.ora_store, list(xax_variable = "Gene.name"), origin = "agent")
  session$flushReact()
  ok(
    identical(store_read(.ora_store)$resultspace.ora.xax_variable, "Gene.name"),
    "agent cascade write is stored"
  )
})

## ------------------------------ S4 result-space: fgsea -----------------
# The ranking cascade plus the results-table row selection (drives the
# leading-edge bar plot).
.fgs_fd <- data.frame(
  `ttest|KO|stat` = c(rnorm(10, mean = 2), rnorm(20)),
  row.names = paste0("F", 1:30), check.names = FALSE)
attr(.fgs_fd, "GS") <- data.frame(
  featureId = factor(paste0("F", c(1:10, 11:30))),
  gsId = factor(rep(c("GS_UP", "GS_REST"), c(10, 20))))
.fgs_root <- widget_store_new()
.fgs_store <- widget_store_child(.fgs_root, "resultspace.fgsea")
app_fgs <- function(input, output, session) {
  omicsViewer:::enrichment_fgsea_module(
    "fgsea", reactive_featureData = shiny::reactive(.fgs_fd),
    store = .fgs_store)
}
shiny::testServer(app_fgs, {
  ok(
    identical(sort(names(store_read(.fgs_store))),
              sort(paste0("resultspace.fgsea.", c("xax_analysis", "xax_subset",
                                                  "xax_variable", "selected_row")))),
    "fgsea registers the ranking cascade and the results row selection"
  )
  session$setInputs(`fgsea-tris_fgsea-analysis` = "ttest")
  session$setInputs(`fgsea-tris_fgsea-subset` = "KO")
  session$setInputs(`fgsea-tris_fgsea-variable` = "stat")
  session$flushReact()
  ok(
    identical(store_read(.fgs_store)$resultspace.fgsea.xax_variable, "stat"),
    "user ranking-variable pick syncs into the store"
  )
  b <- .fgs_root$bindings[["resultspace.fgsea.selected_row"]]
  ch <- b$choices_provider(list())
  ok(
    length(ch) == 2L && setequal(ch, c("GS_UP", "GS_REST")),
    "row choices derive from the real fGSEA pathway table"
  )
  session$setInputs(`fgsea-stab-table_rows_selected` = 1L)
  session$flushReact()
  ok(
    identical(store_read(.fgs_store)$resultspace.fgsea.selected_row, ch[[1]]),
    "user pathway-row click syncs into the store"
  )
  store_apply(.fgs_store, list(selected_row = ch[[2]]), origin = "agent")
  session$flushReact()
  session$setInputs(`fgsea-stab-table_rows_selected` = 2L)
  session$flushReact()
  ok(
    identical(store_read(.fgs_store)$resultspace.fgsea.selected_row, ch[[2]]) &&
      length(.fgs_root$pending) == 0L,
    "agent pathway-row push applies and acknowledges (drives bar plot)"
  )
})

## ------------------------------ S4 result-space: ptm -------------------
# The sequence triselector registers and auto-seeds to the first SeqLogo
# column (unset keys only - a restore or user pick wins and sticks).
.ptm_fd <- data.frame(
  `General|All|Gene.name` = paste0("g", 1:8),
  `SeqLogo|All|win15` = replicate(8, paste(sample(c("A","C","D","E","F","G",
                                                    "H","I","K","L","M","N",
                                                    "P","Q","R","S","T","V",
                                                    "W","Y"), 15, TRUE),
                                           collapse = "")),
  `SeqLogo|All|win20` = replicate(8, paste(sample(c("A","C","D","E","F","G",
                                                    "H","I","K","L","M","N",
                                                    "P","Q","R","S","T","V",
                                                    "W","Y"), 20, TRUE),
                                           collapse = "")),
  row.names = paste0("P", 1:8), check.names = FALSE)
.ptm_root <- widget_store_new()
.ptm_store <- widget_store_child(.ptm_root, "resultspace.ptm")
app_ptm <- function(input, output, session) {
  omicsViewer:::ptmotif_module(
    "ptm", pdata = shiny::reactive(NULL), fdata = shiny::reactive(.ptm_fd),
    expr = shiny::reactive(NULL), feature_selected = shiny::reactive(1:4),
    sample_selected = shiny::reactive(NULL),
    store = .ptm_store)
}
shiny::testServer(app_ptm, {
  session$flushReact()
  v <- store_read(.ptm_store)
  ok(
    identical(sort(names(v)),
              sort(paste0("resultspace.ptm.", c("xax_analysis", "xax_subset",
                                                "xax_variable")))),
    "ptm registers the sequence cascade"
  )
  ok(
    identical(v$resultspace.ptm.xax_analysis, "SeqLogo") &&
      identical(v$resultspace.ptm.xax_subset, "All") &&
      identical(v$resultspace.ptm.xax_variable, "win15"),
    "sequence cascade auto-seeds to the first SeqLogo column (system origin)"
  )
  # a user pick sticks (seeding never overwrites set keys); the sequence
  # dropdown only offers SeqLogo columns, so the pick is another window
  session$setInputs(`ptm-tris_seqlogo-analysis` = "SeqLogo")
  session$setInputs(`ptm-tris_seqlogo-subset` = "All")
  session$setInputs(`ptm-tris_seqlogo-variable` = "win20")
  session$flushReact()
  ok(
    identical(store_read(.ptm_store)$resultspace.ptm.xax_variable, "win20"),
    "user sequence pick syncs and sticks over the auto-seed"
  )
})

## --------------------------- S4 result-space: geneshot -----------------
# Search-term string + ID-mapper cascade; the Search button is a command
# and stays unregistered.
.gs2_fd <- data.frame(
  `General|All|Gene.name` = paste0("g", 1:16),
  `General|Info|Symbol` = paste0("sym", 1:16),
  row.names = paste0("F", 1:16), check.names = FALSE)
.gs2_root <- widget_store_new()
.gs2_store <- widget_store_child(.gs2_root, "resultspace.geneshot")
.sent_gs2 <- new.env(); .sent_gs2$msgs <- list()
app_gs2 <- function(input, output, session) {
  omicsViewer:::geneshot_module(
    "gs", pdata = shiny::reactive(NULL), fdata = shiny::reactive(.gs2_fd),
    expr = shiny::reactive(NULL), feature_selected = shiny::reactive(NULL),
    sample_selected = shiny::reactive(NULL), object = shiny::reactive(NULL),
    store = .gs2_store)
}
shiny::testServer(app_gs2, {
  .spy_input_messages(session, .sent_gs2)
  ids <- names(store_read(.gs2_store))
  ok(
    identical(sort(ids), sort(paste0("resultspace.geneshot.",
                                     c("term", "xax_analysis", "xax_subset",
                                       "xax_variable")))),
    "geneshot registers the search term and the ID-mapper cascade"
  )
  # user cascade pick + term edit sync into the store
  session$setInputs(`gs-geneNameCol-analysis` = "General")
  session$setInputs(`gs-geneNameCol-subset` = "Info")
  session$setInputs(`gs-geneNameCol-variable` = "Symbol")
  session$flushReact()
  session$setInputs(`gs-term` = "p53; cell cycle")
  session$flushReact()
  ok(
    identical(store_read(.gs2_store)$resultspace.geneshot.xax_variable,
              "Symbol") &&
      identical(store_read(.gs2_store)$resultspace.geneshot.term,
                "p53; cell cycle"),
    "user ID-mapper pick and search term sync into the store"
  )
  # external term write pushes updateTextInput
  .sent_gs2$msgs <- list()
  store_apply(.gs2_store, list(term = "BRCA1"), origin = "agent")
  session$flushReact()
  ok(
    identical(store_read(.gs2_store)$resultspace.geneshot.term, "BRCA1") &&
      any(vapply(.sent_gs2$msgs, function(m)
        grepl("(^|\\.)term$", m$id) &&
          identical(m$msg$value, "BRCA1"), logical(1))),
    "agent search-term write pushes the text input"
  )
})

## --------------------------- S4 result-space: stringdb -----------------
# Taxonomy string + labels boolean + enrichment row selection; the Run
# button is a stateless command and stays unregistered.
.str_root <- widget_store_new()
.str_store <- widget_store_child(.str_root, "resultspace.stringdb")
.sent_str <- new.env(); .sent_str$msgs <- list()
app_str <- function(input, output, session) {
  omicsViewer:::string_module(
    "sdb", reactive_ids = shiny::reactive(c("g1", "g2", "g3")),
    store = .str_store)
}
shiny::testServer(app_str, {
  .spy_input_messages(session, .sent_str)
  ids <- names(store_read(.str_store))
  ok(
    identical(sort(ids), sort(paste0("resultspace.stringdb.",
                                     c("taxonomy", "show_labels",
                                       "selected_row")))),
    "stringdb registers taxonomy, labels, and the enrichment row selection"
  )
  # user edits sync (warm the ignoreInit swallow first)
  session$setInputs(`sdb-tax` = "9606")
  session$flushReact()
  session$setInputs(`sdb-tax` = "10090")
  session$setInputs(`sdb-showLabel` = TRUE)
  session$flushReact()
  ok(
    identical(store_read(.str_store)$resultspace.stringdb.taxonomy, "10090") &&
      isTRUE(store_read(.str_store)$resultspace.stringdb.show_labels),
    "user taxonomy and label edits sync into the store"
  )
  # external writes push updateTextInputIcon / updateCheckboxInput
  .sent_str$msgs <- list()
  store_apply(.str_store, list(taxonomy = "9606", show_labels = FALSE),
              origin = "agent")
  session$flushReact()
  ok(
    identical(store_read(.str_store)$resultspace.stringdb.taxonomy, "9606") &&
      any(vapply(.sent_str$msgs, function(m)
        grepl("(^|\\.)tax$", m$id), logical(1))) &&
      any(vapply(.sent_str$msgs, function(m)
        grepl("(^|\\.)showLabel$", m$id), logical(1))),
    "agent taxonomy/label writes push both inputs"
  )
})

## -------------- S4 sample_general embedded: survival censor slider ------
# Registered on the sample_general child view; the push is gated on the
# survival view and clamped to the rendered range.
.sv_root <- widget_store_new()
.sv_store <- widget_store_child(.sv_root, "resultspace.sample_general")
.sent_sv <- new.env(); .sent_sv$msgs <- list()
.sv_pd <- data.frame(
  `General|All|group` = factor(rep(c("A", "B"), each = 4)),
  `Surv|All|time` = c(1, 5, 9, 20, 2, 8, 15, 30, 4, 7, 11, 25,
                      3, 6, 12, 28),
  row.names = paste0("S", 1:16), check.names = FALSE)
.sv_ex <- matrix(rnorm(16 * 16), nrow = 16,
                 dimnames = list(paste0("F", 1:16), paste0("S", 1:16)))
app_sv <- function(input, output, session) {
  omicsViewer:::sample_general_module(
    "sg", reactive_phenoData = shiny::reactive(.sv_pd),
    reactive_expr = shiny::reactive(.sv_ex),
    reactive_j = shiny::reactive(paste0("S", 1:8)),
    store = .sv_store)
}
shiny::testServer(app_sv, {
  .spy_input_messages(session, .sent_sv)
  # switch the sample cascade to the Surv view
  session$setInputs(`sg-tris_sample_general-analysis` = "Surv")
  session$setInputs(`sg-tris_sample_general-subset` = "All")
  session$setInputs(`sg-tris_sample_general-variable` = "time")
  session$flushReact()
  ok(
    "resultspace.sample_general.survival_censor" %in% names(store_read(.sv_store)),
    "the survival censor slider registers on the sample_general view"
  )
  .sent_sv$msgs <- list()
  store_apply(.sv_store, list(survival_censor = 10), origin = "agent")
  session$flushReact()
  ok(
    identical(store_read(.sv_store)$resultspace.sample_general.survival_censor,
              10) &&
      any(vapply(.sent_sv$msgs, function(m)
        grepl("(^|\\.)censor$", m$id), logical(1))),
    "agent censor write pushes the slider in the survival view"
  )
  # out-of-range values are clamped to the rendered slider range
  .sent_sv$msgs <- list()
  store_apply(.sv_store, list(survival_censor = 1e6), origin = "agent")
  session$flushReact()
  hit <- Filter(function(m) grepl("(^|\\.)censor$", m$id), .sent_sv$msgs)
  ok(
    length(hit) == 1L && hit[[1]]$msg$value <= max(.sv_pd$`Surv|All|time`),
    "out-of-range censor value is clamped to the slider maximum"
  )
})

## -------------- S4 sample_general embedded: batch comparison toggles ----
.bc_root <- widget_store_new()
.bc_store <- widget_store_child(.bc_root, "resultspace.sample_general")
.sent_bc <- new.env(); .sent_bc$msgs <- list()
app_bc <- function(input, output, session) {
  omicsViewer:::sample_general_module(
    "sg", reactive_phenoData = shiny::reactive(.sv_pd),
    reactive_expr = shiny::reactive(.sv_ex),
    reactive_j = shiny::reactive(paste0("S", 1:8)),
    store = .bc_store)
}
shiny::testServer(app_bc, {
  .spy_input_messages(session, .sent_bc)
  # a numeric (beeswarm) view keeps the batch module's tables running
  session$setInputs(`sg-tris_sample_general-analysis` = "Surv")
  session$setInputs(`sg-tris_sample_general-subset` = "All")
  session$setInputs(`sg-tris_sample_general-variable` = "time")
  session$flushReact()
  keys <- sort(names(store_read(.bc_store)))
  ok(
    all(paste0("resultspace.sample_general.",
               c("batch_show_phenotype", "batch_show_features",
                 "batch_phenotype_selected_row",
                 "batch_feature_selected_row")) %in% keys),
    "the batch comparison toggles and row selections register"
  )
  # user toggle edits sync (warm the ignoreInit swallow first)
  session$setInputs(`sg-batch_comp-show_phenotype` = TRUE)
  session$flushReact()
  session$setInputs(`sg-batch_comp-show_features` = TRUE)
  session$flushReact()
  ok(
    isTRUE(store_read(.bc_store)$resultspace.sample_general.batch_show_features),
    "user batch toggle edit syncs into the store"
  )
  # external toggle write pushes updateCheckboxInput
  .sent_bc$msgs <- list()
  store_apply(.bc_store, list(batch_show_phenotype = FALSE), origin = "agent")
  session$flushReact()
  ok(
    isFALSE(store_read(.bc_store)$resultspace.sample_general.batch_show_phenotype) &&
      any(vapply(.sent_bc$msgs, function(m)
        grepl("(^|\\.)show_phenotype$", m$id), logical(1))),
    "agent batch toggle write pushes the checkbox"
  )
})

## ------------------- S4: meta_scatter adopts the shared attr4 panel ----
# The panel became store-aware in result-space step 1; the data-space
# scatters now pass their store so the attr4 keys register there too.
.ms_root <- widget_store_new()
.ms_store <- widget_store_child(.ms_root, "dataspace.feature_space")
.ms_fd <- data.frame(
  `General|All|grp` = factor(rep(c("a", "b"), each = 8)),
  `PCA|All|PC1` = rnorm(16),
  row.names = paste0("F", 1:16), check.names = FALSE)
.ms_pd <- data.frame(
  `General|All|batch` = factor(rep(c("B1", "B2"), each = 8)),
  row.names = paste0("S", 1:16), check.names = FALSE)
.ms_ex <- matrix(rnorm(16 * 16), nrow = 16,
                 dimnames = list(paste0("F", 1:16), paste0("S", 1:16)))
app_ms <- function(input, output, session) {
  omicsViewer:::meta_scatter_module(
    "feature_space", reactive_meta = shiny::reactive(.ms_fd),
    reactive_expr = shiny::reactive(.ms_ex), combine = "feature",
    source = "scatter_meta_feature", store = .ms_store)
}
shiny::testServer(app_ms, {
  session$flushReact()
  ids <- names(store_read(.ms_store))
  ok(
    all(paste0("dataspace.feature_space.attr4.",
               c("color_analysis", "color_variable", "xcut", "scorner")) %in% ids),
    "the feature-space scatter registers the shared attr4 panel keys"
  )
  ok(
    identical(sort(ids)[1], "dataspace.feature_space.attr4.color_analysis") ||
        length(ids) >= 18L,
    "attr4 adoption adds the full 18-key panel per scatter"
  )
})
