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
