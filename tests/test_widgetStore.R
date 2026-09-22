# S1 unit suite for the universal widget control plane (auxi_widgetStore.R).
# Pure protocol logic: no browser, no Shiny session required.

library(omicsViewer)
library(unittest, quietly = TRUE)

widget_store_new <- omicsViewer:::widget_store_new
widget_store_child <- omicsViewer:::widget_store_child
widget_binding <- omicsViewer:::widget_binding
store_register <- omicsViewer:::store_register
store_apply <- omicsViewer:::store_apply
store_read <- omicsViewer:::store_read
store_ack <- omicsViewer:::store_ack
store_sync_from_ui <- omicsViewer:::store_sync_from_ui
store_snapshot <- omicsViewer:::store_snapshot
store_restore <- omicsViewer:::store_restore
store_registry_view <- omicsViewer:::store_registry_view
store_describe <- omicsViewer:::store_describe

## ---------------------------------------------------------------- [1] ----
mk <- function() {
  s <- widget_store_new()
  store_register(s,
    widget_binding("app.main_tab", "navbar", label = "Data-space tab",
      help = "Visible data-space tab",
      choices_provider = function(v) c("Feature", "Sample", "Heatmap")),
    widget_binding("app.y_axis", "select_cascaded", label = "Y axis",
      help = "Y-axis annotation variable",
      depends_on = c("app.y_analysis", "app.y_subset"),
      choices_provider = function(v) c("log.fdr", "log.pvalue")),
    widget_binding("app.y_analysis", "select", label = "Y analysis",
      help = "Y-axis analysis",
      choices_provider = function(v) c("ttest", "PCA")),
    widget_binding("app.y_subset", "select", label = "Y subset",
      help = "Y-axis subset", depends_on = "app.y_analysis",
      choices_provider = function(v) c("RE_vs_ME", "All")),
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

ok(
  ut_cmp_error(store_register(mk(), widget_binding("app.y_axis", "select")),
    "already registered"),
  "duplicate ids are rejected"
)
ok(
  ut_cmp_error(
    store_register(widget_store_new(),
      widget_binding("a.b", "select", depends_on = "a.missing")),
    "unregistered id"),
  "unknown dependencies are rejected"
)
ok(
  ut_cmp_error(
    store_register(widget_store_new(),
      widget_binding("a.x", "select", depends_on = "a.y"),
      widget_binding("a.y", "select", depends_on = "a.x")),
    "cycle"),
  "dependency cycles are rejected at registration"
)
ok(
  ut_cmp_error(widget_binding("bad id!", "select"), "dotted identifier"),
  "malformed ids are rejected"
)

## ---------------------------------------------------------------- [2] ----
s <- mk()
ok(
  ut_cmp_error(store_apply(s, list(`app.no_such` = "x")), "Unknown widget id"),
  "unknown patch keys are rejected"
)
ok(
  ut_cmp_error(store_apply(s, list(app.theme = "dark")), "Unknown value"),
  "invalid enum values are rejected with guidance"
)
ok(
  grepl("minimal",
    tryCatch(store_apply(s, list(app.theme = "minmal")), error = function(e) conditionMessage(e))),
  "invalid values carry did-you-mean suggestions"
)
ok(
  ut_cmp_error(store_apply(s, list(app.y_axis = "log.fdrr")),
    "Unknown value"),
  "select values are validated against choices providers"
)
ok(
  ut_cmp_error(store_apply(s, list(app.min_size = 0L)), "must be >="),
  "numeric bounds are enforced"
)
ok(
  ut_cmp_error(store_apply(s, list(app.min_size = 1.5)), "requires an integer"),
  "integer coercion is strict"
)
ok(
  identical(store_apply(s, list(app.show_labels = "true")) |> invisible() |> capture.output(), character()),
  "string booleans are accepted"
)
ok(
  identical(store_read(s, "app.show_labels")$app.show_labels, TRUE),
  "string booleans normalise to logical"
)
ok(
  is.null(store_read(s, "app.y_axis")$app.y_axis),
  "sentinel values register as NULL"
)
r <- store_apply(s, list(app.y_axis = "null", app.theme = "{}"))
ok(
  identical(r$applied, character()),
  "literal null/{} sentinels are treated as omitted (no writes)"
)

## ---------------------------------------------------------------- [3] ----
s <- mk()
store_apply(s, list(app.main_tab = "Sample"))
e1 <- s$epochs$app.main_tab
r <- store_apply(s, list(app.main_tab = "Sample", app.theme = "minimal"))
ok(
  identical(r$applied, "app.theme") && identical(r$skipped, "app.main_tab"),
  "diff-only writes: unchanged keys are skipped"
)
ok(
  identical(s$epochs$app.main_tab, e1),
  "unchanged keys' epochs are not bumped"
)

# dependency ordering: y_axis is written after its dependencies
s <- mk()
r <- store_apply(s, list(app.y_axis = "log.fdr", app.y_subset = "All",
                         app.y_analysis = "PCA"))
ok(
  match("app.y_axis", r$applied) > match("app.y_subset", r$applied) &&
    match("app.y_subset", r$applied) > match("app.y_analysis", r$applied),
  "writes are ordered by the dependency graph"
)

## ---------------------------------------------------------------- [4] ----
s <- mk()
store_apply(s, list(app.theme = "classic"))
ok(
  identical(s$pending$app.theme$value, "classic"),
  "applies record pending (unacknowledged) values"
)
ok(
  identical(store_ack(s, "app.theme", "classic"), FALSE) &&
    is.null(s$pending$app.theme),
  "matching acknowledgement clears the pending entry"
)
s <- mk()
store_apply(s, list(app.theme = "classic"))
ok(
  identical(store_ack(s, "app.theme", "minimal"), TRUE),
  "mismatched widget reports signal a re-assert"
)

# user edits mirror into the store without epoch bumps, and win over pending
s <- mk()
store_apply(s, list(app.theme = "classic"))
e_apply <- s$epochs$app.theme
overridden <- store_sync_from_ui(s, "app.theme", "bw")
ok(
  isTRUE(overridden) && is.null(s$pending$app.theme),
  "manual user edit overrides an in-flight agent write and clears pending"
)
ok(
  identical(s$epochs$app.theme, e_apply) &&
    identical(store_read(s, "app.theme")$app.theme, "bw"),
  "user sync updates the value without bumping epochs"
)
ok(
  identical(s$origins$app.theme, "user"),
  "origins are tracked"
)
ok(
  length(s$override_log) == 1L && identical(s$override_log[[1]]$user, "bw"),
  "user overrides are recorded in the bounded log"
)

## ---------------------------------------------------------------- [5] ----
s <- mk()
store_apply(s, list(app.main_tab = "Heatmap", app.theme = "classic",
                    app.min_size = 15L), origin = "system")
store_apply(s, list(app.internal_counter = 3.5), origin = "restore")
ok(
  ut_cmp_error(store_apply(s, list(app.internal_counter = 9),
    origin = "agent"), "not user-editable"),
  "agent-origin writes cannot touch internal state"
)
snap <- store_snapshot(s)
ok(
  identical(snap$values$app.internal_counter, 3.5),
  "snapshots include internal state"
)
s2 <- mk()
receipt <- store_restore(s2, snap)
ok(
  setequal(receipt$applied,
    c("app.main_tab", "app.theme", "app.min_size", "app.internal_counter")),
  "restore-origin applies reach internal keys"
)
ok(
  identical(store_read(s2)$app.internal_counter, 3.5) &&
    identical(store_read(s2)$app.theme, "classic"),
  "snapshot round-trips exactly"
)
ok(
  identical(store_restore(s2, list(values = list(app.gone = 1)))$unknown_ids,
    "app.gone"),
  "unknown snapshot ids are reported, not applied"
)

## ---------------------------------------------------------------- [6] ----
s <- mk()
view <- store_registry_view(s)
ids <- vapply(view, function(r) r$id, "")
ok(
  "app.internal_counter" %in% ids == FALSE &&
    setequal(ids, setdiff(names(s$bindings), "app.internal_counter")),
  "registry exposes exactly the user-editable widgets"
)
ok(
  identical(store_describe(s, "app.internal_counter"), NULL),
  "internal widgets cannot be described to the agent"
)
d <- store_describe(s, "app.y_axis")
ok(
  identical(d$depends_on, c("app.y_analysis", "app.y_subset")) &&
    "log.fdr" %in% d$allowed_values,
  "descriptions carry dependencies and allowed values"
)

# child views namespace registrations and applies
s <- mk()
child <- widget_store_child(s, "dataspace.feature_space")
store_register(child, widget_binding("x_axis", "select",
  label = "X axis", help = "X-axis variable",
  choices_provider = function(v) c("mean.diff", "log.fdr")))
r <- store_apply(child, list(x_axis = "mean.diff"))
ok(
  identical(store_read(s, "dataspace.feature_space.x_axis")$
              dataspace.feature_space.x_axis, "mean.diff"),
  "child views prefix ids for registration and applies"
)

# global epoch advances once per transaction (future state token)
s <- mk()
g0 <- s$global_epoch
store_apply(s, list(app.theme = "classic", app.min_size = 10L))
ok(
  identical(s$global_epoch, g0 + 1L),
  "the global epoch advances once per transaction"
)

# cascaded patches validate against the in-patch overlay: a jointly valid
# analysis/subset/variable triple installs in ONE transaction even though
# the subset is invalid under the *current* analysis
sc <- widget_store_new()
store_register(sc,
  widget_binding("a.analysis", "select", label = "A", help = "a",
    choices_provider = function(v) c("ttest", "PCA")),
  widget_binding("a.subset", "select", label = "S", help = "s",
    depends_on = "a.analysis",
    choices_provider = function(v)
      if (identical(v[["a.analysis"]], "PCA")) c("All", "removeMissing") else c("RE_vs_ME", "All")),
  widget_binding("a.variable", "select_cascaded", label = "V", help = "v",
    depends_on = c("a.analysis", "a.subset"),
    choices_provider = function(v)
      if (identical(v[["a.analysis"]], "PCA") && identical(v[["a.subset"]], "All"))
        c("PC1", "PC2") else c("mean.diff", "log.fdr")))
r <- store_apply(sc, list(a.analysis = "PCA", a.subset = "All", a.variable = "PC1"))
ok(
  ut_cmp_identical(r$applied, c("a.analysis", "a.subset", "a.variable")),
  "cascaded triples validate jointly and apply in dependency order"
)
ok(
  ut_cmp_error(store_apply(sc, list(a.subset = "removeMissing", a.variable = "PC2")),
    "Unknown value"),
  "cascaded validation still rejects values invalid under current state"
)

# ack-aware user sync: confirming an in-flight write is not an override
s <- mk()
store_apply(s, list(app.theme = "classic"))
ov <- store_sync_from_ui(s, "app.theme", "classic")
ok(
  isFALSE(ov) && is.null(s$pending$app.theme) &&
    identical(s$origins$app.theme, "agent"),
  "widget confirming a pending value acks it without logging an override"
)

## ------------------------------------------------- multi_select kind ----
# S4: multi-selection widgets (heatmap annotations/tooltips). Value is a
# character vector; empty vectors, empty lists, and "" all clear; entries
# must all be allowed; sentinel strings stay omitted optionals.
ms <- widget_store_new()
store_register(ms,
  widget_binding("hm.annot_col", "multi_select", label = "Annotations",
    help = "Annotation columns",
    choices_provider = function(v) c("group", "batch", "stage")),
  widget_binding("hm.annot_row", "multi_select", label = "Row annotations",
    help = "Row annotation columns",
    choices_provider = function(v) c("score", "pathway")))
r <- store_apply(ms, list(hm.annot_col = c("batch", "group")))
ok(
  identical(r$applied, "hm.annot_col") &&
    identical(store_read(ms)$hm.annot_col, c("batch", "group")),
  "multi_select applies a character vector"
)
ok(
  ut_cmp_error(store_apply(ms, list(hm.annot_col = c("group", "celtype"))),
    "Unknown value"),
  "multi_select rejects values outside the allowed set"
)
r <- store_apply(ms, list(hm.annot_col = c("group", "celtype")),
                 strict = FALSE)
ok(
  identical(r$applied, character()) && length(r$rejected) == 1L &&
    grepl("celtype", r$rejected[[1]]$reason),
  "multi_select rejections are per key under strict=FALSE"
)
# JSON arrays arrive as lists (simplifyVector = FALSE); "" clears
r <- store_apply(ms, list(hm.annot_col = list("stage")))
ok(
  identical(store_read(ms)$hm.annot_col, "stage"),
  "multi_select unlists JSON array values"
)
r <- store_apply(ms, list(hm.annot_col = ""))
ok(
  identical(r$applied, "hm.annot_col") &&
    identical(store_read(ms)$hm.annot_col, character(0)),
  "an empty string clears a multi_select"
)
r <- store_apply(ms, list(hm.annot_col = list()))
ok(
  identical(store_read(ms)$hm.annot_col, character(0)),
  "an empty JSON array clears a multi_select"
)
# sentinel strings are omitted optionals, not writes
r <- store_apply(ms, list(hm.annot_col = "[]"))
ok(
  identical(r$applied, character()) &&
    identical(store_read(ms)$hm.annot_col, character(0)),
  "a literal '[]' sentinel is dropped like an omitted optional"
)
# user sync from a cleared multi-select input (character(0))
store_sync_from_ui(ms, "hm.annot_row", c("score", "pathway"))
store_sync_from_ui(ms, "hm.annot_row", character(0))
ok(
  identical(store_read(ms)$hm.annot_row, character(0)),
  "clearing a multi-select input syncs character(0) into the store"
)
# ack semantics on vectors
store_apply(ms, list(hm.annot_row = c("score")))
ok(
  isFALSE(store_ack(ms, "hm.annot_row", c("score"))),
  "acks compare multi_select vectors by identity"
)
# snapshot round-trip keeps vectors (and empties) intact
snap <- store_snapshot(ms)
ms2 <- widget_store_new()
store_register(ms2,
  widget_binding("hm.annot_col", "multi_select", label = "A", help = "a",
    choices_provider = function(v) c("group", "batch", "stage")),
  widget_binding("hm.annot_row", "multi_select", label = "R", help = "r",
    choices_provider = function(v) c("score", "pathway")))
store_restore(ms2, snap)
vals <- store_read(ms2)
ok(
  identical(vals$hm.annot_col, character(0)) &&
    identical(vals$hm.annot_row, "score"),
  "multi_select values survive the snapshot round-trip"
)
# registry view exposes the kind and allowed values
view <- Filter(function(r) identical(r$id, "hm.annot_col"),
               store_registry_view(ms))[[1]]
ok(
  identical(view$kind, "multi_select") &&
    identical(view$allowed_values, c("group", "batch", "stage")),
  "registry view reports multi_select bindings with allowed values"
)
