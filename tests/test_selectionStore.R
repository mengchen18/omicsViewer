library(omicsViewer)
library(unittest, quietly = TRUE)

# ---------------------------------------------------------------------------
# The unified selection bus (auxi_selectionStore.R): one canonical record per
# space; every selection source reports through it, every consumer watches it.
# Core protocol semantics are tested without a browser.
# ---------------------------------------------------------------------------
sb <- omicsViewer:::selection_store_new(c("feature", "sample"))
pf <- omicsViewer:::selection_port(sb, "feature")
ps <- omicsViewer:::selection_port(sb, "sample")

ok(omicsViewer:::selection_read(sb, "feature")$epoch == 0L, "a new store holds an empty record per key")
ok(identical(omicsViewer:::selection_read(sb, "feature")$mirror, TRUE), "the initial mirror is TRUE (no table filter)")
ok(identical(pf$read()$ids, character(0)), "the port reads the same record as the store")

# --- report dedupe: an unchanged report from the same origin is an echo ----
ok(isTRUE(pf$report(origin = "corner", report = list(r = 1), ids = c("a", "b"),
                    anchor = "x1")), "a corner report applies")
ok(identical(pf$read()$ids, c("a", "b")), "the record carries the reported ids")
ok(identical(pf$read()$origin, "corner"), "the record carries the origin")
ok(isFALSE(pf$report(origin = "corner", report = list(r = 1), ids = c("a", "b"),
                     anchor = "x1")), "an identical report from the same origin is skipped")
e0 <- pf$read()$epoch
pf$report(origin = "corner", report = list(r = 1), ids = c("a", "b"), anchor = "x1")
ok(pf$read()$epoch == e0, "an echoed report does not bump the epoch")

# a DIFFERENT origin may report with the same ids and wins (last genuine
# interaction wins - e.g. a table row click after a corner selection)
ok(isTRUE(pf$report(origin = "table", report = c("a"), ids = c("a"), mirror = NULL)),
   "a table report applies")
ok(identical(pf$read()$origin, "table"), "the last genuine report wins the record")
pf$report(origin = "corner", report = list(r = 2), ids = c("a", "b"),
          anchor = "x1", mirror = c("a", "b"))
ok(identical(pf$read()$mirror, c("a", "b")), "a corner report mirrors its ids into the tables")
pf$report(origin = "table", report = c("a"), ids = c("a"), mirror = NULL)
ok(identical(pf$read()$mirror, c("a", "b")), "a table report leaves the mirror unchanged")

# --- the mirror policy ------------------------------------------------------
pf$report(origin = "figure", report = list(sel = c("c")), ids = c("c"),
          anchor = "x2", mirror = c("c"))
pf$report(origin = "clear", report = list(src = "clear"), ids = character(0),
          mirror = TRUE)
ok(identical(pf$read()$mirror, TRUE), "a clear report un-filters the tables")
ok(length(pf$read()$ids) == 0L, "a clear report empties the selection")

# --- apply dedupe and NULL-means-unchanged fields ---------------------------
ok(isTRUE(pf$apply(ids = c("d"), origin = "restore", mirror = c("d"))),
   "a direct apply applies")
ok(isFALSE(pf$apply(ids = c("d"), origin = "restore", mirror = c("d"))),
   "an identical apply is a no-op")
pf$apply(ids = c("d", "e"), origin = "restore")
ok(identical(pf$read()$ids, c("d", "e")) &&
     identical(pf$read()$mirror, c("d")) &&
     identical(pf$read()$origin, "restore"),
   "apply with NULL anchor/mirror/clicked leaves those fields unchanged")

# --- keys are independent ----------------------------------------------------
ps$apply(ids = c("S1"), origin = "heatmap", mirror = c("S1"))
ok(identical(pf$read()$ids, c("d", "e")), "a sample-space write leaves the feature record alone")
ok(identical(ps$read()$ids, "S1"), "the sample record holds its own ids")

# --- snapshot state ----------------------------------------------------------
st <- omicsViewer:::selection_state(sb)
ok(identical(sort(names(st)), c("feature", "sample")), "selection_state covers every key")
ok(identical(st$sample$origin, "heatmap"), "selection_state carries origins")

# --- watch: reactive re-derivation ------------------------------------------
watch_log <- new.env(); watch_log$n <- 0L; watch_log$last <- NULL
shiny::testServer(function(input, output, session) {
  pf$apply(ids = c("z"), origin = "system")
}, {
  w <- pf$watch()
  shiny::observe({
    watch_log$n <- watch_log$n + 1L
    watch_log$last <- w()$ids
  })
  pf$apply(ids = c("f"), origin = "system")
  session$flushReact()
})
ok(watch_log$n >= 1L && identical(watch_log$last, "f"),
   "selection_watch re-derives on an applied change")

# --- origin validation -------------------------------------------------------
res <- tryCatch(pf$apply(ids = "x", origin = "not-an-origin"),
                error = function(e) e)
ok(inherits(res, "error"), "an unknown origin is rejected")

