# WP-C unit suite for the action tier: widget_action / store_register_action /
# store_action_view / store_action_run (auxi_widgetStore.R section [7]),
# agent_action_apply + policy exclusions (auxi_agentActions.R), and the
# action capability records (auxi_agentCapabilities.R).
# Pure protocol logic: no browser, no Shiny session required.

library(omicsViewer)
library(unittest, quietly = TRUE)

widget_store_new <- omicsViewer:::widget_store_new
widget_store_child <- omicsViewer:::widget_store_child
widget_binding <- omicsViewer:::widget_binding
store_register <- omicsViewer:::store_register
store_register_action <- omicsViewer:::store_register_action
widget_action <- omicsViewer:::widget_action
store_action_view <- omicsViewer:::store_action_view
store_action_run <- omicsViewer:::store_action_run
store_snapshot <- omicsViewer:::store_snapshot
store_restore <- omicsViewer:::store_restore
store_apply <- omicsViewer:::store_apply
agent_action_apply <- omicsViewer:::agent_action_apply
agent_capability_records <- omicsViewer:::agent_capability_records

## ---------------------------------------------------------------- [1] ----
# registration + run basics
mk <- function() {
  s <- widget_store_new()
  ran <- 0L
  store_register_action(s,
    widget_action("panel.search", label = "Search", help = "Run the search",
      handler = function(params) ran <<- ran + 1L))
  store_register_action(s,
    widget_action("panel.clear", label = "Clear", help = "Clear the panel",
      handler = function(params) NULL, input_id = "panelUI-clear"))
  list(store = s, ran = function() ran)
}

s1 <- mk()
ok(setequal(names(s1$store$actions), c("panel.search", "panel.clear")),
   "actions registered under canonical ids")
ok(identical(s1$store$action_epochs[["panel.search"]], 0L),
   "action epoch starts at 0")

r1 <- store_action_run(s1$store, "panel.search")
ok(identical(r1$ran, "panel.search"), "run receipt carries the id")
ok(identical(r1$epoch, 1L), "epoch bumps to 1 after first run")
ok(identical(s1$ran(), 1L), "handler body actually ran")
ok(identical(s1$store$action_epochs[["panel.search"]], 1L), "epoch recorded")
ok(identical(s1$store$action_runs[["panel.search"]]$origin, "agent"),
   "default origin is agent")
ok(nzchar(s1$store$action_runs[["panel.search"]]$timestamp),
   "last-run timestamp recorded")
store_action_run(s1$store, "panel.search", origin = "user")
ok(identical(s1$store$action_runs[["panel.search"]]$origin, "user"),
   "origin recorded per run")
ok(identical(s1$store$action_epochs[["panel.search"]], 2L),
   "epoch accumulates")

## ---------------------------------------------------------------- [2] ----
# validation + unknown ids
s2 <- mk()
ok(grepl("Unknown action id", tryCatch(
  store_action_run(s2$store, "panel.nothing"), error = conditionMessage)),
  "unknown action id errors")
ok(grepl("panel.search", tryCatch(
  store_action_run(s2$store, "panel.serch"), error = conditionMessage)),
  "unknown id error carries a closest-match suggestion")
ok(grepl("dotted", tryCatch(
  store_register_action(s2$store,
    widget_action("bad id", handler = function(p) NULL)),
  error = conditionMessage)), "malformed id rejected at construction")
ok(grepl("function", tryCatch(
  store_register_action(s2$store,
    widget_action("panel.x", handler = "nope")),
  error = conditionMessage)), "non-function handler rejected")
ok(grepl("already registered", tryCatch(
  store_register_action(s2$store,
    widget_action("panel.search", handler = function(p) NULL)),
  error = conditionMessage)), "duplicate action id rejected")
s2b <- widget_store_new()
store_register(s2b, widget_binding("panel.clear", "string"))
ok(grepl("widget binding", tryCatch(
  store_register_action(s2b,
    widget_action("panel.clear", handler = function(p) NULL)),
  error = conditionMessage)),
  "an id cannot be both a binding and an action")

## ---------------------------------------------------------------- [3] ----
# destructive + confirm gating + availability
s3 <- widget_store_new()
store_register_action(s3,
  widget_action("panel.wipe", label = "Wipe", handler = function(p) NULL,
    destructive = TRUE),
  widget_action("panel.reset", label = "Reset", handler = function(p) NULL,
    requires_confirm = TRUE),
  widget_action("panel.lazy", label = "Lazy", handler = function(p) NULL,
    availability = function() "needs a dataset"))
ok(grepl("destructive", tryCatch(
  store_action_run(s3, "panel.wipe", origin = "agent"),
  error = conditionMessage)),
  "agent origin cannot run destructive actions")
ok(isTRUE(store_action_run(s3, "panel.wipe", origin = "user")$epoch >= 1L),
  "user origin may run destructive actions (package-author decision)")
ok(grepl("confirm", tryCatch(
  store_action_run(s3, "panel.reset"), error = conditionMessage)),
  "requires_confirm gates an unconfirmed run")
ok(isTRUE(store_action_run(s3, "panel.reset",
  params = list(confirm = TRUE))$epoch >= 1L),
  "confirmed run passes the gate")
ok(grepl("not currently available", tryCatch(
  store_action_run(s3, "panel.lazy"), error = conditionMessage)),
  "unavailable action rejects with the reason")
ok(grepl("needs a dataset", tryCatch(
  store_action_run(s3, "panel.lazy"), error = conditionMessage)),
  "availability reason string is carried into the rejection")

## ---------------------------------------------------------------- [4] ----
# agent-facing view never exposes destructive actions
v3 <- store_action_view(s3)
ok(setequal(vapply(v3, function(a) a$id, character(1)),
            c("panel.reset", "panel.lazy")),
   "store_action_view hides destructive actions")
av <- setNames(vapply(v3, function(a) as.character(a$available), character(1)),
               vapply(v3, function(a) a$id, character(1)))
ok(identical(av[["panel.lazy"]], "needs a dataset"),
   "view carries the availability reason")
ok(identical(av[["panel.reset"]], "TRUE"),
   "view carries available = TRUE")

## ---------------------------------------------------------------- [5] ----
# child-view prefixing
s5root <- widget_store_new()
s5 <- widget_store_child(s5root, "dataspace.tab_pheno")
store_register_action(s5,
  widget_action("clear", handler = function(p) NULL, input_id = "tab_phenoUI-clear"))
ok("dataspace.tab_pheno.clear" %in% names(s5root$actions),
   "child registration prefixes to the canonical namespace")
ok(isTRUE(store_action_run(s5, "clear")$epoch >= 1L),
   "module-local id resolves through the child view")
ok(identical(store_action_view(s5root, "dataspace.tab_pheno")[[1]]$input_id,
             "tab_phenoUI-clear"),
   "input_id survives into the view (dom mapping)")

## ---------------------------------------------------------------- [6] ----
# snapshot carries last-run epochs only; restore never replays
s6 <- mk()
store_action_run(s6$store, "panel.search")
snap <- store_snapshot(s6$store)
ok(identical(snap$actions[["panel.search"]], 1L),
   "snapshot carries the last-run epoch")
ok(is.null(snap$values[["panel.search"]]),
   "snapshot carries no action VALUE (commands have no value)")
before <- s6$ran()
store_restore(s6$store, snap)
ok(identical(s6$ran(), before),
   "restore never replays actions")

## ---------------------------------------------------------------- [7] ----
# agent_action_apply: receipts, policy exclusions, suggestions
s7 <- mk()
r7 <- agent_action_apply(s7$store, "panel.search")
ok(identical(r7$ran, "panel.search") && is.null(r7$error),
   "apply receipt reports a successful run")
ok(identical(r7$epoch, 1L), "apply receipt carries the epoch")
r7b <- agent_action_apply(s7$store, "panel.serch")
ok(length(r7b$ran) == 0L && grepl("panel.search", r7b$error),
   "unknown id returns a suggestion-bearing error receipt")
r7c <- agent_action_apply(s7$store, "action:panel.clear")
ok(identical(r7c$ran, "panel.clear"),
   "the action: record-id prefix is accepted and stripped")

ex <- agent_action_apply(s7$store, "snapshot_save")
ok(length(ex$ran) == 0L && isTRUE(ex$excluded_by_policy),
   "snapshot-lifecycle id rejected as excluded by policy")
ok(grepl("permanently excluded", ex$error),
   "policy rejection names the policy")
for (id in c("app-snapshot_save", "app-selectFile", "snapshot",
             "action:app-snapshot_delete", "dataset_selector"))
  ok(isTRUE(agent_action_apply(s7$store, id)$excluded_by_policy),
     paste("policy exclusion matches normalised form:", id))
ok(!isTRUE(agent_action_apply(s7$store, "panel.search")$excluded_by_policy),
   "a normal action is not flagged as excluded")
ok(grepl("non-empty", agent_action_apply(s7$store, "")$error),
   "empty action_id rejected with a self-describing reason")

## ---------------------------------------------------------------- [8] ----
# capability records gain action records
s8 <- mk()
store_register_action(s8$store,
  widget_action("panel.wipe", handler = function(p) NULL, destructive = TRUE))
recs <- agent_capability_records(s8$store)
acts <- Filter(function(r) identical(r$kind, "action"), recs)
ok(setequal(vapply(acts, function(a) a$id, character(1)),
            c("action:panel.search", "action:panel.clear")),
   "action records use the action: id prefix")
ok(!any(vapply(acts, function(a) a$id, character(1)) == "action:panel.wipe"),
   "destructive actions never enter capability records")
ok(all(vapply(acts, function(a) a$operation, character(1)) == "run_control"),
   "action records carry operation run_control")
ok(identical(acts[[2]]$input_id, "panelUI-clear"),
   "action records carry the dom input mapping")

## ---------------------------------------------------------------- [9] ----
# agent-facing availability is read through the receipt when unavailable
s9 <- widget_store_new()
store_register_action(s9,
  widget_action("panel.lazy", handler = function(p) NULL,
    availability = function() "needs a dataset"))
r9 <- agent_action_apply(s9, "panel.lazy")
ok(length(r9$ran) == 0L && grepl("needs a dataset", r9$error),
   "apply surfaces the availability reason to the model")
