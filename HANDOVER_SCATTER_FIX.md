# HANDOVER — scatter quick-view switch multi-flash fix (DONE, 2026-09-24)

**Status: complete and verified.** Bug 1 (volcano corner on load) stayed
fixed; the quick-view switch multi-render (volcano<->cor AND
volcano->volcano) is fixed; the infinite reactive loop that the previous
attempt left in the working tree is gone. Everything below is for the
record — read `AGENT_ACCURACY_PLAN.md` §6 for the control-plane
architecture and `HANDOVER.md` §"Scatter volcano corner" for the earlier
history.

## Root cause of the previous attempt's infinite loop (diagnosed live)

Reproduced with printf tracing (message() counters) + websocket tracing
(`options(shiny.trace = TRUE)`, see the Shiny "Debugging applications"
article). Two observers disagreed forever about the corner:

- the side-effect observer computed the corner from the pending intent
  (`"volcano"`), while the generic input observer mirrored the raw widget
  report (`input$scorner` = `"None"`) — both wrote the shared
  `.a4_cutoff_key` + `params$cutoff`, each write re-invalidating the other
  (~250 iterations/s);
- `corner_effective` READ `params$cutoff` while the side-effect observer
  WROTE it — a reactiveValues write always re-invalidates readers (no
  diffing), so the graph contained a permanent read/write feedback edge;
- the intended convergence (browser acks the `updateSelectInput` push,
  both observers agree) could never happen: while the flush loop spins,
  Shiny never drains the outgoing message queue (12 s of spinning, 47
  SEND messages, zero scorner updates), so the ack starved — a
  self-sustaining loop.

## The fix (constraints 1–3 of the previous handover, all satisfied)

In `R/module_figureAttr4.R` (attr4selector_module):

- `corner_effective` resolves as: gated intent (`pendingCorner` +
  `corner_apply_gate`) %||% `widget_corner()` — it no longer reads
  `params$cutoff`, so the feedback edge is gone.
- `widget_corner()` is the ack-filtered widget state: while a pushed
  corner value awaits the browser ack, a differing widget report is stale
  and the pushed value wins (same pattern as `pendingAnalysis` in
  module_triselector.R). The choices-rebuild observer respects the pushed
  value too, and the push sends `choices` together with `selected`
  (selectize silently drops a setValue for a value not among the options).
- The side-effect observer is the ONLY writer of `params$cutoff`, the
  store and the scorner widget, content-deduped via `.a4_cutoff_key`
  (the historical `attr(l, "seed") <- Sys.time()` full-redraw seed is
  gone). The generic input observer now only mirrors widget state into
  the status snapshot reactiveVals.
- Any scorner widget report retires the auto-intent (`observeEvent(
  input$scorner)`) — an ack keeps the value via `widget_corner`, a
  user/restore choice always wins over the automatic volcano corner.

In `R/module_meta_scatter.R`:

- `pre_vol` is the pure store-watcher predicate (transitions exactly once
  per view change; never dips between two volcano views — constraint 1).
- `.scatter_axes_converged` (store watchers == displayed v1/v2 triples)
  is passed as `corner_apply_gate`, so the corner applies at display
  convergence (constraint 2).
- `rectval` consumes `params$cutoff_reactive` (`cutoff_effective`), so
  the resolved corner lands in the SAME reactive recompute as the new
  axes — one paint, rects included (constraint 2).
- Side effects are content-deduped (constraint 3).

## Verification (all green)

- Tier A (tests/e2e_agent/tier_a.mjs) x2 full runs: 97/97 each, including
  the tightened/added assertions:
  - load: volcano badge active, scorner=volcano, 2 rects;
  - volcano -> cor: exactly ONE render, corner already cleared (`|0`);
  - cor -> volcano: exactly ONE render with both corner rects (`|2`);
  - NEW volcano -> volcano (RE_vs_ME -> RE_vs_LE): exactly ONE render,
    rects stay 2, scorner carries over.
- 10 ms-resolution frame sampling in a live browser (stricter than the
  40 ms Tier A poller): all four switch directions (volcano->cor,
  cor->volcano, volcano->volcano both ways) = exactly one distinct frame
  transition, no intermediate garbage frames.
- Load no longer hangs: plotly container renders, R process 0.2 %
  instantaneous CPU after idle (the loop pinned ~95 %).
- Unit board green: widgetStore 52, triselectorCascade 6,
  scatterSelection 6, quickViews 17, appState 37, agentAssistant 64,
  agentFigures 67, agentLogging 21, agentWidgets 105, aiAssistantTools
  33, agentLogSummary 32, shinyAuxi 10, tableWidgetState 5.

## Environment reminders (unchanged)

- ports 7775–7799: kill orphans before/after browser runs
  (`ss -tlnp | grep -E ':77[0-9][0-9]'`); an orphan serving an OLD build
  poisons verification runs.
- Chrome headless flags: `--disable-webgl --disable-webgl2` (mandatory).
- Reinstall after every source change:
  `Rscript -e "install.packages("<repo>", repos = NULL, type = "source")"`.
