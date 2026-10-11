# WP-B finding — the table `columns` write "anomaly" (F4) is not a bug

**Outcome: NO product code change.** All three legs of the F4 evidence
dissolve under direct verification; the write path is fully functional.
One WP-0-gate harness bug was found and fixed along the way (read-back
field), because it *produced* one leg of the original evidence.

## What F4 claimed (ANALYSIS.md, spotcheck2_results.json, sweep receipts)

1. `set_widgets` on `dataspace.tab_*.columns` reports `applied`, but the
   applied value is recorded as a **scalar string** (the sent 1-element
   array "mangled").
2. The store read-back immediately after is **NULL** — "the write appears
   to be immediately overwritten" (suspected push/ack echo race).
3. The sweep recorded **no-ui-effect** for the columns cases.

## What is actually true (evidence per claim)

### Claim 1 — "scalar mangle": DISPLAY-ONLY, no store corruption

- `agent_args_sanitize` (R/auxi_agentTools.R:84-100) canonicalizes a
  homogeneous JSON array to an R **atomic vector** — by design, and for a
  1-element array `["X"]` the result is `"X"`, which in R *is* the
  ordinary length-1 character vector (R has no separate scalar type).
  Nothing is lost; `c("X")` and `"X"` are the same object.
- The multi_select validator already accepts both wire shapes (a scalar
  string validates as a 1-element selection; arrays are unlist'ed). The
  plan's "make the validator accept both" was already the behavior.
- The scalar *display* comes from the test-hook receipt serialization:
  `jsonlite::toJSON(..., auto_unbox = TRUE)` (R/module_agentTestHooks.R:200)
  renders length-1 vectors as bare JSON strings. The store holds the
  correct value — verified browser-side (wp_b_results.json): a
  `['General|All|Origin']` write reads back as exactly that value and the
  table renders exactly one column, `Origin`.

### Claim 2 — "read-back NULL": HARNESS BUG (now fixed)

- `sweep_harness.mjs` read `snap.values` from `runHook`'s `{result, ms}`
  wrapper — `(snap.values || {})` is always `{}`, so **every** applied
  case in every run recorded `store_readback: null` (159/159 applied
  cases in the 2026-10-12 WP-A full-run scratch). The "immediately
  overwritten" evidence was universal harness breakage, not a table
  anomaly. Fixed in the harness (read `snap.result.values`); no product
  change involved.

### Claim 3 — "no-ui-effect": OBSERVATION LIMIT of the input-diff probe

- A columns write acts through `scn(okcols)` → a **table re-render** (an
  output-level effect; R/module_dataTable.R push observer). The sweep's
  ui-effect probe diffs `Shiny.shinyapp.$inputValues` — outputs are
  invisible to it, so a *working* columns write registers zero input
  changes by construction. Page/filter pushes go through the DT proxy
  and tick `input$table_state`, which is why those keys are observable
  at all — and even those register nothing while the owning tab is
  hidden (deterministic: 0 non-hook input changes for columns/page/
  column_filters in every isolation run, while the same page push
  verifiably lands on a visible table).
- spotcheck2's "no effect after columns write" was a **rejected** write:
  its receipt shows `applied: []` — it sent `'featureId'`, which is not
  a `tab_feature` column (the rejection reason channel wasn't captured
  by that probe).

## The suspected push/ack echo race: NOT CONFIRMED (disproven)

The plan suspected the load-time sync observer echoes a stale `scn()`
over the fresh store value before the pushed re-render acknowledges.
Browser verification (verify_wp_b_table_columns.mjs, 16/16, two runs)
shows the opposite: after every write the store read-back equals the
written set once the echo settles — including after a re-render, after a
page push that follows a columns-driven re-render, and after opening a
previously-unopened tab whose table then renders the stored set. The
`store_sync_from_ui(store, "columns", scn())` observer acks the armed
pending entry (identical value) exactly as designed; no override log
entries, no clobber.

## Reproduction-rate evidence (instruction step 1)

Isolation re-runs of the `dataspace.tab_*` cases (`--only dataspace.tab_`,
full inventory + 12 cases each) — raw status and non-hook input changes:

| run | columns ×3 raw | non-hook Δ | page ×3 raw | non-hook Δ | gate |
|---|---|---|---|---|---|
| WP-A full run (q8K9oj) | success | 0 | success | 0 (expr: 2) | PASS |
| isolation 1 (mHTtoh)   | success | 0 | success | 0 | PASS |
| isolation 2 (IYcw3B)   | success | 0 | success | 0 | PASS |
| isolation 3 (4OE3LI)   | success | 0 | success | 0 | PASS |

(columns: 12/12 observations, zero non-hook input changes, zero failures;
probe runs ×2: 32/32 end-to-end checks green. Isolation 3 ran after the
harness read-back fix and reads back real store values, confirming it.)

The cases do NOT flap and do NOT regress: they deterministically produce
zero *input-level* changes while the write itself demonstrably works
(16/16 browser proof). This is neither the "congestion artifact" nor the
"reliable failure" branch — it is a measurement limitation of the harness
status rule (`applied && any-input-diff`), which the original curation
pass already corrected case-by-case (hook-input churn is not a UI effect;
the three columns cases were curated to no-ui-effect on exactly this
basis). The curated baseline therefore remains correct as-is: the three
`dataspace.tab_*.columns` no-ui-effect statuses are the documented
residue of an output-level effect that an input-diff probe cannot see.

## Consequences for the fix plan

- WP-B code change: **not warranted** (nothing to fix; both plan
  hypotheses — scalar mangling in the store, echo-race clobber — are
  disproven by direct evidence).
- WP-B2 (expr-heatmap congestion, F6): same caution applies — statuses
  derived from input-diff polling under main-thread congestion are
  unreliable evidence; re-verify live before any code change.
- Sweep methodology note: a binding write that acts through a re-render
  needs a DOM/output-level probe (as in verify_wp_b_table_columns.mjs)
  to be graded; the input-diff rule alone will always grade it
  no-ui-effect. Future WPs should treat no-ui-effect on re-render-type
  bindings as "verify with the DOM probe first".

## Artifacts

- `tests/e2e_agent/verify_wp_b_table_columns.mjs` — the deterministic
  reproduction/verification (16 checks), runnable any time:
  `node verify_wp_b_table_columns.mjs`
- `tests/e2e_agent/widget_sweep/wp_b_results.json` — its evidence output
- harness read-back fix: `tests/e2e_agent/widget_sweep/sweep_harness.mjs`
  (`snap.result.values`), benefits all future sweeps
