# HANDOVER — agent accuracy & widget control plane

Written 2026-09-23 at the end of the S4 steps-2/3 session (dataTableDownload
decision settled; ora/fgsea/ptm/geneshot/string/survival/batch migrated;
meta_scatter attr4 adopted). Read this first in a fresh context, then
`AGENT_ACCURACY_PLAN.md` (canonical plan + status table, new row `2`) and
`AGENTS.md` (environment + commands). Delete or trim this file once absorbed.

## Where we are / how to resume

Branch `agent-driven-exploration`, working tree clean through this
session's commit. DESCRIPTION 2.1.1. Every user-editable widget in the app
is now on the store except the documented exclusions below.

| Change | What |
|---|---|
| this session | dataTableDownload row-selection decision + shared `store`/`store_key` support; ORA + fGSEA + PTMotif + Geneshot + String + survival censor + batch comparison migrations; meta_scatter attr4 adoption; Tier A drift-race harness fix |

**Next session should close S4**: the only remaining S4 deliverable is the
**snapshot save/restore re-route through the store** (plan §6.4 row S4:
"snapshot round-trip equals store state exactly"). Every module currently
keeps dual, idempotent status paths (legacy tab_status/observeEvent +
store). Plan of attack: make `app_state`/`.ESS` snapshot assembly read
`store_snapshot()` as the single source, keep per-key-resilient restores
(strict=FALSE), then verify a save→reload round-trip equals the store
state exactly (new test; watch for the observer-GC and seeding-race traps
below). After that: WP1 (`sections` param on `get_omics_viewer_state`),
WP3 (figure spec round-trip), WP4 (log summarizer).

### The settled dataTableDownload decision (do not re-litigate)

Row selection registers ONLY where the click drives a downstream view with
no other agent interface: ORA results→overlap table, fGSEA results→
leading-edge barplot, STRING enrichment→network highlight, batch tables→
sample link (all done). gsList deliberately OFF: its click writes the
app-wide feature selection that `set_omics_viewer_state` already controls
(circular choices — the table lists the CURRENT selection's memberships);
snapshot persistence stays on its existing tab_status path. No-effect
tables (ORA overlapTab, doseResponse ×2, PTMotif ×4, geneshot autorif,
feature/sample general meta tables) stay OFF — cosmetic row highlight is
zero capability value. Action buttons (Run, Search) are commands, not
widget values — never registered. Full text in plan §6.2.

Documented no-widget no-ops: barplotGsea (helper module, never
instantiated in the app), doseResponse, factorIndependency
(contTableStats), plot_roc_pr — display-only, nothing to register.

## What landed this session (architecture notes)

- `R/module_dataTableDownload.R`: `store`/`store_key`/`store_label`/
  `store_help` params; single `select` binding with choices =
  `reactive_row_ids()` (call sites now pass semantic ids: ORA/fgsea
  pathway, string term, batch variable/feature names). **New push
  pattern**: a desired-row-id reactiveVal feeds `formatTab`'s
  `selection$selected`, so pushes re-render the table with the row
  preselected; DT's report acks via the normal sync observer. Guard: a
  push-time snapshot of `input$table_rows_selected` marks the stale
  pre-render report as not-a-user-override; any real input CHANGE wins
  (override log). Id↔row mapping goes through `tabsort()$index` in both
  directions (ORA display order ≠ ids order: tabsort re-sorts by p.value).
- `R/module_ora.R` / `R/module_fgsea.R`: `store` param; xax cascade via
  store_watch + store_epoch (meta_scatter pattern); `selected_row` on the
  results table; status restores routed through store_apply(origin=
  "restore") with the legacy path kept for store-less callers.
- `R/module_PTMotif.R`: cascade + seeding that writes the store ONCE per
  unset key (the old observer reset the selection on every fdata change —
  now a restore/user pick sticks).
- `R/module_geneshot.R`: `term` (string; empty never seeded/pushed) +
  ID-mapper cascade; Search button unregistered.
- `R/module_string.R`: taxonomy (string) + show_labels (boolean) +
  `selected_row` (term ids; strtab_df factored out of the eventReactive).
- `R/module_survival.R`: censor slider registered on the sample_general
  child view directly (key `survival_censor` — no child-of-child needed);
  push gated on the survival checkpoint and clamped to the rendered range.
- `R/module_batch_comparison.R`: two boolean toggles + two row-selection
  keys, all on the sample_general child (`batch_*`).
- `R/L1_module_result_space.R`: child stores created + passed for
  ora/fgsea/stringdb/ptm/geneshot.
- `R/module_meta_scatter.R`: attr4selector call now passes `store` (both
  data-space scatters register the 18-key panel under
  `dataspace.{feature,sample}_space.attr4.*`).

## Rules discovered this session (do not re-learn)

1. **Selectize choices land after the value**: an updateSelectInput
   `selected` message can be processed before the `choices` message, so a
   browser test that selectize.setValue()s right after the input map shows
   the value races an empty option list (reproduced on clean HEAD —
   machine timing, not S4 code). Fix pattern: waitForFunction that the
   selectize has >0 options before driving it (tier_a 3b).
2. **testServer runs the INSTALLED package** — stashing source changes
   does not change app behavior; always `install.packages(..., type=
   'source')` before drawing conclusions from app-level tests. Stale R
   processes on the port answer with the OLD app (cost a confusing
   debug: "welcome screen" was a previous run's app).
3. **Receipt `skipped` is serialized as `unchanged`** by
   agent_widget_apply — browser no-op assertions must check
   `receipt.unchanged`, not `receipt.skipped`. Diff-only receipts also
   mean "applied lists 3 keys" assertions are wrong when a key was
   already equal — re-apply the patch and assert all keys are unchanged.
4. vectORATall needs BOTH GS columns as factors and ≥ minOverlap
   enriched features — synthetic test gene sets must be non-degenerate
   (identical sets filter to zero rows via the p/OR cutoff).
5. DT push flow (new pattern, browser-verified end to end): agent patch →
   epoch observer sets desired id → re-render preselects → DT reports →
   sync observer acks. A user click during the push window overrides and
   logs — do not "fix" that away; it is the store's user-wins contract.

## Validation status (all green — as of this session's commit)

- Unit board: widgetStore 52, agentWidgets **105** (17 new: dtd 6, ora 4,
  fgsea 4, ptm 3, geneshot 2, string 2, survival 2, batch 2, attr4 2),
  aiAssistantTools 19, appState 30, tableWidgetState 5, agentAssistant 38,
  agentFigures 26, agentLogging 21, triselectorCascade 6, quickViews 17,
  scatterSelection 6, shinyAuxi 10, ora 23, motif 6, dose_response 60
- Tier A browser: **83/83 ×2 consecutive runs** (new section 4f: fGSEA
  ranking cascade, user row click → store no-op receipt, agent pathway
  push → table re-render, unknown-pathway rejection; plus the 3b
  selectize-options race fix). ORA specifics are unit-covered (live ORA
  row-driving needs ≥4 selected features and demo enrichment luck —
  fGSEA covers the identical dataTableDownload machinery in-browser).
- NOT yet run: live Tier B (tasks 15–16 in tests/e2e_agent/tier_b_tasks.md
  cover the S4 surface) — run when a provider is configured.

## Migration recipe (unchanged, proven; row selection variant in dataTableDownload)

1. `store = NULL` param + roxygen `@param store`.
2. `store_register` every user-editable widget (kinds: string/numeric/
   integer/boolean/enum/select/select_cascaded/multi_select/tabset/
   navbar/checkbox/slider).
3. Choices providers must never `req()`; read module reactives through
   defensive no-req helpers.
4. UI→store sync per widget (`store_sync_from_ui`, ack-aware); module
   state (reactiveVals) watches the val itself.
5. Seeding: one observer, `origin="system", strict=FALSE`, skip held keys.
6. Store→UI push: ONE epoch observer; cascaded triselectors need no
   direct push (store_watch selectors drive the triselector module);
   DT row selection pushes through the desired-id reactiveVal (re-render).
7. Observer retention mandatory (`.xx_keep(...)` lists).
8. Keep legacy status() paths (dual, idempotent) until the final re-route.
9. Tests: extend tests/test_agentWidgets.R + a tier_a section per module
   (visible DOM effect + rejection); run the board + Tier A ×2.

## Session workflow reminders

- Reinstall after every source change
  (`Rscript -e "install.packages('<repo>', repos=NULL, type='source')"`)
- Kill orphaned R processes on ports 7775–7786 (and 7790+) before browser
  runs — they answer with the stale app
- Tier A convention is ×2 consecutive runs; provider.env is gitignored
- Log timestamps are UTC (user is UTC+2)
