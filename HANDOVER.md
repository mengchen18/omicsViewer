# HANDOVER — agent accuracy & widget control plane

Written 2026-09-23 at the end of the S4-3 result-space step-1 session
(feature_general + sample_general + analyst navbar + shared attr4 panel).
Read this first in a fresh context, then `AGENT_ACCURACY_PLAN.md`
(canonical plan + status table, new row `1⅞+⅞`) and `AGENTS.md`
(environment + commands). Delete or trim this file once absorbed.

## Where we are / how to resume

Branch `agent-driven-exploration`, working tree clean through this
session's commit. DESCRIPTION 2.1.1. Data-space is COMPLETE (S4-1/S4-2);
result-space step 1 is COMPLETE:

| Change | What |
|---|---|
| this session | result-space step 1: `resultspace.analyst_tab` navbar + `resultspace.feature_general.*` (23 keys) + `resultspace.sample_general.*` (21 keys) + store-aware shared attr4 panel (`<prefix>.attr4.*`, 18 keys per instance) |

**Next session continues S4 result-space** in this suggested order:
gsList → ora → fgsea → barplotGsea → string → survival → PTMotif →
doseResponse → batch → contTableStats → roc_pr → geneshot → attr4
(figure tab). meta_scatter's embedded attr4 can now simply pass its
child store to `attr4selector_module` (the panel became store-aware this
session; meta_scatter still runs it store-less — migrating it is a
two-line change plus test updates). Embedded dataTableDownload instances
need nothing (buttons only).

## What landed this session (architecture notes)

- `R/module_figureAttr4.R`: `store` param; registers 5 cascades
  (color/shape/size/tooltip/search) + xcut/ycut strings + scorner
  (select; choices derive from the patch's own effective cutoffs —
  volcano needs both). Status restore routes through the store when
  store-backed; legacy s1/s2/s3 path kept for store-less callers.
  The hidden `searchon` select stays UNregistered (dormant wiring:
  searchValue is never fed from input$searchon).
- `R/module_feature_general.R`: `store` param; xax cascade (Surv
  excluded), plot_type enum, regression_line mirrored through the
  showRegLine reactiveVal (pushes re-render regTickBox). Status xax
  restore goes through the store; plot_type/attr4 legacy paths kept.
- `R/module_sample_general.R`: twin; cascade keeps Surv (survival
  view); batch-comparison link routes through
  `store_apply(origin="system")` via `.sg_apply_triple`.
- `R/L1_module_result_space.R`: `store` param; child stores
  `resultspace.{feature_general,sample_general}` + navbar binding
  `resultspace.analyst_tab` (kind navbar; choices provider mirrors the
  renderUI tab logic incl. GS/ResponseCurve/StringDB/SeqLogo and
  additionalTabs).
- `R/L0_module_app.R`: passes `store = app_store` to the result space.

## Rules discovered this session (do not re-learn)

1. **`c()` splices list ARGUMENTS one level**: `c(binds, widget_binding(...))`
   flattens the trailing record into its fields and store_register dies
   with `$ operator is invalid for atomic vectors`. Wrap appended
   records: `c(binds, list(widget_binding(...)))`.
2. **Child-of-child store views are unsupported** (`.widget_store_key`
   and `store_register` assume one prefix level). attr4 therefore builds
   its `<prefix>.attr4.*` namespace directly off the ROOT via
   `widget_store_child(root, paste0(prefix, ".attr4"))`.
3. **Pre-existing state-bridge gap (NOT a regression, confirmed on
   stashed code)**: `apply_agent_state` feature patches are transient —
   the `esv_status(NULL→full_state)` restore roundtrip re-derives `ri`
   from the data-space module return and empties it. Analysis-panel
   content therefore needs selections made through the real UI. Candidate
   future fix (out of scope): patch `panels$data_space` table selections
   inside apply_agent_state, or defer the v1()-watch while restoring.
4. **DT table interaction quirks** (Tier A): os-select style — plain
   clicks accumulate, ctrl-click REPLACES; the DT redraw on tab return
   can swallow one click (self-correcting loop in tier_a 4e); tab
   re-renders reset the multi-selection switch to its default (DOM
   bootstrap-switch class is the reliable state source, not
   `$inputValues` right after a tab switch).
5. Demo data: `General|All|MDR` is NUMERIC (beeswarm, not contingency —
   use TP53.Status/Origin for categorical views); raw demo.RDS lacks the
   GS attr but `readESVObj → tallGS` attaches it, so ORA/fGSEA tabs DO
   exist in the live app (don't use ORA as a "dataset-absent tab"
   negative — Response/SeqLogo work).

## Validation status (all green)

- Unit board: widgetStore 52, agentWidgets **71** (20 new result-space
  tests), aiAssistantTools 19, appState 30, tableWidgetState 5,
  agentAssistant 38, agentFigures 26, agentLogging 21,
  triselectorCascade 6, quickViews 17, scatterSelection 6, shinyAuxi 10
- Tier A browser: **74/74 × 2 consecutive runs** (new section 4e:
  result-space navbar, sample cascade Surv→survival / Origin→contingency,
  feature cascade + Curve→ROC/PR, attr4 color cascade, three rejection
  cases; features selected through the real Feature table per rule 3)
- NOT yet run: live Tier B (tasks 15–16 in tests/e2e_agent/tier_b_tasks.md
  cover the S4 surface incl. these keys) — run when a provider is
  configured; no new task added for result-space step 1 (covered by 15–16)

## Migration recipe (unchanged, proven five times now)

1. `store = NULL` param + roxygen `@param store`.
2. `store_register` every user-editable widget (kinds: string/numeric/
   integer/boolean/enum/select/select_cascaded/multi_select/tabset/
   navbar/checkbox/slider).
3. Choices providers must never `req()`; read module reactives through
   defensive no-req helpers (see `.a4_ts`, `.fg_ts`, `.sg_ts`).
4. UI→store sync per widget (`store_sync_from_ui`, ack-aware); module
   state (reactiveVals) watches the val itself.
5. Seeding: one observer, `origin="system", strict=FALSE`, skip held keys.
6. Store→UI push: ONE epoch observer; cascaded triselectors need no
   direct push (store_watch selectors drive the triselector module).
7. Observer retention mandatory (`.xx_keep(...)` lists).
8. Keep legacy status() paths (dual, idempotent) until the final re-route.
9. Tests: extend tests/test_agentWidgets.R + a tier_a section per module
   (visible DOM effect + rejection); run the board + Tier A ×2.

## Session workflow reminders

- Reinstall after every source change
  (`Rscript -e "install.packages('<repo>', repos=NULL, type='source')"`)
- Kill orphaned R processes on ports 7775–7786 before browser runs
- Tier A convention is ×2 consecutive runs; provider.env is gitignored
- Log timestamps are UTC (user is UTC+2)
