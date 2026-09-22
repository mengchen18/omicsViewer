# HANDOVER — agent accuracy & widget control plane

Written 2026-09-22 at the end of the S4-1/S4-2 session (heatmap
completion + data-space tables). Read this first in a fresh context,
then `AGENT_ACCURACY_PLAN.md` (canonical plan + status table) and
`AGENTS.md` (environment + commands). Delete or trim this file once
absorbed.

## Where we are / how to resume

Branch `agent-driven-exploration`, working tree clean through commit
`b076fb6`. DESCRIPTION 2.1.1. Data-space migration is COMPLETE:

| Change | What |
|---|---|
| S4-2 `b076fb6` | dataTable store migration (multi_selection + columns, min=1 bound) |
| S4-1 `af2bdbe` | multi_select kind + heatmap sorting/clustering/annotation registration |

**Next session starts with: S4 step 3 — result-space modules.** In
`R/L1_module_result_space.R`, create child stores per module (naming:
`resultspace.<module_id>`, mirroring the `dataspace.*` convention) and
migrate module by module in this suggested order (simplest surface
first): feature_general → sample_general → gsList → ora → fgsea →
barplotGsea → string → survival → PTMotif → doseResponse → batch →
contTableStats → roc_pr → attr4 (figure tab). The `dataTableDownload`
instances embedded in those modules need nothing (buttons only).

Migration recipe (proven three times now — meta_scatter, heatmap,
tables; copy the heatmap block in `R/heatmapshinyApp.R` as the
template):

1. `store = NULL` param + roxygen `@param store`.
2. `store_register` every user-editable widget. Kinds available:
   string/numeric/integer/boolean/enum/select/select_cascaded/
   multi_select/tabset/navbar/checkbox/slider.
3. **Choices providers must never `req()`** — a shiny.validation
   condition inside a provider aborts the whole `store_apply`
   transaction (`tryCatch(error=)` cannot catch it). Read module
   reactives through defensive no-req helpers that mirror exactly what
   the UI choices-update observers offer (see `.heatmap_hcl_names`,
   `.dt_col_choices`).
4. UI→store sync: `observeEvent(input$<id>, store_sync_from_ui(...),
   ignoreInit = TRUE)` per widget. For module-state (non-input)
   values like the table column set, watch the reactiveVal itself
   (ack-aware, so store pushes ack through the same observer).
5. Seeding: one observer, fires once when inputs exist, applies
   defaults with `origin = "system", strict = FALSE`, skips keys the
   store already holds (restore-first-wins).
6. Store→UI push: ONE observer watching `store_epoch(store)` (created
   ONCE), reads pending entries, dispatches the right updater
   (`updateSelectizeInput` for server-side selectize ids, else
   `updateSelectInput`/`updateSliderInput`/`updateSwitchInput`/direct
   reactiveVal writes).
7. **Observer retention is mandatory**: keep every store-glue observer
   in a module-level list (`.xx_keep(...)` pattern) — unreferenced
   observers are garbage collected between flushes and later pushes
   silently never fire.
8. Keep the legacy `status()`/`tab_status` restore path as-is (dual
   path, idempotent — legacy snapshots need it). The full re-route is
   the LAST S4 step, after every module is migrated.
9. Tests: extend `tests/test_agentWidgets.R` (registration count,
   seed/sync/push via testServer — see the dataTable section for the
   `.spy_input_messages` + status-contract observation pattern, and
   note MockShinySession only RELAYS update messages, so assert on
   spied messages, not input values); extend `tests/e2e_agent/
   tier_a.mjs` with a section per module (visible DOM effect +
   rejection case); run the board + Tier A ×2.

testServer quirks that WILL bite (documented in AGENTS.md, verified
again this session): `ignoreInit = TRUE` swallows the first value
change (warm observers with the seed value); update ids relay
unprefixed (match `(^|\.)id$`); testServer-local assignments don't
leak (use a sink env); choices providers reading module reactives
must be invoked inside the live session (destroyed-session reactives
error).

## Validation status (all green — do not re-litigate without evidence)

- Unit board: widgetStore 52, agentWidgets 51, aiAssistantTools 19,
  appState 30, tableWidgetState 5, agentAssistant 38, agentFigures 26,
  agentLogging 21, triselectorCascade 6, quickViews 17,
  scatterSelection 6, shinyAuxi 10
- Tier A browser: **55/55 × 2 consecutive runs** (sections 4c heatmap
  completion, 4d tables)
- Manual browser smoke (user-verified): user-edit drift detection,
  snapshot save/restore round-trip, registry surface (59 widgets in
  `dataspace` after S4: 26 at S3 + 27 heatmap keys + 6 table keys)
- NOT yet run: live Tier B on the S4 surface — tasks 15–16 were added
  to `tests/e2e_agent/tier_b_tasks.md` for this; run them (plus the
  smoke subset) when a provider is configured, ideally before starting
  result-space work so the data-space baseline is model-verified too

## Carry-forward invariants (accumulated S1–S4)

- Governing principle: agent-controllability mirrors user-controllability
  exactly; only user-editable widgets are registered/agent-visible.
- `multi_select` kind: character vector; empty vector / empty JSON
  array / `""` clear; min/max bound ENTRY COUNTS (tables use min=1);
  literal `"[]"`/`"null"` sentinel strings stay omitted optionals
  (shared AGENT_SENTINEL_STRINGS, normalized at every boundary).
- Sentinel sweeps in test_aiAssistantTools/test_agentAssistant/
  test_agentFigures must stay green for any new tool surface.
- Diff-only writes; restores are per-key resilient (strict=FALSE);
  seeding is restore-first-wins; ack-aware sync (a widget confirming a
  pushed value acks it, a differing value is a logged user override).
- Root-store + full canonical ids is the generic-tier wire contract;
  patches are JSON-object strings (dynamic keys); ellmer tibble/array
  coercions normalized in `.agent_widget_normalize_patch`.
- Store-backed vs store-less triselector regimes both exist
  (`reactive_selector1() %||% input$analysis` in module_triselector);
  test_triselectorCascade guards both — don't break the store-less one.
- Demo data quirks: truncated sample default axes (`PCA|All|PC1(`)
  cleanly rejected per key; demo pdata/fdata columns for test
  fixtures: `General|All|Cell.line`, `General|All|MDR`,
  `General|All|Gene.name`, `mean|Origin|RE`.
- Open threads (harmless): pending entries for system-origin seeds
  never acked (multi_select empty seeds hit this regularly); provider
  settings re-save kills in-flight ellmer turn (cosmetic).

## Session workflow reminders

- Reinstall after every source change
  (`Rscript -e "install.packages('<repo>', repos=NULL, type='source')"`)
- Kill orphaned R processes on ports 7775–7786 before browser runs
- NOT_CRAN=true needed for the browser-effect test; Tier A convention
  is ×2 consecutive runs
- provider.env is gitignored and holds the live key — never commit,
  never echo; artifacts under tests/e2e_agent/artifacts/ are gitignored
- Log timestamps are UTC (user is UTC+2)
