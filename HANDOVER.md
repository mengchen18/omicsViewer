# HANDOVER — agent accuracy & widget control plane

Written 2026-09-23 at the end of the S4-completion session + decision
review (snapshot save/restore re-routed through the store; round-trip
equality verified; the WP1/WP7 open decisions settled with evidence).
**S4 is COMPLETE — the whole §6 control plane (S1–S4) is done and the
Phase-1 design decisions are settled.** Read this first in a fresh
context, then `AGENT_ACCURACY_PLAN.md` (status table rows through `2½`;
WP1 decision + disclosure ladder in §3-WP1; WP7/WP10 updated) and
`AGENTS.md` (environment + commands). Delete or trim this file once
absorbed.

**Update 2026-09-24: WP1 is COMPLETE** (see plan §3-WP1 STATUS) —
overview/sections ladder landed, overview 1,185 B vs 18,693 B full on
demo.RDS (15.8×), unit board green (agentAssistant 64, aiAssistantTools
22), Tier A 93/93 ×2. `agent_state` is now an on-demand builder function
(`state(sections)`), not a reactive — keep that contract when touching
L0_module_app.R / ai_assistant_module / agent_test_hooks_module.
**Next: WP3** (figure spec round-trip, hard prerequisite per WP7), then
WP4 (log summarizer).

## Where we are / how to resume

Branch `agent-driven-exploration`, working tree clean through this
session's commits. DESCRIPTION 2.1.1.

| Change | What |
|---|---|
| last session | S4 completion: store re-routes, `dataspace.active_tab`, round-trip tests (unit appState + tier_a 4g) |
| decision review | WP1 overview settled (scatter_view on the ground floor, read from the store); WP7 settled (update_figure confirmed core surface → WP3 hard prerequisite); WP10 closed by construction |

**Next session goes straight to implementation** (all groundwork +
decisions done):

1. **WP1 — `sections` parameter on `get_omics_viewer_state`** — **DONE
   2026-09-24** (decision was settled; implemented as designed):
   overview = dataset, active tabs, available tabs, selection counts +
   ≤20 example IDs, quick-view id+label lists, **`scatter_view`** (x/y
   axis triples + axis_mode of BOTH data-space scatters via
   `store_read`, keys
   `dataspace.{feature,sample}_space.{x,y}_{analysis,subset,variable}` +
   `.axis_mode`), `available_sections`, `state_policy`. Requested
   sections return today's full payloads. Touch points landed:
   `auxi_agentAssistant.R` (`agent_compact_state` +
   `agent_scatter_view_from_store` + `.agent_normalize_state_sections`),
   `module_aiAssistant.R` (tool schema + description), `L0_module_app.R`
   (`agent_state` is now an on-demand isolated builder),
   `module_agentTestHooks.R` (`overview` op forwards `payload.sections`),
   `constants.R` (`AGENT_STATE_SECTIONS`). Tests: test_agentAssistant.R
   (+26), test_aiAssistantTools.R (+3), tier_a.mjs (+4).
2. **WP3 — figure spec round-trip** (HARD PREREQUISITE — user-confirmed
   that `update_figure` is core surface): `render_assistant_figure()`
   includes the full normalized spec in the tool result; store normalized
   specs in the session `figures()` registry. Tests:
   `test_aiAssistantTools.R` (create → result contains spec; spec
   re-submitted as update normalizes identically).
3. **WP4 — log summarizer** `agent_summarize_log(path)` in
   `auxi_agentLogging.R`: event counts, per-tool call counts, error
   taxonomy (regex classes: unknown_id/column/tab, invalid_figure_spec,
   invalid_argument, no_dataset, request_limit, provider_failure),
   first-attempt success rate, retry outcomes, most-rejected arguments.
   New `tests/test_agentLogSummary.R` with a synthetic JSONL fixture
   (REAL fixtures now archived at
   `tests/e2e_agent/artifacts/decision-evidence-20260923/` — 9 JSONL
   logs; 3 of 5 real get_state calls wanted the current view).
4. Then WP5 (benchmark tasks — `tests/e2e_agent/tier_b_tasks.md`, tasks
   15–16 cover the S4 surface; optionally add a task 17 for the snapshot
   round-trip), WP6 figure templates, and WP7 patch-mode `update_figure`
   (gated on tasks 10–11 revision failures after WP1–WP3).

Also open (task, not decision): live Tier B baseline — the app prints
`Using model = "gpt-5.6-terra"` (env-configured), but Tier B needs
`tests/e2e_agent/provider.env` (gitignored) + `npm install` there once.
Check that first for before/after numbers around WP1.

## What landed this session (architecture notes)

- **Re-routes** (all conditional on a store being present; legacy paths
  kept for store-less callers — the standalone heatmap app and old
  snapshots):
  - `R/heatmapshinyApp.R`: the 13-param status observer now first applies
    a translated patch via `store_apply(origin="restore", strict=FALSE)`
    (`.heatmap_status_patch` via the `.heatmap_store_keys` map; multi_
    select keys may be empty), then the direct `update*Input` calls run
    as before (idempotent; also covers values the store rejects).
  - `R/module_dataTable.R`: `showColumns`/`multiSelection` status →
    `columns`/`multi_selection` store patch (the multi_select validator
    intersects with live colnames — safer than the old raw `scn(i)`).
  - `R/module_feature_general.R`: `plotType`/`showRegLine` status →
    store patches.
  - `R/L1_module_result_space.R`: `analyst_active_tab` status →
    `store_apply` on `store_rs`.
  - `R/L1_module_data_space.R`: NEW binding `dataspace.active_tab`
    (kind navbar, static 9-tab values: Feature, Feature table, Sample,
    Sample table, Cor, Heatmap, Dynamic heatmap, Expression, GSList) +
    standard sync/seed/push glue + status re-route. This was the last
    user-editable widget not on the store; the semantic
    `data_space_tab` path in `set_omics_viewer_state` remains the
    Tier-1 interface.
- `R/L0_module_app.R`: `app_module(store = NULL)` — callers may inject a
  store (embedding contexts, tests read it directly); default unchanged.
- `R/module_agentTestHooks.R`: new `store` op returning
  `isolate(store_snapshot(store))` (JSON via the existing renderPrint).
- **Round-trip tests**:
  - `tests/test_appState.R`: full `app_module` under testServer —
    agent-path widget writes, real snapshot save (.ESS to disk with
    embedded `widget_store`), drift, restore via the savedSS cell
    selection, `store_snapshot` before/after IDENTICAL across all keys.
  - `tests/e2e_agent/tier_a.mjs` section 4g: same through the real
    snapshot modal in the browser; drifted heatmap palette visibly
    reverts; deep-equality of all store values (canon-normalized for
    auto_unbox'd single-element arrays); cleans its .ESS out of
    inst/extdata.

## Rules discovered this session (do not re-learn)

1. **MockShinySession output poisoning**: any observer whose event
   expression ERRORS during a flush makes every subsequent
   `output[[...]]` read fail with "unexpected error resolving its
   promise". The app's `v1()` chain errors while `input$eset` is unset
   under mock (NULL → `if (sta$eset_active_tab != ...)` length-zero).
   App-level testServers must `setInputs("app-dataspace-eset", ...)`
   BEFORE the first flush. (Output reads then work fine.)
2. **Mock input maps don't follow update messages**: `updateNavbarPage`
   relays but the mock `$inputValues`-equivalent stays stale, so a saved
   panel-status payload can disagree with the widget_store in testServer
   (in real browsers they always agree — both derive from live bindings).
   When a round-trip test needs them consistent, drive the widget through
   its INPUT like a browser would.
3. **Status-path-after-store_restore ordering**: on .ESS restore,
   `store_restore` runs synchronously and the module status observers
   apply at the next flush — the status path can overwrite widget_store
   values when the two disagree (they normally don't; see rule 2). If a
   future cross-session restore mismatch appears, check this ordering.
4. **Content-dependent keys** (dataTableDownload `selected_row`, table
   rows in tab_status) restore exactly in-session (same-session save →
   restore revalidates against the unchanged live table); cross-session
   reloads re-validate against live choices and drop invalid keys
   per-key-resilient (documented, accepted).
5. `names(output)` on a mock session lists `impl|ns`, NOT registered
   outputs — probe output existence by reading, and expect rule 1 to
   mask the real cause.
6. Tier A JS: top-level `const s1` is already taken (session 1); name
   round-trip variables distinctly (`rtSaved`/`rtAfter`...).

## Validation status (all green — as of this session's commit)

- Unit board: widgetStore 52, agentWidgets 105, aiAssistantTools 19,
  appState **37** (7 new round-trip assertions), tableWidgetState 5,
  agentAssistant 38, agentFigures 26, agentLogging 21,
  triselectorCascade 6, quickViews 17, scatterSelection 6, shinyAuxi 10,
  ora 23, motif 6, dose_response 60, stats 7
- Tier A browser: **89/89 ×2 consecutive runs** (new section 4g:
  browser-level snapshot round-trip through the real modal + full store
  deep-equality). Note: `inst/extdata/ESVSnapshot_*_test{0,6}.ESS` are
  pre-existing gitignored leftovers from older runs — harmless.
- NOT run: live Tier B (needs provider.env).

## The store surface at a glance (for WP1 context)

~200 keys across `dataspace.*` (active_tab, feature/sample_space axes +
mode + attr4 panels, three heatmaps × 13 keys, three tables × 2 keys),
`resultspace.*` (analyst_tab, feature/sample_general + attr4 + survival
censor + 4 batch keys, ora/fgsea/ptm/geneshot/stringdb cascades and
selected rows). Every user-editable widget in the app is registered;
gsList row clicks, no-effect table selections, and action buttons are
deliberately OUT (plan §6.2 decision). The snapshot modal save/restore is
the canonical persistence path (widget_store embedded in .ESS).

## Session workflow reminders

- Reinstall after every source change
  (`Rscript -e "install.packages('<repo>', repos=NULL, type='source')"`)
  — testServer/browser harnesses run the INSTALLED package.
- Kill orphaned R processes on ports 7775–7795 before browser runs.
- Tier A convention is ×2 consecutive runs; provider.env is gitignored.
- Log timestamps are UTC (user is UTC+2).
