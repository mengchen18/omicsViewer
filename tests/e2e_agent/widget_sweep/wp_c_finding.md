# WP-C finding — action tier + run_control (F1 core)

**Outcome: implemented as planned; one pre-existing issue surfaced by the
acceptance run and filed (not WP-C-caused, A/B-proven).**

## What landed (commit: see git log)

1. **Store (additive section [7], `R/auxi_widgetStore.R`)**:
   `widget_action()` / `store_register_action()` / `store_action_view()` /
   `store_action_run()`; `store$actions` + `store$action_epochs` +
   `store$action_runs`; snapshot carries last-run epochs only
   (`store_snapshot()$actions`); restore never replays. Ids share the
   binding namespace; destructive actions are registrable
   (package-author decision) but never agent-runnable and never in
   agent-facing views.
2. **Agent surface (`R/auxi_agentActions.R` new + registry)**:
   `agent_action_apply()` — the run_control server pathway (mirrors
   `agent_widget_apply`): settled-policy exclusions enforced BY NAME
   (`AGENT_EXCLUDED_CONTROLS`: snapshot lifecycle, dataset selector —
   never registered anywhere, rejected as "permanently excluded" with
   `excluded_by_policy = TRUE`), closest-match suggestions for unknown
   ids, availability + confirm gating surfaced as self-describing
   receipts. One new tool `run_control(action_id, params, _intent)`
   through the `agent_tool()` seam.
3. **Discovery**: capability records gain `kind = "action"` records
   (`id = "action:<canonical>"`, `operation = "run_control"`,
   `input_id` dom mapping, `requires_confirm`, `available`);
   `store_action_view` skips destructive actions.
4. **Module glue (10 widget instances, 5 files)**: geneshot submit,
   stringdb run, dataTable clear ×3, iheatmapClear ×3 instances,
   scatter clear ×2. Each observeEvent body became a named local
   function; the button observer and the action handler call the same
   function (`withReactiveDomain(session, ...)`). stringdb's Run drives
   TWO consumers (network observer + enrichment eventReactive), so both
   the button and the action bump one `.run_signal` reactiveVal.
5. **Hard exclusions honored**: snapshot save/restore/delete and the
   dataset selector were NOT touched, NOT registered, NOT stubbed.
   Downloads/export: nothing built (cancelled by Chen).

## Acceptance evidence (`verify_wp_c_actions.mjs`, 23/23)

- **Criterion 1** (buttons operable): all 10 actions run via run_control
  with epoch bumps 0→1 (epoch recorded only after the handler body
  returns — store-level proof, immune to the input-diff blindness);
  `feature_space.clear` additionally shows the real side effect
  (selection 3 → 0 through the selection bus). 10 action capability
  records discoverable via find_controls.
- **Criterion 2** (geneshot chain): term set → `run_control
  action:resultspace.geneshot.submit` → **result table populated with
  519 real entries** (live AutoRIF network query).
- **Criterion 3**: unknown id (`dataspace.tab_pheno.columns`) →
  `"Unknown action id ... Closest matches: dataspace.tab_pheno.clear."`
- **Criterion 4**: `snapshot_save` / `app-snapshot_delete` /
  `app-selectFile` → `"permanently excluded from the assistant surface
  by policy"` with `excluded_by_policy = true`.

## Filed: pre-existing stringdb enrichment blockage (NOT WP-C)

The stringdb chain's last leg — enrichment table populates →
`resultspace.stringdb.selected_row` allowed values materialize — does
NOT complete in this environment. **A/B proof it is pre-existing**: with
`R/module_string.R` stashed at the pre-WP-C commit (real button clicked
in the UI, STRING tab open, default selection), the network renders but
the enrichment table stays empty for 112+ s (no error toast — consistent
with the `req(!overflow())` silent-return path when the selection mirror
reads as TRUE ⇒ 2702 genes > `STRING_MAX_GENES` = 300). With the WP-C
action (same conditions), identical outcome (network renders; the
isolated `reactive_ids()` read can also return empty at evaluation time
⇒ "No genes provided" toast). The enrichment REST endpoint itself works
(curl: string-db.org/api/tsv/enrichment returns data for demo genes).
The action tier is therefore **behaviorally equivalent to the button**
— the `*_selected_row` materialization investigation is exactly WP-E's
"batch *_selected_row" scope in the fix plan.

Also noted: agent-set selections carry `mirror = NULL` by selection-bus
design — data tables do not filter on them, so the dataTable-clear
side-effect demo via `set_omics_viewer_state` is not applicable; its
body-ran proof is the epoch bump (handler = `selectedRowOrCol(TRUE)`,
structurally identical to the button observer's call).

## Guard board (before → after, identical unless noted)

widgetStore 69→69, agentWidgets 104→104, agentCapabilities 46→46,
agentToolSeam 48→48, agentAssistant 99→99, agentContracts 7→7,
aiAssistantTools 48→48, tableWidgetState 10→10, snapshotRoundTrip
63→63, renderStability 46→46, agentUiEffects 113 + Tier A 112/112
(browser, before and after), **agentActions NEW 49/49**.

Golden fixtures updated deliberately (exactly two): the assistant tool
allowlists in `tests/test_agentCapabilities.R` and
`tests/test_aiAssistantTools.R` gained `"run_control"` — the deliberate
WP-C tool-surface change, nothing else moved.

## Sweep integration

`sweep_harness.mjs`: new `action` case category (every registered action
run through the `runcontrol` hook == the exact tool pathway; status
decided on the store-level epoch bump, per the WP-B/WP-B2 methodology),
and dom-only enumeration skips buttons mapped to registered actions
(`input_id`). `check_regression.mjs`: `action` count guard = error /
timeout / unchanged (availability-`rejected` is a legitimate,
self-describing outcome, not a failure).
