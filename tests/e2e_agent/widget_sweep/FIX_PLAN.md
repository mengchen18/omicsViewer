# Fix plan — systematic remediation of the widget-sweep findings

Companion to `ANALYSIS.md` (same directory). Every work package (WP) names:
root cause, change, files, blast-radius class, the guard suites that must be
green BEFORE merge (per AGENTS.md decision rule), and an acceptance criterion
expressed in sweep-matrix terms so `sweep_harness.mjs` is the regression gate.

Baseline matrix to beat (see `summary.json`):

| category | success | no-ui-effect | timeout | no-probe | unaddressable |
|---|---|---|---|---|---|
| binding (170) | 141 | 11 | 2 | 16 (8 benign) | — |
| dom-only (122) | — | — | — | — | 122 |

## Policy decisions (settled 2026-10-11)

1. **Destructive actions stay out permanently.** Snapshot save, restore,
   delete and their confirm dialogs are never exposed to the LLM. No
   action tier, no `run_control` surface, no capability records for them.
   They remain on the permanent exclusion list below.
2. **No dataset switching.** `app-selectFile` stays permanently outside
   the agent's reach (it is also the path-traversal-sensitive surface at
   `R/L0_module_app.R:367`). The agent operates ONE dataset per session,
   the one the user loaded.

Everything excluded by these decisions is counted as **by-design** in the
acceptance matrix, not as residual debt.

---

## Phase 0 — turn the sweep into the regression gate (prereq, ~0.5 d)

**WP-0 Sweep as acceptance suite. — DONE (2026-10-12)**
- Commit `widget_sweep/` (harness + baseline artifacts); add a runner
  `tests/e2e_agent/widget_sweep/run_regression.sh` that boots the app,
  runs `sweep_harness.mjs`, and diffs the status matrix against
  `baseline_summary.json` (committed), failing on any regression
  (`timeout`/`no-ui-effect`/`unaddressable` count increase or
  success→non-success transition).
- Gate behind `NOT_CRAN=true` exactly like `tests/test_agentUiEffects.R`.
- After each WP below: re-run, refresh baseline only for the WP's
  intended improvements, review the diff in the PR.
- Guard: none (it IS the guard). Acceptance: one command, deterministic
  pass on unchanged master.

Implemented as: `run_regression.sh` (orchestrator: NOT_CRAN gate, env
checks with clean skip, port-7794 hygiene, optional `--install`, auto
retry of regressed cases in isolation), `check_regression.mjs` (diff
engine + `--selftest`, failure rules above plus count guards for
`no-probe`/`rejected`/`unchanged`/`error`), `refresh_baseline.mjs`
(deliberate baseline refresh from a reviewed results dir), and
`baseline_summary.json` (337 per-case statuses from the curated run).
Harness gained two additive knobs: `WIDGET_SWEEP_OUT` (scratch output so
gate runs never clobber committed evidence) and `WIDGET_SWEEP_ONLY`
(caseId filter for retries/spot runs). Workflow per WP: run gate → review
transitions → `refresh_baseline.mjs --results <scratch>/results` only for
the WP's intended improvements → commit baseline with the code change.

---

## Phase 1 — correctness holes (bug fixes; land on master, cherry-pick to RELEASE_3_22)

### WP-A Close the dyn-heatmap validation/discovery hole (F3) — ~1 d — DONE (2026-10-12)

Root cause: `heatmapshinyApp` choices providers for the dynamic-heatmap
instance read module reactives that are `req()`-gated on widget
initialization, so while the tab is inactive (a) `store_registry_view`
records `allowed_values: NULL` (model cannot discover values) and (b)
`.widget_store_validate_value` treats NULL allowed as "no constraint"
(garbage like `GARBAGE|X|Y` is accepted, then diverges from the UI).

Change, two independent layers:
1. **Validator hardening (store core — cluster edit!)**:
   in `.widget_store_validate_value` (`R/auxi_widgetStore.R:321`), for
   kinds that derive allowed values from a `choices_provider`, a NULL/empty
   provider result must REJECT with a self-describing reason ("choices
   currently unavailable for <id>; open/activate the owning panel or query
   again") instead of accepting anything. Static `values` enums keep
   current behavior. Add opt-out only if a legitimate "free" multi-select
   exists (none known today — assert in review).
2. **Provider fix (module level — leaf edit)**:
   `heatmapshinyApp`'s annot/tooltip/sort providers read the dataset
   reactives (`pd()`, `fdColWithGS()`) rather than instance-gated state, so
   choices are activation-independent. Mirror of what the expr/cor
   instances already effectively do (their records had allowed values
   while inactive).

Blast radius: (1) touches the widget-store cluster → full unit board +
Tier A BEFORE the change; (2) is module-local.
Guards (run first): `test_widgetStore`, `test_agentWidgets`,
`test_snapshotRoundTrip`, `test_renderStability`, `test_tableWidgetState`,
`test_agentUiEffects` (browser).
Acceptance (sweep): the 5 dyn-heatmap `no-probe` cases become probeable
(`success` or `rejected` with the new reason while inactive + success after
provider fix); the garbage write from `supplementary_dyn_heatmap.json` now
rejects; zero regressions elsewhere.

Implemented as: (1) `.widget_store_allowed_state` in
`R/auxi_widgetStore.R` replaces `.widget_store_allowed_values` at the
three validation sites (select/select-ish, multi_select, mapping): a
declared-but-empty provider result now rejects with a self-describing
reason; static `values` enums and provider-less bindings unchanged;
empty-vector multi_select clears still accepted (min/max bounds kept);
review assertion held — every multi_select/mapping binding in the package
declares a provider. (2) `iheatmapModule` gained `choices_pd`/
`choices_fd` (default: the module's own pd/fd, so expr/cor are
byte-identical); the dyn instance passes the dataset-level `pdata`/
`fdata`; choices-update observers, store providers and the M9 push all
derive from them; `.heatmap_hcl_names` reads the display matrix
defensively (cold matrix ⇒ no HCL entries, not an abort). Third finding
while verifying: choice-set changes on tab activation re-defaulted the
selectize and the UI→store sync reverted a landed write (M8 clobber
class) — every choices push now preserves the current selection while it
stays valid (`.hmap_preserve`). Evidence:
`tests/e2e_agent/verify_wp_a_dyn_heatmap.mjs` (26/26: discovery healed,
garbage rejects, valid writes apply + reach the UI while inactive,
values survive activation). Baseline refreshed surgically (the 5 cases
no-probe→success; the 4 dyn dom-only inputs deliberately KEPT at
unaddressable as congestion canaries — see baseline note). Full sweep
gate PASS, zero regressions.

### WP-B Fix the table `columns` write anomaly (F4) — ~1 d — DONE (2026-10-12, NO product code change)

Root cause (hypothesis to confirm first, ~2 h): the sweep sent a 1-element
array; `agent_args_sanitize` canonicalizes it to a scalar (by design); the
receipt then records `applied_values` as a scalar, and the store read-back
is NULL — the write appears to be immediately overwritten. Prime suspect:
the columns push/ack echo in `module_dataTable.R` (~lines 362-395: the
push arms `scn(okcols)` but the load-time sync observer
(`store_sync_from_ui(store, "columns", scn())`, line ~285) may echo a
stale `scn()` back over the fresh store value before the pushed table
re-render acknowledges — the pending/ack window is not armed for columns
the way page/filters are via the proxy.

**Outcome (see `wp_b_finding.md` for the full evidence): all three F4
claims dissolved under direct verification; the write path is fully
functional; the plan's two hypotheses are disproven.** (1) The "scalar
mangle" is display-only: sanitize canonicalizes homogeneous arrays to
atomic vectors — for 1-element arrays the result IS the ordinary length-1
R character vector (no scalar type exists) — and the hook serializes
receipts with toJSON(auto_unbox=TRUE), rendering length-1 vectors as bare
JSON strings; the validator already accepts both wire shapes; zero store
corruption (browser-verified: read-back exact, table renders exactly the
written set). (2) The "read-back NULL" was a HARNESS bug — sweep_harness.mjs
read `snap.values` from runHook's `{result, ms}` wrapper, so every applied
case in every run read back null (159/159 in the WP-A full run); fixed in
the harness. (3) The "no-ui-effect" is the input-diff probe's observation
limit: columns writes act through a table RE-RENDER (output-level),
invisible to $inputValues diffing — and spotcheck2's "no effect" was a
rejected write (`'featureId'` is not a tab_feature column). The suspected
push/ack echo race is disproven: after every write the store equals the
written set once the echo settles (16/16 checks ×2 runs,
`verify_wp_b_table_columns.mjs`), including after re-renders, page pushes
that follow columns re-renders, and deferred renders on tab open.
Isolation re-runs ×3 + full run: columns cases deterministically show
zero non-hook input changes (not flapping, not failing) — the curated
no-ui-effect statuses remain correct as the documented residue. Sweep
methodology note recorded for future WPs: re-render-type bindings need a
DOM/output-level probe before believing no-ui-effect.

Change:
- Reproduce deterministically (extend `spotcheck2.mjs` into a unit test
  with `testServer` or keep browser-level): write columns while table
  initialized; assert store read-back equals the written set after the
  echo settles.
- Arm the columns push with the pending/ack protocol (mark pending like
  the other keys; clear only on the settled `scn()` echo), and make the
  multi-select validator accept both scalar and array for 1-element
  writes if the receipt mangling is in the diff display only (inspect
  `store_apply` diff construction for length-1 atomic vectors).
Blast radius: dataTable module glue + possibly `store_apply` diff
(avoid core if possible).
Guards (first): `test_tableWidgetState`, `test_widgetStore`,
`test_snapshotRoundTrip`, `test_agentUiEffects`.
Acceptance: sweep cases `dataspace.tab_*.columns` → `success` with
non-noise UI diff; read-back equals written value.

### WP-B2 expr-heatmap congestion follow-up (F6) — ~0.5 d — DONE (2026-10-12, NO code change: pass branch)

- Re-run cases 022/026 (`row_sort_by`, `annot_row`) live in isolation
  (spotcheck3). If they pass: mark harness congestion artifact in the
  baseline notes, no code change. If they reproduce: same provider/push
  analysis as WP-A/WP-B applied to the expr-heatmap instance.
Acceptance: both cases green or a filed bug with reproduction.

Outcome: both cases GREEN in every isolation run — 6/6 spotcheck3 runs
(built for this WP; the plan's named tool did not exist) plus an
in-harness isolation run (13/13 expr_heatmap success, 022 at 3.0 s and
026 at 0.3 s case wall time). Effects land within the hook roundtrip
(<~150 ms); the 2702×60 re-render stalls the main thread only ~1.2 s;
store read-backs exact; no echo clobber. Root cause of the original
no-ui-effect grades: the harness's waitUiEffect returns on the FIRST
input diff — always the hook's own controls (~20 ms) — so the real
observation window was ~20 ms, not 30 s (committed evidence:
waited_ms 21-22). Curated statuses kept as documented residue; full
evidence + methodology consequences in wp_b2_finding.md; baseline note
field updated. Unit suites recorded: renderStability 46/46,
snapshotRoundTrip 63/63.

---

## Phase 2 — the action tier: make every button operable (F1 core)

### WP-C `action` bindings + `run_control` tool — ~3-4 d — DONE (2026-10-12)

Root cause: the store is value-based by design (13 kinds); Shiny
`actionButton`s are event sources with no canonical representation, so 10
module action buttons + snapshot modal controls are structurally
unreachable.

Design (additive, mirrors the binding pattern):
1. **Store**: new registry section `store$actions` with
   `widget_action(id, label, help, handler, destructive = FALSE,
   requires_confirm = FALSE, availability = function() TRUE)` +
   `store_register_action()`, `store_action_run(store, id, origin)`,
   actions included in snapshot as last-run epoch only (no values).
   Ids share the binding namespace convention
   (`resultspace.geneshot.submit`). No changes to value validation or
   epochs — purely additive sections on the root store.
2. **Module glue (leaf edits per button)**: refactor each
   `observeEvent(input$<btn>, {...})` body into a named local function;
   the observer calls it, and the registered action handler calls it
   through `withReactiveDomain(session_domain, ...)` (same discipline as
   the existing apply callbacks). Phase-2 button inventory (demo, all
   non-destructive): `geneshot submit`, `stringdb run`,
   `tab_{pheno,feature,expr} clear`, heatmap `iheatmapClear`,
   `feature/sample space clear` — 8 buttons. Snapshot save/restore/delete
   are PERMANENTLY excluded (settled policy 1).
3. **Agent surface**: one new tool `run_control(action_id, params = NULL,
   _intent)` registered via the existing `agent_tool()` seam (append-only
   in `agent_tool_registry`); per-action `destructive_hint` annotations
   from the action record; rejection with closest-match suggestions for
   unknown ids (reuse `.agent_suggest_text`).
   **Policy guard**: `widget_action()` gains a `destructive` flag that is
   checked at registration — registering a destructive action is a
   package-author decision, and the assistant-facing registry simply never
   exposes destructive actions (settled policy 1). Snapshot lifecycle
   buttons are therefore not refactored onto the action tier at all.
4. **Discovery**: capability records gain action records
   (`id = "action:<canonical>", kind = "action", operation =
   "run_control"`) generated from `store$actions` in
   `agent_capability_records` — find_controls sees them automatically.

Blast radius: store core (additive sections) + registry + one new tool +
~8 leaf module edits. The three clusters are touched but only additively;
state which suites run BEFORE starting:
`test_widgetStore`, `test_agentWidgets`, `test_agentCapabilities`,
`test_agentToolSeam`, `test_agentAssistant`, `test_agentContracts`,
`test_aiAssistantTools`, browser `test_agentUiEffects`, plus new unit
tests (`test_agentActions.R`: register/run/reject/confirm-gating/
snapshot round-trip of action registry).
Acceptance (sweep): the 8 phase-2 buttons move from `unaddressable` to a
new `action` case category with `success` (button effect observed via the
existing per-button observable: geneshot search runs → table populates →
`selected_row` choices materialize, i.e. F5's stringdb chain unblocks);
`run_control` on unknown id → suggestion-bearing rejection (new semantic
case); `run_control` on a snapshot-lifecycle id → rejected as
"permanently excluded" (policy assertion case).

Implemented as planned (evidence: wp_c_finding.md +
verify_wp_c_actions.mjs 23/23 + test_agentActions.R 49/49). All 10
widget instances (8 button kinds × module instantiation) run through
run_control with handler-ran epoch proof; geneshot chain fully green
(519 real entries); unknown-id and policy rejections as specified.
Sweep: new action category (9 success, 1 availability-rejected in sweep
order — stringdb.run with the selection cleared by feature_space.clear
earlier in the loop; fully green with features selected), the 10 buttons
leave dom-only (122→112), gate PASS zero regressions, baseline refreshed
surgically. ONE pre-existing issue surfaced and FILED (not WP-C-caused,
A/B-proven against the stashed pre-WP-C module): the stringdb ENRICHMENT
table does not populate in this environment even from the real button
(network renders; enrichment empty), so stringdb selected_row choices do
not materialize — that investigation is exactly WP-E's `*_selected_row`
scope. Golden fixtures updated deliberately (exactly two allowlists
gained "run_control").

### WP-D Bind the four `searchon` gene-highlight pickers — ~1 d — DONE (2026-10-12)

The flagship demo gap (`verify_search_highlight.mjs` B4).
- Register `search_value` as a `multi_select` binding under each attr4
  namespace (`dataspace.feature_space.attr4.search_value`, sample space,
  resultspace feature/sample general) with
  `choices_provider = function(v) vv()` (values of the selected Search
  column; bounded by the registry's head(50) for discovery, full
  membership check at validation).
- Push glue in `module_figureAttr4.R`: pending-marked writes call
  `updateSelectInput("searchon", selected=…)`; `input$searchon` syncs via
  `store_sync_from_ui` (replacing the bare `searchValue` reactiveVal as
  the source of truth; keep `params$status$searchValue` seeded from the
  store for snapshot compat — the legacy path is documented at
  `R/module_figureAttr4.R:123-129`).
- Migration note: snapshot format gains the new keys; `store_restore`
  already merges unknown keys safely — old snapshots keep working.
Blast radius: attr4 module + store registration (leaf); observer-GC rule
applies (keep new observers in a module-level list like `.a4_keep`).
Guards: `test_widgetStore`, `test_snapshotRoundTrip`,
`test_scatterSelection`, `test_renderStability`, browser
`verify_search_highlight.mjs` extended with a phase C (agent sets
search_value through set_widgets → ring trace appears).
Acceptance: sweep adds 4 `binding` cases for `attr4.search_value`, all
`success`; `verify_search_highlight` 16/16.

Implemented as planned, plus two joint-patch discoveries: (1) the choices
provider derives from the EFFECTIVE view when the same transaction
installs the Search column (the triselector recommits asynchronously,
one observer flush later — varSelector is a pure function of the triple
and the data, so `depends_on = search_variable` + an effective-view
fallback makes a single 4-key patch valid); (2) a late-bind re-assert
re-applies a search_value push whose updateSelectInput fired before the
conditional searchon input existed (dropped by the client). The choices
observer gained selection precedence (legacy restore > store value >
preserved current selection — the WP-A re-default lesson). Legacy
status path (params$status$searchValue → pre_search) kept: old
snapshots migrate into the store key through the input round trip.
Evidence: verify_search_highlight 16/16 (per-check output incl. joint
patch phase C: 4/4 applied, read-back c("CNOT1"), ring trace at the
gene coordinates); sweep gate PASS zero regressions with the 4 new
binding cases success (binding 146→150 curated; raw 163), baseline
refreshed surgically to 341 cases; guards identical (widgetStore 69,
snapshotRoundTrip 63, scatterSelection 40, renderStability 46, plus
attr4SearchHighlight 10, attr4TooltipDefault 23, triselectorCascade 10).

### WP-E Dead-end chains & conditional widgets (F5) — ~1 d

- `selected_row` (stringdb): unblocked by WP-C (choices materialize after
  `run`). Add a sweep case: run → then set selected_row → success.
- Batch `*_selected_row`: investigate why rows were not materialized
  (activation-gated provider? same class as WP-A; apply the same fix).
- `survival_censor`: capability record gains a `precondition` note
  ("requires a survival-type sample variable selected") via the binding
  `help` text (cheap) + verify the push lands when the slider renders
  (the module already handles late-render pushes — add a testServer case:
  set xax_variable to a `surv` column, then write censor).
Acceptance: stringdb selected_row end-to-end success; censor write
succeeds once precondition met (new sweep/semantic case); remaining
no-probe count = 8 benign.

---

## Phase 3 — performance guardrail (F2)

### WP-F Categorical scatter-axis freeze — ~2 d

Root cause: a character axis with 2 702 unique values forces plotly to lay
out 2 702 category labels; software rendering blocks the main thread
~200 s. This also freezes the UI for human users picking such a column by
hand — an app bug the agent merely exposed.

Two layers:
1. **App render fix** (`module_scatter.R` render path, medium radius):
   when a character x/y exceeds ~500 unique values, render as a numeric
   index axis with hover tooltips carrying the real labels (and note the
   substitution in the axis title), instead of `type: 'category'`.
   Configurable cutoff via constant.
2. **Agent-side early guard** (leaf): in
   `agent_normalize_scatter_view` / the axis bindings' validators
   (`apply_agent_scatter_view` already re-checks numerics at
   `R/L0_module_app.R:886`), reject an axis whose cardinality exceeds the
   render cutoff with a message naming the count and suggesting numeric
   columns — the model gets an immediate, correctable error instead of a
   200 s freeze.
Guards: `test_scatterSelection`, `test_renderStability`,
`test_quickViews`, `test_scatterDefaults`, `test_agentUiEffects`,
`verify_search_highlight.mjs`.
Acceptance: sweep cases 028/030 become `rejected` (agent guard) in <5 s;
the manual UI path renders the bounded index axis; no timeout remains in
the matrix.

---

## Phase 4 — understanding (the first original complaint)

### WP-G1 Capability enrichment + semantic search — ~2 d

- Records gain: `tab` (owning navbar tab label — derive from the
  namespace map already implied by `.agent_widget_panel_label`,
  `R/auxi_agentCapabilities.R:186`; keep one static table in
  `constants.R`), `effect` (one sentence "what changes when you set
  this"), and `group` (attr4 group / cascade role).
- Surface the existing sr-only module prose (11 modules, e.g. STRING,
  Geneshot, Dynamic heatmap) as `panel` capability records so
  find_controls("highlight") / ("gene") return the Search cascade and the
  gene picker with real explanations (this directly addresses the
  discovery failure recorded in `verify_search_highlight` B3).
- Upgrade `agent_capability_search` from literal `grepl` to word-boundary
  + `agrep` fuzzy fallback (bounded: fuzzy only when zero literal hits;
  max_results unchanged 50).
Guards: `test_agentCapabilities`, `test_agentAssistant`,
`test_aiAssistantTools` (golden fixtures updated deliberately, not
blindly).
Acceptance: new sweep discovery cases: find_controls("highlight") returns
the search cascade + picker; find_controls("gene name") non-zero; all
records carry tab/effect fields (assert in a unit test).

### WP-G2 Downloads only (destructive actions and dataset switching are settled: permanently out)

- Downloads (`export_data` tool returning a file into the chat via
  shinychat attachment / download link message): optional, deferred, and
  to be built behind `OMICSVIEWER_LLM_ENABLE_EXPORT` (default off) if/when
  wanted. Not required for the program's acceptance matrix.
- The snapshot lifecycle and the dataset selector move to the sweep's
  **permanent exclusion list**: their `dom-only` cases get an explicit
  `excluded-by-policy` status in the harness (still enumerated, never
  tested) so the acceptance matrix counts them as by-design.

---

## Sequencing & effort

| order | WP | depends on | effort | risk class |
|---|---|---|---|---|
| 0 | WP-0 sweep gate | — | 0.5 d | none |
| 1 | WP-A validator+providers | WP-0 | 1 d | store cluster (guarded) |
| 1 | WP-B table columns echo | WP-0 | 1 d | dataTable glue |
| 2 | WP-C action tier + run_control | WP-0 | 3-4 d | store+registry+tool (additive) |
| 2 | WP-D searchon bindings | WP-0 | 1 d | attr4 leaf |
| 3 | WP-E chains/conditional | WP-C, WP-A | 1 d | leaf |
| 3 | WP-F axis freeze | WP-0 | 2 d | scatter render (medium) |
| 4 | WP-G1 enrichment/search | WP-0 | 2 d | capabilities leaf |
| 4 | WP-G2 downloads (optional) | decision made: not required | — | — |

Total ~11-13 working days. Bug fixes (WP-A/B) are master-first +
RELEASE_3_22 cherry-picks; the feature WPs land on `extend-llm-tools`
(and master after review).

## Target matrix (acceptance for the whole program)

| category | baseline | target |
|---|---|---|
| binding success | 141/170 | ≥ 165/170 (+action-adjacent chains) |
| binding timeout | 2 | 0 |
| binding no-ui-effect | 11 | ≤ 2 (documented) |
| binding no-probe | 16 (8 benign) | 8 benign |
| action cases (new) | — | 8 success (phase-2 buttons) |
| dom-only unaddressable | 122 | ≈ 15 by-design, 0 unaccounted |

By-design accounting of the ≈15 (each carries `excluded-by-policy` or
`covered-semantically` in its result file): assistant own UI (5,
intentional), snapshot modal controls (4, settled policy 1), dataset
selector (1, settled policy 2), downloads (8 controls → the remaining
download cases stay `unaddressable` until/unless WP-G2 is ever built —
counted here as the only optional residue), quick-view badges (7,
permanently covered semantically by `set_scatter_view`), add-column
triselector inputs + attr4 component selects (documented as
covered-semantically by their parent bindings), DT pagination/sort
controls (covered by page/columns bindings).

Every number above is produced mechanically by re-running
`sweep_harness.mjs` — no manual claims.
