# Agent Tool-Call Accuracy & Efficiency Plan

Branch: `agent-driven-exploration`
Scope: the optional ellmer-backed assistant (`R/module_aiAssistant.R`,
`R/auxi_agentAssistant.R`, `R/auxi_agentFigures.R`, `R/auxi_agentLogging.R`,
wiring in `R/L0_module_app.R`).

Status: DRAFT — for review. No code has been changed.

---

## 1. Problem statement

The assistant currently exposes 7 tools:

| Tool | Kind | Notes |
|------|------|-------|
| `get_omics_viewer_state` | read | returns dataset, tabs, selection, **full annotation catalog**, **100 quick-view records**, panels, **figure grammar** |
| `search_annotations` | read | IDs / column names / bounded value matches |
| `summarize_annotation` | read | one-column summary |
| `set_omics_viewer_state` | write | tabs + semantic selections |
| `set_scatter_view` | write | quick view or custom axes |
| `create_figure` | write | declarative ggplot spec |
| `update_figure` | write | full spec resend, parented |

Identified failure/efficiency risks, in priority order:

1. **State payload bomb.** The system prompt mandates calling
   `get_omics_viewer_state` before answering anything, and that call returns
   the full annotation catalog (up to 200 columns × 2 spaces with per-column
   stats), up to 100 quick-view *full records*, the complete figure grammar,
   and panel state. For real datasets this is tens of KB of JSON on every
   conversational turn — token cost, latency, and distraction-driven accuracy
   loss. This is the single biggest lever.
2. **No self-correction affordances.** Validation errors are precise
   ("Unknown feature-space X-axis annotation: p value") but offer no
   suggestions, so a retry often repeats the same guess. The model must
   discover `search_annotations` on its own.
3. **Figure spec invisibility.** `create_figure`/`update_figure` return
   *metadata* only (geoms + mappings; no params, transforms, theme, palette,
   facet). `update_figure` requires the **complete** spec — which the model
   never received. It must reconstruct its previous spec from memory. This
   makes even simple revisions ("make the points red", "add a title")
   error-prone.
4. **No measurement loop.** Logging is comprehensive (JSONL: user messages,
   tool requests/results with errors, stream status/failures) but nothing
   consumes it. We have no first-attempt-success numbers, no failure
   taxonomy, and therefore no way to verify that any change helps.
5. **Complex figure composition.** The ggplot grammar is powerful but the
   model must compose nested specs (layers/mappings/params) for common plots
   (volcano, boxplot) that are conceptually one template.

Non-problems (already handled well, no action needed):

- Tool count and granularity (7 grouped semantic tools; not per-widget).
- Validation against *fresh* state — `apply_agent_state` /
  `apply_agent_scatter_view` validate proposals against current IDs, tabs,
  columns at apply time, so stale values are rejected.
- Security posture — allowlists everywhere, credentials excluded from state,
  logs redacted.

---

## 2. Design principles

1. **Progressive disclosure, not prompt bloat.** The model gets a compact
   overview plus the ability to request detail. No giant static UI dump in
   the system prompt; no tool-per-widget.
2. **Errors as a control surface.** Every validation failure should (a) name
   the exact problem, (b) offer closest valid matches, (c) point at the tool
   that finds more matches. Suggestions only — never silent fuzzy acceptance.
3. **Round-trip state.** Anything the model may be asked to revise (figure
   specs) must be returned to it in full, normalized form.
4. **Measure before/after.** A scripted benchmark task set + log summarizer
   gates each phase. Changes that don't move first-attempt success get
   revisited, not stacked.
5. **No new dependencies, no changes to non-agent code paths.** Everything
   lands inside the existing agent files (plus optional test-only helpers).

---

## 3. Phase 1 — predict-driven fixes (highest confidence)

### WP0: headless validation harness (two tiers)

**STATUS: Tier A complete (26/26 assertions green via `Rscript tests/test_agentUiEffects.R`).**

**First catch:** Tier A's initial run exposed a real product bug in the live
agent-restore path — `apply_agent_scatter_view` (and any future live
snapshot restore) left the feature scatter on stale axes. Root cause
(`R/module_triselector.R`): after the external axis request observer sends
`updateSelectInput("analysis", selected = "PCA")`, the init/self-heal
observer (`observeEvent(list(names(input), validated_x()))`) re-fired while
the browser ack was still in flight, preferred the stale `input$analysis`,
and reverted the select ~150 ms later. Fixed by tracking the last
server-sent-but-unacknowledged value (`pendingAnalysis`) so the self-heal
observer treats an in-flight restore as current. Verified: full cascade
`analysis → subset → variable` now applies and plotly titles equal the
requested axes; quickViews/appState/agentAssistant/agentFigures suites all
still pass.

**Environment notes** (verified the hard way):

- shinytest2/chromote is unusable on this machine: Chrome 147 + CDP hangs on
  any `http://` `Page.navigate` (`data:` URLs work). **Playwright + system
  Chrome is the driver** (`tests/e2e_agent/node_modules`, no browser
  download).
- Chrome must launch with `--disable-webgl --disable-webgl2`: headless
  SwiftShader reports `webgl_supported = true`, sending the scatter down the
  `toWebGL()` path, which stalls on software GL. Disabling WebGL makes the
  app's own detection pick the SVG path — matching the documented no-GPU
  desktop environment.
- The app selectizes every select; drive them via the `selectize` JS API.
- `#app-dataspace-eset` is a navbarPage binding, not a `<select>` — read tab
  state from `Shiny.shinyapp.$inputValues`.
- In plotly.js 2.x the graph container is the output div itself
  (`.js-plotly-plot` holds `_fullLayout`; there is no inner
  `.plotly-graph-div`).
- Spawned app processes must be killed on every exit path; orphans hold
  the port and poison subsequent runs with stale sessions.

**Verified environment**: headless Chrome works on this machine; `shinytest2` +
`chromote` installed; Playwright 1.63 via `npx`. No GPU/WebGL — software
rendering + the app's WebGL fallback cover DOM-level assertions.

**Tier A — UI-effect regression (no LLM, deterministic).**
`tests/test_agentUiEffects.R` (shinytest2 `AppDriver`):

- launch app on port 7775 with `inst/extdata` demo dataset;
- exercise the same state-transition semantics the agent tools invoke
  (`apply_agent_state` / `apply_scatter_view`: snapshot state → mutate →
  `esv_status()` round-trip) by driving the equivalent UI inputs;
- assert DOM outcomes: active tab class, plotly container presence, axis
  label text, selection table row counts;
- **cross-session isolation** (periodic, not per-commit): two concurrent
  AppDrivers — separate assistant sessions, separate diagnostic-log files,
  no config/state leakage (skill checklist item 1);
- runs in every test cycle; no API key, no network.

**Tier B — end-to-end agent validation (real provider, env-gated).**
`tests/e2e_agent/run_e2e.mjs` (Playwright, Node) + task list from WP5:

- spawns the Shiny app (R subprocess, port 7775, log `/tmp/omicsviewer_shiny.log`),
  with `OMICSVIEWER_LLM_LOG=true`, `OMICSVIEWER_LLM_LOG_DIR` pointed at a run
  directory, and a low `OMICSVIEWER_LLM_MAX_REQUESTS` (e.g. 8) for cost control;
- provider/key come from the environment (`agent_environment_config` picks
  them up automatically — no settings-modal interaction needed);
- opens the assistant drawer, types each benchmark prompt into the shinychat
  composer (web component → CSS selector + keyboard `Enter`; chromote/JS-inject
  fallback if flaky), fresh chat per task;
- waits generously for effects (LLM latency; 60–120 s timeouts), then asserts
  DOM outcomes the same way Tier A does (tab switched, scatter axes changed,
  `omicsviewer-ai-figure` img appeared, assistant status line healthy);
- saves per-step **screenshots as artifacts for human review** (the current
  model cannot view images; DOM/text assertions are the automated ground
  truth);
- on completion, hands the run's JSONL logs to WP4's summarizer for scoring
  (first-attempt success, retries, failure taxonomy);
- **stream-cancel recovery** task: cancel a prompt mid-stream, verify status
  returns to idle, the UI stays consistent, and the next submission works
  (skill checklist item 5);
- **negative-path run** (consumes no provider quota): deliberately invalid
  key / base URL → graceful error surfaced in chat or status line, app fully
  usable afterwards, diagnostic log records the failure (item 7);
- gating: runs only when `OMICSVIEWER_E2E=1` **and** a provider key is present;
  never in ordinary test cycles; benchmark pass = majority over 3 runs
  (LLM nondeterminism is expected and measured, not hidden).

The Tier A/B assertion lists are crosswalked against the review checklist in
`bioSkills/agent-development/shiny-ellmer-agents/usage-guide.md` §12 so the
harness inherits battle-tested failure modes instead of only the ones we
predicted.

Both tiers land **before** WP1–WP3 so every later WP gets UI-level regression
cover and the baseline (WP5a) is machine-scored rather than hand-waved.

### WP1: `sections` parameter on `get_omics_viewer_state`

**Change.** `agent_compact_state()` gains a `sections` argument.
`get_omics_viewer_state(sections = [...], _intent)`.

- **Default (`"overview"`)** returns: dataset, active tabs, available tabs,
  selection **counts + up to 20 example IDs** (currently 100), quick-view
  **id + label lists only** (currently full records), an
  `available_sections` field (`annotations`, `quick_views`, `panels`,
  `figure_grammar`), and `state_policy`.
- Requested sections return today's full payloads:
  - `annotations` — the column catalog;
  - `quick_views` — full records (id, label, x, y, description, source);
  - `panels` — bounded panel state;
  - `figure_grammar` — the allowlisted grammar description.

**Tool schema.** `sections`: optional array of enum
`["annotations","quick_views","panels","figure_grammar"]`, default
overview-only. Description updated: *"Call with no sections for a compact
overview; request sections only when needed."*

**Touch points.** `auxi_agentAssistant.R` (`agent_compact_state`),
`module_aiAssistant.R` (tool wrapper + description), `L0_module_app.R`
(pass `sections` through — note `agent_state` reactive currently builds the
whole thing eagerly; the tool wrapper should build on demand from parts).
Tests: extend `tests/test_agentAssistant.R`.

**Expected effect.** Turn-1 payload shrinks from ~tens of KB to ~1–2 KB for
large datasets; detail is one cheap local call away (no provider round-trip
for the tool itself — only the continuation request).

**Decision point for review:** does overview include panel state summary
(e.g. current scatter axes) or keep `panels` fully opt-in? I lean: include
*current scatter x/y axes* in overview (tiny, frequently asked), everything
else opt-in.

### WP2: did-you-mean suggestions in all validation errors

**STATUS: complete.** `.agent_suggest()` (substring → prefix → containment →
segment-aware edit distance, adist skipped for candidate sets > 5,000) wired
into tab/ID/quick-view/axis/column/facet/mapping rejections and zero-hit
searches. Also normalized the literal-`"null"`-string artifact (glm flash
serializes omitted optional strings as `"null"`; observed live when
`set_scatter_view` failed 3× with `Unknown quick view: null`). Tier B run
confirmed the model now self-corrects typos via `search_annotations` +
suggestions without user intervention.

**Change.** New internal helper in `auxi_agentAssistant.R`:

```r
.agent_suggest(value, candidates, max = 5)
```

Scoring: exact case-insensitive → substring containment (both directions) →
prefix → `adist` similarity, top-N. **Performance guard:** for candidate sets
> 5,000 (feature IDs on real datasets), restrict to substring/prefix
matching only (no `adist`); cap scan length.

Wire into:

| Site | Error today | Error after |
|------|-------------|-------------|
| `agent_normalize_state_update` (tabs) | "Unknown data-space tab: X" | + "Closest: Feature, Feature table" |
| `agent_normalize_state_update` (IDs) | "Unknown feature ID(s): X" | + up to 3 closest + "Use search_annotations(query=\"X\") for exact IDs" |
| `agent_normalize_scatter_view` (quick view) | "Unknown feature quick view: X" | + closest ids/labels |
| `agent_normalize_scatter_view` (axes) | "Unknown feature X-axis annotation: X" | + closest columns + naming-convention reminder |
| `agent_summarize_annotation` (column) | "Unknown feature annotation column: X" | + closest columns |
| `agent_normalize_figure_spec` (mappings/facet) | "...unavailable column(s): X" | + closest columns for the relevant data_source |
| `agent_search_annotations` (0 hits) | (empty result) | + `suggestions` field with closest column names (not an error) |

**Tests.** Extend `test_agentAssistant.R` / `test_agentFigures.R` with
typo cases (`"p value"` → `ttest|A_vs_B|log.fdr`; `"featrue"` tab;
`Gene1 ` whitespace; case variants).

### WP3: figure spec round-trip

**Change.** `render_assistant_figure()` in `module_aiAssistant.R` includes
the **full normalized spec** (`rendered$normalized`) in the tool result value
(alongside metadata). It is already JSON-safe (that is what
`agent_normalize_figure_spec` produces). `update_figure` keeps full-resend
semantics — now viable because the model holds its own last spec.

**Also.** Store the normalized spec in the session `figures()` registry
(needed later for patch mode, WP7).

**Tests.** Extend `test_aiAssistantTools.R`: create → result contains spec;
spec re-submitted as update normalizes identically.

### WP4: log summarizer (`agent_summarize_log`)

**Change.** New function in `auxi_agentLogging.R`:

```r
agent_summarize_log(path)   # or agent_summarize_logs(directory)
```

Parses one JSONL log; returns:

- counts by event type; session duration;
- tool call counts per tool;
- **error taxonomy** per failed `tool_result` (classified by regex on
  message + tool name): `unknown_id`, `unknown_column`, `unknown_tab`,
  `invalid_figure_spec`, `invalid_argument`, `no_dataset`,
  `request_limit`, `provider_failure`;
- **first-attempt success rate**: per `tool_call_id`, did the matching
  `tool_result` error?
- **retry outcome**: same tool name re-invoked within the following K events
  and succeeded;
- most-rejected arguments (argument key frequencies in errors).

Console-friendly `print` method. No UI — this is a developer tool run via
`Rscript -e "omicsViewer:::agent_summarize_log(...)"`.

**Tests.** New `tests/test_agentLogSummary.R` with a synthetic JSONL fixture
covering each taxonomy class + a successful retry sequence.

### WP5: benchmark task set + prompt workflows

**Change A — benchmark.** `tests/agent_benchmark_tasks.md` (or inst/):
~12 scripted tasks against `demo.RDS`:

1. "What dataset is loaded?" (state orientation)
2. "How many genes are selected?"
3. "Switch to the Sample tab"
4. "Make a volcano plot of the current results" (scatter quick view)
5. "Plot log FDR vs mean difference with custom axes"
6. "Find genes matching 'kinase' and select the first five"
7. "Summarize the sample group annotation"
8. "Boxplot of expression for the selected genes grouped by sample group"
9. "Create a histogram of the t-test log FDR"
10. "Change the last figure's title to X" (revision round-trip)
11. "Color the volcano points by category" (revision)
12. "What can you change in the app?" (capability disclosure)

Protocol: executed by WP0 Tier B (not manually): scripted Playwright run,
fresh chat per task, logging on, scored via WP4. Run **before** any code
change to establish the machine-scored baseline.

**Change B — system prompt.** Append 4 concise workflows (volcano; find +
select genes; expression boxplot; revise last figure) and the exact-ID
contract ("never guess IDs; use values from tool results; on error, use the
suggested matches or search_annotations"). Keep the whole prompt ≤ ~40 lines.

---

## 4. Phase 2 — composition accuracy & new capabilities

### WP6: figure templates

`create_figure(template = "...", ...)` where template is an enum
(`volcano`, `scatter`, `boxplot`, `barplot`, `histogram`, `density`,
`line`) and a few well-named arguments (`x`, `y`, `color`, `label_top_n`).
Server-side `agent_figure_template_spec()` expands to a full validated spec
(also returned via WP3 round-trip, so follow-up customization flows through
`update_figure`). Generic grammar stays as the advanced path. `pca`
template optional — needs server-side projection; decide later.

### WP7: patch-mode `update_figure`

`update_figure(figure_id, changes = list(labels = ..., theme = ...))` where
`changes` has the same shape as a spec but all fields optional. Server
merges onto the stored normalized spec (WP3 registry), validates the merged
result, renders. Full-resend remains supported (spec present ⇒ patch
ignored). This is mostly convenience once WP3 exists; do it if benchmark
task 10/11 show revision failures persisting.

### WP8: first new capability tools (enrichment / table view)

`set_enrichment_parameters` (ORA/fGSEA panel), `set_table_view`
(feature/sample table page + filter). Each ships with its help text stored
in a **shared metadata structure** (see WP9) from day one, not retrofitted.
This is the point where the capability registry starts earning its keep.

---

## 5. Phase 3 — architecture for scale (log-gated)

### WP9: UI capability registry + discovery tools

Only when the writable surface grows past ~4 capability tools / the tab
count grows. A server-side list of capability records
(`id`, `panel`, `label`, `description`, `writable`, `allowed_values`,
`operation`, `help_text`) generated from the same metadata used for UI
tooltips (one source of truth). Tools: `search_ui_capabilities(query)`,
`get_ui_capability(id)`. `get_omics_viewer_state` overview then lists only
capability counts, not contents.

### WP10: state token (optional)

`get_*` results carry `state_version`; `set_*` accept an optional echoed
token; mismatch ⇒ instructive error. **Trigger:** only if WP4 logs show
stale-state failures (e.g. model overwrites a user's mid-conversation
manual selection). Current validation already rejects stale IDs/tabs, so
expected benefit is narrow; friction is real.

### WP11: conversation-in-snapshot (history persistence)

Enable `chat_server(history = TRUE)` with the shinychat `on_save()` /
`on_restore()` auxiliary-state pattern (usage-guide §8): persist the chat
transcript **together with** the figure registry and current selections
inside the existing `.ESS` snapshot machinery (pillar 1 ↔ pillar 3 bridge).
Restoring a snapshot then revives both the widgets and the conversation
context. Guardrails: no credentials or API keys ever in history values;
restored tool evidence must not restore any execution authority (re-validate
on restore); opt-in per snapshot since transcripts may contain sensitive
dataset content.

### WP12: governance upgrade (cost ceiling / deputy evaluation)

Two cheap-to-expensive steps, log-gated:

1. **Session cost ceiling now**: we already log per-turn `tokens`/`cost`
   (`agent_log_turn`); extend the existing `on_request_start` request-limit
   hook with an optional cumulative cost/token budget
   (`OMICSVIEWER_LLM_MAX_COST_USD`), mirroring deputy's `UsageLimits`
   semantics without adopting the package.
2. **deputy migration evaluation later**: when `deputy::Agent` stabilizes,
   evaluate wrapping our ellmer client in it to gain permissions,
   allowlists, approvals, and budget enforcement out of the box
   (skill: raw ellmer is acceptable for benign custom tools — which we have —
   but deputy is the intended end-state for governed agents).

---

## 6. Universal widget control plane (NEW PILLAR — designed 2026-09-21)

**Goal:** the agent (and snapshots, and any future automation) can flexibly
and reliably read and control *every individual widget/component*, with no
unrequested side effects — via architecture, not per-bug patching.

**Evidence this is needed:** five defects found in a single day, all in one
module's state machinery, each a distinct symptom of the same gap — the app
has no single source of truth and no synchronization protocol between
widgets, module-internal models, snapshot state, and agent writes:

| Defect (date) | Symptom | Missing mechanism |
|---|---|---|
| stale-input revert | in-flight restore reverted ~150 ms later | acknowledgement protocol |
| mid-restore `req()` abort | restore event silently consumed | transactional applies |
| internal xax/yax drift | applies no-op after manual edits | UI→store sync (canonical model) |
| axisMode side effect | mode tab flips when only axes requested | diff-based (minimal) writes |
| `"null"`/`"{}"` args | spurious validation errors | per-kind boundary validation |

### 6.1 Architecture

Unidirectional data flow with a command/ack protocol (Elm/Redux lineage,
adapted to Shiny; extends the app-level canonical-state bridge from the
bioSkills skill down to widget level):

```text
                      ┌───────────────────────────────────────────┐
   user edits ───────▶│  canonical widget store (per-app env of    │
   (input bindings,   │  reactiveVals keyed by canonical IDs,      │
    origin="user")    │  each with origin + epoch metadata)        │
                      └──────┬──────────────┬──────────────┬──────┘
                 projections │              │              │
              (renderPlotly, │              │              │
               tables, …)    │              │              │
                            ▼              ▼              ▼
                     status collectors  snapshot      agent tools
                     (read store,       (store →     (validate →
                      never raw input)    .ESS)        state_apply)

   external writers ──▶ state_apply(patch) ──▶ diff ──▶ ordered setters
   (agent / snapshot   (transactional function, never an observer cascade)
    restore / badge)
                            │
                            └─▶ epoch bump + re-assert loop (bounded)
                                until widget ack == store value
```

The five mechanisms:

1. **Canonical store.** Every controllable widget is *registered* with a
   binding: canonical id, kind (`select`, `selectize_server`, `radio`,
   `tabset`, `navbar`, `checkbox`, `slider`, …), value getter, setter
   (correct `update*Input` + cascade position), choices provider, upstream
   dependencies (variable depends on analysis+subset), and owning module.
   The store is the single source of truth for desired state; widgets are a
   view.
2. **Transactional applies.** All external writes go through one function
   `state_apply(patch, origin=)` that (a) validates against the registry,
   (b) computes the diff *vs the store* (not vs widget inputs — kills the
   drift class), (c) writes store keys, (d) applies setters in dependency
   order, (e) bumps an epoch for touched keys. It is a plain function, not
   an observer — transient `req()` failures cannot consume it.
3. **Diff-only writes.** Untouched keys are never written: the axisMode
   side effect disappears by construction (custom-axes applies never touch
   the mode; quick-view applies do, because the badge state is part of the
   requested view).
4. **Ack + re-assert.** Setters mark in-flight values (generalized
   `pendingAnalysis`); a bounded verify loop re-sends until acknowledged.
   Status collectors and `get_omics_viewer_state` read the *store*, so
   state reads are never stale mid-flight.
5. **UI→store sync.** Each binding subscribes to its input and mirrors user
   edits into the store (origin="user", no re-assert). Internal models can
   no longer drift from reality.

### 6.2 Registration contract (ergonomics)

Hand-maintaining a registry would rot. Bindings are declared co-located
with the UI, via wrappers that return the tag unchanged:

```r
# module UI code
agent_input(
  selectInput(ns("variable"), NULL, choices = NULL, selectize = TRUE),
  id = "dataspace.feature_space.y_axis",
  kind = "select_cascaded",
  depends_on = c("dataspace.feature_space.y_analysis",
                 "dataspace.feature_space.y_subset"),
  help = "Y-axis variable; choices depend on the analysis and subset above"
)
```

The WP9 capability registry (labels, help text, allowed values, writability)
is **generated from these declarations** — one source of truth serves
humans (tooltips), the model (registry/discovery), and the apply protocol
(validation). New widgets become agent-controllable by registering, not by
writing a tool.

**Governing principle (settled):** agent-controllability mirrors
user-controllability *exactly*. Every widget a user can change is
agent-settable and agent-discoverable; anything a user cannot change is
registered (at most) as internal state for snapshot purposes and is never
exposed to the agent — discovery, descriptions, and tools stay free of
components the model cannot meaningfully use.

### 6.3 Agent tool surface (three tiers)

- **Tier 1 — semantic capability tools** (unchanged strategy):
  `set_scatter_view`, future `set_enrichment_parameters`, … implemented as
  thin validate + `state_apply` wrappers. Best accuracy for common intents.
- **Tier 2 — generic widget tools** (now safe because the substrate is
  reliable): `list_widgets(section)`, `get_widget(id)`,
  `set_widgets(patch)` — one batched setter, registry-driven validation,
  WP2-style suggestions. This is what "flexibly control every widget"
  literally means; it was previously rejected because the substrate could
  not guarantee it.
- **Tier 3 — discovery**: `search_ui_capabilities` from WP9, generated from
  the registry. Progressive disclosure stays (WP1).

Prompt contract: prefer Tier 1 for known intents; Tier 2 only for what
Tier 1 doesn't cover; the registry tells the model what exists.

### 6.4 Migration phases (each independently valuable)

| Phase | Deliverable | Acceptance criteria |
|---|---|---|
| S1 | Store + binding spec + `state_apply` + ack/epoch machinery, pure logic, unit-tested standalone (`R/auxi_widgetStore.R`) | unit suite: diff, ordering, ack, re-assert, origin tracking — no browser needed |
| S2 | Migrate meta_scatter (x/y axes, axisMode, attr4): store-backed; **delete** `pendingAnalysis`, the hand-rolled cascade ordering, and the restore observer's bespoke versioning; mode side effect gone | today's four repro scenarios + manual-edit stickiness, as generated Tier A cases; no observer-cascade races reproducible under stress (10 rapid interleaved edits/applies) |
| S3 | Tier-2 generic tools + registry generation from bindings; WP9 re-scoped to consume the registry | agent changes any registered widget (incl. one never exposed before, e.g. heatmap params) in one Tier B run each |
| S4 | Migrate remaining data-space modules (tables, heatmaps), then result-space; snapshot save/restore re-routed through the store | snapshot round-trip equals store state exactly; Tier A green across modules |

**Known issues assigned here** (explicitly *not* individually patched):
- axisMode flips on axis-only applies → fixed by S2 diff-only writes
  (interim: tool description tells the model to prefer explicit axes and
  not to switch modes unprompted — a prompt line, not a code patch)
- provider `"null"`/`"{}"` artifacts → S1 per-kind validators normalize
  known sentinels at the boundary (generalizes the current ad-hoc
  `.agent_nullable_scalar`)

### 6.5 Risks & mitigations

- **Dynamic UIs re-bind** (renderUI outputs): bindings keyed by canonical
  id, re-registration idempotent; epoch detects stale bindings.
- **Server-side selectize**: choices live server-side; setters go through
  the same `updateSelectizeInput(server=TRUE)` channel — protocol unchanged,
  but S2 must include a server-selectize widget as an acceptance case.
- **Feedback loops**: UI→store sync marks origin="user" and never re-asserts;
  store→UI writes only for diffed keys; ack loop is bounded (3 retries).
- **Performance**: store is O(registered widgets), re-assert only on diff;
  stress test in S2 acceptance.
- **Scope creep**: S1/S2 are self-contained; if S3+ stall, the curated
  tools still benefit (they move onto the store in S2).

### 6.6 Relationship to existing plan

- WP1–WP7 unchanged (accuracy of curated tools; they gain reliability by
  being re-implemented on the store in S2).
- WP8 (new capability tools) becomes cheap: register widgets + thin tool.
- WP9 (capability registry) is re-scoped: **generated from bindings**
  instead of hand-maintained.
- WP10 (state token) becomes trivial: the store epoch *is* the token.
- The "no tool per widget / no generic set_input" non-goal is amended:
  generic widget *tools* remain rejected for model-facing primary use, but
  a registry-validated generic tier becomes permitted once S1/S2 land.

---

## 7. Sequencing & effort

| Order | WP | Effort | Files |
|-------|----|--------|-------|
| 0a | WP0 Tier A UI-effect harness | **DONE** | tests/test_agentUiEffects.R, tests/e2e_agent/tier_a.mjs + R/module_agentTestHooks.R (env-gated hooks) |
| 0a′ | **Unplanned: live-restore revert fix** | **DONE** | R/module_triselector.R (found by Tier A) |
| 0a″ | **Unplanned: figure-spec limit raise 8→12 + logging sanitizer fix** | **DONE** | R/auxi_agentFigures.R, R/module_aiAssistant.R, R/auxi_agentLogging.R (unnamed-list crash silently killed `tool_request` logging for figure specs; drops are now never silent) |
| 0b′ | **Unplanned: ellmer tibble-coercion fix in figure specs** | **DONE** | R/auxi_agentFigures.R — ellmer converts `type_array(type_object)` args into tibbles, so `length(spec$layers)` counted 15 columns, not layers: **every chat-path `create_figure` failed the layer cap regardless of count** (the model was innocent; it even said "empty objects misparsed"). Layers/params now coerced to row-lists before counting; JSON-null→NA params fall back to defaults. Verified live end-to-end (Tier B driver, glm-5.3-flash): volcano prompt → get_state → create_figure (4 layers) → fig_1 rendered, 2702 rows, 416 KB PNG |
| 0b | WP0 Tier B e2e driver | **first run done** (`tests/e2e_agent/tier_b.mjs`, `node tier_b.mjs "prompt"`) | screenshots + archived log under tests/e2e_agent/artifacts/ |
| 1 | WP2 suggestions | **DONE** | auxi_agentAssistant.R, auxi_agentFigures.R, tests |
| 1½ | **§6 control plane S1: widget store + state_apply** | **DONE (34/34)** | R/auxi_widgetStore.R (single file, comment-sectioned: registration/validation/apply/ack/snapshot/registry) + tests/test_widgetStore.R |
| 1¾ | **§6 control plane S2: migrate meta_scatter onto store** | **DONE (Tier A 33/33 incl. stress + right-panel guard)** | meta_scatter x/y axes + axis_mode store-backed; xax/yax/axisRequest machinery deleted; apply_agent_scatter_view never writes axisMode (side effect fixed); cascades derive choices from the store, not in-flight inputs; restores apply per-key-resilient (strict=FALSE); placeholder states never enter the store. Found en route: demo.RDS stores truncated sample default axes (PCA\|All\|PC1( ) — old code silently no-opped them, now cleanly rejected per key
| 1⅞ | **§6 control plane S3: generic widget tier + heatmap registration + snapshot bridge** | **DONE (unit 31+16+37; Tier A 40/40 ×2; live Tier B: model chains list_widgets → set_widgets on heatmap colors, zero rejections)** | R/auxi_agentWidgets.R (agent_widget_list/describe/apply: registry-driven, per-key rejections with WP2 suggestions, JSON-object-string patch + tibble/record normalization); ellmer tools list_widgets/get_widget/set_widgets registered when a store is passed (L0 wiring); heatmap colors/scale/margins registered on all three heatmap instances (dataspace.{cor,expr,dyn}_heatmap.*) — the never-exposed-before acceptance case; prompt-contract line (prefer semantic tools); store_snapshot embedded in .ESS saves + store_restore on load (S4 start); test hooks gained a widgets op. **Two latent bugs found and fixed en route**: (a) store_read on a child view with ids=NULL returned empty (names() read on the child env); (b) **observer GC** — store-glue observers whose only dependencies are weakly held by the reactive graph are garbage collected between flushes, silently killing later store→UI pushes (reproduced deterministically with gc(); explains the sporadic Tier A stress flakiness). Fix: create epoch reactives once and keep observer references in module-level lists (heatmap + meta_scatter glue); select/enum validators now append allowed values alongside suggestions; store_restore is per-key resilient (strict=FALSE) |
| 1⅞+½ | **§6 control plane S4 (heatmap completion): multi_select kind + sorting/clustering/annotation widgets** | **DONE (unit 48+44+19; Tier A 48/48 ×2)** | Store gained the `multi_select` kind (character vector; empty vector/empty JSON array/`""` all clear; entries validated against choices; literal `[]`/`null` sentinels stay omitted optionals). iheatmapModule now registers the **full 13-widget parameter panel** on all three instances (dataspace.{cor,expr,dyn}_heatmap.*): palette/scale/margins (S3) + col/row sorting (dynamic choices mirroring the UI, incl. HCL dendrograms and GS| gene sets), clustering distance/linkage (static), and annot_col/annot_row/tooltip_info multi_selects. UI→store sync for all 13 (multi-select clearing syncs character(0)), seeding covers the new keys, store→UI push routes server-side selectize ids (rowSortBy/annotRow) through updateSelectizeInput. **Robustness rule discovered**: choices providers must never `req()` — a shiny.validation condition raised inside choices_provider aborts the whole store_apply transaction (tryCatch(error=) can't catch it); heatmap providers read module reactives through defensive no-req helpers. Status-path snapshot restore for heatmap keys intentionally kept (legacy pre-S3 snapshots + seeding race) — both paths are idempotent; full re-route stays deferred until all modules are migrated |

**S2 regression post-mortem (user-reported):** the initial S2 cascade
required reactive_selector1/2 unconditionally, but five modules drive
triselectors WITHOUT the store (feature_general, fgsea, geneshot, tables,
attr4 pass restore-only selectors) — their cascades never fired and the
analysis panel stayed blank after selections. Fixed by
`reactive_selector1() %||% input$analysis` fallbacks. Test debt repaid:
new tests/test_triselectorCascade.R spies sendInputMessage to assert
cascades fire in BOTH regimes (validated: fails on the broken code), and
Tier A gained the missing right-panel guard (analysis-panel populate +
render after selection — the exact missed regression) |
| 1″ | **Unplanned: stale-internal-axes fix (user-reported, 19:49 session)** | **DONE** | meta_scatter's internal xax/yax never track manual triselector edits, so a restore targeting values the internal model already holds changed no reactive and never touched the widgets (quick-badge path was immune via its axisRequest bump). The restore path now bumps axisRequest too, and triselector_module accepts reactive_axis_request so analysis/subset/variable observers re-assert on version bumps. Reproduced: manual y=log.pvalue drift + volcano quick-view apply previously a silent no-op; now corrects. Manual edits still stick (no bump on user input) |
| 1′ | **Unplanned: mid-restore req-abort fix** | **DONE** | R/module_meta_scatter.R — `current_axes <- .scatter_axis_signature(isolate(v1()), isolate(v2()))` ran before the axis assignment; v1()/v2() are req(input$variable)-guarded, so during an in-flight triselector cascade the req silently aborted the restore observer and the requested axes were lost (reproduced: apply during init lands on defaults). current_axes now tryCatch-guarded; verified racy and settled apply paths |
| 2 | WP1 sections | M | auxi_agentAssistant.R, module_aiAssistant.R, L0 wiring, tests |
| 3 | WP3 spec round-trip | S | module_aiAssistant.R, tests |
| 4 | WP4 log summarizer | M | auxi_agentLogging.R, new test file |
| 5 | WP5b prompt workflows | S | module_aiAssistant.R |
| 6 | re-run benchmark (Tier B), compare | S | tests/e2e_agent/ |
| 7+ | WP6, WP7, WP8 (Phase 2) | M each | figures/app modules |
| last | WP9–WP12 (Phase 3) | M/L | new file(s) |

WP2 before WP1 deliberately: suggestions are self-contained and safe;
`sections` changes the tool contract the model sees and benefits from being
validated against a logged baseline that already exists.

Each WP is an independent commit + test run; nothing here touches non-agent
code paths or adds dependencies.

---

## 8. Explicit non-goals (rejected approaches)

- Tool per widget; generic `set_input`; generic R/JS/CSS generation.
- Model reading DOM/HTML/tooltips directly.
- Giant static UI description in the system prompt.
- Silent fuzzy-match acceptance of IDs/columns (suggestions only).
- `validate_figure_spec` as a separate tool (costs the same provider
  request as a failed create; WP2/WP6 attack the root cause instead).

---

## 9. Decisions (settled 2026-09-21)

| # | Topic | Decision |
|---|-------|----------|
| 12 | Store scope (Q1) | Single per-app store with hierarchical canonical ids; child views (`widget_store_child`) provide NS()-like namespace prefixes. Implementation lives in ONE file (`R/auxi_widgetStore.R`) organized by numbered comment sections — no file sprawl |
| 13 | Agent visibility (Q2) | Agent-controllability mirrors user-controllability exactly: every user-editable widget is agent-settable and discoverable; anything the user cannot change is at most internal snapshot state and is never exposed to the agent (registry views and tools exclude it — context hygiene). Semantic tools remain as convenience, never as access gates |

All eleven open questions resolved as follows; these are binding for
implementation:

| Q | Topic | Decision |
|---|-------|----------|
| 1 | WP1 overview | Include current scatter space + x/y axes; selection counts + ≤ 20 example IDs; no separate selection section |
| 2 | WP1 mechanism | Single `sections` array argument (no new tools); `available_sections` advertised in every overview; sections enumerated in tool description |
| 3 | WP4 exposure | Internal (`omicsViewer:::agent_summarize_log()`); consider export after log format stabilizes post-Phase-1 |
| 4 | WP6 templates | Ship `volcano`, `boxplot`, `histogram`, `scatter`; `density`/`barplot` log-gated; `pca` deferred |
| 5 | WP7 gate | Build patch-mode only if tasks 10–11 first-attempt success < 2/3 after WP1–WP3 |
| 6 | Selection semantics | Keep replace-absolute; document the get-state → union → set pattern in WP5b workflows |
| 7 | E2E provider | Pin one model across baseline and post-change runs; exported per-run via env; mid-tier model (measures tool-use competence, not model floor) |
| 8 | Flakiness | Majority-of-3 + triage: identical repeated failures = harness bug (fix); varying failures = model nondeterminism (outvote, archive artifacts) |
| 9 | Runtime budget | 5-task smoke (tasks 1, 4, 6, 8, 10) × 1 run per change; full 12 × 3 at phase gates; serial execution |
| 10 | WP11 scope | Transcript + figure registry + selections together; opt-in per snapshot; byte-cap transcript |
| 11 | Cost ceiling | Default unlimited-but-logged; `OMICSVIEWER_LLM_MAX_COST_USD` opt-in for admins |

---

## References

- Skill: `bioSkills/agent-development/shiny-ellmer-agents/` — SKILL.md
  (implementation rules, architecture choice), usage-guide.md §3 (canonical
  state bridge), §8 (history + auxiliary state), §11 (security checklist),
  §12 (review checklist → WP0 assertions).
- Current implementation: `R/module_aiAssistant.R`, `R/auxi_agentAssistant.R`,
  `R/auxi_agentFigures.R`, `R/auxi_agentLogging.R`, wiring in
  `R/L0_module_app.R`.
