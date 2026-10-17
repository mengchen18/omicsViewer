# UI-mimic findings — remediation plan (F2–F8)

> **STATUS 2026-10-18 (WP-9 close-out): COMPLETE.** All WPs shipped on
> branch `llm-ui-hardening` — per-WP resolution records with live
> validation numbers live in `ui_mimic_findings.md` (F1–F8 all
> resolved); full replay 23/23 turns ok, sweep gate PASS. See the
> "Program close-out" section there for the final guard board and
> deviations.

Source of truth for the findings: `tests/e2e_agent/ui_mimic_findings.md`
(2026-10-17 evaluation; evidence in `~/.omicsviewer-llm-eval/`). F1 is
already fixed (commit `2d4601f`, isolate in `apply_agent_scatter_view`).
This plan dispositions F2–F8. Style follows `widget_sweep/FIX_PLAN.md`:
decision-complete work packages, guards named up front, one commit per WP.

**Program rules** (carried over from the widget-sweep program):
- Branch `llm-ui-hardening` off `master`; ONE commit per WP; no push, no
  merge, no RELEASE_3_22 cherry-picks without Chen's word.
- Guards run via plain `Rscript tests/test_X.R`; quote the real plan line.
  Browser suites with `NOT_CRAN=true`. Reinstall the package after every
  R change.
- Every WP ends with a **live mimic validation** turn through the real
  chat UI (`ui_mimic.mjs` + a turns file committed with the WP) — the
  findings came from real sessions, so the fixes are proven there, not
  only in unit tests.
- Blast-radius rule (AGENTS.md): prefer call-site patches, new functions,
  defaulted parameters. The clusters to respect: widget-store control
  plane (`auxi_widgetStore.R`), agent sanitizer seam
  (`auxi_agentAssistant.R` validators), triselector cascade internals.
  Where a WP must touch one, the guard budget is stated BEFORE the edit.

Known-good guard counts (2026-10-17): agentAssistant 99, aiAssistantTools
48, agentCapabilities 47, agentContracts 7, agentWidgets 104, agentFigures
133, agentToolSeam 48, agentActions 49, agentUiEffects 113 + Tier A
112/112, shinyAuxi 10, scatterSelection 40, quickViews 24, scatterDefaults
22, widgetStore 69, snapshotRoundTrip 63, renderStability 46, sweep gate
PASS (344 cases). Any count change is justified in the WP commit.

---

## Phase 1 — model-quality quick wins (independent, low risk)

### WP-1 (F3): annotation catalog — `search_annotations` enumeration path
**Problem.** "Colored by group"-style prompts burn 5–8 blind synonym
searches (observed in all 4 sessions; batch1 t2 never recovered):
`search_annotations` matches literal substrings only, and the
`annotations` state section is hard-truncated at the 8,192-byte tool
budget (19,811 bytes wanted), so the model has no way to SEE the column
vocabulary.

**Goal.** One call that answers "what columns exist", bounded to the
tool budget, so enumeration takes ≤2 tool calls in a live session.

**Design.**
1. New leaf helper `agent_annotation_catalog(feature_data, sample_data)`
   (in `R/auxi_agentAssistant.R` next to `agent_search_annotations`,
   ~line 666): for each space, group columns by `Category|Subcategory`
   prefix → `n_columns`, up to 8 example column full-names, and a type
   hint (numeric/categorical) per example. Hard byte cap (~4 KB total,
   well under the 8,192 tool budget); over-budget categories list counts
   only, with a hint to query the category name for details.
2. Entry branch in the `search_annotations` tool handler: empty/whitespace
   `query` (post-sanitizer it arrives as `""`) returns the catalog instead
   of erroring. Literal-search behavior for non-empty queries is
   byte-identical — the branch sits BEFORE the existing `!nzchar(query)`
   stop (auxi_agentAssistant.R:670), replacing it.
3. Tool description + capability record (`tool:search_annotations`) gain
   one sentence: "Pass an empty query to list the annotation catalog by
   category before searching."
4. `search_hint` plumbing (auxi_agentAssistant.R:840) unchanged.

**Files.** `R/auxi_agentAssistant.R` (helper + branch), tool registry
description (`R/auxi_agentToolRegistry.R`), tests.

**Blast radius.** Additive branch at one validator entry + new function —
no cluster internals. The sanitizer seam tests (agentToolSeam 48) must
stay green unchanged.

**Guards.** Before/after: agentAssistant, aiAssistantTools, agentToolSeam,
agentCapabilities. New tests: empty query → catalog with both spaces,
byte-bounded, no literal-path change for "group"/"EGFR" style queries.

**Acceptance.** (a) Unit: catalog ≤4 KB for demo.RDS, both spaces present.
(b) Live: fresh mimic session, "Show a PCA scatter plot of the samples,
colored by group." → enumeration completes in ≤2 annotation calls (was
5–8), turn settles <120 s, and the scatter actually switches.

**Effort.** ~0.5 d. **Risk.** Low — worst case the model ignores the
catalog; the live turn decides.

### WP-2 (F5): `find_controls` empty query → group digest
**Problem.** The model's opening move ("enumerate the controls") is
rejected: `A non-empty capability search query is required`
(auxi_agentCapabilities.R:329) — a wasted call in every session.

**Design.** Empty query returns a bounded digest built from the WP-G1
records: per owning tab → control count + action count + cascade groups
present, plus one line of usage hints ("query by widget meaning, e.g.
find_controls('heatmap clustering')"). Non-empty literal/fuzzy paths
unchanged.

**Files.** `R/auxi_agentCapabilities.R` (validator + digest builder),
tests.

**Blast radius.** Same class as WP-1: additive branch; capability-count
fixtures may need deliberate updates (panel/tool counts unchanged — the
digest is a new return shape for one input, not new records).

**Guards.** agentCapabilities, aiAssistantTools, agentAssistant.

**Acceptance.** Unit: empty query returns digest listing ≥ the 6 top tabs;
live: a discovery turn uses ONE find_controls call where sessions
previously opened with reject-then-retry.

**Effort.** ~0.3 d.

### WP-3 (F6): figure filter errors carry the column vocabulary
**Problem.** `update_figure` failed with `Figure filter column is
unavailable: featureId.` (auxi_agentFigures.R:551). The site already
appends a closest-match suggestion (`.agent_suggest_text`) and the model
recovered in one retry — but the message still hides the vocabulary.

**Design.**
1. At the error site, replace the bare suggestion with
   `.agent_suggest_text(column, colnames(data))` PLUS a bounded
   available-column note (first 12 columns + "+N more"). Build it with a
   new tiny helper `.agent_available_columns_note(cols, max = 12L)` —
   error-string construction only, no validator logic change.
2. Optional (decide during implementation): accept `featureId` as an
   alias for the feature-id column when the figure's data source is
   `feature_annotation` — figures lose nothing, the alias is the model's
   natural vocabulary. If it complicates the data-source plumbing, drop
   it; the note alone halves the cost of the miss.

**Files.** `R/auxi_agentFigures.R` (error site + helper), tests.

**Blast radius.** `.agent_figure_*` validators have fan-in 22–26 — but the
edit is at ONE stop() call's message; the agentFigures battery (133 tests)
runs before/after; any golden error-string fixtures are updated
deliberately with per-change justification.

**Guards.** agentFigures, agentToolSeam, agentAssistant.

**Acceptance.** Unit: the error message for a bogus column lists real
columns; agentFigures count unchanged or fixture diff justified. Live:
the highlight-genes turn completes with ≤1 failed filter call (observed:
1).

**Effort.** ~0.3 d.

---

## Phase 2 — transport reliability (F2, the top bug; F8)

### WP-4 (F2/F8): repro + instrumentation (no fix yet)
**Problem.** After WP13 janitor `history_stubbed` events, the next user
message is silently dropped in 3 of 4 sessions: bubble renders
client-side, server never receives it (no `user_message`, no
`provider_request_start`), websocket open, no disconnect overlay. NOT
deterministic yet: batch3 survived 3 stubs + 4 compactions — stubbing is
necessary-but-not-sufficient. F8 is the milder sibling: the FIRST
submission of a fresh session also lagged >25 s once (minirepro).

**Goal.** A deterministic ≤5-minute repro, plus enough diagnostic signal
to discriminate the hypotheses below.

**Hypotheses & discriminating experiments.**
| # | Hypothesis | Discriminator |
|---|------------|---------------|
| H-a | Janitor's `client$set_turns()` at idle leaves shinychat's ExtendedTask believing a stream is pending → composer submissions ignored | Probe `chat_object$status()` + shinychat busy state at the client right after each stub event |
| H-b | Stubbing DURING an unsettled flush (vs true idle) corrupts the pending-turn queue | Force stub only after N s of verified idle; if drop disappears, timing is the trigger |
| H-c | Summary/compaction requests (coro) interleaving with stubbing wedges the task queue | batch3 survived compaction+stub; reproduce with compaction forced ON/OFF |
| H-d | Client-side composer queue disables submit after specific server-published state; bubble is purely optimistic | Instrument submit handler: log client submit attempts + what shinychat emits on the websocket for a dropped message |

**Design.**
1. `tests/e2e_agent/ui_probe_stub_drop.mjs`: fresh app (hooks + logging),
   scripted heavy turns (`create_figure` ×2–3 with wide data, one ORA
   turn) to force ≥1 stub, then N=5 probe messages ("reply with the word
   X") at 30 s spacing; per message: server receipt yes/no, ws frame
   observation (Playwright `page.on('websocket')` framereceived),
   cancel-control visibility, `stream_status` tail of the jsonl. Loop
   until a drop is observed; print the exact event timeline around it.
2. Server diagnostics (small, temporary-but-kept): new log events
   `janitor_run` (stub details + post-action `chat_object$status()`) and
   `chat_submit_seen` at the earliest server entry for user input (the
   existing `input$chat_user_input` observer already sits there — add
   the event, rate-limited). These events stay: they are the transport's
   flight recorder.
3. Re-run the three mimic batches' turn sequences against the
   instrumented build; correlate every drop with the new events.

**Files.** `R/module_aiAssistant.R` (two log events), probe script,
`ui_mimic_findings.md` updated with the root-cause verdict.

**Guards.** agentAssistant, agentContracts (log-event surface), agentUiEffects.

**Acceptance.** Either (a) a repro script that drops ≥1 of 5 probe
messages in ≤5 min, run 3× with ≥2 drops total, or (b) a documented
non-repro PLUS a discriminating observation that kills ≥2 hypotheses.
No fix ships from this WP.

**Effort.** ~1 d. **Risk.** The bug may need provider-specific timing —
the probe must use the real provider (provider.env), so runs cost tokens.

### WP-5 (F2): the fix (decision-gated on WP-4)
**Goal.** Post-stub (and generally: always) a submitted message either
reaches the server or the user sees why it didn't — no silent loss.

**Candidate fixes, in preference order (pick by WP-4 verdict):**
1. *If H-a/H-b (janitor timing):* move the stub/compact action from the
   idle observer into a guarded transaction: after any janitor action,
   verify liveness by evaluating the chat input path once (a zero-cost
   sentinel: `isolate(chat_object$status())` + a `chat_heartbeat` input
   echo the client answers); if the heartbeat doesn't round-trip within
   3 s, mark the chat object for rebuild-on-next-submit and log
   `chat_recovered_by_rebuild`. The rebuild re-uses the archived context
   (WP13 already keeps stubs + archive ids).
2. *If H-d (client queue):* server-side nudge — publish an explicit
   `stream_status: idle` (or shinychat's equivalent) after janitor
   actions so the client re-enables submission.
3. *Regardless of verdict — belt and braces:* hard server-side watchdog:
   the `chat_submit_seen` flight recorder gains a timer; a submission
   seen with no `provider_request_start` within 15 s logs
   `chat_stuck_detected` and triggers candidate 1's rebuild path.

**Files.** `R/module_aiAssistant.R` (janitor section, ~800–1000), maybe a
small client JS nudge; tests.

**Blast radius.** module_aiAssistant janitor is self-contained (WP13
seam) — NOT one of the three AGENTS.md clusters; still, the observer-GC
rule applies (any new observer joins `.context_keep`).

**Guards.** agentAssistant, agentContracts, agentUiEffects + Tier A,
agentToolSeam; the WP-4 probe becomes the regression gate (0 drops over
3 runs × 5 messages); full sweep gate if any store-adjacent line moved
(not expected).

**Acceptance.** WP-4 probe: 0 dropped messages across 3 consecutive
runs; a 5-heavy-turn mimic session delivers ALL subsequent prompts
(batch1/2/4 replay); no behavior change for short sessions (probe1
6/6 baseline unchanged).

**Effort.** ~1.5–2 d. **Risk.** Highest of the program — mitigated by
the gate; if neither candidate lands cleanly, ship candidate 3 alone
(detect + recover) and keep the silent-drop class extinct by
construction.

### WP-6 (F8): delivery feedback in the composer (ships regardless)
**Problem.** Independent of F2's root cause, the composer gives zero
feedback when a submission doesn't start a stream (first-submission lag
observed; post-stub drops pre-WP-5).

**Design.** Small client-side watchdog injected with the assistant UI:
on submit, a 10 s timer; if no streaming indicator and no new assistant
content, show an inline notice under the composer — "Message not
delivered — the assistant may still be finishing. [Retry]" — where
Retry re-fires the composer submission (and, post-WP-5, the server
watchdog has usually already recovered). DOM-only addition + one
`shiny:conditionalPanel`-free custom message binding; no shinychat
internals forked.

**Files.** `R/module_aiAssistant.R` (UI + JS), tests.

**Guards.** agentUiEffects (+ new case: notice appears when the server
is stalled via the test hooks), agentAssistant.

**Acceptance.** Browser-verified: stalled server (hook-simulated) →
notice within ~10 s; Retry delivers; normal turns never show it.

**Effort.** ~0.5 d.

---

## Phase 3 — capability

### WP-7 (F4): `panel_results` — readable enrichment tables
**Problem.** After a successful fGSEA run the assistant cannot report
top-3 NES pathways — no tool exposes panel result tables (it checked
every channel, then refused honestly). DT column sort is also
gesture-only.

**Design.**
1. New accessor `read_enrichment_results(method, n)` implemented in the
   result-space module (next to the enrichment reactives it already
   owns): returns the head(n, default 15) of the ORA overlap table /
   fGSEA table with key columns only (pathway/id, size, NES or p/q
   values as applicable), values rounded, byte-capped (~3 KB), plus
   `columns_available` when truncation drops columns.
2. Expose through progressive disclosure, NOT a new tool:
   `AGENT_STATE_SECTIONS` (constants.R:199) gains `"panel_results"`;
   `get_omics_viewer_state(sections="panel_results")` calls the accessor
   through the existing app-module wiring (the same callback pattern as
   `apply_enrichment`). Empty/uncomputed results return an explicit
   "not run yet — see set_enrichment_parameters" note.
3. Capability prose (WP-G1 `AGENT_PANEL_PROSE` ORA/fGSEA records) gains
   one sentence: results are readable via the `panel_results` state
   section; sorting is click-only.
4. STRING enrichment table joins the accessor only if trivial; otherwise
   explicitly out of scope for this WP.

**Files.** `R/constants.R`, `R/auxi_agentAssistant.R` (state builder),
`R/L1_module_result_space.R` + `R/L0_module_app.R` (accessor wiring),
tests.

**Blast radius.** State builder + L0/L1 wiring are call-site level; no
store/triselector internals. agentCapabilities fixtures get deliberate
updates (section list + prose strings).

**Guards.** agentCapabilities, agentAssistant, aiAssistantTools,
agentContracts; new unit tests for the accessor (computed vs uncomputed
panels); sweep gate untouched (no store change).

**Acceptance.** Live: "Run fGSEA on RE vs ME and tell me the top 3
pathways by NES" returns three REAL pathway rows matching the app's
table (verified against the browser table via hooks); the honest-refusal
path only remains for genuinely unrun panels.

**Effort.** ~1 d. **Risk.** ORA/fGSEA result reactives must be reachable
from the app module — if the module scoping fights back, fall back to a
bounded `results` field inside `apply_enrichment`'s receipt (one step,
same information).

---

## Phase 4 — perceived latency (F7)

### WP-8 (F7): visible progress + soft tool-call budget
**Problem.** Median provider request 12.6 s (max 42 s); one user turn
chained 33 requests over 4+ minutes with nothing user-visible between
tool calls.

**Design.**
1. Status line under the composer ("Working — step 3/…:
    set_enrichment_parameters"), driven by the existing
   `provider_request_start`/`tool_request` events already logged —
   publish a small reactive to the client; DOM-only, no shinychat
   fork. Step counter resets per user turn.
2. System-prompt guidance (constants where the prompt is assembled):
   "Prefer at most ~6 tool calls per user request. If a task needs
   more, complete the core step, report progress, and confirm before
   long continuations." Soft only — no hard cap in this WP (a hard cap
   risks chopping legitimate ORA-style pipelines; revisit only if the
   soft rule under-delivers).
3. Optional stretch (decide at implementation): interim assistant
   commentary between chained calls (ellmer stream interleave) — only
   if it falls out naturally; otherwise skip.

**Files.** `R/module_aiAssistant.R` (status line + event publish),
prompt-assembly constants, tests.

**Guards.** agentUiEffects, agentAssistant, agentCapabilities (prompt
text may appear in fixtures — deliberate).

**Acceptance.** Live ORA/fgsea turn: longest user-invisible quiet
stretch ≤15 s (measured from the session jsonl timestamps + the status
line's DOM assertions in the driver); turn tool-call counts visibly
drop for open-ended prompts (catalog from WP-1 also helps here).

**Effort.** ~0.5–1 d. **Risk.** Prompt-text changes can shift model
behavior broadly — validate with the full batch1–4 replay, not just one
turn.

---

## Phase 5 — program close-out

### WP-9: full validation + documentation refresh
1. Reinstall; run the ENTIRE known-good guard board (all suites listed
   above) + the full sweep gate (`NOT_CRAN=true
   tests/e2e_agent/widget_sweep/run_regression.sh`) — any baseline drift
   justified per case, surgically refreshed.
2. Replay the batch1–4 turn files through `ui_mimic.mjs` against the
   fixed build; regenerate the per-turn outcome table; diff against the
   2026-10-17 baseline.
3. Update `ui_mimic_findings.md` with resolution status per finding
   (F2–F8), including measured before/after numbers (tool calls to
   enumerate, quiet-stretch length, dropped-message counts).
4. NEWS entry; version decision (patch bump vs 2.2.0) left to Chen.

**Effort.** ~0.5–1 d (much of it gate runtime).

---

## Sequencing & dependencies

```
Phase 1  WP-1 (F3)  WP-2 (F5)  WP-3 (F6)     [independent, any order]
Phase 2  WP-4 (repro) ──> WP-5 (F2 fix) ──> WP-6 (F8 watchdog)
                               (WP-6 independent of WP-5's verdict)
Phase 3  WP-7 (F4)                            [after Phase 1 ideally, else anytime]
Phase 4  WP-8 (F7)                            [last; benefits from all above]
Phase 5  WP-9                                  [close-out]
```

Total estimate: ~6–8 working days including gate runtime and live
validation. Phase 2 is the risk concentrator; Phases 1/3 are safe to
interleave if provider access is unavailable on a given day (WP-4/5
live runs need `provider.env`).

## Standing non-goals
- No downloads/export surface (cancelled by Chen, 2026-10-12).
- Snapshot lifecycle + dataset switching stay permanently excluded.
- No hard per-turn tool-call cap in this program (soft guidance only).
- No changes to the widget-store control plane or triselector cascade
  internals — every WP above avoids them by design.
