# UI-mimic LLM evaluation — findings & improvement list (2026-10-17)

Method: a Playwright driver (`ui_mimic.mjs`) plays a realistic user through
the real chat UI (one app process per session, real provider glm-5.3-flash,
`OMICSVIEWER_LLM_LOG=true` into a durable dir). Sessions + full diagnostic
logs: `/home/chen/.omicsviewer-llm-eval/{batch1..4,minirepro,minirepro_fixed,probe1..4}`.
Driver: `node ui_mimic.mjs <turns.json> <outName>`; turns files:
`ui_mimic_turns1..3.json`. Probes: `ui_probe_transport.mjs`,
`ui_probe_fresh_scatter.mjs`, `ui_probe_scatter_seam.mjs`,
`ui_probe_scatter_timing.mjs`.

Evidence base: 4 multi-turn sessions (23 prompts attempted, 18 delivered),
2 minimal repros, 4 instrumented probes, 154 provider requests, 139 tool
calls logged.

## Fixed during the evaluation

### F1 — set_scatter_view dead on fresh sessions (FIXED, commit 2d4601f)
Every fresh-session "switch the scatter view" request failed 3/3 with the
misleading rejection *"Choices are currently unavailable for
dataspace.feature_space.x_analysis"*; the model then fell back to raw
set_widgets (worked) or gave up (batch1/batch2 PCA turns, 118s+ flailing).
Root cause: `apply_agent_scatter_view` ran `store_apply()` bare — in
ellmer's async tool task no reactive context is active, so the choices
providers threw and surfaced as "provider unavailable". Every sibling
path already isolates. Fix: wrap both applies in `shiny::isolate()`
(call-site only). Verified: same one-turn repro returns a clean receipt;
batch4's PCA request (switch to PCA quick view + color) completed with
real UI effects (quick-view badge action + `attr4.color_variable`).

## Open findings (priority order)

### F2 — silent message drop (top bug) — CONTAINED (WP-4/5/6: root cause + detect/notify + visible retry)
In 3 of 4 multi-turn sessions the chat wedges exactly after the WP13
janitor's `history_stubbed` events: the user's next message renders a
bubble client-side (shinychat optimistic render) but NEVER reaches the
server (no `user_message`, no `provider_request_start`, no error, no
disconnect overlay; websocket open).

**Root-cause verdict 2026-10-17 (WP-4 probe: `ui_probe_stub_drop.mjs`,
modes idle/race/midstream, 6 instrumented runs + 2 batch replays; logs
in ~/.omicsviewer-llm-eval/stub_probe_*):**

1. **Primary mechanism (client-side, deterministic, 4/4 runs):** while
   a stream is ACTIVE, shinychat's composer accepts Enter, renders the
   optimistic user bubble, and NEVER SENDS the input update — the probe
   observed `bubbleDelta≈300` (bubble in DOM) with `inputFrames=0` (no
   websocket frame) and no server event, ever. There is no queue and no
   drain: the message is silently lost. This fires on EVERY active
   stream (even with the cancel control invisible — `busy=false` at
   submit — so users cannot tell). The stubbing janitor was pure
   correlation (it runs after every settle); midstream runs with zero
   stubs and zero compactions still drop 2/2 queued messages.
2. **Amplifier (server-side, history-size-dependent):** when a long
   chain aborts at the session request limit with a LARGE history
   (~80 turns; batch2 + wp4_replay2 tails), the stream NEVER settles:
   the log ends at `request_limit_reached → history_stubbed` with NO
   `stream_failure`/`stream_status: idle` — the composer stays busy
   forever and ALL subsequent submissions are discarded (batch2's
   permanent death incl. the failed retry). With small histories
   (~12 turns; stub_probe_mid1–3) the same abort settles cleanly
   (`stream_failure` logged, status idle) and the composer recovers.
3. Hypotheses killed: janitor-leaves-bad-status (all `janitor_run`
   passes report `status_after=idle`, ≤153 ms), stub-racing (15 race
   probes clean), compaction interleaving (drops occur with zero
   compactions; survivors compacted most).

New flight-recorder events (kept): `chat_submit_seen` (rate-limited,
earliest server receipt marker) and `janitor_run` (on.exit post-action
status + turn count). Repro gate: `cd tests/e2e_agent && PROBE_MODE=
midstream OMICSVIEWER_LLM_MAX_REQUESTS=6 node ui_probe_stub_drop.mjs <n>`
→ exactly 2 of 5 drops in ≤5 min, 4/4 runs. Fix directions land in
WP-5 (make the abort always settle; visible busy state) + WP-6
(composer delivery watchdog covers the client-side discard class).

**WP-5 shipped (detect + notify; silent → loud):** submit watchdog —
15 s after `chat_submit_seen` with no provider request ⇒
`chat_stuck_detected` logged; 30 s ⇒ `chat_stuck_confirmed` + error
notification (reload guidance). Abort-settle check — 5 s after a
request-limit/budget abort with the stream still `streaming` ⇒
`chat_abort_stuck` + notification (the batch2 wedge shape). No
automatic re-drive (shinychat exposes no server-side cancel for a hung
stream; double-submit risk) — recovery stays user-initiated, and
WP-6's composer Retry covers it. Validated live: healthy turns and 4/4
limit aborts produced zero false stuck events; the 2 client-side
drops per repro run are unchanged BY DESIGN until WP-6 makes them
visible + retryable.

**Program close-out (WP-9 replay, 2026-10-18):** the full batch1–4
replays (replay_batch1..4 in ~/.omicsviewer-llm-eval/) show the F2
class contained: 0 `chat_stuck_detected` events, 0 permanently-wedged
sessions, every session ran to its last turn (the baselines died
mid-session in batch2/batch4), and the single driver-flagged
submission (replay_batch1 t1 — the driver's known false negative on
slow first tokens, handover quirk 3) round-tripped on its one retry.
23/23 turns delivered, 0 timeouts (baseline: t7-figure TIMEOUT,
b2-down-enriched DROP). Known residue: the client-side discard while
a stream is active still needs the user to notice the WP-6 notice and
press Retry (by design — no auto re-drive); provider-side 500s
(bigmodel 网络错误) end a turn visibly with the error text in the
transcript.

### F3 — no cheap column enumeration ("colored by group" flail) — FIXED (WP-1)
Every session burned 5–8 tool calls guessing synonyms (group, genotype,
condition, class, cell, type...) because `search_annotations` matches
only literal substrings and `get_omics_viewer_state(annotations)` is
hard-truncated at 8,192 bytes (19,811 wanted). The model recovers only
by luck (state catalog excerpt) — batch1 t2 never recovered and the PCA
never happened. Fix: allow an empty/broad `search_annotations` query to
return a grouped catalog (category → count + example columns), bounded
to fit the tool budget; mention it in the tool help.

**Resolution 2026-10-17 (WP-1).** `agent_annotation_category_catalog()`
groups columns by `Category|Subcategory` prefix (unprefixed columns
stand as their own single-column group), bounded to ~4 KB by degrading
tail categories to counts-only; empty/whitespace query on
`search_annotations` returns it (tool description + capability record
updated). Live validation (`wp1_validation`, `wp1_validation_b` in
~/.omicsviewer-llm-eval/): both sessions enumerated in 2
search_annotations calls (1 blind synonym miss → empty-query catalog →
straight to `General|All|Origin`); the catalog digest is 3,489 bytes
for demo.RDS (21 feature + 3 sample categories); PCA scatter turns
settled in 98 s / 54 s and the app scatter switched via
`set_scatter_view(sample, pca)`.

### F4 — panel result tables are unreadable — FIXED (WP-7)
After a successful fGSEA run the assistant cannot report the top-3 NES
pathways: no tool exposes enrichment/table results (it checked state,
widgets, controls, actions — then gave an exemplary honest refusal:
"I won't guess or reconstruct them"). Same class: DT column sort is
gesture-only.

**Resolution 2026-10-17 (WP-7).** New `panel_results` section for
`get_omics_viewer_state`: the ORA/fGSEA modules expose their current
result tables as attributes on their module returns (the data-space
quickViews pattern; fgsea behind the same R-H3 visibility gate the
outputs use, so the computation never becomes eager), the L1
result-space module builds a bounded `read_enrichment_results()`
reader, and `agent_enrichment_table_excerpt()` byte-shapes the payload
(key columns first — pathway, NES/p-value — values rounded, list
columns dropped, ~3 KB cap with `columns_available` when truncation
drops columns). Uncomputed panels return explicit `not_run` notes
pointing at `set_enrichment_parameters`; snapshot/state surfaces never
see raw tables (attributes are dropped by `as.list()`/sanitize).
System prompt, tool description, capability records and the ORA/fGSEA
panel prose now teach the path; column sorting stays click-only.
Guards: new tests/test_agentPanelResults.R (26), Tier A 115/115 with
panel_results cases (excerpt rows match the browser's own fGSEA
table), agentAssistant/aiAssistantTools/agentCapabilities/
agentContracts/agentToolSeam/agentWidgets/ora/snapshotRoundTrip green.
Live (`wp7_validation`): the fGSEA top-3-by-NES turn = 4 tool calls,
76 s, three REAL rows (cell.cell.adhesion NES 2.11 / size 72,
cell.junction.assembly, angiogenesis) matching the app's table
(verified against an independent fgsea run — fgseaMultilevel
stochasticity explains ±0.03 NES between runs); the follow-up quoted
real padj/size numbers "directly from the fGSEA result table". The
honest-refusal class is extinct for enrichment tables; replay batch4
b4-fgsea-read: 246 s → 87 s with the correct search →
set_enrichment_parameters → panel_results chain (one provider-side
500 cut that live run short — provider flake, same prompt clean in
wp7_validation).

### F5 — find_controls rejects empty queries — DIGEST SHIPPED (WP-2)
The model's opening move ("enumerate the controls") is rejected with
"A non-empty capability search query is required" in every session.
Allow empty query → top-level groups (WP-G1 records already carry them).

**Resolution 2026-10-17 (WP-2).** `agent_capability_digest()` groups
all records by owning tab with control/action/panel counts, cascade
groups, and covering tools; empty/whitespace/NULL query on
`find_controls`/`agent_capability_search` returns it (tool description,
query-arg description, and capability help_text updated). The sentinel
consequence changed deliberately: sanitized NULL now enumerates instead
of erroring (bounded, read-only). Live validation (`wp2_validation`):
the enumeration opener used ONE find_controls call (empty query →
digest) and produced an accurate panel tour in 53 s.

**Still open (amendment, needs Chen's call — plan freezes non-empty
behavior):** fuzzy matching drowns precise queries — live evidence
(wp1_validation_b, wp2_validation t2): `find_controls("heatmap
clustering")` → 38 KB truncated result; `"dataspace.sample_space.color"`
→ 207/214 fuzzy matches with the 3 literal matches buried; drill-downs
burn 5–10 retries. Candidate: literal matches rank before fuzzy ones
and truncation keeps the literal head; OR fuzzy only when the query has
no id-separator characters.

**Amendment shipped 2026-10-17 (Chen-approved).** Matching semantics
unchanged; results now (1) rank fuzzy hits by matched-word count
(the 3 true color widgets surface first instead of past the cut),
(2) byte-shape the payload — full records while they fit, then slim
id/panel/label/operation/kind/tab records (describe_control(id)
returns full), stepping the count down to a floor of 10 so the codec's
opaque 8–38 KB truncation marker never fires again, and (3) id-shaped
fuzzy misses get a hint pointing at prefix/meaning queries. Live
re-validation (wp2a_validation): "sample scatter color" → 1 call,
color cascade ranked first, correct disambiguation question (Origin
vs TP53.Status — "group" is genuinely ambiguous in demo.RDS), and the
follow-up "use Origin" turn SET the color (`color_variable:
--select-- → Origin` via set_widgets) — the exact outcome the baseline
never reached; "heatmap cluster" → 1 call, 63 s (was ~10 calls,
38 KB truncation, 120 s). Sweep gate PASS (find_controls cases
assert ids/labels, which slim records keep).

### F6 — figure filter errors don't list usable columns — FIXED (WP-3)
`update_figure` failed with *"Figure filter column is unavailable:
featureId"*; the model self-recovered by retrying with Gene.name, but
the error should name the figure's available columns (same
self-describing-rejection policy as WP-A).

**Resolution 2026-10-17 (WP-3).** The filter-miss error now appends a
bounded available-column note (first 12 columns + "+N more") via
`.agent_available_columns_note()`; and `featureId`/`sampleId` resolve
to the reserved `__feature_id__`/`__sample_id__` columns when present
(the annotation catalog labels id columns featureId/sampleId, the
figure grammar reserves the dunder spelling — both now work, the
reserved name stays canonical). Live validation (`wp3_validation`):
the highlight-genes turn used the canonical `__feature_id__` spelling
from the figure grammar with ZERO failed filter calls (was 1); the
vocabulary note + alias are defense in depth for models that guess.

### F7 — long silent chains, no progress feedback — FIXED (WP-8)
Median provider request 12.6 s (max 42 s); one user turn chained 33
requests over 4+ minutes with nothing user-visible between tool calls
(ORA pipeline). Recommendation: interim streaming progress notes or a
soft per-turn tool-call budget with a checkpoint question.

**Resolution 2026-10-17 (WP-8).** Two shipped pieces: (1) a DOM-only
progress line under the composer, driven by
`omicsviewerAssistantProgress` custom messages published from the
events the module already logs — submission ("Working — waiting for
the model..."), provider request start ("Working — the model is
responding..."), every tool call ("Working — step N: <tool>", per-turn
counter reset on each user submission), and stream settle (line
hides, including governed aborts). No shinychat internals touched.
(2) a soft budget in the system prompt ([PACE]: prefer about 6 tool
calls per user request; report progress and confirm before long
continuations) — deliberately no hard cap (ORA-class pipelines are
legitimate long chains). Browser-verified 10/10
(`verify_wp8_progress.mjs`: submit→model→idle message sequence
against a dead endpoint, step rendering, hook-driven idle clear);
Tier A 120/120. Live (`wp8_validation`): open-ended ORA turn 17 tool
calls / 238 s with the line visible 48 of 79 poll rounds ('step 11:
find_controls' sampled mid-chain) and every ≥15 s gap carrying a
visible Working notice; fGSEA follow-up 4 calls / 64 s. Replay
tool-call counts on the baseline's heaviest turns: batch1 t5-ora
26→5, batch4 pca-color 31→7, batch1 t2-pca 12→1 (soft budget +
WP-1/2 discovery cuts the flail, not a cap — chains that legitimately
need more calls still run).

### F8 — composer can swallow the FIRST submission too — FIXED (WP-6)
minirepro turn 1 waited >25 s before the message registered (retry
saved it); same silent-swallow class as F2 but pre-stub. Delivery
feedback in the composer is needed regardless.

**WP-4 note:** same client-side discard class as F2 mechanism 1 —
any submission racing an active stream is dropped without a frame;
the watchdog (WP-6) makes it visible and retryable.

**Resolution 2026-10-17 (WP-6).** A composer delivery watchdog ships
with the assistant UI: every Enter arms a 10 s timer; the server acks
every received submission (`omicsviewerAssistantDeliveryAck` custom
message fired by the `chat_submit_seen` path), and the ack clears the
watchdog. No ack within 10 s ⇒ inline notice — "Message not delivered
- the assistant may still be finishing. Retry when it is idle." —
whose Retry posts the stashed text to the server, which re-drives the
submission through shinychat's own `update_user_input(submit = TRUE)`
(owned APIs only; no synthetic events; a busy stream defers the retry
with a status message instead of double-submitting). Browser-verified
7/7 (`verify_wp6_delivery_notice.mjs`: notice at 10 s, retry clears +
delivers, normal submits never show it). Live end-to-end (wp6_validate
midstream probe): a queued-during-chain discard surfaced the notice and
Retry delivered it — `drops 1 of 5, notice=true, retryRecovered=true`.
Known limitation: the watchdog tracks the most recent armed
submission, so several messages queued during one long chain surface
one notice/pending at a time; the rest are re-armed on their next
Enter.

## What worked well (keep)
- Refusal quality: CSV export → clean upfront refusal listing what it
  CAN do (excluded-by-policy surface holding up in real use).
- WP-C action tier: `run_control` built the STRING network end-to-end.
- WP-G1 capability records: "The StringDB tab has four changeable
  things plus one button" — discovery reads accurately.
- WP-A self-describing rejections: wrong id `feature_table.columns`
  → closest match suggested, model recovered in one step.
- Full pipelines complete: volcano figures with highlights (fig layer
  retry), ORA (area-select → enrichment → row select), heatmap
  clustering, PCA+coloring (post-F1).
- Sanitizer seam: glm's literal `"null"`/`"{}"` strings never reached a
  validator (set_scatter_view "null" axes sanitized correctly).

## Program close-out (WP-9, 2026-10-18)

All findings F1–F8 resolved. Branch `llm-ui-hardening`, one commit per
WP: WP-1 57fafeb (F3 catalog), WP-2 8b807bf + amendment 0bafe8c +
follow-up e8554e7 (F5 digest/ranking/byte-shaping), WP-3 593e338
(F6 column vocabulary), WP-4 afc6226 (F2 repro + flight recorders),
WP-5 680bcdf (F2 detect + notify), WP-6 54ef9b6 (F8 composer
watchdog), WP-7 75f153a (F4 panel_results), WP-8 87445ec (F7 progress
line + soft budget), WP-9 5efa0f3 (sweep-gate repair) + this
validation.

**Full replay outcome (replay_batch1..4, 23 turns):** ok 23/23
(baseline: 18 delivered, 1 TIMEOUT, 1 permanent DROP, 2 sessions
died mid-run). Duration per turn comparable or better (b4-fgsea-read
246→87 s, t6-heatmap 210→57 s, b3-fgsea 135→78 s); tool-call counts
on the baseline's heaviest turns drop 2–5× (26→5, 31→7, 12→1) with no hard cap; the progress line was
live during chains (ui_mimic now records progress.visible_rounds and
samples per turn).

**Guard board at close (all green, counts vs 2026-10-17 plan):**
agentAssistant 111 (99+), aiAssistantTools 53 (48+),
agentCapabilities 62 (47+), agentContracts 7, agentWidgets 104,
agentFigures 137 (133+, WP-3), agentToolSeam 48, agentActions 49,
agentUiEffects/Tier A 120/120 (113/112+), new test_agentPanelResults
26, widgetStore 69, snapshotRoundTrip 63, renderStability 46,
scatterSelection 40, quickViews 24, scatterDefaults 22, shinyAuxi 10,
appState 51, ora 37; **sweep gate PASS (344 cases)** after the WP-9
gene-picker vocabulary repair (pre-existing e8554e7 regression,
bisected with a two-build probe; content-only fix, no baseline
refresh).

**Deviations & notes for review:**
- The [PACE] soft budget under-delivers on exploratory turns that
drill find_controls (wp8_validation t1: 17 calls — 7 of them
find_controls hunting a selection control). No hard cap per plan;
revisit only if live use shows users bothered by it.
- The F2 client-side discard class is contained, not eliminated:
during an active stream an Enter still needs the user to notice the
WP-6 notice + Retry. A queue/drain would need shinychat changes
upstream.
- Provider-side 500s (bigmodel) end turns visibly; the model does not
retry them (it reports the error) — acceptable.
- ui_mimic driver's 25 s bubble check false-negatives on slow first
tokens and re-delivers; both deliveries reach the server when idle
(the duplicate is absorbed or runs twice — harmless for evaluation,
noted in handover quirk 3).
- Version/NEWS decision left to Chen (suggested: patch bump within
  1.15.x or 2.2.0 given the feature additions).
