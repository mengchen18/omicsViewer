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

### F2 — silent message drop after context stubbing (top bug)
In 3 of 4 multi-turn sessions the chat wedges exactly after the WP13
janitor's `history_stubbed` events: the user's next message renders a
bubble client-side (shinychat optimistic render) but NEVER reaches the
server (no `user_message`, no `provider_request_start`, no error, no
disconnect overlay; websocket open). Repro: fresh app, 3–4 heavy tool
turns (figures/enrichment) until `history_stubbed` appears in the jsonl,
then send anything. Candidates: the janitor's `client$set_turns()` on
idle breaks shinychat's next submission (the module's own comments
document ellmer 0.5.0 set_turns desyncs); fix by verifying liveness
after stubbing (probe input path), deferring stubbing to the next
request start, and adding a client-side "message not delivered" notice
when no stream starts within ~10 s.

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

### F4 — panel result tables are unreadable
After a successful fGSEA run the assistant cannot report the top-3 NES
pathways: no tool exposes enrichment/table results (it checked state,
widgets, controls, actions — then gave an exemplary honest refusal:
"I won't guess or reconstruct them"). Same class: DT column sort is
gesture-only. Fix: bounded `panel_results` section for
get_omics_viewer_state (top rows of ORA/fGSEA tables) or a read tool.

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

### F7 — long silent chains, no progress feedback
Median provider request 12.6 s (max 42 s); one user turn chained 33
requests over 4+ minutes with nothing user-visible between tool calls
(ORA pipeline). Recommendation: interim streaming progress notes or a
soft per-turn tool-call budget with a checkpoint question.

### F8 — composer can swallow the FIRST submission too
minirepro turn 1 waited >25 s before the message registered (retry
saved it); same silent-swallow class as F2 but pre-stub. Delivery
feedback in the composer is needed regardless.

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
