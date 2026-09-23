# HANDOVER — agent accuracy & widget control plane

Written 2026-09-24 at the close of the WP7-gate session. Read this first
in a fresh context, then `AGENT_ACCURACY_PLAN.md` (§4 WP7 gate evidence +
benchmark findings are the source of truth) and `AGENTS.md` (environment,
commands, quirks). Delete or trim this file once absorbed.

## Where we are

Branch `agent-driven-exploration`, DESCRIPTION 2.1.1. Everything through
S1–S4 + WP1–WP6 is done (see git history). This session:

| What landed | Evidence |
|---|---|
| **WP7 gate run: patch-mode update_figure SHELVED (decision 5)** | full ×3 Tier B benchmark (12 core tasks, glm-5.3-flash): tasks 10+11 combined first-attempt **5/6 = 83% ≥ 2/3**; every completed revision used full-spec `update_figure` (the WP3 round-trip), zero spec failures; the single miss was a mid-turn provider stream hang (infra), not composition error |
| **WP6b: create_figure surface fixes** | (a) spec-present ⇒ spec authoritative, template args ignored + warning in result — the old "not both" hard error was never recovered from (27 rejected calls, up to 16 identical retries in one benchmark task); (b) templates accept `features`/`samples` ID subsets — the task-8 need (107 selected > 50-feature cap) was previously inexpressible on the template path and pushed models into the collision; volcano `label_top_n` ranks within an explicit subset; expression boxplot works with NO selection (features list supplied) |
| **Agent selection revert FIXED** | `apply_agent_state` patched `selection$features` but not the panel-status mirror; the data-space restore observer pushed the OLD selection back through `v1()` and overwrote `ri`/`rh` within one flush (deterministic: apply reports 5, overview reverts to the initial 107 in <2 s, on every tab). The apply now patches every mirror (`eset_selected_features/samples` + nulls stale table `rows_selected`). Verified live: selection sticks on Feature / Feature table / Sample tabs; MAPK task-106 first-attempt **3/3** (was impossible). This resolves the S4 "transient feature patches" note |
| **Benchmark infrastructure** | `tests/e2e_agent/tier_b_full.mjs` (multi-prompt tasks, fresh app+chat per task, request cap 8, log-aware settle ≈ 6 s) + `score_tier_b_gate.R`; scoring = log ground truth (tool errors per prompt span) + DOM observations |

MAPK wording variants (tasks 106/108) exist because "kinase" matches
nothing in demo.RDS (search returns 0 hits; honest models ask for
clarification — 0/3 on the canonical wording is a dataset artifact).
Task 108 is 2/3 outcome — the single miss was a clarification question
("sample group" is ambiguous: Cell.line/Gender/Origin all qualify), not
a tool failure.

## Validation snapshot (all green at HEAD)

- Unit: agentFigures **67** (7 new WP6b), aiAssistantTools 33 (precedence
  test re-pointed), agentAssistant 64, agentWidgets 105, agentLogSummary 32,
  agentLogging 21, appState 37, quickViews 17, scatterSelection 6,
  triselectorCascade 6, tableWidgetState 5, widgetStore 52, shinyAuxi 10,
  stats 7, ora 23
- Tier A browser: 93/93 ×2 (one intermediate 92/93 = the documented
  sporadic stress flake, clean twice after)
- Live Tier B: gate ×3 (artifacts/gate: logs + records + scores committed;
  bulk PNGs kept local only), MAPK ×3 (artifacts/mapk)

## Flakes / rules discovered (cumulative; do not re-learn)

- **Settle on LOG truth, not UI probes**: shinychat's cancel-control probe
  both misses silent streams AND lingers after idle; the diagnostic log's
  `assistant_response` followed by `stream_status: idle` (after the last
  user_message) is the authoritative turn-complete signal. A backward scan
  for this is subtly wrong (request_start precedes response in the normal
  ending — forward-scan the last indexes instead); the fixed predicate is
  in tier_b_full.mjs `logScan`/`logTurnComplete`, unit-checked against
  real logs (stream-hang log correctly stays in-flight).
- **Request cap is load-bearing**: `OMICSVIEWER_LLM_MAX_REQUESTS=8` in the
  spawned app stops runaway retry loops (16 identical create_figure
  rejections in one task otherwise burn 5 minutes and tokens).
- Provider stream hangs happen (task 11 run 1: preamble text, then silence
  for 5+ min; log ends at `provider_request_start`). The benchmark scores
  them as misses but the taxonomy should distinguish them (infra vs model)
  — noted for WP4 follow-up if it recurs.
- **Selection state has (had) multiple mirrors**: app `selection$features`,
  `panels$data_space$eset_selected_features`, and per-table
  `rows_selected`. Any external write must patch ALL of them before the
  `esv_status(NULL); esv_status(full_state)` transaction or the module
  restore reverts the patch. Keep this in mind for future apply paths.
- Empty-object + literal-string sentinels on optional tool args: every
  optional reader treats length-0, NA, and AGENT_SENTINEL_STRINGS like
  omitted; the new id-array param (`.agent_figure_ids_param`) follows the
  same rule and also normalizes tibbles/record lists.
- Tier A stress flake, MockShinySession quirks, observer GC, WebGL:
  AGENTS.md (unchanged).

## Next session: WP8

WP7 is closed (gate NO). Proceed to **WP8 — first new capability tools**
(`set_enrichment_parameters` for the ORA/fGSEA panel, `set_table_view`
for feature/sample table page+filter), each shipping help text in the
shared metadata structure from day one (WP9's one-source-of-truth
principle). The widget store already registers most of these widgets
(S4: resultspace.ora.*, resultspace.fgsea.*, dataspace.tab_*.*); the WP8
work is the thin semantic tier on top plus discovery help.

Housekeeping before starting: the stateful UI gap docs
(STATEFUL_UI_GAPS.md) and any stale plan sections (§7 table) were
updated this session — keep them that way per change.

Session workflow reminders (unchanged): reinstall after every source
change; roxygen after touching roxygen blocks (none touched this session
beyond params — check man/ diff); kill orphaned R processes on 7775–7795
before browser runs; Tier A ×2 convention; Tier B smoke subsets per
change, full ×3 at gates; log timestamps are UTC (user is UTC+2).
