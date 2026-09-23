# HANDOVER — agent accuracy & widget control plane

Written 2026-09-24 at the close of the WP6 session (Phase 1 complete since
the previous handover; WP6 now landed on top, one commit). Read this first
in a fresh context, then `AGENT_ACCURACY_PLAN.md` (§3 statuses + §4 WP6
notes are the source of truth; §6 control-plane background) and `AGENTS.md`
(environment, commands, quirks). Delete or trim this file once absorbed.

## Where we are

Branch `agent-driven-exploration`, DESCRIPTION 2.1.1. **Done and
validated:** S1–S4 control plane, WP1–WP5 (previous sessions; see the WP5
handover in git history for details), and now:

| WP | What landed | Notes |
|---|---|---|
| WP6 | Figure templates: `create_figure(template=, x=, y=, color=, label_top_n=, title=, space=)` with enum `volcano`/`scatter`/`boxplot`/`histogram` (decision 4: density/barplot log-gated, pca deferred). Server-side `agent_figure_template_spec()` in `R/auxi_agentFigures.R` expands to a full validated spec and flows through the exact WP3 render/round-trip path — template figures carry the echoed `spec` and stay revisable via `update_figure`. Templates advertised in the `figure_grammar` state section (`agent_figure_templates()`); prompt workflows updated (volcano + common-figures lines) | unit test_agentFigures 60, test_aiAssistantTools 33; Tier A 93/93 ×2; live glm-5.3-flash smoke: "volcano, label top 5" → ONE get_state + ONE template create_figure, first-attempt success |

Template semantics worth remembering:
- `space='feature'|'sample'` disambiguates column names present in BOTH
  annotation spaces (PCA columns are the classic collision); omitting it on
  an ambiguous column is a hard error naming both options.
- volcano `label_top_n` REORDERS `spec$features` by significance (y desc —
  the app's log.pvalue/log.fdr convention, mirroring meta_scatter's volcano
  detection) so the capped `label` layer (head-N semantics) marks the most
  significant features; scatter `label_top_n` labels the first N rows in
  data order (no ranking semantics — documented in the grammar).
- boxplot WITHOUT `y` = expression mode: `x` must be a sample column,
  features come from the current selection (empty selection → clear error).
- `label_top_n` on boxplot/histogram errors honestly (unsupported).
- New sentinel rule absorbed: `.agent_figure_{numeric,integer}_param` now
  treat literal sentinel STRINGS ("null" etc.) like omitted too (previously
  only length-0/NA); `.agent_figure_choice` errors now carry
  closest-match suggestions.

## Flakes / rules discovered (cumulative; do not re-learn)

- **Tier A "rapid successive applies converge" is a known sporadic flake**
  (45 s settle race in the 6-apply stress step): reproduced on PRE-WP6
  code (1/3 runs) and with WP6 (2/6 runs then 93/93 ×2). Do not chase it as
  a regression; only investigate if it fails with a NEW error signature or
  fails 3+ consecutive runs.
- Empty-object + literal-string sentinels on optional tool args: every
  optional reader treats length-0, NA, and AGENT_SENTINEL_STRINGS like
  omitted. Check this class FIRST when adding optional tool arguments
  (WP6 template args were built this way from the start — zero live
  sentinel failures).
- Tier B settle heuristics (shinychat cancel probe misses silent streams;
  log-aware activity; 18 s quiet after a figure) — all in tier_b.mjs.
- MockShinySession quirks, observer GC, WebGL: see AGENTS.md.

## Validation snapshot (all green at the WP6 commit)

- Unit board: agentFigures 60, aiAssistantTools 33, widgetStore 52,
  agentWidgets 105, agentAssistant 64, agentLogSummary 32, agentLogging 21,
  appState 37, quickViews 17, scatterSelection 6, triselectorCascade 6,
  tableWidgetState 5, shinyAuxi 10, stats 7, ora 23
- Tier A browser: 93/93 ×2 consecutive (plus the known sporadic stress
  flake documented above)
- Live Tier B: WP6 smoke only (volcano via template, first-attempt
  success; log + screenshot archived under tests/e2e_agent/artifacts/)

## Next session: WP7 gate decision, then WP7 or WP8

**WP7 is GATED** (decision 5): before building patch-mode `update_figure`,
run the full ×3 Tier B benchmark (12 tasks) and score tasks 10–11
("change the last figure's title", "color the volcano points by category")
first-attempt success. Build patch mode only if < 2/3. The merge base is
the WP3 registry (`figures()[[id]]$spec`, normalized shape) and the
registry now also records `template` (WP6) — a template-created figure can
be patched exactly like a spec-created one.

If the gate says NO: proceed to **WP8** (`set_enrichment_parameters` ORA/
fGSEA panel + `set_table_view` feature/sample table page+filter), each
shipping help text in the shared metadata structure from day one.

Session workflow reminders (unchanged): reinstall after every source
change; roxygen after touching roxygen blocks; kill orphaned R processes on
7775–7795 before browser runs; Tier A ×2 convention; Tier B smoke subsets
per change, full ×3 at gates; log timestamps are UTC (user is UTC+2).
