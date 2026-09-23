# HANDOVER — agent accuracy & widget control plane

Written 2026-09-24 at the close of the WP1–WP5 session (Phase 1 + WP5
complete; five commits, 9a452a9..c5cdbfc). Read this first in a fresh
context, then `AGENT_ACCURACY_PLAN.md` (§3 statuses are the source of
truth; §6 control plane background) and `AGENTS.md` (environment,
commands, quirks). Delete or trim this file once absorbed.

## Where we are

Branch `agent-driven-exploration`, working tree clean, DESCRIPTION 2.1.1.
**Done and validated:** S1–S4 control plane (previous sessions), and now

| WP | What landed | Commit |
|---|---|---|
| WP1 | `sections` + compact overview on `get_omics_viewer_state`; overview 1,185 B vs 18,693 B full on demo.RDS (15.8×); `agent_state` in L0 is an on-demand builder `state(sections)` — NOT a reactive (keep that contract for `ai_assistant_module` / `agent_test_hooks_module`) | 9a452a9 |
| WP3 | Figure spec round-trip: tool results carry the spec in the **echo shape** (`agent_figure_spec_echo()`, flat layer aesthetics — ellmer's `convert_from_type` drops schema-foreign keys like `mappings` from echoed args); `figures()` registry stores the **normalized** spec (WP7 merge base); normalizer is idempotent (accepts both shapes) and exempts a verbatim full sample set from the 200-sample cap | 3b0963f |
| WP4 | `agent_summarize_log(path)` / `agent_summarize_logs(dir)` + registered `print.omicsViewerAgentLogSummary`; 8-class error taxonomy over failed tool results + stream_failures; first-attempt success, retry recovery, most-rejected args; `tests/test_agentLogSummary.R` (32, incl. real archived-fixture parse) | acf5a6e |
| WP5/5b | Task 17 (snapshot round-trip) in tier_b_tasks.md; system prompt workflows + exact-ID contract; tier_b.mjs settle loop log-aware + quiet-after-figure; live glm-5.3-flash smoke validated WP1 overview + WP3 round-trip | c5cdbfc |

Live-smoke evidence (glm-5.3-flash, archived under
`tests/e2e_agent/artifacts/tier_b_*.png` + logs): orientation task = ONE
`get_omics_viewer_state` call returning the 1.6 KB overview incl.
`scatter_view`; create→update revision round-trip reuses the echoed spec
and self-corrects through validation errors.

## New rule discovered this session (do not re-learn)

- **Empty-object sentinels on the round-trip path**: glm flash echoes
  omitted optionals as `{}` (empty JSON object), not just the literal
  `"null"`/`"{}"` strings already in `AGENT_SENTINEL_STRINGS`. Observed
  live as `facet_ncol = {}` failing twice with "must be an integer
  between 1 and 6". Fix pattern: every param/optional reader must treat
  length-0 values (empty list, `integer(0)`) exactly like omitted —
  `.agent_figure_{numeric,integer}_param` and `params$se` now do; check
  this class first when adding any new optional tool argument. Live
  before/after: 3 spec failures → 0.
- **Tier B settle heuristics**: the shinychat cancel-control probe misses
  silent glm-flash streams; activity must ALSO be read from the
  diagnostic log (newest file ends with `provider_request_start` ⇒ still
  working), and after a figure renders, wait for 18 s of quiet so a
  pending revision (second figure) is captured. Both are in tier_b.mjs
  now.
- (Still true from before: MockShinySession output poisoning, mock input
  maps not following update messages, observer GC — see AGENTS.md
  "Environment quirks" and the S4 notes below.)

## Validation snapshot (all green at c5cdbfc)

- Unit board: widgetStore 52, agentWidgets 105, agentAssistant 64,
  agentFigures 37, aiAssistantTools 25, agentLogSummary 32,
  agentLogging 21, appState 37, quickViews 17, scatterSelection 6,
  triselectorCascade 6, tableWidgetState 5, shinyAuxi 10
- Tier A browser: 93/93 (×2 of the last 3 runs; one single-assertion
  timing flake in between — the known sporadic class)
- Live Tier B: smoke subset only (orientation + create/revise); full ×3
  benchmark re-run is a phase-gate activity (before WP7's gate decision)

## Next session: WP6 (figure templates)

Plan §4-WP6: `create_figure(template = "...", ...)` with template enum
(`volcano`, `scatter`, `boxplot`, `barplot`, `histogram`, `density`,
`line`) + a few well-named args (`x`, `y`, `color`, `label_top_n`);
server-side `agent_figure_template_spec()` in `auxi_agentFigures.R`
expands to a full validated spec **and returns it via the WP3 round-trip**
so follow-up customization flows through `update_figure`. Generic grammar
stays the advanced path. `pca` template optional (needs server-side
projection — decide later). Tests: unit expansion cases + tool-level in
`test_aiAssistantTools.R`; one Tier B smoke ("make a volcano plot" via
template).

After WP6: WP7 patch-mode `update_figure` — **GATED** (decision 5): build
only if tasks 10–11 first-attempt success < 2/3 at the full ×3 phase-gate
run; the WP3 registry (`figures()[[id]]$spec`) is the merge base. Then
WP8 (enrichment/table-view tools with shared help metadata).

## Session workflow reminders

- Reinstall after every source change
  (`Rscript -e "install.packages('<repo>', repos=NULL, type='source')"`)
  — testServer/browser harnesses run the INSTALLED package; roxygen
  (`roxygen2::roxygenise('.')`) after touching roxygen blocks.
- Kill orphaned R processes on ports 7775–7795 before browser runs.
- Tier A convention ×2 consecutive runs; Tier B needs
  `tests/e2e_agent/provider.env` (present; glm-5.3-flash via
  open.bigmodel.cn) and consumes provider quota — smoke subsets only per
  change, full ×3 at gates.
- Log timestamps are UTC (user is UTC+2).
