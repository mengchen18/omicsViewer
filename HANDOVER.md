# HANDOVER — agent accuracy & widget control plane

Written 2026-09-21 at the end of a long working session. Read this first in a
fresh context, then `AGENT_ACCURACY_PLAN.md` (canonical plan + status) and
`AGENTS.md` (environment + commands). Delete or trim this file once absorbed.

## Where things stand

Branch `agent-driven-exploration`, everything committed through `dd7029d`
(working tree clean). DESCRIPTION is 2.1.1. The assistant works end-to-end
with a real provider (figures, scatter views, selections, self-correcting
typos).

Commit map (most recent first):

| Commit | What |
|---|---|
| `dd7029d` | fix: store-less triselector cascades (analysis-panel regression) |
| `29c640a` | S2: meta_scatter store migration |
| `08dadf2` | S1: canonical widget store + transactional apply |
| `7be1777` | docs: §6 control-plane design |
| `61a0d13` | fix: mid-restore req-abort |
| `aac602f` | fix: stale-internal-axes (manual-edit drift) |
| `43afb53` | WP2: self-correcting validation errors |
| `0cdf159` | assistant + harness milestone (46 files) — includes scrubbing a vim swap file that briefly contained the live API key (history verified clean; never pushed) |

## Verified working (do not re-litigate without evidence)

- Figures through chat: volcano prompt → `get_state` → `create_figure` →
  rendered PNG + download (live-verified twice with glm-5.3-flash)
- Scatter applies in all regimes: fresh/init-racing, settled, manual-drift
  correction, rapid succession (Tier A 33/33)
- Axis changes never touch the quick/custom display mode (user requirement)
- Right panel reacts to selections (regression fixed + guarded)
- Test suites: 6+37+6+17+30+35+22+21 unit, Tier A 33/33 browser

## The one architectural rule to remember

`R/module_triselector.R` cascades use
`reactive_selector1() %||% input$analysis`. Both regimes matter:
store-backed (meta_scatter — requested state wins over in-flight inputs)
and store-less (feature_general, fgsea, geneshot, dataTable, attr4 —
selectors are NULL until a restore; inputs drive).
`tests/test_triselectorCascade.R` guards both; it was red/green validated.

## Known open threads (none blocking S3)

1. Re-saving provider settings mid-stream kills the in-flight ellmer turn
   (`AssistantTurn` validation error) — cosmetic, backlog
2. WP4 log summarizer not built yet — logs exist and are clean; classification
   helper pending
3. WP1 (`sections` param) + WP3 (figure spec round-trip) pending
4. Snapshot restore through the store (S4) pending; currently only
   meta_scatter routes restores through `store_apply`
5. Demo's truncated sample defaults (`PCA|All|PC1(`) — cleanly rejected;
   regenerating demo.RDS defaults is optional polish

## Next stage: S3 (per plan §6.3/§6.4)

Registry-generated generic agent tools. Concrete steps:

1. Expose the store to the assistant module: `app_module` already creates
   `app_store`; pass it (or a projection) into `ai_assistant_module`
2. New ellmer tools (thin, registry-driven):
   - `list_widgets(section?)` → `store_registry_view()` (user-editable only)
   - `get_widget(id)` → `store_describe()` + current value
   - `set_widgets(patch)` → `store_apply(origin="agent")`; surface
     `receipt$rejected` per key with WP2 suggestions
3. Wire snapshot save/restore of store state (`store_snapshot`/
   `store_restore`) for the registered axes/mode keys (start of S4)
4. Tier B benchmark additions: one generic-tier task on a widget never
   exposed before (e.g. heatmap colors)
5. Prompt contract line: prefer semantic tools; generic tier only for what
   they don't cover

Design decisions already settled (plan §9): single store + namespace-prefixed
children; agent visibility == user editability, exactly; no per-widget tools.

## Session workflow reminders

- Reinstall after every source change
  (`Rscript -e "install.packages('<repo>', repos=NULL, type='source')"`)
- Kill orphaned R processes on ports 7775–7786 before browser runs
- `NOT_CRAN=true` needed for shinytest2-based tests
- provider.env is gitignored and holds the live key — never commit, never
  echo; screenshots/artifacts under tests/e2e_agent/artifacts/ are gitignored
- Log timestamps are UTC (user is UTC+2)
