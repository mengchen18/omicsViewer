# HANDOVER — LLM-assistant UI-mimic evaluation & hardening program

Written 2026-10-17 at the end of the evaluation session. Read together
with `tests/e2e_agent/ui_mimic_findings.md` (evidence) and
`tests/e2e_agent/ui_mimic_FIX_PLAN.md` (the remediation plan to execute).

## Where things stand

- Branch `master`, tree clean. **3 local unpushed commits** (Chen pushes):
  - `2d4601f` fix: isolate the two store_apply calls in
    apply_agent_scatter_view (finding F1 — set_scatter_view was dead on
    every fresh session; verified by minimal repro before/after)
  - `a7f09e5` eval: UI-mimic harness + findings (ui_mimic.mjs, turn
    files, 4 probe scripts, ui_mimic_findings.md)
  - `f4d91b1` plan: remediation plan for findings F2–F8 (9 WPs, 5 phases)
- Before those, `extend-llm-tools` was merged to master and pushed
  (widget-sweep program WP-0..WP-G2, version 2.1.79 = `6cc0952`).

## The next task

Execute `tests/e2e_agent/ui_mimic_FIX_PLAN.md` — start with Phase 1
(WP-1 annotation catalog, WP-2 find_controls empty-query digest, WP-3
figure-filter error vocabulary; independent, any order), then Phase 2
(WP-4 repro probe for the F2 silent message drop → WP-5 fix → WP-6
composer watchdog), Phase 3 (WP-7 panel_results), Phase 4 (WP-8 status
line + soft tool budget), Phase 5 (WP-9 close-out). Plan rule: branch
`llm-ui-hardening` off master, ONE commit per WP, no push/merge without
Chen, guards via plain `Rscript tests/test_X.R` (never test_file), live
mimic validation turn per WP.

Finding status: F1 FIXED (2d4601f). F2–F8 open, each with evidence,
design, guards, acceptance in the plan.

## How to run the evaluation harness

```bash
cd tests/e2e_agent
node ui_mimic.mjs <turns.json> <outName>     # real chat UI, real provider
# turns files: ui_mimic_turns1..3.json; write new ones per validation
```
- Needs `provider.env` (gitignored, present on this machine:
  openai/glm-5.3-flash at bigmodel). Costs real tokens.
- Logs land in `~/.omicsviewer-llm-eval/<outName>/` (durable, outside
  the repo): `session_report.json` (per-turn receipts, ui diffs, texts),
  `llm-logs/*.jsonl` (the assistant diagnostic log), `artifacts/*.png`.
- Kill stray app processes after runs: ports 7785–7789
  (`kill $(ss -tlnp | grep :778x | grep -oP 'pid=\K[0-9]+')`).
- Package must be reinstalled after R edits:
  `Rscript -e "install.packages('/media/chen/Elements/github/omicsViewer', repos = NULL, type = 'source')"`.

## Session gotchas (hard-won, all cost debug time)

1. **ESM `require()` trap**: any `require('node:fs')` inside a `.mjs`
   throws ReferenceError; if wrapped in try/catch it fails SILENTLY
   (produced two false-negative "no receipts" episodes). Import fs
   functions at the top of .mjs files.
2. **Duplicated `const f = () => { const f = ... }`** parses fine (inner
   shadowing) and makes the outer call a silent no-op — cost a full
   probe run of empty output.
3. **shinychat renders the user bubble client-side (optimistic)**: a
   DOM-text-growth check proves NOTHING about delivery. Server receipt
   = a new `user_message` event in the session jsonl. The mimic driver
   waits-for-idle (no cancel control AND log's last event is not
   `provider_request_start`) before each turn and retries a swallowed
   submission once.
4. **First-turn "52 ui changes" = app warm-up churn** (cascade seeding),
   not agent actions — the driver's turn-1 diff is meaningless; ignore
   or add a settle wait.
5. `Shiny.shinyapp.$inputValues` is a plain object (not iterable); use
   `Object.keys()`. Exclude `.clientdata` keys.
6. Test-hook bridge needs `OMICSVIEWER_TEST_HOOKS=true` in the app env
   and `display:block` on `#app-agentTestHooks-container` before
   fill/click (it is display:none by design).
7. Long tool chains false-settle a text-stability poll: always ALSO
   check the log's in-flight state (see 3) — 24 s of stable text during
   a 12–42 s provider request is normal.
8. Session death signature (F2): last jsonl events end at
   `assistant_response`/`history_stubbed`; subsequent prompts produce
   client bubbles but no `user_message`. Repro recipe in the plan WP-4.

## Evidence pointers (~/.omicsviewer-llm-eval/)

- `batch1..4/` — the four multi-turn sessions (reports + logs)
- `minirepro/` — F1 failing repro (pre-fix); `minirepro_fixed/` — post-fix
- `probe1/` — transport probe (websocket lifecycle, 6/6 trivial turns)
- `probe2..4/` — store/scatter-seam probes that isolated F1
- Key numbers: 154 provider requests, 139 tool calls, median request
  12.6 s (max 42 s), one turn chained 33 requests / 4+ min.

## Standing constraints (unchanged)

- No push / no merge to master / no RELEASE_3_22 cherry-picks without
  Chen's explicit instruction.
- Downloads/export: CANCELLED by Chen — build nothing.
- Snapshot lifecycle + dataset switching: permanently outside the agent
  surface (action tier, tools, capability records).
- Sweep gate stays green: `NOT_CRAN=true tests/e2e_agent/widget_sweep/
  run_regression.sh` (full run ~15–25 min; baseline 344 cases).
- WP-E deviation still awaiting Chen's accept/redirect (no-probe 11 =
  8 benign + 3 documented-conditional vs the plan's original "8") — see
  `widget_sweep/wp_e_finding.md`.
