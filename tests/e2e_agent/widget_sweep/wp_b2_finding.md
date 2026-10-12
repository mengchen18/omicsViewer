# WP-B2 finding — expr-heatmap congestion (F6): harness artifact, not a bug

**Outcome: NO product code change (the plan's expected "pass" branch).**
Cases 022 (`dataspace.expr_heatmap.row_sort_by`) and 026
(`dataspace.expr_heatmap.annot_row`) are green in every isolation run;
the original no-ui-effect grades are a measurement artifact of the
sweep harness, now understood precisely.

## Evidence (all from live runs, 2026-10-12)

`spotcheck3.mjs` — built for this WP (the plan named it; it did not
exist). Per run: fresh app boot (OMICSVIEWER_TEST_HOOKS), write through
the exact tool path, poll the SPECIFIC input key (hook controls
excluded), measure main-thread stalls via a 50 ms heartbeat, verify the
store read-back (post-WP-B fixed path).

| run | mode | 022 applied | 022 time-to-effect | 022 store | stalls | 026 applied | 026 effect | 026 store |
|---|---|---|---|---|---|---|---|---|
| 1 | tab inactive | yes | ≤ roundtrip* | "none" | 1×1.19 s | yes | 65 ms | correct |
| 2 | tab inactive | yes | ≤ roundtrip* | "none" | 1×1.19 s | yes | 52 ms | correct |
| 3 | tab inactive | yes | ≤ roundtrip* | "none" | 1×1.18 s | yes | 65 ms | correct |
| 4 | Heatmap tab active | unchanged** | ≤ roundtrip* | "none" | none | yes | 60 ms | correct |
| 5 | tab inactive, hardened probe*** | yes | ≤ roundtrip* | ok | 1×1.16 s | yes | 61 ms | correct |
| 6 | tab inactive, hardened probe*** | yes | ≤ roundtrip* | ok | 1×1.20 s | yes | 64 ms | correct |

\*  the push lands DURING the hook roundtrip (<~150 ms) — `input_before`
   was `""`/`[]` and the written value is present at the first post-hook
   check, so the effect is faster than any poll interval can resolve.
\** run 4 wrote 'none' which the store already held (converged state).
\*** runs 5-6 wrote `General|All|Gene.name` (≠ any default) for 022 to
   exclude the "already there" ambiguity; input transitioned
   `""` → `General|All|Gene.name`.

In-harness isolation run (`run_regression.sh --only dataspace.expr_heatmap`,
gate PASS): all 13 expr_heatmap cases success — row_sort_by 3033 ms,
annot_row 335 ms (case wall times, well inside the 30 s budget).

Unit suites (unchanged code, counts recorded): test_renderStability 46/46,
test_snapshotRoundTrip 63/63.

## Root cause — sharper than "congestion"

The F6 note suspected the 2702×60 re-render "blocks the browser during
the poll window". Measured: the re-render blocks the main thread for
**~1.2 s** (022's heartbeat stall), which is nowhere near the 30 s
budget. The actual mechanism is in the harness's `waitUiEffect`:

- it polls `Shiny.shinyapp.$inputValues` and **returns on the FIRST
  input diff** — but the hook's own op/payload/run controls always
  differ within ~20 ms of the click, so the effective observation
  window was ~20 ms, not 30 s (committed evidence: `waited_ms: 21-22`
  on both curated cases);
- whether the real input change was seen depended on whether the
  store→UI push (plus, for 022, the ~1.2 s re-render sharing the same
  flush) completed within that ~20-250 ms first-poll window — hence the
  flapping between the curated run (missed) and the WP-A full run
  (caught).

This is the same class as the WP-B finding: the raw harness grades
`applied && any-input-diff` where the diff is hook churn, and the
curation pass corrected case-by-case. The curated no-ui-effect statuses
for 022/026 therefore remain as documented residue (measurement
artifacts on a working write path), consistent with the WP-B treatment
of `tab_*.columns`.

## Recorded consequences

- No expr-heatmap provider/push defect exists: providers are warm
  (dataset-level reactives, the WP-A pattern already), pushes land
  within the hook roundtrip, read-backs are exact, the pending/ack echo
  never clobbers (store equals written set in all runs).
- Future WPs: treat `no-ui-effect` grades from the input-diff probe as
  WEAK evidence (carried from WP-B); re-verify with a targeted probe
  (spotcheck3 pattern) before planning any code change.
- A systematic harness fix (wait for a non-hook diff, or a DOM-level
  effect check) would re-grade large parts of the matrix and require a
  deliberate re-baseline — out of WP-B2 scope, noted for the gate's
  backlog.

## Artifacts

- `spotcheck3.mjs` (this directory) — runnable: `node spotcheck3.mjs`
  (options: `--tabbed` activates the Heatmap tab first; `ROW_WRITE=...`
  overrides the 022 probe value)
- `spotcheck3_results.json` — run 6 evidence (hardened probe)
- baseline_summary.json note field references this finding
