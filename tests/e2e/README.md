# omicsViewer — end-to-end UI tests (Playwright)

User-mimicking, end-to-end tests for the omicsViewer Shiny application,
driven through the browser with Playwright — the same interaction surface a
real user has (no test hooks, no privileged access). Every check clicks,
types, drags and reads like a person would, then asserts the app reacted
correctly.

**Dual-backend + parity.** The identical journey runs against BOTH dataset
formats and the user-visible results are compared against each other:

- `rds` — `inst/extdata/demo.RDS` copied into a temp data dir
- `sqlite` — the SAME `demo.RDS`, converted with the app's own save path
  (`saveOmicsViewerDb`) into `demo.sqlite` in its own temp dir

The parity block compares dataset dimensions, selection counts, quick-view
badge ids, analysis table row counts and first-row statistics, the survival
censor range and expression cell spot-checks — with numeric tolerance
rel 1e-9 / abs 1e-12 (SQLite REAL stores exact IEEE-754 doubles; the
tolerance absorbs last-bit display drift only). Plotly internals and DOM
structure are deliberately NOT compared. Mouse-geometry results (box
selections) are per-arm checks only. The sqlite arm additionally verifies
the connection lifecycle: after the journey no file descriptor may remain
open on the database (app process still alive), an external
`saveOmicsViewerDb` overwrite must succeed, and no `-wal`/`-shm` journal
files may be left behind.

> **Dataset saving:** the UI has no dataset-object save control — only an
> xlsx table export and app-state snapshots (`.ESS`). The snapshot
> save/restore journey is the supported save mechanism and is asserted
> state-equal (pre-save vs post-restore).

## What it covers

One continuous user journey against `inst/extdata/demo.RDS` (2,702
features × 60 samples, gene sets, t-tests, PCA, survival, drug response,
seq-logo and STRING annotations):

| # | Area | Checks |
|---|------|--------|
| 1 | Boot & dataset load | page load, dataset selector, real disk-load of `demo.RDS`, all 9 data tabs + all 8 analysis tabs render, defaults |
| 2 | Feature scatter | default volcano (axes, corner rects), corner auto-selection = 107 features mirrored into the table, quick-view badge switches (RE_vs_LE = 239, cor_MDR), custom triselector axes (PCA PC1/PC2), **real box-selection drag**, clear-selection button |
| 3 | Feature table | pagination, row click selects, column sort, multi-selection toggle + ctrl-click |
| 4 | Sample space | PCA scatter with all 60 points, custom axes (MDR vs doubling time), box-selection of samples |
| 5 | Sample table | all 60 rows, row click selects |
| 6 | Heatmaps | correlation heatmap and expression heatmap render |
| 7 | Expression & GSList | expression matrix table, gene-set membership table, gene-set row click selects its member feature |
| 8 | Analysis ▸ Feature | boxplot for 107 selected features, stats table, single feature + categorical → beeswarm, ROC/PR curves for a binary variable |
| 9 | Analysis ▸ Sample | batch comparison table, beeswarm (continuous), contingency table (categorical), Kaplan–Meier survival |
| 10 | Analysis ▸ ORA | enrichment table for the 107 corner features, overlap table on row click |
| 11 | Analysis ▸ fGSEA | ranked enrichment results, ranking-line plot, leading-edge bars on pathway row click |
| 12 | Analysis ▸ Response | dose-response curve + EC50/Hill parameter table (single feature) |
| 13 | Analysis ▸ SeqLogo | PTM motif seq-logo + sequence table |
| 14 | Dynamic heatmap | renders, cell click selects exactly one feature (+sample), Clear button |
| 15 | Analysis ▸ Geneshot | UI renders (live AutoRIF API deliberately not called) |
| 16 | Analysis ▸ StringDB | UI renders, Run with empty selection → graceful error notification |
| 17 | Snapshots | modal, save (notification), restore round-trip (tab state), delete |
| 18 | Export | xlsx download (5.7 MB workbook) |
| 19 | AI assistant | panel opens with content (no LLM call) |
| 20 | Session health | zero JavaScript page errors across the whole journey |

## Running

```bash
cd tests/e2e
node e2e.spec.mjs            # both arms + parity (default)
node e2e.spec.mjs rds        # RDS arm only
node e2e.spec.mjs sqlite     # SQLite arm only
# or: npm test
```

Requirements:

- Google Chrome at `/usr/bin/google-chrome` (same as the existing
  `tests/e2e_agent` harness — no browser download).
- `playwright` is resolved through the `node_modules` symlink to
  `tests/e2e_agent/node_modules`; if that is missing, run
  `npm install` in `tests/e2e_agent` first.
- The R package must be loadable (`omicsViewer::omicsViewer`); the suite
  spawns `Rscript` itself.

Environment variables:

- `E2E_PORT` — port for the rds arm (default 7780; the sqlite arm uses
  port + 1)
- `E2E_BACKEND` — default backend mode when no CLI arg is given
  (`rds` | `sqlite` | `both`)
- `E2E_HEADLESS=0` — watch the journey in a visible browser
- `E2E_KEEP=1` — keep the temp data dirs under `artifacts/data-<backend>`
  for debugging

Artifacts:

- `e2e_results.json` — machine-readable report: per-arm checks, the parity
  comparison table and combined totals
- `artifacts/<backend>/fail-*.png` — screenshot at the moment of each
  failure
- console log — live PASS/FAIL lines with timings

## Design notes

- The app runs against a **temp copy** of `demo.RDS` so the dataset
  selector flow is exercised and snapshot files never touch the repo.
- Headless Chrome runs with WebGL disabled (the documented no-GPU scatter
  path, matching `tests/e2e_agent/tier_a.mjs`).
- Selection assertions read the DataTable info text (“… of N entries”) —
  the same surface a user sees — rather than privileged app state. Note
  the analysis stats table intentionally renders blank while no features
  are selected (`req()` on the selection), so “empty” is checked via the
  feature-table mirror instead.
- If a check fails, the suite records a screenshot and continues; a
  stuck session state is recovered (scatter clear / fresh session) so
  later checks stay meaningful.
