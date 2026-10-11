// =============================================================================
// refresh_baseline.mjs — (re)build baseline_summary.json from a sweep results
// directory.
// =============================================================================
// The regression gate (run_regression.sh) diffs every fresh sweep against the
// committed baseline_summary.json. This script builds that baseline from the
// per-case result files of a run you have REVIEWED (the initial baseline comes
// from the curated 2026-10-11 sweep: statuses reconciled with the live
// spot-checks recorded in spotcheck_results.json / spotcheck2_results.json).
//
// Workflow per FIX_PLAN.md: after a work package lands, re-run the sweep,
// review the diff (check_regression.mjs shows every transition), then — and
// ONLY for the WP's intended improvements — point this script at the fresh
// results dir and commit the new baseline together with the code change.
//
// Usage:
//   node refresh_baseline.mjs [--results <dir>] [--out <file>] [--note "..."]
//     --results  default: ./results  (the curated evidence in this directory)
//     --out      default: ./baseline_summary.json
// =============================================================================

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));

const args = process.argv.slice(2);
const flag = (name, dflt) => {
  const i = args.indexOf(name);
  return i >= 0 && i + 1 < args.length ? args[i + 1] : dflt;
};
const resultsDir = path.resolve(flag('--results', path.join(here, 'results')));
const outFile = path.resolve(flag('--out', path.join(here, 'baseline_summary.json')));
const note = flag('--note', '');

const files = fs.readdirSync(resultsDir).filter(f => f.endsWith('.json')).sort();
if (!files.length) {
  console.error(`no result files in ${resultsDir}`);
  process.exit(1);
}

const cases = {};
const byCat = {};
let bad = 0;
for (const f of files) {
  let r;
  try { r = JSON.parse(fs.readFileSync(path.join(resultsDir, f), 'utf8')); }
  catch (e) { console.error(`unparseable: ${f}: ${e.message}`); bad++; continue; }
  if (!r.caseId || !r.status || !r.category) {
    console.error(`missing caseId/status/category: ${f}`); bad++; continue;
  }
  if (cases[r.caseId]) { console.error(`duplicate caseId: ${r.caseId} (${f})`); bad++; continue; }
  cases[r.caseId] = { status: r.status, case: r.case, category: r.category };
  byCat[r.category] = byCat[r.category] || {};
  byCat[r.category][r.status] = (byCat[r.category][r.status] || 0) + 1;
}
if (bad) process.exit(1);

const baseline = {
  generated: new Date().toISOString(),
  note: note || 'baseline; see README.md (Regression gate) and FIX_PLAN.md WP-0',
  source_results: resultsDir,
  totals: { cases: files.length },
  by_category_status: byCat,
  cases
};
fs.writeFileSync(outFile, JSON.stringify(baseline, null, 2) + '\n');

console.log(`baseline written: ${outFile}`);
console.log(`  ${files.length} cases from ${resultsDir}`);
for (const [cat, st] of Object.entries(byCat))
  for (const [s, n] of Object.entries(st)) console.log(`  ${cat.padEnd(14)} ${s.padEnd(20)} ${n}`);
