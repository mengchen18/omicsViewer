// =============================================================================
// check_regression.mjs — diff a fresh widget sweep against the committed
// baseline (baseline_summary.json). Pure logic, no browser: run_regression.sh
// feeds it result directories.
// =============================================================================
// FAILURE RULES (FIX_PLAN.md WP-0: "failing on any regression"):
//   1. count increase in a failure status:
//        binding:  timeout / no-ui-effect / no-probe / rejected / unchanged / error
//        dom-only: unaddressable / no-effect / unchanged / error
//        semantic: no-ui-effect / error
//   2. per-case success -> non-success transition (any category).
// Never a failure (informational): improvements, cases that disappear
// (widget removed -> refresh baseline deliberately), new cases,
// dom-only unexpected-success (a finding, reported loudly).
//
// Multiple --results dirs: later directories OVERRIDE earlier ones on a
// per-case basis (caseId keyed). This is how retries work: the gate re-runs
// only the regressed caseIds via WIDGET_SWEEP_ONLY and overlays the verdicts.
//
// Usage:
//   node check_regression.mjs [--baseline <file>] --results <dir> [--results <dir> ...]
//                             [--write-regressions <file>] [--selftest]
// Exit code: 0 = no regression, 1 = regression (or selftest failure), 2 = usage.
// =============================================================================

import fs from 'node:fs';
import path from 'node:path';
import os from 'node:os';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));

const args = process.argv.slice(2);
const resultsDirs = [];
for (let i = 0; i < args.length; i++) {
  if (args[i] === '--results') resultsDirs.push(args[++i]);
  else if (args[i] === '--baseline') i; // handled by flag() below
}
const flag = (name, dflt) => {
  const i = args.indexOf(name);
  return i >= 0 && i + 1 < args.length ? args[i + 1] : dflt;
};
const baselineFile = flag('--baseline', path.join(here, 'baseline_summary.json'));
const regressionsFile = flag('--write-regressions', null);

// statuses whose count must never exceed the baseline, per category
const COUNT_GUARD = {
  binding: ['timeout', 'no-ui-effect', 'no-probe', 'rejected', 'unchanged', 'error'],
  // WP-G2: excluded-by-policy is the settled-status home of the snapshot
  // lifecycle + dataset selector cases; count increases beyond the
  // baseline fail (nothing may join the exclusion list silently)
  'dom-only': ['unaddressable', 'no-effect', 'unchanged', 'error', 'excluded-by-policy'],
  'semantic-tool': ['no-ui-effect', 'error'],
  // WP-C action cases: 'rejected' is a legitimate outcome (availability
  // gating with a self-describing reason, e.g. 'set a search term first');
  // only harness/infra failures guard the count
  action: ['error', 'timeout', 'unchanged']
};

// Runtime-generated ids (Shiny auto ids like btn-394616653 on id-less
// buttons) churn on every app boot: the same structural widget gets a new
// id each run, so per-case matching and count guards would be nondeterministic.
// Such cases are excluded from both rules and reported as a separate bucket;
// only a NON-unaddressable verdict on a volatile id fails the gate (it would
// mean the store accepted a garbage id).
const VOLATILE_ID = /^(?:dom:)?(?:button:)?btn-\d+$/;
const isVolatileId = id => VOLATILE_ID.test(id);

const loadBaseline = file => JSON.parse(fs.readFileSync(file, 'utf8'));

// read one results dir into { caseId: {status, case, category} }
function loadResults(dir) {
  const out = {};
  const files = fs.readdirSync(dir).filter(f => f.endsWith('.json')).sort();
  if (!files.length) throw new Error(`no result files in ${dir}`);
  for (const f of files) {
    const r = JSON.parse(fs.readFileSync(path.join(dir, f), 'utf8'));
    if (!r.caseId || !r.status) throw new Error(`${f}: missing caseId/status`);
    out[r.caseId] = { status: r.status, case: r.case || f, category: r.category || 'unknown' };
  }
  return out;
}

function compare(baseline, freshDirs) {
  // overlay: later dirs win
  const freshAll = {};
  for (const dir of freshDirs) Object.assign(freshAll, loadResults(dir));
  const split = map => {
    const stable = {}, volatile = {};
    for (const [id, v] of Object.entries(map)) (isVolatileId(id) ? volatile : stable)[id] = v;
    return { stable, volatile };
  };
  const base = split(baseline.cases);
  const fresh = split(freshAll);

  const counts = src => {
    const c = {};
    for (const v of Object.values(src)) {
      c[v.category] = c[v.category] || {};
      c[v.category][v.status] = (c[v.category][v.status] || 0) + 1;
    }
    return c;
  };
  const bc = counts(base.stable);
  const fc = counts(fresh.stable);

  const regressions = [];   // {kind: 'count'|'case', ...}
  const improvements = [];
  const warnings = [];

  // rule 1: guarded-status count increases
  const cats = new Set([...Object.keys(bc), ...Object.keys(fc)]);
  for (const cat of cats) {
    const guard = COUNT_GUARD[cat] || [];
    const statuses = new Set([...Object.keys(bc[cat] || {}), ...Object.keys(fc[cat] || {})]);
    for (const st of statuses) {
      const b = (bc[cat] || {})[st] || 0;
      const f = (fc[cat] || {})[st] || 0;
      if (f > b && guard.includes(st))
        regressions.push({ kind: 'count', category: cat, status: st, before: b, after: f });
      else if (f !== b && (st === 'success' || st === 'out-of-scope'))
        improvements.push(`${cat}/${st}: ${b} -> ${f}`);
    }
    if ((fc[cat] || {}).hasOwnProperty('unexpected-success'))
      warnings.push(`${cat}: unexpected-success x${fc[cat]['unexpected-success']} (agent could set a widget believed unaddressable — investigate)`);
  }

  // rule 2: per-case success -> non-success (stable ids only)
  const baseIds = Object.keys(base.stable);
  for (const id of baseIds) {
    const b = base.stable[id];
    const f = fresh.stable[id];
    if (!f) {
      if (b.status === 'success') warnings.push(`case vanished (widget removed?): ${id} [${b.case}]`);
      else improvements.push(`case resolved/removed: ${id} ${b.status} -> (absent)`);
      continue;
    }
    if (b.status === 'success' && f.status !== 'success')
      regressions.push({ kind: 'case', caseId: id, before: b.status, after: f.status, file: f.case });
    else if (b.status !== 'success' && f.status === 'success')
      improvements.push(`${id}: ${b.status} -> success`);
    else if (b.status !== f.status)
      improvements.push(`${id}: ${b.status} -> ${f.status}`);
  }
  for (const id of Object.keys(fresh.stable)) {
    if (!base.stable[id]) warnings.push(`new case not in baseline: ${id} [${fresh.stable[id].status}] (${fresh.stable[id].case})`);
  }

  // volatile auto-id bucket: never guarded by counts; churn is expected
  const bv = Object.values(base.volatile).length;
  const fv = Object.values(fresh.volatile).length;
  warnings.push(`volatile auto-id cases (btn-<random>, excluded from rules): baseline ${bv}, fresh ${fv}${bv === fv ? '' : ' — CHURN'}`);
  for (const [id, v] of Object.entries(fresh.volatile)) {
    if (v.status !== 'unaddressable')
      regressions.push({ kind: 'case', caseId: id, before: 'unaddressable (expected)', after: v.status, file: v.case });
  }

  return { regressions, improvements, warnings, baseCounts: bc, freshCounts: fc,
           nBase: baseIds.length, nFresh: Object.keys(fresh.stable).length };
}

function report(res, freshDirs) {
  const line = (bc, fc, cat, st) =>
    `  ${cat.padEnd(14)} ${st.padEnd(20)} ${String((bc[cat] || {})[st] || 0).padStart(4)}  -> ${String((fc[cat] || {})[st] || 0).padStart(4)}`;
  console.log('== widget sweep regression check ==');
  console.log(`baseline: ${res.nBase} cases | fresh: ${res.nFresh} cases from ${freshDirs.length} results dir(s)`);
  console.log('status matrix (baseline -> fresh):');
  const cats = [...new Set([...Object.keys(res.baseCounts), ...Object.keys(res.freshCounts)])].sort();
  for (const cat of cats) {
    const sts = [...new Set([...Object.keys(res.baseCounts[cat] || {}), ...Object.keys(res.freshCounts[cat] || {})])].sort();
    for (const st of sts) console.log(line(res.baseCounts, res.freshCounts, cat, st));
  }
  if (res.improvements.length) {
    console.log(`\nimprovements (${res.improvements.length}):`);
    for (const i of res.improvements) console.log(`  + ${i}`);
  }
  if (res.warnings.length) {
    console.log(`\nwarnings (${res.warnings.length}):`);
    for (const w of res.warnings) console.log(`  ! ${w}`);
  }
  if (res.regressions.length) {
    console.log(`\nREGRESSIONS (${res.regressions.length}):`);
    for (const r of res.regressions) {
      if (r.kind === 'count') console.log(`  x [count] ${r.category} ${r.status}: ${r.before} -> ${r.after}`);
      else console.log(`  x [case] ${r.caseId}: ${r.before} -> ${r.after}`);
    }
  } else {
    console.log('\nNO REGRESSION — gate PASS');
  }
  return res.regressions.length === 0;
}

// ---------------------------------------------------------------------------
// selftest: synthetic fixtures exercising every rule (no browser needed)
// ---------------------------------------------------------------------------
function selftest() {
  const mk = (id, category, status, file) => ({ [id]: { status, case: file || id, category } });
  const base = {
    cases: {
      ...mk('b.ok1', 'binding', 'success', '001-b_ok1'),
      ...mk('b.ok2', 'binding', 'success', '002-b_ok2'),
      ...mk('b.bad1', 'binding', 'no-ui-effect', '003-b_bad1'),
      ...mk('b.tout', 'binding', 'timeout', '004-b_tout'),
      ...mk('b.noprobe', 'binding', 'no-probe', '005-b_noprobe'),
      ...mk('s.ok1', 'semantic-tool', 'success', '006-s_ok1'),
      ...mk('d.u1', 'dom-only', 'unaddressable', '007-d_u1'),
      ...mk('d.u2', 'dom-only', 'unaddressable', '008-d_u2'),
      ...mk('o.g1', 'out-of-scope', 'out-of-scope', '009-o_g1')
    }
  };
  const write = (dir, entries) => {
    fs.mkdirSync(dir, { recursive: true });
    for (const [id, v] of Object.entries(entries))
      fs.writeFileSync(path.join(dir, `${v.case}.json`),
        JSON.stringify({ caseId: id, ...v }));
  };

  const tmp = fs.mkdtempSync(path.join(os.tmpdir(), 'sw-check-selftest-'));
  const dir = n => path.join(tmp, n);

  const cases = [];
  // NOTE: compare() loads results dirs from disk; write fixtures, then compare.
  const run = (name, baselineObj, dirs, expectPass) => {
    const res = compare(baselineObj, dirs);
    const ok = (res.regressions.length === 0) === expectPass;
    cases.push([name, ok]);
    console.log(`  ${ok ? 'ok  ' : 'FAIL'} ${name}${ok ? '' : ` (regressions: ${JSON.stringify(res.regressions)})`}`);
  };

  // 1. identical fresh run -> pass
  write(dir('r1'), {
    'b.ok1': { status: 'success', case: '001-b_ok1', category: 'binding' },
    'b.ok2': { status: 'success', case: '002-b_ok2', category: 'binding' },
    'b.bad1': { status: 'no-ui-effect', case: '003-b_bad1', category: 'binding' },
    'b.tout': { status: 'timeout', case: '004-b_tout', category: 'binding' },
    'b.noprobe': { status: 'no-probe', case: '005-b_noprobe', category: 'binding' },
    's.ok1': { status: 'success', case: '006-s_ok1', category: 'semantic-tool' },
    'd.u1': { status: 'unaddressable', case: '007-d_u1', category: 'dom-only' },
    'd.u2': { status: 'unaddressable', case: '008-d_u2', category: 'dom-only' },
    'o.g1': { status: 'out-of-scope', case: '009-o_g1', category: 'out-of-scope' }
  });
  run('identical matrix passes', base, [dir('r1')], true);

  // 2. success -> no-ui-effect -> fail
  write(dir('r2'), {
    'b.ok1': { status: 'no-ui-effect', case: '001-b_ok1', category: 'binding' },
    'b.ok2': { status: 'success', case: '002-b_ok2', category: 'binding' },
    'b.bad1': { status: 'no-ui-effect', case: '003-b_bad1', category: 'binding' },
    'b.tout': { status: 'timeout', case: '004-b_tout', category: 'binding' },
    'b.noprobe': { status: 'no-probe', case: '005-b_noprobe', category: 'binding' },
    's.ok1': { status: 'success', case: '006-s_ok1', category: 'semantic-tool' },
    'd.u1': { status: 'unaddressable', case: '007-d_u1', category: 'dom-only' },
    'd.u2': { status: 'unaddressable', case: '008-d_u2', category: 'dom-only' },
    'o.g1': { status: 'out-of-scope', case: '009-o_g1', category: 'out-of-scope' }
  });
  run('success -> no-ui-effect fails', base, [dir('r2')], false);

  // 3. improvement only (bad1 fixed) -> pass
  write(dir('r3'), {
    'b.ok1': { status: 'success', case: '001-b_ok1', category: 'binding' },
    'b.ok2': { status: 'success', case: '002-b_ok2', category: 'binding' },
    'b.bad1': { status: 'success', case: '003-b_bad1', category: 'binding' },
    'b.tout': { status: 'success', case: '004-b_tout', category: 'binding' },
    'b.noprobe': { status: 'no-probe', case: '005-b_noprobe', category: 'binding' },
    's.ok1': { status: 'success', case: '006-s_ok1', category: 'semantic-tool' },
    'd.u1': { status: 'unaddressable', case: '007-d_u1', category: 'dom-only' },
    'd.u2': { status: 'unaddressable', case: '008-d_u2', category: 'dom-only' },
    'o.g1': { status: 'out-of-scope', case: '009-o_g1', category: 'out-of-scope' }
  });
  run('pure improvements pass', base, [dir('r3')], true);

  // 4. NEW case with a guarded bad status (count increase without transition) -> fail
  write(dir('r4'), {
    'b.ok1': { status: 'success', case: '001-b_ok1', category: 'binding' },
    'b.ok2': { status: 'success', case: '002-b_ok2', category: 'binding' },
    'b.bad1': { status: 'no-ui-effect', case: '003-b_bad1', category: 'binding' },
    'b.tout': { status: 'timeout', case: '004-b_tout', category: 'binding' },
    'b.noprobe': { status: 'no-probe', case: '005-b_noprobe', category: 'binding' },
    'b.newbad': { status: 'no-ui-effect', case: '010-b_newbad', category: 'binding' },
    's.ok1': { status: 'success', case: '006-s_ok1', category: 'semantic-tool' },
    'd.u1': { status: 'unaddressable', case: '007-d_u1', category: 'dom-only' },
    'd.u2': { status: 'unaddressable', case: '008-d_u2', category: 'dom-only' },
    'o.g1': { status: 'out-of-scope', case: '009-o_g1', category: 'out-of-scope' }
  });
  run('new no-ui-effect case fails (count guard)', base, [dir('r4')], false);

  // 5. retry overlay: r2 regression fixed by overlay dir -> pass
  write(dir('r5fix'), {
    'b.ok1': { status: 'success', case: '101-b_ok1', category: 'binding' }
  });
  run('retry overlay clears regression', base, [dir('r2'), dir('r5fix')], true);

  // 6. baseline case missing entirely -> pass with warning (widget removed)
  write(dir('r6'), {
    'b.ok2': { status: 'success', case: '002-b_ok2', category: 'binding' },
    'b.bad1': { status: 'no-ui-effect', case: '003-b_bad1', category: 'binding' },
    'b.tout': { status: 'timeout', case: '004-b_tout', category: 'binding' },
    'b.noprobe': { status: 'no-probe', case: '005-b_noprobe', category: 'binding' },
    's.ok1': { status: 'success', case: '006-s_ok1', category: 'semantic-tool' },
    'd.u1': { status: 'unaddressable', case: '007-d_u1', category: 'dom-only' },
    'd.u2': { status: 'unaddressable', case: '008-d_u2', category: 'dom-only' },
    'o.g1': { status: 'out-of-scope', case: '009-o_g1', category: 'out-of-scope' }
  });
  run('removed success case warns but passes', base, [dir('r6')], true);

  // 7. dom-only unaddressable count increase -> fail
  write(dir('r7'), {
    'b.ok1': { status: 'success', case: '001-b_ok1', category: 'binding' },
    'b.ok2': { status: 'success', case: '002-b_ok2', category: 'binding' },
    'b.bad1': { status: 'no-ui-effect', case: '003-b_bad1', category: 'binding' },
    'b.tout': { status: 'timeout', case: '004-b_tout', category: 'binding' },
    'b.noprobe': { status: 'no-probe', case: '005-b_noprobe', category: 'binding' },
    's.ok1': { status: 'success', case: '006-s_ok1', category: 'semantic-tool' },
    'd.u1': { status: 'unaddressable', case: '007-d_u1', category: 'dom-only' },
    'd.u2': { status: 'unaddressable', case: '008-d_u2', category: 'dom-only' },
    'd.u3': { status: 'unaddressable', case: '011-d_u3', category: 'dom-only' },
    'o.g1': { status: 'out-of-scope', case: '009-o_g1', category: 'out-of-scope' }
  });
  run('dom-only unaddressable increase fails', base, [dir('r7')], false);

  // 8. fresh error status -> fail even where baseline had same-count worse statuses
  write(dir('r8'), {
    'b.ok1': { status: 'error', case: '001-b_ok1', category: 'binding' },
    'b.ok2': { status: 'success', case: '002-b_ok2', category: 'binding' },
    'b.bad1': { status: 'no-ui-effect', case: '003-b_bad1', category: 'binding' },
    'b.tout': { status: 'timeout', case: '004-b_tout', category: 'binding' },
    'b.noprobe': { status: 'no-probe', case: '005-b_noprobe', category: 'binding' },
    's.ok1': { status: 'success', case: '006-s_ok1', category: 'semantic-tool' },
    'd.u1': { status: 'unaddressable', case: '007-d_u1', category: 'dom-only' },
    'd.u2': { status: 'unaddressable', case: '008-d_u2', category: 'dom-only' },
    'o.g1': { status: 'out-of-scope', case: '009-o_g1', category: 'out-of-scope' }
  });
  run('error status fails', base, [dir('r8')], false);

  // 9. no-probe count increase via transition of a success case -> fail
  write(dir('r9'), {
    'b.ok1': { status: 'success', case: '001-b_ok1', category: 'binding' },
    'b.ok2': { status: 'no-probe', case: '002-b_ok2', category: 'binding' },
    'b.bad1': { status: 'no-ui-effect', case: '003-b_bad1', category: 'binding' },
    'b.tout': { status: 'timeout', case: '004-b_tout', category: 'binding' },
    'b.noprobe': { status: 'no-probe', case: '005-b_noprobe', category: 'binding' },
    's.ok1': { status: 'success', case: '006-s_ok1', category: 'semantic-tool' },
    'd.u1': { status: 'unaddressable', case: '007-d_u1', category: 'dom-only' },
    'd.u2': { status: 'unaddressable', case: '008-d_u2', category: 'dom-only' },
    'o.g1': { status: 'out-of-scope', case: '009-o_g1', category: 'out-of-scope' }
  });
  run('success -> no-probe fails', base, [dir('r9')], false);

  // 10. volatile auto-id churn (btn-<random>): churn never fails; only a
  //     non-unaddressable volatile verdict does
  write(dir('r10'), {
    'b.ok1': { status: 'success', case: '001-b_ok1', category: 'binding' },
    'b.ok2': { status: 'success', case: '002-b_ok2', category: 'binding' },
    'b.bad1': { status: 'no-ui-effect', case: '003-b_bad1', category: 'binding' },
    'b.tout': { status: 'timeout', case: '004-b_tout', category: 'binding' },
    'b.noprobe': { status: 'no-probe', case: '005-b_noprobe', category: 'binding' },
    's.ok1': { status: 'success', case: '006-s_ok1', category: 'semantic-tool' },
    'd.u1': { status: 'unaddressable', case: '007-d_u1', category: 'dom-only' },
    'd.u2': { status: 'unaddressable', case: '008-d_u2', category: 'dom-only' },
    'o.g1': { status: 'out-of-scope', case: '009-o_g1', category: 'out-of-scope' },
    'dom:btn-999111': { status: 'unaddressable', case: '020-dom_btn_999111', category: 'dom-only' },
    'dom:btn-999222': { status: 'unaddressable', case: '021-dom_btn_999222', category: 'dom-only' },
    'dom:btn-999333': { status: 'unaddressable', case: '022-dom_btn_999333', category: 'dom-only' }
  });
  const baseVol = { cases: { ...base.cases,
    ...mk('dom:btn-111222', 'dom-only', 'unaddressable', '030-dom_btn_111222'),
    ...mk('dom:btn-333444', 'dom-only', 'unaddressable', '031-dom_btn_333444') } };
  run('volatile id churn passes (count 2 -> 3)', baseVol, [dir('r10')], true);

  write(dir('r11'), {
    'b.ok1': { status: 'success', case: '001-b_ok1', category: 'binding' },
    'b.ok2': { status: 'success', case: '002-b_ok2', category: 'binding' },
    'b.bad1': { status: 'no-ui-effect', case: '003-b_bad1', category: 'binding' },
    'b.tout': { status: 'timeout', case: '004-b_tout', category: 'binding' },
    'b.noprobe': { status: 'no-probe', case: '005-b_noprobe', category: 'binding' },
    's.ok1': { status: 'success', case: '006-s_ok1', category: 'semantic-tool' },
    'd.u1': { status: 'unaddressable', case: '007-d_u1', category: 'dom-only' },
    'd.u2': { status: 'unaddressable', case: '008-d_u2', category: 'dom-only' },
    'o.g1': { status: 'out-of-scope', case: '009-o_g1', category: 'out-of-scope' },
    'dom:btn-555666': { status: 'unexpected-success', case: '023-dom_btn_555666', category: 'dom-only' }
  });
  run('volatile id with unexpected-success fails', baseVol, [dir('r11')], false);

  fs.rmSync(tmp, { recursive: true, force: true });
  const failed = cases.filter(([, ok]) => !ok);
  console.log(`selftest: ${cases.length - failed.length}/${cases.length} ok`);
  return failed.length === 0;
}

// ---------------------------------------------------------------------------
if (args.includes('--selftest')) {
  process.exit(selftest() ? 0 : 1);
}
if (!resultsDirs.length) {
  console.error('usage: node check_regression.mjs --results <dir> [--results <dir> ...] [--baseline <file>] [--write-regressions <file>] | --selftest');
  process.exit(2);
}

let ok;
try {
  const baseline = loadBaseline(baselineFile);
  const res = compare(baseline, resultsDirs.map(d => path.resolve(d)));
  ok = report(res, resultsDirs);
  if (regressionsFile) {
    const ids = res.regressions.map(r => r.caseId).filter(Boolean);
    fs.writeFileSync(regressionsFile, JSON.stringify({ caseIds: ids }, null, 2) + '\n');
  }
} catch (e) {
  console.error(`check_regression: ${e.message}`);
  process.exit(2);
}
process.exit(ok ? 0 : 1);
