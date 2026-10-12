// =============================================================================
// EXHAUSTIVE WIDGET SWEEP — every user-facing interactive widget, driven
// through the LLM assistant tool surface, validated on the UI side.
// =============================================================================
// Methodology (see README.md):
//   * The app runs with OMICSVIEWER_TEST_HOOKS=true. The hidden hook module
//     exposes the EXACT server-side pathway of the assistant tools:
//       op 'widgets'      == set_widgets              (agent_widget_apply)
//       op 'scatter'      == set_scatter_view         (apply_agent_scatter_view)
//       op 'state'        == set_omics_viewer_state   (apply_agent_state)
//       op 'enrichment'   == set_enrichment_parameters
//       op 'tableview'    == set_table_view
//       op 'overview'     == get_omics_viewer_state
//       op 'capabilities' == find_controls / describe_control
//     So every operation below is what an LLM tool call would do.
//   * Enumeration is UI-authoritative: every Shiny input binding in the DOM
//     across all tabs, dropdown panels and the snapshot modal, plus DT global
//     search boxes, download controls and (listed, NOT tested) plotly /
//     network gesture surfaces (out of scope by instruction).
//   * Every widget gets a case file (cases/) and a result file (results/).
//   * A widget operation whose observable UI effect takes > 30 s is recorded
//     as failed (status 'timeout').
//
// Run:  node sweep_harness.mjs          (from tests/e2e_agent/widget_sweep/)
// =============================================================================

import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import path from 'node:path';
import fs from 'node:fs';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../../..');
const EXTDATA = path.join(REPO, 'inst/extdata');
const PORT = 7794;
// Output redirection: the regression gate (run_regression.sh) sets
// WIDGET_SWEEP_OUT so a gate run writes to a scratch directory instead of
// clobbering the committed evidence artifacts (cases/, results/,
// summary.json, ...). Default: this directory (full-evidence mode).
const OUT = process.env.WIDGET_SWEEP_OUT || here;
const CASES = path.join(OUT, 'cases');
const RESULTS = path.join(OUT, 'results');
const ART = path.join(OUT, 'artifacts');
for (const d of [CASES, RESULTS, ART]) fs.mkdirSync(d, { recursive: true });

const EFFECT_WAIT_MS = 30000;   // > 30 s => failed (user rule)
const POLL_MS = 250;
// Optional case filter: comma-separated substrings matched against caseId
// (e.g. WIDGET_SWEEP_ONLY=dataspace.expr_heatmap.row_sort_by,030-dataspace).
// Used by run_regression.sh for targeted retries and by developers for
// spot re-runs. Unset / empty = full sweep (the regression-gate default).
const ONLY = (process.env.WIDGET_SWEEP_ONLY || '').split(',').map(s => s.trim()).filter(Boolean);

// ground-truth ids from demo.RDS (for selection probes)
const FEATURE_IDS = ['X32.TYW5', 'X52.CNOT1', 'X55.MAP1LC3B.MAP1LC3B2'];
const SAMPLE_IDS = ['X786O_NCI60', 'A498_NCI60', 'A549_NCI60'];

const r = spawn('Rscript', ['-e', `
  options(shiny.port = ${PORT}, shiny.host = '127.0.0.1')
  eset <- readRDS(file.path('${EXTDATA.replace(/'/g, "\\'")}', 'demo.RDS'))
  omicsViewer::omicsViewer(dir = '${EXTDATA.replace(/'/g, "\\'")}', ESVObj = eset)
`], { cwd: REPO, env: { ...process.env, OMICSVIEWER_TEST_HOOKS: 'true' }, stdio: ['ignore', 'pipe', 'pipe'] });
process.on('exit', () => { try { r.kill('SIGKILL'); } catch {} });
process.on('SIGINT', () => { try { r.kill('SIGKILL'); } catch {}; process.exit(130); });
let rlog = '';
r.stderr.on('data', d => { rlog += d; });
r.stdout.on('data', d => { rlog += d; });

const waitPort = async (port, ms = 120000) => {
  const t0 = Date.now();
  while (Date.now() - t0 < ms) {
    try { const res = await fetch(`http://127.0.0.1:${port}/`); if (res.status > 0) return true; } catch {}
    await new Promise(s => setTimeout(s, 500));
  }
  return false;
};

const HOOK = {
  box: '#app-agentTestHooks-container',
  op: '#app-agentTestHooks-op',
  payload: '#app-agentTestHooks-payload',
  run: '#app-agentTestHooks-run',
  result: '#app-agentTestHooks-result'
};

// ---------------------------------------------------------------------------
// hook bridge == the LLM tool calls
// ---------------------------------------------------------------------------
async function runHook(page, op, payload) {
  const t0 = Date.now();
  const before = (await page.innerText(HOOK.result)).trim();
  await page.evaluate(({ sel, value }) => {
    const el = document.querySelector(sel);
    if (el.selectize) el.selectize.setValue(value);
    else { el.value = value; el.dispatchEvent(new Event('change', { bubbles: true })); }
  }, { sel: HOOK.op, value: op });
  await page.fill(HOOK.payload, JSON.stringify(payload ?? {}));
  const runCount = await page.evaluate(sel => {
    const el = document.querySelector(sel);
    const n = (parseInt(el.dataset.runs || '0', 10) + 1);
    el.dataset.runs = String(n);
    return n;
  }, HOOK.run);
  await page.click(HOOK.run);
  // 35 s primary wait, then a recovery poll: a slow reactive flush (e.g.
  // software-rendered scatter redraws) can delay the hook result far
  // beyond the 30 s budget; we still want the receipt (and sequencing
  // safety for the next case) rather than corrupting the sweep
  let done = false;
  try {
    await page.waitForFunction(({ sel, before, n }) => {
      const t = (document.querySelector(sel)?.innerText || '').trim();
      return t !== before && t.includes(`"hook_run": ${n}`);
    }, { sel: HOOK.result, before, n: runCount }, { timeout: 35000 });
    done = true;
  } catch { /* fall through to recovery */ }
  if (!done) {
    const t0r = Date.now();
    while (Date.now() - t0r < 240000) {
      const txt = (await page.innerText(HOOK.result).catch(() => '')).trim();
      if (txt.includes(`"hook_run": ${runCount}`) && txt !== before) { done = true; break; }
      await new Promise(s => setTimeout(s, 1000));
    }
  }
  if (!done) return { result: { hook_error: 'no hook response within recovery window' }, ms: Date.now() - t0 };
  let parsed = null;
  try { parsed = JSON.parse((await page.innerText(HOOK.result)).trim()); }
  catch (e) { parsed = { hook_error: 'unparseable hook result: ' + e.message }; }
  return { result: parsed, ms: Date.now() - t0 };
}

// ---------------------------------------------------------------------------
// UI-side observation
// ---------------------------------------------------------------------------
const snapshotInputs = (page) => page.evaluate(() => {
  const out = {};
  const iv = Shiny.shinyapp.$inputValues || {};
  for (const k of Object.keys(iv)) {
    if (k.startsWith('.clientdata')) continue;
    let v = iv[k];
    if (v !== null && typeof v === 'object') {
      try { v = JSON.parse(JSON.stringify(v)); } catch { v = String(v); }
    }
    out[k] = v === undefined ? null : v;
  }
  // selected table rows per visible table (DT proxy effects)
  out.__selectedRows = Array.from(document.querySelectorAll('table.dataTable'))
    .filter(t => t.offsetParent !== null)
    .map(t => ({ id: (t.closest('[id]') || { id: '?' }).id, n: t.querySelectorAll('tr.selected').length }));
  return out;
});

const diffInputs = (a, b) => {
  const keys = new Set([...Object.keys(a || {}), ...Object.keys(b || {})]);
  const changes = [];
  for (const k of keys) {
    const ja = JSON.stringify(a?.[k] ?? null), jb = JSON.stringify(b?.[k] ?? null);
    if (ja !== jb) changes.push({ input: k, before: a?.[k] ?? null, after: b?.[k] ?? null });
  }
  return changes;
};

const axisTitles = (page) => page.evaluate(() => {
  const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
    .filter(e => e.offsetParent !== null && e._fullLayout);
  if (!els.length) return null;
  const l = els[0]._fullLayout;
  const t = a => String((a && a.title && (a.title.text || a.title)) || '');
  return { x: t(l.xaxis), y: t(l.yaxis) };
});

// wait until the input map changes (any non-clientdata key) or predicate ok
async function waitUiEffect(page, before, predicate, budget = EFFECT_WAIT_MS) {
  const t0 = Date.now();
  let last = [];
  while (Date.now() - t0 < budget) {
    const now = await snapshotInputs(page);
    last = diffInputs(before, now);
    if (last.length || (predicate && await predicate())) return { changed: last, ms: Date.now() - t0, timedOut: false };
    await new Promise(s => setTimeout(s, POLL_MS));
  }
  return { changed: last, ms: Date.now() - t0, timedOut: true };
}

// ---------------------------------------------------------------------------
// DOM enumeration (UI-authoritative inventory)
// ---------------------------------------------------------------------------
const enumerateDom = (page) => page.evaluate(() => {
  const iv = Shiny.shinyapp.$inputValues || {};
  const vis = el => !!(el && el.offsetParent !== null);
  // DT per-column filter inputs (user-facing, not Shiny inputs)
  const dtColFilters = Array.from(document.querySelectorAll('.dataTable th input'))
    .map(el => ({ placeholder: el.placeholder || null, closestId: (el.closest('[id]') || { id: null }).id, visible: vis(el) }));
  const shiny = [];
  document.querySelectorAll('.shiny-bound-input').forEach(el => {
    const id = el.id;
    if (!id) return;
    let options = null, optionCount = null;
    if (el.tagName === 'SELECT') {
      if (el.selectize) {
        optionCount = Object.keys(el.selectize.options).length;
        options = Object.keys(el.selectize.options).slice(0, 8);
      } else {
        optionCount = el.options.length;
        options = Array.from(el.options).slice(0, 8).map(o => o.value);
      }
    }
    let v = iv[id];
    if (v !== null && typeof v === 'object') { try { v = JSON.parse(JSON.stringify(v)); } catch { v = String(v); } }
    shiny.push({
      id, tag: el.tagName, inputType: el.type || null,
      classes: (el.className || '').toString().slice(0, 80),
      visible: vis(el), currentValue: v === undefined ? null : v,
      optionCount, sampleOptions: options
    });
  });
  const dtSearch = Array.from(document.querySelectorAll('.dataTables_filter input'))
    .map(el => ({ id: el.id || null, placeholder: el.placeholder || null, visible: vis(el),
      closestId: (el.closest('[id]') || { id: null }).id }));
  const downloads = Array.from(document.querySelectorAll(
    '.shiny-download-link, a.downloadButton, button.downloadButton, [data-testid*="download"]'))
    .map(el => ({ id: el.id || null, tag: el.tagName, text: (el.innerText || '').trim().slice(0, 40),
      testid: el.getAttribute('data-testid'), visible: vis(el) }));
  // gesture-only surfaces (out of scope by instruction, listed for completeness)
  const gestures = Array.from(document.querySelectorAll('.js-plotly-plot, .forceNetwork, .html-widget'))
    .map(el => ({ id: el.id || null, cls: (el.className || '').toString().slice(0, 60), visible: vis(el) }));
  // plain visible buttons that are NOT shiny inputs (aux UI)
  const plainButtons = Array.from(document.querySelectorAll('button'))
    .filter(b => !b.classList.contains('shiny-bound-input') && !b.closest('.selectize-control, .dropdown-menu') && vis(b))
    .map(b => ({ id: b.id || null, text: (b.innerText || '').trim().slice(0, 40), visible: true }))
    .filter(b => b.text);
  return { shiny, dtSearch, dtColFilters, downloads, gestures, plainButtons };
});

async function exploreUi(page, log) {
  const rounds = [];
  const seenIds = new Set();
  const absorb = (tag, inv) => {
    let added = 0;
    for (const e of inv.shiny) if (!seenIds.has(e.id)) { seenIds.add(e.id); added++; }
    rounds.push({ step: tag, newInputs: added, totals: { inputs: inv.shiny.length, dtSearch: inv.dtSearch.length, downloads: inv.downloads.length, gestures: inv.gestures.length } });
    return inv;
  };
  let merged = null;
  const m = inv => {
    if (!merged) merged = inv;
    else for (const k of Object.keys(inv)) {
      const have = new Set((merged[k] || []).map(x => x.id || x.closestId || x.text));
      for (const e of inv[k]) if (!have.has(e.id || e.closestId || e.text)) merged[k].push(e);
    }
  };
  m(absorb('base', await enumerateDom(page)));

  // every bootstrap tab in both navbars + nested tabsets, multi-pass:
  // nested tabsets (heatmap Parameters/Legend, modal tabsets) become
  // visible only after their parent tab is activated
  const visitedTabs = new Set();
  for (let pass = 0; pass < 6; pass++) {
    let clicked = 0;
    const targets = await page.evaluate(() =>
      Array.from(document.querySelectorAll('a[data-toggle="tab"]'))
        .filter(a => a.offsetParent !== null && a.innerText.trim() && !a.dataset.sweepTab)
        .map(a => ({ text: a.innerText.trim() })));
    for (const t of targets) {
      try {
        await page.evaluate(txt => {
          const a = Array.from(document.querySelectorAll('a[data-toggle="tab"]'))
            .find(x => x.offsetParent !== null && x.innerText.trim() === txt && !x.dataset.sweepTab);
          if (a) { a.dataset.sweepTab = '1'; a.click(); }
        }, t.text);
        clicked++;
        await page.waitForTimeout(900);
        m(absorb(`tab:${t.text}`, await enumerateDom(page)));
      } catch (e) { log(`tab click failed: ${t.text}: ${e.message.split('\n')[0]}`); }
    }
    visitedTabs.add(`pass${pass}:${clicked}`);
    if (!clicked) break;
  }
  // dropdown panels (gear boxes, Add column, Map ID, ...)
  const openedDrop = new Set();
  for (let i = 0; i < 15; i++) {
    const target = await page.evaluate(() => {
      const els = Array.from(document.querySelectorAll('.dropdown-toggle'))
        .filter(e => e.offsetParent !== null && !e.closest('.selectize-control'));
      const el = els.find(e => !e.dataset.sweepOpened);
      if (el) { el.dataset.sweepOpened = '1'; return { text: (el.innerText || '').trim().slice(0, 30), idx: Array.prototype.indexOf.call(document.querySelectorAll('.dropdown-toggle'), el) }; }
      return null;
    });
    if (!target) break;
    try {
      await page.evaluate(() => {
        const el = Array.from(document.querySelectorAll('.dropdown-toggle'))
          .filter(e => e.offsetParent !== null && e.dataset.sweepOpened === '1' &&
                       !e.dataset.sweepClicked)[0];
        if (el) { el.dataset.sweepClicked = '1'; el.click(); }
      });
      await page.waitForTimeout(600);
      m(absorb(`dropdown:${target.text}`, await enumerateDom(page)));
    } catch (e) { log(`dropdown open failed: ${target.text}: ${e.message.split('\n')[0]}`); }
  }
  // snapshot modal
  try {
    await page.click('[data-testid="app-snapshot-button"]', { timeout: 5000 });
    await page.waitForTimeout(900);
    m(absorb('modal:snapshot', await enumerateDom(page)));
    await page.keyboard.press('Escape');
    await page.waitForTimeout(400);
  } catch (e) { log(`snapshot modal failed: ${e.message.split('\n')[0]}`); }
  // assistant settings modal (best effort)
  try {
    const opened = await page.evaluate(() => {
      const b = document.querySelector('#app-assistant-settings');
      if (b && b.offsetParent !== null) { b.click(); return true; }
      return false;
    });
    if (opened) {
      await page.waitForTimeout(700);
      m(absorb('modal:assistant-settings', await enumerateDom(page)));
      const closed = await page.evaluate(() => {
        const b = document.querySelector('#app-assistant-settings_cancel');
        if (b) { b.click(); return true; } return false;
      });
      if (closed) await page.waitForTimeout(300);
    }
  } catch (e) { log(`assistant settings failed: ${e.message.split('\n')[0]}`); }

  return { merged, rounds };
}

// ---------------------------------------------------------------------------
// store / capability enumeration (the agent-addressable surface)
// ---------------------------------------------------------------------------
async function fetchStore(page) {
  const snapshot = (await runHook(page, 'store', {})).result;
  const bindingIds = Object.keys(snapshot.values || {});
  const prefixes = new Set();
  for (const id of bindingIds) {
    const parts = id.split('.');
    if (parts.length >= 3) prefixes.add(parts.slice(0, 3).join('.'));
    else if (parts.length === 2) prefixes.add(parts.slice(0, 2).join('.'));
  }
  prefixes.add('tool:');
  const records = [];
  for (const q of prefixes) {
    const res = (await runHook(page, 'capabilities', { query: q, max_results: 50 })).result;
    if (res.hook_error) { log(`capabilities query ${q} failed: ${res.hook_error}`); continue; }
    records.push(...(res.capabilities || []));
  }
  // dedupe by id
  const byId = new Map();
  for (const rec of records) byId.set(rec.id, rec);
  return { snapshot, records: [...byId.values()],
           actions: snapshot.action_registry || [] };
}

// topological order: dependencies first
function topoSort(bindings) {
  const byId = new Map(bindings.map(b => [b.id, b]));
  const done = new Set(), out = [];
  const visit = (b, chain) => {
    if (done.has(b.id)) return;
    if (chain.has(b.id)) return; // cycle guard
    chain.add(b.id);
    for (const d of (b.depends_on || [])) if (byId.has(d)) visit(byId.get(d), chain);
    done.add(b.id); out.push(b);
  };
  for (const b of bindings) visit(b, new Set());
  return out;
}

// fresh allowed values straight before a case (cascade state dependent)
async function freshRecord(page, id) {
  const res = (await runHook(page, 'capabilities', { query: id, max_results: 50 })).result;
  return (res.capabilities || []).find(c => c.id === id) || null;
}

function normalizeAllowed(a) {
  if (a == null) return [];
  if (Array.isArray(a)) return a.map(String);
  if (typeof a === 'object') return Object.keys(a).map(String);
  return [String(a)];
}

function pickProbe(binding, current, marker) {
  const kind = binding.kind;
  const allowed = normalizeAllowed(binding.allowed_values);
  switch (kind) {
    case 'boolean': case 'checkbox':
      return { value: !(current === true), why: 'boolean flip' };
    case 'enum': case 'select': case 'select_cascaded': case 'navbar': case 'tabset': {
      const opts = (allowed || []).filter(v => String(v) !== String(current ?? ''));
      if (!opts.length) return null;
      return { value: opts[0], why: `first allowed != current (${opts[0]})` };
    }
    case 'multi_select': {
      const opts = allowed || [];
      if (!opts.length) return null; // no discoverable choices -> no-probe
      const cur = Array.isArray(current) ? current : [];
      const pick = opts.find(o => !cur.includes(o)) ?? opts[0];
      return { value: [pick], why: `single allowed entry (${pick})` };
    }
    case 'mapping': {
      const keys = allowed || [];
      if (!keys.length) return { value: {}, why: 'no allowed keys: clear mapping' };
      const k = keys.find(x => !String(x).startsWith('GS|')) || keys[0];
      return { value: { [k]: 'AGENTPROBE' }, why: `one key probe (${k})` };
    }
    case 'numeric': case 'integer': case 'slider': {
      const base = (typeof current === 'number' && isFinite(current)) ? current : 0;
      return { value: base + 1, why: `current+1 (${base + 1}); range rejection auto-retried with current-1`, retryDelta: -2 };
    }
    case 'string':
      if (/xcut|ycut/.test(binding.id)) return { value: '1', why: 'cutoff probe (valid numeric string)' };
      return { value: `AGENTPROBE-${marker}`, why: 'distinctive string marker' };
    default:
      return null;
  }
}

// ---------------------------------------------------------------------------
// case runner
// ---------------------------------------------------------------------------
let caseCounter = 0;
const slug = s => String(s).replace(/[^a-zA-Z0-9]+/g, '_').replace(/^_+|_+$/g, '').slice(0, 80);
async function runCase(page, caseSpec, consoleBuf) {
  if (ONLY.length && !ONLY.some(s => String(caseSpec.caseId).includes(s))) return null;
  caseCounter += 1;
  const num = String(caseCounter).padStart(3, '0');
  const name = `${num}-${slug(caseSpec.caseId)}`;
  const file = path.join(CASES, `${name}.json`);
  fs.writeFileSync(file, JSON.stringify(caseSpec, null, 2));

  const result = {
    case: name, caseId: caseSpec.caseId, category: caseSpec.category,
    bindingId: caseSpec.bindingId || null,
    tool: caseSpec.tool, widget: caseSpec.widget || null,
    request: caseSpec.request || null,
    status: null, duration_ms: null, note: '',
    receipt: null, attempts: [], ui_effect: null, store_readback: null,
    console_errors: null
  };
  if (caseSpec.category === 'out-of-scope') {
    result.status = 'out-of-scope';
    result.note = caseSpec.note || 'gesture-only surface; excluded from testing by instruction';
    fs.writeFileSync(path.join(RESULTS, `${name}.json`), JSON.stringify(result, null, 2));
    return result;
  }
  const t0 = Date.now();
  try {
    const before = await snapshotInputs(page);
    const axesBefore = caseSpec.checkAxes ? await axisTitles(page) : null;
    consoleBuf.length = 0;
    let receipt, hookMs = 0;
    const call = async (op, payload) => {
      const h = await runHook(page, op, payload);
      hookMs = Math.max(hookMs, h.ms);
      return h.result;
    };
    if (caseSpec.category === 'dom-only') {
      // attempt every plausible canonical id through set_widgets
      for (const guess of caseSpec.guesses) {
        receipt = await call('widgets', { patch: { [guess]: caseSpec.probeValue } });
        result.attempts.push({ guess, receipt });
        if ((receipt.applied || []).length) break;
      }
    } else {
      receipt = await call(caseSpec.hookOp, caseSpec.request);
      result.attempts.push({ op: caseSpec.hookOp, receipt });
      // no-op retry: the write was skipped (same value); try an alternative
      if (caseSpec.altRequest && (receipt.unchanged || []).length && !receipt.hook_error) {
        receipt = await call(caseSpec.hookOp, caseSpec.altRequest());
        result.attempts.push({ op: 'alt-retry', receipt });
      }
      // numeric range retry
      if (caseSpec.retry && (receipt.rejected || []).length && !receipt.hook_error) {
        const alt = caseSpec.retry();
        if (alt != null) {
          receipt = await call(caseSpec.hookOp, alt.request);
          result.attempts.push({ op: 'range-retry', receipt });
        }
      }
    }
    result.receipt = receipt;
    result.hook_roundtrip_ms = hookMs;
    const eff = await waitUiEffect(page, before, caseSpec.predicate || null);
    const axesAfter = caseSpec.checkAxes ? await axisTitles(page) : null;
    result.ui_effect = {
      changed_inputs: eff.changed.map(c => ({ input: c.input, before: c.before, after: c.after })),
      waited_ms: eff.ms
    };
    if (caseSpec.checkAxes) result.ui_effect.axis_titles = { before: axesBefore, after: axesAfter };

    // store read-back for binding cases. WP-B/F4 note: runHook returns
    // {result, ms}; reading snap.values here (not snap.result.values) made
    // EVERY read-back NULL - the "store read-back is NULL" evidence in the
    // original F4 analysis was this harness bug, not a store anomaly.
    if (caseSpec.category === 'binding' && (receipt.applied || []).length) {
      const snap = await runHook(page, 'store', {});
      result.store_readback = { [caseSpec.bindingId]: ((snap.result || {}).values || {})[caseSpec.bindingId] ?? null };
    }
    result.console_errors = consoleBuf.slice(0, 5);

    // status decision
    const applied = (receipt.applied || []).length > 0;
    const rejected = (receipt.rejected || []).length > 0 || !!receipt.hook_error;
    const uiChanged = eff.changed.length > 0 ||
      (caseSpec.predicate ? await caseSpec.predicate() : false) ||
      (caseSpec.checkAxes && axesBefore && axesAfter && JSON.stringify(axesBefore) !== JSON.stringify(axesAfter));
    if (caseSpec.category === 'action') {
      // WP-C: primary evidence is store-level - the per-action epoch bump
      // is recorded only after the handler body returned - never the
      // input diff (hook churn / re-render blindness do not apply)
      const snapAfter = await runHook(page, 'store', {});
      const after = ((snapAfter.result || {}).actions || {})[caseSpec.actionId] ?? -1;
      result.store_readback = { [caseSpec.actionId]: after };
      if (receipt.error) {
        result.status = 'rejected';
        result.note = String(receipt.error).slice(0, 200);
      } else if ((receipt.ran || []).length && after > (caseSpec.epochBefore ?? -1)) {
        result.status = 'success';
        result.note = `handler ran (epoch ${caseSpec.epochBefore ?? -1} -> ${after})`;
      } else {
        result.status = 'error';
        result.note = 'action did not run (no epoch bump)';
      }
    } else if (eff.timedOut && !uiChanged) result.status = 'timeout';
    else if (hookMs > EFFECT_WAIT_MS) { result.status = 'timeout'; result.note = `tool round-trip itself took ${hookMs}ms (> 30 s budget)`; }
    else if (caseSpec.category === 'semantic-tool') {
      if (receipt.hook_error) { result.status = 'rejected'; result.note = receipt.hook_error; }
      else if ((receipt.rejected || []).length && !(receipt.applied || []).length) { result.status = 'rejected'; result.note = (receipt.rejected[0] || {}).reason || ''; }
      else if (uiChanged) result.status = 'success';
      else if (Object.keys(receipt).filter(k => !['hook_run'].includes(k)).length) result.status = 'success';
      else { result.status = 'no-ui-effect'; result.note = 'tool returned but no browser input/plot change observed'; }
    }
    else if (caseSpec.category === 'dom-only') {
      result.status = applied ? 'unexpected-success' : (rejected ? 'unaddressable' : 'no-effect');
      result.note = caseSpec.note || '';
    }
    else if (applied && uiChanged) result.status = 'success';
    else if (applied && !uiChanged) { result.status = 'no-ui-effect'; result.note = 'store write applied but no browser input/plot change observed'; }
    else if (rejected) { result.status = 'rejected'; result.note = (receipt.rejected || [])[0]?.reason || receipt.hook_error || ''; }
    else if (caseSpec.skipReason) { result.status = 'no-probe'; result.note = 'no meaningful probe value derivable (kind/choices) in this session state'; }
    else result.status = 'unchanged';
  } catch (e) {
    result.status = 'error';
    result.note = String(e && e.message || e).split('\n')[0];
  }
  result.duration_ms = Date.now() - t0;
  fs.writeFileSync(path.join(RESULTS, `${name}.json`), JSON.stringify(result, null, 2));
  console.log(`  ${name} [${result.status}] ${result.duration_ms}ms`);
  return result;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
const logLines = [];
const log = m => { logLines.push(`[${new Date().toISOString()}] ${m}`); console.log(m); };

let exitCode = 0;
try {
  if (!(await waitPort(PORT))) throw new Error('app did not start on port ' + PORT);
  const browser = await chromium.launch({
    executablePath: '/usr/bin/google-chrome', headless: true,
    args: ['--no-sandbox', '--disable-gpu', '--disable-dev-shm-usage', '--disable-webgl', '--disable-webgl2']
  });
  const ctx = await browser.newContext();
  const page = await ctx.newPage();
  const consoleBuf = [];
  page.on('console', m => { if (m.type() === 'error') consoleBuf.push(m.text().slice(0, 200)); });
  page.on('pageerror', e => consoleBuf.push(String(e).slice(0, 200)));
  await page.goto(`http://127.0.0.1:${PORT}/`, { timeout: 90000 });
  await page.waitForFunction(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-eset'] &&
    document.querySelector('#app-contents').offsetParent !== null, null, { timeout: 90000 });
  await page.evaluate(sel => { document.querySelector(sel).style.display = 'block'; }, HOOK.box);
  log('app ready');

  // ================= 1. inventory =================
  const { merged: dom, rounds } = await exploreUi(page, log);
  fs.writeFileSync(path.join(ART, 'dom_inventory.json'), JSON.stringify(dom, null, 2));
  fs.writeFileSync(path.join(ART, 'exploration_log.json'), JSON.stringify(rounds, null, 2));
  log(`DOM inventory: ${dom.shiny.length} shiny inputs, ${dom.dtSearch.length} DT search boxes, ${dom.downloads.length} download controls, ${dom.gestures.length} gesture surfaces`);

  const store = await fetchStore(page);
  fs.writeFileSync(path.join(ART, 'store_snapshot.json'), JSON.stringify(store.snapshot, null, 2));
  fs.writeFileSync(path.join(ART, 'capability_registry.json'), JSON.stringify(store.records, null, 2));
  const widgetRecords = store.records.filter(c => !c.id.startsWith('tool:'));
  const internalIds = Object.keys(store.snapshot.values || {}).filter(id => !widgetRecords.some(w => w.id === id));
  fs.writeFileSync(path.join(ART, 'internal_store_ids.json'), JSON.stringify(internalIds, null, 2));
  log(`store: ${Object.keys(store.snapshot.values || {}).length} ids total, ${widgetRecords.length} agent-writable widgets, ${internalIds.length} internal`);

  const overviewHook = await runHook(page, 'overview', {});
  const overview = overviewHook.result;
  fs.writeFileSync(path.join(ART, 'overview_state.json'), JSON.stringify(overview, null, 2));

  // ================= 2. binding cases (every store widget) =================
  const ordered = topoSort(widgetRecords);
  const bindingResults = [];
  let marker = 0;
  const currentFrom = (vals, id) => {
    let v = vals?.[id] ?? null;
    return v === undefined ? null : v;
  };
  for (const b of ordered) {
    marker += 1;
    const snap = await runHook(page, 'store', {});
    const current = currentFrom(snap.result.values, b.id);
    let fresh = (await freshRecord(page, b.id)) || b;
    let probe = pickProbe(fresh, current, marker);
    for (let attempt = 0; attempt < 2 && !probe; attempt++) {
      // choices may be settling (cascade parents, lazy providers); re-probe
      await new Promise(s => setTimeout(s, 3000));
      fresh = (await freshRecord(page, b.id)) || fresh;
      probe = pickProbe(fresh, current, marker);
    }
    const altProbe = pickProbe({ ...fresh, allowed_values: normalizeAllowed(fresh.allowed_values).slice(1) }, current, marker + 1000);
    const spec = {
      caseId: b.id, category: 'binding', tool: 'set_widgets',
      bindingId: b.id,
      widget: { binding_id: b.id, kind: fresh.kind, label: fresh.label, help: fresh.help_text, depends_on: fresh.depends_on, operation: fresh.operation },
      hookOp: 'widgets',
      request: { patch: { [b.id]: probe ? probe.value : current } },
      fresh_current: current,
      probeNote: probe ? probe.why : 'no probe derivable (kind/choices)'
    };
    if (!probe) { spec.skipReason = 'no probe'; }
    else if (altProbe) spec.altRequest = () => ({ patch: { [b.id]: altProbe.value } });
    if (probe && probe.retryDelta !== undefined) {
      const base = (typeof current === 'number' && isFinite(current)) ? current : 0;
      spec.retry = () => ({ patch: { [b.id]: Math.max(0, base + probe.retryDelta) } });
    }
    const bres = await runCase(page, spec, consoleBuf);
    if (bres) bindingResults.push(bres);
  }
  log(`binding cases: ${bindingResults.length} done`);

  // ================= 2b. action cases (WP-C run_control tier) =================
  // Every registered action runs through the runcontrol hook == the
  // exact run_control tool pathway. Status is decided on the store-level
  // epoch bump (proof the handler body ran), never on input diffs.
  const actionResults = [];
  for (const a of (store.actions || [])) {
    const snap = await runHook(page, 'store', {});
    const epochBefore = ((snap.result || {}).actions || {})[a.id] ?? -1;
    const ares = await runCase(page, {
      caseId: a.id, category: 'action', tool: 'run_control',
      actionId: a.id,
      widget: { action_id: a.id, label: a.label, help: a.help,
                input_id: a.input_id, requires_confirm: a.requires_confirm,
                available: a.available },
      hookOp: 'runcontrol',
      request: { action_id: a.id },
      epochBefore
    }, consoleBuf);
    if (ares) actionResults.push(ares);
  }
  log(`action cases: ${actionResults.length} done`);

  // ================= 3. semantic tool cases =================
  const semantic = [];
  const push = async s => { const res = await runCase(page, s, consoleBuf); if (res) semantic.push(res); };
  await push({ caseId: 'tool:get_omics_viewer_state', category: 'semantic-tool', tool: 'get_omics_viewer_state',
    hookOp: 'overview', request: {}, predicate: null });
  await push({ caseId: 'tool:find_controls(color)', category: 'semantic-tool', tool: 'find_controls',
    hookOp: 'capabilities', request: { query: 'color', max_results: 5 } });
  await push({ caseId: 'tool:describe_control(axis_mode)', category: 'semantic-tool', tool: 'describe_control',
    hookOp: 'capabilities', request: { id: 'dataspace.feature_space.axis_mode' } });
  await push({ caseId: 'tool:set_scatter_view(feature,quick)', category: 'semantic-tool', tool: 'set_scatter_view',
    hookOp: 'scatter', request: { space: 'feature', quick_view_id: 'volcano_RE_vs_ME' },
    widget: { semantic: 'feature-space scatter quick view' } });
  await push({ caseId: 'tool:set_scatter_view(feature,axes)', category: 'semantic-tool', tool: 'set_scatter_view',
    hookOp: 'scatter', request: { space: 'feature', x_axis: 'ttest|RE_vs_ME|mean.diff', y_axis: 'ttest|RE_vs_ME|pvalue' },
    checkAxes: true, widget: { semantic: 'feature-space scatter custom axes' } });
  const sampleQv = (overview.quick_views?.sample || [])[0]?.id;
  await push({ caseId: 'tool:set_scatter_view(sample,quick)', category: 'semantic-tool', tool: 'set_scatter_view',
    hookOp: 'scatter', request: { space: 'sample', quick_view_id: sampleQv || 'pca' },
    widget: { semantic: 'sample-space scatter quick view' } });
  await push({ caseId: 'tool:set_scatter_view(sample,axes)', category: 'semantic-tool', tool: 'set_scatter_view',
    hookOp: 'scatter', request: { space: 'sample', x_axis: 'General|All|MDR', y_axis: 'General|All|Doubleing.Time' },
    checkAxes: true, widget: { semantic: 'sample-space scatter custom axes' } });
  await push({ caseId: 'tool:set_omics_viewer_state(data_space_tab)', category: 'semantic-tool', tool: 'set_omics_viewer_state',
    hookOp: 'state', request: { data_space_tab: 'Heatmap' } });
  await push({ caseId: 'tool:set_omics_viewer_state(analysis_space_tab)', category: 'semantic-tool', tool: 'set_omics_viewer_state',
    hookOp: 'state', request: { analysis_space_tab: 'Geneshot' } });
  await push({ caseId: 'tool:set_omics_viewer_state(features)', category: 'semantic-tool', tool: 'set_omics_viewer_state',
    hookOp: 'state', request: { features: FEATURE_IDS } });
  await push({ caseId: 'tool:set_omics_viewer_state(samples)', category: 'semantic-tool', tool: 'set_omics_viewer_state',
    hookOp: 'state', request: { samples: SAMPLE_IDS } });
  await push({ caseId: 'tool:set_omics_viewer_state(features,clear)', category: 'semantic-tool', tool: 'set_omics_viewer_state',
    hookOp: 'state', request: { features: [] } });
  await push({ caseId: 'tool:set_enrichment_parameters(ora)', category: 'semantic-tool', tool: 'set_enrichment_parameters',
    hookOp: 'enrichment', request: { method: 'ora', collapse: 'General|All|Gene.name' } });
  await push({ caseId: 'tool:set_enrichment_parameters(fgsea)', category: 'semantic-tool', tool: 'set_enrichment_parameters',
    hookOp: 'enrichment', request: { method: 'fgsea', collapse: 'ttest|RE_vs_ME|mean.diff' } });
  await push({ caseId: 'tool:set_enrichment_parameters(bogus_pathway)', category: 'semantic-tool', tool: 'set_enrichment_parameters',
    hookOp: 'enrichment', request: { method: 'ora', selected_pathway: 'AGENTPROBE_NOT_A_PATHWAY' } });
  await push({ caseId: 'tool:set_table_view(feature_table)', category: 'semantic-tool', tool: 'set_table_view',
    hookOp: 'tableview', request: { table: 'feature_table', multi_selection: true, column_filters: { 'General|All|Gene.name': 'CNOT' }, page: 2 } });
  await push({ caseId: 'tool:set_table_view(sample_table,clear)', category: 'semantic-tool', tool: 'set_table_view',
    hookOp: 'tableview', request: { table: 'sample_table', clear_filters: true } });
  await push({ caseId: 'tool:set_table_view(expression_table)', category: 'semantic-tool', tool: 'set_table_view',
    hookOp: 'tableview', request: { table: 'expression_table', multi_selection: false } });

  // ============ 3b. WP-E chained semantic case: stringdb run -> select row
  // (F5 dead-end chain): resultspace.stringdb.selected_row is only
  // operable AFTER a successful run whose gene set has STRING
  // interactions. The generic binding case cannot probe it (choices are
  // empty pre-run); this case proves the WHOLE chain with real
  // preconditions: STRING-connected features selected, tab opened (the
  // enrichment computes lazily and the table must initialize), run via
  // run_control, choices materialized, selected_row written through
  // set_widgets, and the table row actually selected (DOM-level proof).
  {
    const CHAIN_FEATURES = ['X1452.EGFR', 'X1477.MYC', 'X1570.TP53',
                           'X2033.BRAF', 'X2371.CCND1'];
    const SR = 'resultspace.stringdb.selected_row';
    const t0c = Date.now();
    const cres = {
      case: '', caseId: 'tool:stringdb_chain_run_then_select',
      category: 'semantic-tool',
      tool: 'set_omics_viewer_state + run_control + set_widgets',
      widget: { chain: 'select STRING-connected features -> run -> set selected_row' },
      request: { features: CHAIN_FEATURES, run: 'resultspace.stringdb.run', set: SR },
      status: null, duration_ms: null, note: '', receipt: null, attempts: [],
      ui_effect: null, store_readback: null, console_errors: [] };
    const step = async (op, payload) => (await runHook(page, op, payload)).result;
    try {
      await page.evaluate(() => {
        // the result-space navbar (and its StringDB tab) materializes
        // after the dataset settles; poll briefly for the anchor
        const a = Array.from(document.querySelectorAll('a[data-toggle="tab"]'))
          .find(x => x.offsetParent !== null && /StringDB/i.test(x.innerText));
        if (a) a.click();
      });
      await new Promise(s => setTimeout(s, 1500));
      const st = await step('state', { features: CHAIN_FEATURES });
      // earlier sweep cases probe the taxonomy binding with a marker
      // string (AGENTPROBE-*); restore a valid NCBI taxid before the
      // run or the STRING API answers for an unknown species (empty
      // results) - exactly what a well-behaved model would do after
      // reading the current taxonomy value
      const tax = await step('widgets', {
        patch: { 'resultspace.stringdb.taxonomy': '9606' } });
      const run = await step('runcontrol', { action_id: 'resultspace.stringdb.run' });
      cres.attempts.push({
        feature_count: st?.feature_count ?? st?.selection?.features?.count ?? null,
        taxonomy: tax?.applied ?? null,
        run: { ran: run?.ran, epoch: run?.epoch, error: run?.error } });
      if (!run || !run.ran || run.error)
        throw new Error(`run_control failed: ${JSON.stringify(run?.error || run)}`);
      let terms = null;
      // two run attempts x 30 s: the STRING enrichment endpoint is a real
      // external API and can transiently fail/timeout; a second signal
      // bump re-fires the lazy eventReactive chain
      for (let attempt = 0; attempt < 2 && !(terms && terms.length); attempt++) {
        if (attempt > 0) {
          const rerun = await step('runcontrol', { action_id: 'resultspace.stringdb.run' });
          cres.attempts.push({ rerun: { ran: rerun?.ran, epoch: rerun?.epoch } });
        }
        for (let i = 0; i < 15; i++) {
          await new Promise(s => setTimeout(s, 2000));
          const rec = await step('capabilities', { query: SR, max_results: 3 });
          const r = (rec.capabilities || []).find(c => c.id === SR);
          terms = r?.allowed_values || [];
          if (terms.length) break;
        }
      }
      if (!terms.length) {
        const notif = await page.evaluate(() =>
          Array.from(document.querySelectorAll('.shiny-notification'))
            .map(n => n.innerText.slice(0, 80)));
        throw new Error(`enrichment terms did not materialize (60 s x 2 runs);` +
          ` notifications: ${JSON.stringify(notif)}`);
      }
      const term = terms[0];
      const set = await step('widgets', { patch: { [SR]: [term] } });
      cres.receipt = set;
      cres.store_readback = { [SR]: term };
      const appliedOK = set.applied === SR ||
        (Array.isArray(set.applied) ? set.applied.includes(SR)
                                    : String(set.applied || '').includes(SR));
      let rowSelected = false;
      for (let i = 0; i < 15 && !rowSelected; i++) {
        await new Promise(s => setTimeout(s, 1000));
        rowSelected = await page.evaluate(() => {
          for (const t of Array.from(document.querySelectorAll('table.dataTable'))) {
            let el = t;
            while (el && !(el.id && /stringdb/.test(el.id))) el = el.parentElement;
            if (el && t.offsetParent !== null &&
                t.querySelectorAll('tbody tr.selected').length) return true;
          }
          return false;
        });
      }
      cres.ui_effect = { strtab_row_selected: rowSelected, term };
      cres.status = appliedOK && rowSelected ? 'success' : 'error';
      if (cres.status !== 'success')
        cres.note = `applied=${JSON.stringify(set.applied)} rowSelected=${rowSelected}`;
    } catch (e) {
      cres.status = 'error';
      cres.note = String(e.message || e).slice(0, 200);
    }
    cres.duration_ms = Date.now() - t0c;
    caseCounter += 1;
    cres.case = `${String(caseCounter).padStart(3, '0')}-${slug(cres.caseId)}`;
    fs.writeFileSync(path.join(RESULTS, `${cres.case}.json`), JSON.stringify(cres, null, 2));
    console.log(`  ${cres.case} [${cres.status}] ${cres.duration_ms}ms`);
    semantic.push(cres);
  }

  // ============ 3c. WP-G1 discovery cases: find_controls surface.
  // Discovery failures were the original complaint (B3/B3b of
  // verify_search_highlight): find_controls("highlight") must return the
  // Search cascade + the gene picker; find_controls("gene name") must
  // return something. Custom status logic - the generic semantic rule
  // would grade any non-error receipt as success regardless of matches.
  const discoveryCase = async (caseId, query, okIds, okLabelRe) => {
    const t0d = Date.now();
    const res = (await runHook(page, 'capabilities', { query, max_results: 50 })).result;
    const caps = res.capabilities || [];
    const ids = caps.map(c => c.id);
    const idOK = okIds.every(re => ids.some(id => re.test(id)));
    const labelOK = !okLabelRe || caps.some(c => okLabelRe.test(c.label || ''));
    const dres = {
      case: '', caseId, category: 'semantic-tool', tool: 'find_controls',
      widget: { query }, request: { query, max_results: 50 },
      status: (!res.hook_error && caps.length > 0 && idOK && labelOK)
        ? 'success' : 'error',
      duration_ms: Date.now() - t0d,
      note: `matches=${caps.length} first=${ids[0] || '(none)'}` +
        (idOK ? '' : ` MISSING expected ${okIds.map(r => r.source).join(',')}`),
      receipt: { match_count: res.match_count, capability_count: res.capability_count,
                 fuzzy: res.fuzzy, matched_ids: ids.slice(0, 12) },
      attempts: [], ui_effect: null, store_readback: null, console_errors: [] };
    caseCounter += 1;
    dres.case = `${String(caseCounter).padStart(3, '0')}-${slug(caseId)}`;
    fs.writeFileSync(path.join(RESULTS, `${dres.case}.json`), JSON.stringify(dres, null, 2));
    console.log(`  ${dres.case} [${dres.status}] ${dres.duration_ms}ms`);
    semantic.push(dres);
  };
  await discoveryCase('tool:find_controls(highlight)', 'highlight',
    [/attr4\.search_variable$/, /attr4\.search_value$/], /Search/);
  await discoveryCase('tool:find_controls(gene name)', 'gene name',
    [/attr4\.search_value$/], null);

  log(`semantic tool cases: ${semantic.length} done`);

  // ================= 4. empirical binding -> DOM input mapping =================
  const mapping = {};
  for (const res of bindingResults) {
    if (res.status === 'success' || res.status === 'no-ui-effect')
      mapping[res.caseId] = (res.ui_effect?.changed_inputs || []).map(c => c.input);
  }
  for (const res of semantic) {
    if (res.ui_effect) mapping[`semantic:${res.caseId}`] = (res.ui_effect.changed_inputs || []).map(c => c.input);
  }
  fs.writeFileSync(path.join(ART, 'mapping_findings.json'), JSON.stringify(mapping, null, 2));
  const mappedInputs = new Set(Object.values(mapping).flat());

  // ================= 5. dom-only cases (user-facing widgets with NO store binding) =================
  const guessIds = domId => {
    const base = domId.replace(/^app-/, '');
    const dotted = base.replace(/-/g, '.');
    const out = new Set([dotted]);
    const m = base.match(/^(dataspace|resultspace)[.-](.*)$/);
    if (m) { out.add(`${m[1]}.${m[2].replace(/-/g, '.')}`); }
    return [...out].slice(0, 3);
  };
  const domOnly = [];
  // WP-C: buttons that back registered actions are covered by the action
  // tier - they leave the dom-only inventory (mapped through input_id,
  // which carries the full dom id including any app- prefix)
  const actionDomIds = new Set((store.actions || [])
    .filter(a => a.input_id)
    .map(a => String(a.input_id).startsWith('app-') ? a.input_id : 'app-' + a.input_id));
  const addDomOnly = (ident, desc, probeValue, note) => domOnly.push({ ident, desc, probeValue, note });
  for (const e of dom.shiny) {
    if (mappedInputs.has(e.id)) continue;
    if (actionDomIds.has(e.id)) continue;
    if (/agentTestHooks/.test(e.id)) continue; // test infrastructure, not user-facing
    const covered = /-select-(analysis|subset|variable)$/.test(e.id) && /tab_(pheno|feature|expr)/.test(e.id);
    const badge = /quickViews-badge_/.test(e.id);
    addDomOnly(e.id,
      `${e.tag}${e.inputType ? ':' + e.inputType : ''} visible=${e.visible}`,
      e.tag === 'SELECT' ? (e.sampleOptions?.[0] ?? 'x') : (typeof e.currentValue === 'boolean' ? !e.currentValue : 'AGENTPROBE'),
      covered ? 'Add-column triselector: content covered semantically by the table columns widget' :
      badge ? 'Quick-view badge: covered semantically by set_scatter_view(quick_view_id)' : '');
  }
  for (const e of dom.dtSearch) addDomOnly(e.id || e.closestId || 'dt-global-search', 'DT global search input', 'AGENTPROBE',
    'DataTables global search box - not a Shiny input binding');
  const seenColFilter = new Set();
  for (const e of (dom.dtColFilters || [])) {
    const key = e.closestId || 'dt';
    if (seenColFilter.has(key)) continue; seenColFilter.add(key);
    addDomOnly(`dtcolfilter:${key}`, 'DT per-column filter input', 'AGENTPROBE',
      'per-column filter box - covered semantically by the column_filters store binding (set_table_view / set_widgets)');
  }
  const seenBtn = new Set();
  for (const e of (dom.plainButtons || [])) {
    const key = e.id || e.text;
    if (seenBtn.has(key)) continue; seenBtn.add(key);
    addDomOnly(`button:${key}`, `plain button (${e.text})`, 'AGENTPROBE',
      /previous|next|First|Last/i.test(e.text) ? 'DT pagination button - covered by the page store binding' :
      /Sort/i.test(e.text) ? 'DT sort control - column ordering not exposed to the agent' : '');
  }
  for (const e of dom.downloads) addDomOnly(e.id || e.testid || 'download', `download control (${e.text || e.testid})`, 'AGENTPROBE',
    'download output - no value semantics');
  const domResults = [];
  marker = 0;
  for (const w of domOnly) {
    marker += 1;
    const guesses = w.ident.startsWith('app-') ? guessIds(w.ident) : ['app', 'dataspace.global_search', 'resultspace.global_search'];
    const dres = await runCase(page, {
      caseId: `dom:${w.ident}`, category: 'dom-only', tool: 'set_widgets(attempt)',
      widget: { dom_id: w.ident, description: w.desc },
      guesses, probeValue: w.probeValue, note: w.note
    }, consoleBuf);
    if (dres) domResults.push(dres);
  }
  log(`dom-only cases: ${domResults.length} done`);

  // ================= 6. out-of-scope gesture surfaces =================
  const oosResults = [];
  for (const g of dom.gestures) {
    const ores = await runCase(page, {
      caseId: `gesture:${g.id || g.cls.split(' ')[0]}`, category: 'out-of-scope',
      widget: g, note: 'plotly/network interactive gesture surface - excluded from testing by instruction (no mouse-gesture mimicry)'
    }, consoleBuf);
    if (ores) oosResults.push(ores);
  }

  // ================= 7. summary =================
  const all = [...bindingResults, ...actionResults, ...semantic, ...domResults, ...oosResults];
  const byCat = {};
  for (const res of all) {
    byCat[res.category] = byCat[res.category] || {};
    byCat[res.category][res.status] = (byCat[res.category][res.status] || 0) + 1;
  }
  const summary = {
    generated: new Date().toISOString(),
    totals: { cases: all.length, case_files: all.length, result_files: all.length },
    by_category_status: byCat,
    binding_failures: bindingResults.filter(x => x.status !== 'success').map(x => ({ case: x.case, status: x.status, note: x.note, binding: x.bindingId })),
    semantic_failures: semantic.filter(x => x.status !== 'success' && x.status !== 'rejected').map(x => ({ case: x.case, status: x.status, note: x.note })),
    dom_only_addressed: domResults.filter(x => x.status === 'unexpected-success').map(x => x.case),
    unaddressable_count: domResults.filter(x => x.status === 'unaddressable').length,
    slowest: [...all].sort((a, b) => (b.duration_ms || 0) - (a.duration_ms || 0)).slice(0, 10)
      .map(x => ({ case: x.case, ms: x.duration_ms, status: x.status }))
  };
  fs.writeFileSync(path.join(OUT, 'summary.json'), JSON.stringify(summary, null, 2));

  const md = [
    `# Widget sweep summary`, ``, `Generated: ${summary.generated}`, ``,
    `| category | status | count |`, `|---|---|---|`];
  for (const [cat, st] of Object.entries(byCat)) for (const [s, n] of Object.entries(st)) md.push(`| ${cat} | ${s} | ${n} |`);
  md.push('', '## Binding (store widget) non-success', '');
  for (const f of summary.binding_failures) md.push(`- **${f.binding}** [${f.status}] ${f.note || ''}`);
  md.push('', '## Semantic-tool non-success', '');
  for (const f of summary.semantic_failures) md.push(`- **${f.case}** [${f.status}] ${f.note || ''}`);
  md.push('', '## DOM-only widgets the agent unexpectedly COULD set', '');
  md.push(...(summary.dom_only_addressed.map(c => `- ${c}`)));
  fs.writeFileSync(path.join(OUT, 'summary.md'), md.join('\n'));

  fs.writeFileSync(path.join(OUT, 'app_log_tail.txt'), rlog.slice(-4000));
  await browser.close();
  console.log(`\nSWEEP COMPLETE: ${all.length} cases, files in cases/ and results/`);
  console.log(JSON.stringify(summary.by_category_status, null, 2));
} catch (e) {
  console.error('FATAL', e);
  fs.writeFileSync(path.join(OUT, 'fatal.log'), String(e && e.stack || e));
  exitCode = 1;
}
process.exit(exitCode);
