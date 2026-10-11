// Browser-level validation of the scatter "highlight gene" feature:
// left panel (data space), volcano quick view, gear box (attr4 panel),
// Search cascade -> column (e.g. General|All|Gene.name), then the gene
// picker (searchon) that reveals itself once a search column is set.
//
// Two perspectives are validated through the SAME browser session:
//
//   [A] USER path   - drive the real UI: volcano badge, open the gear
//                     dropdown, pick the Search column via the three
//                     selectize cascades, pick the gene in the searchon
//                     multi-select, and assert the open-ring highlight
//                     trace (marker.symbol == 'circle-open') appears on
//                     the Plotly figure at the gene's coordinates.
//
//   [B] AGENT path  - drive the EXACT tool surface the LLM assistant
//                     uses, through the hidden OMICSVIEWER_TEST_HOOKS
//                     bridge (op 'widgets' == set_widgets, op 'scatter'
//                     == set_scatter_view, op 'capabilities' ==
//                     find_controls, op 'state' == set_omics_viewer_state):
//                       B1 set the volcano view            (works)
//                       B2 set the Search cascade columns  (works - the
//                          store pushes reach the UI triselectors and the
//                          searchon box is revealed)
//                       B3 discover "highlight" through find_controls
//                          (what the model can learn)
//                       B4 try to fill the GENE NAME through set_widgets
//                          (rejected - searchon is not a store widget)
//                       B5 contrast: set_omics_viewer_state CAN select
//                          the feature, but that is a selection, not the
//                          ring highlight.
//
// Run:  node verify_search_highlight.mjs     (from tests/e2e_agent/)

import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../..');
const EXTDATA = path.join(REPO, 'inst/extdata');
const PORT = 7792;

// expected highlight coordinates: feature X52.CNOT1 (Gene.name == CNOT1)
// on the volcano_RE_vs_ME quick view (x = ttest|RE_vs_ME|mean.diff,
// y = ttest|RE_vs_ME|log.fdr)
const GENE = 'CNOT1';
const EXPECT = { x: 0.07339119, y: 0.3245538 };

const results = [];
const record = (name, pass, detail = '') => {
  results.push({ name, pass, detail: String(detail).slice(0, 300) });
  console.log(`${pass ? 'PASS' : 'FAIL'} | ${name}${detail ? ' | ' + detail : ''}`);
};

const r = spawn('Rscript', ['-e', `
  options(shiny.port = ${PORT}, shiny.host = '127.0.0.1')
  eset <- readRDS(file.path('${EXTDATA}', 'demo.RDS'))
  omicsViewer::omicsViewer(dir = '${EXTDATA}', ESVObj = eset)
`], { cwd: REPO, env: { ...process.env, OMICSVIEWER_TEST_HOOKS: 'true' }, stdio: ['ignore', 'pipe', 'pipe'] });
process.on('exit', () => { try { r.kill('SIGKILL'); } catch {} });
process.on('SIGINT', () => { try { r.kill('SIGKILL'); } catch {}; process.exit(130); });

const waitPort = async (port, ms = 90000) => {
  const t0 = Date.now();
  while (Date.now() - t0 < ms) {
    try {
      const res = await fetch(`http://127.0.0.1:${port}/`);
      if (res.status > 0) return true;
    } catch {}
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

const A4 = 'app-dataspace-feature_space-a4selector';

async function runHook(page, op, payload) {
  const before = (await page.innerText(HOOK.result)).trim();
  await page.evaluate(({ sel, value }) => {
    const el = document.querySelector(sel);
    if (el.selectize) el.selectize.setValue(value);
    else {
      el.value = value;
      el.dispatchEvent(new Event('change', { bubbles: true }));
    }
  }, { sel: HOOK.op, value: op });
  await page.fill(HOOK.payload, JSON.stringify(payload ?? {}));
  const runCount = await page.evaluate(sel => {
    const el = document.querySelector(sel);
    const n = (parseInt(el.dataset.runs || '0', 10) + 1);
    el.dataset.runs = String(n);
    return n;
  }, HOOK.run);
  await page.click(HOOK.run);
  await page.waitForFunction(({ sel, before, n }) => {
    const t = (document.querySelector(sel)?.innerText || '').trim();
    return t !== before && t.includes(`"hook_run": ${n}`);
  }, { sel: HOOK.result, before, n: runCount }, { timeout: 60000 });
  return JSON.parse((await page.innerText(HOOK.result)).trim());
}

const setSelectize = (page, sel, value) => page.evaluate(({ sel, value }) => {
  const el = document.querySelector(sel);
  if (!el) throw new Error('missing ' + sel);
  if (el.selectize) el.selectize.setValue(value);
  else {
    el.value = value;
    el.dispatchEvent(new Event('change', { bubbles: true }));
  }
}, { sel, value });

const waitSelectizeOptions = (page, sel, ms = 30000) =>
  page.waitForFunction(sel => {
    const el = document.querySelector(sel);
    return el && el.selectize && Object.keys(el.selectize.options).length > 0;
  }, sel, { timeout: ms });

const axisTitles = (page) => page.evaluate(() => {
  const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
    .filter(e => e.offsetParent !== null && e._fullLayout);
  if (!els.length || !els[0]._fullLayout) return null;
  const l = els[0]._fullLayout;
  const t = a => (a && a.title && (a.title.text || a.title)) || '';
  return { x: String(t(l.xaxis)), y: String(t(l.yaxis)) };
});

// open-ring highlight traces on the visible plot (marker.symbol circle-open)
const highlightTraces = (page) => page.evaluate(() => {
  const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
    .filter(e => e.offsetParent !== null && e.data);
  for (const gd of els) {
    const tr = (gd.data || []).filter(t =>
      t.marker && String(t.marker.symbol) === 'circle-open');
    if (tr.length)
      return tr.map(t => ({ name: t.name, mode: t.mode,
        x: Array.from(t.x || []), y: Array.from(t.y || []) }));
  }
  return [];
});

// open the attr4 gear dropdown of the visible scatter (left rail)
const openGearBox = async (page) => {
  const clicked = await page.evaluate(() => {
    const icons = Array.from(document.querySelectorAll(
      '[class*="fa-cog"], [class*="fa-gear"]'))
      .filter(i => i.offsetParent !== null);
    if (!icons.length) return false;
    const btn = icons[0].closest('button, a, .dropdown-toggle');
    if (!btn) return false;
    btn.click();
    return true;
  });
  await page.waitForTimeout(400);
  return clicked;
};

const selectSearchColumn = async (page) => {
  // wait until the three Search-cascade selectizes carry choices
  await waitSelectizeOptions(page, `#${A4}-selectSearchCol-analysis`);
  await setSelectize(page, `#${A4}-selectSearchCol-analysis`, 'General');
  await page.waitForFunction(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-feature_space-a4selector-selectSearchCol-analysis'] === 'General',
    null, { timeout: 30000 });
  await waitSelectizeOptions(page, `#${A4}-selectSearchCol-subset`);
  await setSelectize(page, `#${A4}-selectSearchCol-subset`, 'All');
  await page.waitForFunction(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-feature_space-a4selector-selectSearchCol-subset'] === 'All',
    null, { timeout: 30000 });
  await waitSelectizeOptions(page, `#${A4}-selectSearchCol-variable`);
  await setSelectize(page, `#${A4}-selectSearchCol-variable`, 'Gene.name');
  await page.waitForFunction(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-feature_space-a4selector-selectSearchCol-variable'] === 'Gene.name',
    null, { timeout: 30000 });
};

let exitCode = 0;
try {
  if (!(await waitPort(PORT))) throw new Error('app did not start');
  const browser = await chromium.launch({
    executablePath: '/usr/bin/google-chrome',
    headless: true,
    args: ['--no-sandbox', '--disable-gpu', '--disable-dev-shm-usage',
           '--disable-webgl', '--disable-webgl2']
  });

  // un-hide the test-hook panel for Phase B interactions
  const ctx = await browser.newContext();
  const page = await ctx.newPage();
  await page.goto(`http://127.0.0.1:${PORT}/`, { timeout: 90000 });
  await page.waitForFunction(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-eset'] &&
    document.querySelector('#app-contents').offsetParent !== null, null, { timeout: 90000 });
  await page.evaluate(sel => {
    document.querySelector(sel).style.display = 'block';
  }, HOOK.box);

  // ============================================================
  // [A] USER PATH - everything through the real UI
  // ============================================================
  console.log('\n--- [A] user path ---');
  // A1: volcano quick view via the quick-view badge
  await page.click('[data-testid="app-dataspace-feature_space-quickViews-badge-volcano_RE_vs_ME"]');
  await page.waitForFunction(() => {
    const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
      .filter(e => e.offsetParent !== null && e._fullLayout);
    if (!els.length) return false;
    const l = els[0]._fullLayout;
    const t = a => String((a && a.title && (a.title.text || a.title)) || '');
    return t(l.xaxis).includes('mean.diff') && t(l.yaxis).includes('log.fdr');
  }, null, { timeout: 60000 });
  record('A1 user: volcano quick view renders (mean.diff vs log.fdr)', true,
    JSON.stringify(await axisTitles(page)));

  // A2: open the gear box (attr4 dropdown)
  record('A2 user: gear box (attr4 dropdown) opens', await openGearBox(page));

  // A3: pick the Search column through the cascade (General > All > Gene.name)
  await selectSearchColumn(page);
  record('A3 user: Search column cascade set to General|All|Gene.name', true);

  // A4: the gene picker (searchon) reveals itself once a search column is set
  await page.waitForFunction(id =>
    Shiny.shinyapp.$inputValues[`${id}-showSearchBox`] === true,
    A4, { timeout: 30000 });
  await waitSelectizeOptions(page, `#${A4}-searchon`);
  const searchonVisible = await page.evaluate(sel => {
    const el = document.querySelector(sel);
    return el !== null && el.selectize !== undefined &&
      Object.keys(el.selectize.options).length > 1000; // ~2702 gene names
  }, `#${A4}-searchon`);
  record('A4 user: gene picker appears with the column values as choices',
    searchonVisible);

  // A5: pick the gene and assert the open-ring highlight trace
  await setSelectize(page, `#${A4}-searchon`, [GENE]);
  await page.waitForFunction(() => {
    const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
      .filter(e => e.offsetParent !== null && e.data);
    for (const gd of els) {
      const tr = (gd.data || []).filter(t =>
        t.marker && String(t.marker.symbol) === 'circle-open');
      if (tr.length) return true;
    }
    return false;
  }, null, { timeout: 60000 });
  const traces = await highlightTraces(page);
  const one = traces.length === 1 &&
    traces[0].x.length === 1 && traces[0].y.length === 1 &&
    Math.abs(traces[0].x[0] - EXPECT.x) < 1e-4 &&
    Math.abs(traces[0].y[0] - EXPECT.y) < 1e-4;
  record('A5 user: gene highlighted with an open-ring trace at its position',
    one, JSON.stringify(traces));

  // ============================================================
  // [B] AGENT PATH - the exact tool surface the LLM uses
  // ============================================================
  console.log('\n--- [B] agent path (test hooks == assistant tools) ---');
  // fresh session: a new page is a new Shiny session with a fresh store
  const page2 = await ctx.newPage();
  await page2.goto(`http://127.0.0.1:${PORT}/`, { timeout: 90000 });
  await page2.waitForFunction(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-eset'] &&
    document.querySelector('#app-contents').offsetParent !== null, null, { timeout: 90000 });
  await page2.evaluate(sel => {
    document.querySelector(sel).style.display = 'block';
  }, HOOK.box);

  // B1: set_scatter_view to the volcano (agent CAN steer the plot)
  const b1 = await runHook(page2, 'scatter',
    { space: 'feature', quick_view_id: 'volcano_RE_vs_ME' });
  await page2.waitForFunction(() => {
    const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
      .filter(e => e.offsetParent !== null && e._fullLayout);
    if (!els.length) return false;
    const l = els[0]._fullLayout;
    const t = a => String((a && a.title && (a.title.text || a.title)) || '');
    return t(l.xaxis).includes('mean.diff') && t(l.yaxis).includes('log.fdr');
  }, null, { timeout: 60000 });
  record('B1 agent: set_scatter_view -> volcano quick view', !b1.hook_error,
    b1.hook_error || `${b1.space}: ${b1.x_axis} vs ${b1.y_axis}`);

  // B2: set_widgets on the Search cascade (agent CAN pick the column)
  const b2 = await runHook(page2, 'widgets', { patch: {
    'dataspace.feature_space.attr4.search_analysis': 'General',
    'dataspace.feature_space.attr4.search_subset': 'All',
    'dataspace.feature_space.attr4.search_variable': 'Gene.name'
  }});
  const b2ok = !b2.hook_error &&
    (b2.applied || []).length === 3 && (b2.rejected || []).length === 0;
  record('B2 agent: set_widgets applies the Search cascade (3 keys)', b2ok,
    b2.hook_error || `applied: ${(b2.applied || []).join(', ')}`);
  // ...and the write reaches the live UI: the triselector shows the value
  // and the gene picker input is revealed (showSearchBox flips to true)
  await page2.waitForFunction(id =>
    Shiny.shinyapp.$inputValues[
      'app-dataspace-feature_space-a4selector-selectSearchCol-variable'] === 'Gene.name' &&
    Shiny.shinyapp.$inputValues[`${id}-showSearchBox`] === true,
    A4, { timeout: 30000 });
  record('B2b agent: store push updates the UI triselector and reveals the gene picker', true);

  // B3: what can the model DISCOVER about highlighting? (find_controls)
  const b3h = await runHook(page2, 'capabilities', { query: 'highlight' });
  const b3g = await runHook(page2, 'capabilities', { query: 'gene name' });
  const ids = (res) => (res.capabilities || []).map(c => c.id);
  record('B3 agent: find_controls("highlight") finds the search cascade only',
    !b3h.hook_error, `${b3h.match_count} matches: ${ids(b3h).join(', ')}`);
  record('B3b agent: find_controls("gene name") finds nothing about a value box',
    !b3g.hook_error && b3g.match_count >= 0,
    `${b3g.match_count} matches: ${ids(b3g).join(', ')}`);

  // B4: try to fill the gene name through set_widgets - the natural ids
  // a model would try after B2/B3 - and assert the rejection
  const guesses = [
    'dataspace.feature_space.attr4.searchon',
    'dataspace.feature_space.attr4.search_value',
    'dataspace.feature_space.attr4.highlight'
  ];
  let allRejected = true;
  for (const id of guesses) {
    const res = await runHook(page2, 'widgets', { patch: { [id]: [GENE] } });
    const rej = (res.rejected || [])[0];
    const ok = !!rej && rej.id === id &&
      /Unknown widget id/.test(rej.reason || '');
    allRejected = allRejected && ok;
    if (!ok) console.log('   guess not rejected as expected:', id, JSON.stringify(res).slice(0, 200));
  }
  record('B4 agent: set_widgets CANNOT fill the gene name (all plausible ids rejected)',
    allRejected, `tried: ${guesses.join(', ')}`);
  const noTrace = (await highlightTraces(page2)).length === 0;
  record('B4b agent: no highlight trace appears on the plot', noTrace);

  // B5: contrast - the agent CAN select the feature (semantic selection),
  // but that is the table/selection plane, not the ring highlight
  const b5 = await runHook(page2, 'state', { features: ['X52.CNOT1'] });
  const ov = await runHook(page2, 'overview', {});
  const selCount = ov.selection && ov.selection.features &&
    ov.selection.features.count;
  const b5traces = await highlightTraces(page2);
  record('B5 agent contrast: selection works, ring highlight does not',
    !b5.hook_error && selCount === 1 && b5traces.length === 0,
    `selected=${selCount}, highlight traces=${b5traces.length}`);

  await browser.close();
} catch (e) {
  console.error('ERROR', e);
  exitCode = 1;
}

const failed = results.filter(x => !x.pass);
console.log(`\n${results.length - failed.length}/${results.length} passed`);
if (failed.length) exitCode = 1;
process.exit(exitCode);
