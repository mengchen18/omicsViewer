// WP0 Tier A — UI-effect regression for the agent state bridge.
//
// Spawns the omicsViewer Shiny app with demo.RDS preloaded (ESVObj path,
// which bypasses the file-selector dropdown), then drives the exact
// apply_agent_state / apply_scatter_view callbacks used by the ellmer
// assistant tools via the env-gated test hooks, and asserts visible UI
// outcomes. No LLM provider is required.
//
// Run standalone:   node tier_a.mjs          (from tests/e2e_agent/)
// Run via repo:     Rscript tests/test_agentUiEffects.R
//
// Writes tier_a_results.json and screenshots under artifacts/.

import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import { mkdirSync, writeFileSync, readdirSync, unlinkSync, existsSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../..');
const EXTDATA = path.join(REPO, 'inst/extdata');
const ARTIFACTS = path.join(here, 'artifacts');
mkdirSync(ARTIFACTS, { recursive: true });

const PORT = 7778;
const results = [];
const record = (name, pass, detail = '') => {
  results.push({ name, pass, detail: String(detail).slice(0, 300) });
  console.log(`${pass ? 'PASS' : 'FAIL'} | ${name}${detail ? ' | ' + detail : ''}`);
};

// ------------------------------------------------------------------- app
const rCode = `
  options(shiny.port = ${PORT}, shiny.host = '127.0.0.1')
  eset <- readRDS(file.path('${EXTDATA}', 'demo.RDS'))
  omicsViewer::omicsViewer(dir = '${EXTDATA}', ESVObj = eset)
`;
const r = spawn('Rscript', ['-e', rCode], {
  cwd: REPO,
  env: { ...process.env, OMICSVIEWER_TEST_HOOKS: 'true' },
  stdio: ['ignore', 'pipe', 'pipe']
});
let rLog = '';
r.stdout.on('data', d => { rLog += d; });
r.stderr.on('data', d => { rLog += d; });
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
  throw new Error(`app did not start on :${port}\n${rLog}`);
};

// ---------------------------------------------------------------- helpers
const TAB_ID = 'app-dataspace-eset';
const getTab = (page) => page.evaluate(
  (id) => (window.Shiny && Shiny.shinyapp && Shiny.shinyapp.$inputValues[id]) || null, TAB_ID);
const waitTab = (page, value, timeout = 30000) => page.waitForFunction(
  ({ id, v }) => Shiny.shinyapp.$inputValues[id] === v, { id: TAB_ID, v: value }, { timeout });
const HOOK = {
  box: '#app-agentTestHooks-container',
  op: '#app-agentTestHooks-op',
  payload: '#app-agentTestHooks-payload',
  run: '#app-agentTestHooks-run',
  result: '#app-agentTestHooks-result'
};

async function openSession(browser) {
  const ctx = await browser.newContext();
  const page = await ctx.newPage();
  const pageErrors = [];
  page.on('pageerror', e => pageErrors.push(String(e).slice(0, 200)));
  await page.goto(`http://127.0.0.1:${PORT}/`, { timeout: 60000 });
  // dataset is preloaded via ESVObj: wait for the contents panel + data-space tab.
  // NOTE: #app-dataspace-eset is a navbarPage binding (not a <select>), so the
  // reliable source of truth is Shiny's client-side input value map.
  await page.waitForFunction(() => {
    const tab = Shiny.shinyapp.$inputValues['app-dataspace-eset'];
    const box = document.querySelector('#app-contents');
    return !!tab && !!box && box.offsetParent !== null;
  }, null, { timeout: 90000 });
  // un-hide the test-hook panel so ordinary Playwright actions work on it
  await page.evaluate(sel => {
    document.querySelector(sel).style.display = 'block';
  }, HOOK.box);
  return { ctx, page, pageErrors };
}

async function runHook(page, op, payload) {
  const before = (await page.innerText(HOOK.result)).trim();
  // The app selectizes every select, including the hidden test-hook panel's;
  // drive it through the selectize API (with a plain-select fallback).
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
  }, { sel: HOOK.result, before, n: runCount }, { timeout: 30000 });
  const txt = (await page.innerText(HOOK.result)).trim();
  return JSON.parse(txt);
}

const axisTitles = (page) => page.evaluate(() => {
  const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
    .filter(e => e.offsetParent !== null && e._fullLayout);
  if (!els.length || !els[0]._fullLayout) return null;
  const t = (ax) => {
    const v = ax && ax.title ? (ax.title.text || ax.title) : '';
    return typeof v === 'string' ? v : '';
  };
  return [t(els[0]._fullLayout.xaxis), t(els[0]._fullLayout.yaxis)];
});

const waitAxisContains = (page, needle, timeout = 45000) =>
  page.waitForFunction((nd) => {
    const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
      .filter(e => e.offsetParent !== null && e._fullLayout);
    if (!els.length) return false;
    const t = (ax) => {
      const v = ax && ax.title ? (ax.title.text || ax.title) : '';
      return typeof v === 'string' ? v : '';
    };
    const tt = [t(els[0]._fullLayout.xaxis), t(els[0]._fullLayout.yaxis)];
    return tt.some(x => x.includes(nd));
  }, needle, { timeout });

// ------------------------------------------------------------------- run
let exitCode = 0;
try {
  await waitPort(PORT);
  const browser = await chromium.launch({
    executablePath: '/usr/bin/google-chrome',
    headless: true,
    // --disable-webgl mirrors the documented no-GPU desktop environment: the
    // app's WebGL detection then picks the SVG scatter path, which is the
    // configuration users of this machine actually experience.
    args: ['--no-sandbox', '--disable-gpu', '--disable-dev-shm-usage',
           '--disable-webgl', '--disable-webgl2']
  });

  // ---- session 1: load & baseline -----------------------------------
  const s1 = await openSession(browser);
  const p1 = s1.page;
  record('app loads with preloaded dataset', true);
  const initialTab = await getTab(p1);
  record('initial data-space tab is Feature', initialTab === 'Feature', initialTab);

  await p1.waitForFunction(() =>
    Array.from(document.querySelectorAll('.js-plotly-plot'))
      .some(e => e.offsetParent !== null && e._fullLayout), null, { timeout: 90000 });
  record('feature scatter renders a visible plotly container', true);
  record('test hooks are rendered', await p1.locator(HOOK.box).count() === 1);

  // ---- 0. volcano corner regression (load state) --------------------
  // demo.RDS defaults to the RE_vs_ME volcano: the volcano quick view must
  // be the active badge AND the corner auto-selection must have applied -
  // scorner=volcano with both top-corner rectangles drawn and the features
  // inside them selected. Shipped broken: the seed's unacknowledgeable
  // pending entry re-asserted scorner=None over the volcano auto-selection.
  {
    const t0 = Date.now();
    let st = null;
    while (Date.now() - t0 < 30000) {
      st = await p1.evaluate(() => {
        const p = Array.from(document.querySelectorAll('.js-plotly-plot'))
          .filter(e => e.offsetParent !== null && e._fullLayout)[0];
        const scorner = document.querySelector('select[id$="feature_space-a4selector-scorner"]');
        const badge = document.querySelector('[data-quick-view-id="volcano_RE_vs_ME"]');
        return {
          shapes: p && p._fullLayout ? (p._fullLayout.shapes || []).length : -1,
          scorner: scorner ? scorner.value : null,
          badgeActive: badge ? badge.classList.contains('btn-primary') : false,
          y: p && p._fullLayout ? (p._fullLayout.yaxis.title.text || '') : ''
        };
      });
      if (st.shapes === 2 && st.scorner === 'volcano') break;
      await new Promise(s => setTimeout(s, 1000));
    }
    record('load: volcano quick view active with both top-corner rects selected',
      st && st.shapes === 2 && st.scorner === 'volcano' && st.badgeActive &&
        /log\.fdr|log\.pvalue/.test(st.y || ''),
      st ? JSON.stringify(st) : 'no state');
  }

  // ---- 0b. quick-view switch renders (no multi-flash) ----------------
  // cor -> volcano must be a single render (axes + rects together);
  // volcano -> cor may clear the outgoing corner first, so at most two
  // coherent renders. The pre-fix behavior flickered the volcano rects
  // 2 -> 0 -> 2 and drew them over the outgoing correlation plot.
  const startPlotWatch = (page) => page.evaluate(() => {
    window.__plotStates = [];
    window.__lastPlot = null;
    window.__plotTimer = setInterval(() => {
      const p = Array.from(document.querySelectorAll('.js-plotly-plot'))
        .filter(e => e.offsetParent !== null && e._fullLayout)[0];
      if (!p || !p._fullLayout) return;
      const s = (p._fullLayout.xaxis.title.text || '') + '|' +
        (p._fullLayout.yaxis.title.text || '') + '|' +
        (p._fullLayout.shapes || []).length;
      if (s !== window.__lastPlot) { window.__lastPlot = s; window.__plotStates.push(s); }
    }, 40);
  });
  const readPlotStates = (page) => page.evaluate(() => {
    clearInterval(window.__plotTimer);
    return window.__plotStates;
  });
  {
    await startPlotWatch(p1);
    await p1.click('[data-quick-view-id="cor_MDR"]');
    await new Promise(s => setTimeout(s, 4500));
    const states = await readPlotStates(p1);
    const changes = states.slice(1);
    const final = states[states.length - 1] || '';
    record('volcano -> cor quick view: at most two coherent renders, no leftover rects',
      changes.length <= 2 && /Cor\|MDR/.test(final) && /\|0$/.test(final),
      states.join(' ; '));

    await startPlotWatch(p1);
    await p1.click('[data-quick-view-id="volcano_RE_vs_ME"]');
    await new Promise(s => setTimeout(s, 4500));
    const states2 = await readPlotStates(p1);
    const changes2 = states2.slice(1);
    const final2 = states2[states2.length - 1] || '';
    record('cor -> volcano quick view: exactly one render with both corner rects',
      changes2.length === 1 && /log\.(fdr|pvalue)/.test(final2) && /\|2$/.test(final2),
      states2.join(' ; '));
  }

  // ---- 1. state: tab + selection ------------------------------------
  const res1 = await runHook(p1, 'state', {
    data_space_tab: 'Sample',
    features: ['X32.TYW5', 'X52.CNOT1', 'X55.MAP1LC3B.MAP1LC3B2'],
    samples: ['X786O_NCI60', 'A498_NCI60']
  });
  record('state update applies without error', !res1.hook_error, res1.hook_error || '');
  record('hook reports Sample tab', res1.data_space_tab === 'Sample');
  record('hook reports 3 selected features', res1.feature_count === 3);
  await waitTab(p1, 'Sample');
  record('visible data-space tab switched to Sample', true);

  // ---- right-panel regression guard (analysis space must react to a
  // selection): the feature_general cascade must populate and its content
  // render. This exact scenario shipped broken once - the analysis panel
  // stayed blank because non-store cascades never fired.
  const rightPanelState = () => p1.evaluate(() => {
    const iv = Shiny.shinyapp.$inputValues;
    return {
      sub: iv['app-resultspace-feature_general-tris_feature_general-subset'],
      vr: iv['app-resultspace-feature_general-tris_feature_general-variable'],
      content: !!document.querySelector('[id*="feature_general-boxplotly"], [id*="feature_general"] .plotly')
    };
  });
  {
    const t0 = Date.now();
    let st = null;
    while (Date.now() - t0 < 30000) {
      st = await rightPanelState();
      if (st.sub && st.vr && st.vr !== '' && st.content) break;
      await new Promise(s => setTimeout(s, 1000));
    }
    record('analysis-panel cascade populates and renders after a selection',
      !!(st && st.sub && st.vr && st.vr !== '' && st.content),
      st ? JSON.stringify(st) : 'no state');
  }

  // ---- 2. scatter: custom axes ---------------------------------------
  const res2 = await runHook(p1, 'scatter', {
    space: 'feature',
    x_axis: 'PCA|All|PC1(10.5%)',
    y_axis: 'PCA|All|PC2(7.2%)'
  });
  record('custom scatter view applies without error', !res2.hook_error, res2.hook_error || '');
  await waitAxisContains(p1, 'PC1');
  record('plotly x-axis reflects requested PC1 column', true);
  const res3 = await runHook(p1, 'scatter', {
    space: 'feature',
    x_axis: 'PCA|All|PC1(10.5%)',
    y_axis: 'PCA|All|PC3(5.4%)'
  });
  record('second custom scatter view applies without error', !res3.hook_error, res3.hook_error || '');
  await waitAxisContains(p1, 'PC3');
  record('plotly y-axis updates to PC3', true);

  // ---- 3. scatter: quick view ----------------------------------------
  // WP1 progressive disclosure: the no-sections overview is compact (no
  // annotation catalog / panels / figure grammar) but still carries the
  // id+label quick-view menu and the store-backed scatter_view anchors;
  // full quick-view records need the quick_views section.
  const ovCompact = await runHook(p1, 'overview', {});
  record('WP1 overview omits full-detail sections',
    ovCompact.annotations === undefined && ovCompact.panels === undefined &&
      ovCompact.figure_grammar === undefined &&
      Array.isArray(ovCompact.available_sections) &&
      ovCompact.available_sections.includes('quick_views'),
    JSON.stringify(Object.keys(ovCompact)));
  const sv = ovCompact.scatter_view || {};
  record('WP1 overview carries store-backed scatter_view',
    !!(sv.feature && sv.feature.x && sv.feature.x.name && sv.feature.axis_mode) &&
      sv.feature.x.name.startsWith('PCA|All|PC1'),
    JSON.stringify(sv.feature || null));
  const qvMenu = (ovCompact.quick_views && ovCompact.quick_views.feature) || [];
  record('WP1 overview quick views are id+label only',
    qvMenu.length > 0 && qvMenu.every(v => v.id && v.label && v.x === undefined),
    qvMenu.slice(0, 3).map(v => v.id).join(','));
  const ovFull = await runHook(p1, 'overview', { sections: ['quick_views'] });
  const fviews = (ovFull.quick_views && ovFull.quick_views.feature) || [];
  record('WP1 quick_views section returns full records',
    fviews.length > 0 && fviews.every(v => typeof v.x === 'string' && typeof v.y === 'string'),
    fviews.length + ' records');
  record('runtime feature quick views available', fviews.length > 0,
    fviews.map(v => v.id).slice(0, 6).join(','));
  if (fviews.length > 0) {
    const qv = fviews[0];
    const res4 = await runHook(p1, 'scatter', { space: 'feature', quick_view_id: qv.id });
    record(`quick view '${qv.id}' applies without error`, !res4.hook_error, res4.hook_error || '');
    record('applied quick view axes match runtime definition',
      res4.x_axis === qv.x && res4.y_axis === qv.y, `${res4.x_axis} / ${res4.y_axis}`);
    const yVar = qv.y.split('|').pop();
    await waitAxisContains(p1, yVar);
    record('plotly axes reflect the applied quick view', true);
  }

  // ---- 3b. S2 acceptance: mode preservation + manual drift -------------
  const getMode = () => p1.evaluate(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-feature_space-axisMode']);
  const yVar = () => p1.evaluate(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-feature_space-tris_main_scatter2-variable']);
  const modeBefore = await getMode();
  record('custom-axes apply leaves the display mode untouched',
    (await getMode()) === modeBefore, `${modeBefore} -> ${await getMode()}`);

  // return to the volcano axes so log.pvalue is a valid y choice, then
  // manually drift the y variable (user behavior)
  await runHook(p1, 'scatter', { space: 'feature', quick_view_id: 'volcano_RE_vs_ME' });
  await p1.waitForFunction(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-feature_space-tris_main_scatter2-variable'] === 'log.fdr',
    null, { timeout: 30000 });
  // wait for the variable selectize to actually carry its choices before
  // drifting (the selected value lands via one message, the options via
  // another; a fast value-wait can win the race on a slow machine)
  await p1.waitForFunction(() => {
    const el = document.querySelector('#app-dataspace-feature_space-tris_main_scatter2-variable');
    return el && el.selectize && Object.keys(el.selectize.options).length > 0;
  }, null, { timeout: 30000 });
  await p1.evaluate(() => {
    const el = document.querySelector('#app-dataspace-feature_space-tris_main_scatter2-variable');
    if (el && el.selectize) el.selectize.setValue('log.pvalue');
  });
  await p1.waitForFunction(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-feature_space-tris_main_scatter2-variable'] === 'log.pvalue',
    null, { timeout: 15000 });
  record('manual variable edit sticks in the widget', (await yVar()) === 'log.pvalue');
  await p1.waitForTimeout(1500);

  // agent apply targeting the pre-edit value (the drift scenario): the
  // quick view's axes must correct the manually drifted widget
  const resDrift = await runHook(p1, 'scatter', { space: 'feature', quick_view_id: 'volcano_RE_vs_ME' });
  record('drift-correcting apply succeeds', !resDrift.hook_error, resDrift.hook_error || '');
  let corrected = false;
  for (let attempt = 0; attempt < 3 && !corrected; attempt++) {
    // settle first: the UI->store sync of the manual edit must land before
    // the apply diffs against it, otherwise the apply is a legitimate no-op
    await p1.waitForTimeout(2500);
    corrected = await p1.evaluate(() =>
      Shiny.shinyapp.$inputValues['app-dataspace-feature_space-tris_main_scatter2-variable'] === 'log.fdr');
    if (!corrected && attempt < 2)
      await runHook(p1, 'scatter', { space: 'feature', quick_view_id: 'volcano_RE_vs_ME' });
  }
  record('apply corrects the manually drifted axis to the quick view', corrected,
    corrected ? '' : 'y stayed ' + (await yVar()));
  record('drift-correcting apply still leaves the mode untouched',
    (await getMode()) === modeBefore);

  // stress: 6 rapid interleaved applies; final state must equal the last
  const targets = [
    ['PCA|All|PC1(10.5%)', 'PCA|All|PC2(7.2%)'],
    ['ttest|RE_vs_ME|mean.diff', 'ttest|RE_vs_ME|log.fdr'],
    ['PCA|All|PC3(5.4%)', 'PCA|All|PC4(4.8%)'],
  ];
  let stressOk = true, lastErr = '';
  for (const [xx, yy] of targets) {
    const rr = await runHook(p1, 'scatter', { space: 'feature', x_axis: xx, y_axis: yy });
    if (rr.hook_error) { stressOk = false; lastErr = rr.hook_error; break; }
  }
  await p1.waitForFunction(() => {
    const iv = Shiny.shinyapp.$inputValues;
    return iv['app-dataspace-feature_space-tris_main_scatter1-variable'] === 'PC3(5.4%)' &&
           iv['app-dataspace-feature_space-tris_main_scatter2-variable'] === 'PC4(4.8%)';
  }, null, { timeout: 45000 }).catch(() => { stressOk = false; lastErr = 'final axes never landed'; });
  record('rapid successive applies converge on the last request', stressOk, lastErr);

  // ---- 4. validation negatives ---------------------------------------
  const n1 = await runHook(p1, 'state', { data_space_tab: 'Samples-typo' });
  record('invalid tab rejected with descriptive error',
    !!n1.hook_error && /Unknown data-space tab/.test(n1.hook_error), n1.hook_error || '');
  const n2 = await runHook(p1, 'scatter', { space: 'feature', x_axis: 'p value', y_axis: 'log.fdr' });
  record('invalid custom axis rejected with descriptive error',
    !!n2.hook_error && /X-axis annotation/.test(n2.hook_error), n2.hook_error || '');
  const n3 = await runHook(p1, 'scatter', { space: 'feature', quick_view_id: 'nope' });
  record('unknown quick view rejected with descriptive error',
    !!n3.hook_error && /quick view/i.test(n3.hook_error), n3.hook_error || '');
  record('failed updates leave visible tab unchanged',
    (await getTab(p1)) === 'Feature');
  const outErrs = await p1.evaluate(() =>
    document.querySelectorAll('.shiny-output-error').length);
  record('no shiny output errors after negatives', outErrs === 0, String(outErrs));

  await p1.screenshot({ path: path.join(ARTIFACTS, 'tier_a_session1_final.png'), fullPage: false });

  // ---- 4b. S3 acceptance: generic widget tier -------------------------
  // Drive set_widgets through the test hooks (the exact store path the
  // ellmer tool uses) against a widget family never exposed to the agent
  // before: heatmap parameters. Verifies visible UI effect + per-key
  // self-correction feedback.
  await runHook(p1, 'state', { data_space_tab: 'Heatmap' });
  await waitTab(p1, 'Heatmap');
  const hmColor = () => p1.evaluate(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-heatmapViewer-heatmapColors']);
  const hmMargin = () => p1.evaluate(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-heatmapViewer-marginBottom']);
  const w1 = await runHook(p1, 'widgets', {
    patch: {
      'dataspace.expr_heatmap.heatmap_colors': 'RdGy',
      'dataspace.expr_heatmap.margin_bottom': 9
    }
  });
  record('generic widget apply succeeds', !w1.hook_error, w1.hook_error || '');
  record('receipt lists both keys applied',
    (w1.applied || []).includes('dataspace.expr_heatmap.heatmap_colors') &&
    (w1.applied || []).includes('dataspace.expr_heatmap.margin_bottom'),
    JSON.stringify(w1.applied || []));
  await p1.waitForFunction(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-heatmapViewer-heatmapColors'] === 'RdGy',
    null, { timeout: 30000 });
  record('heatmap color select reflects agent-set palette', (await hmColor()) === 'RdGy');
  await p1.waitForFunction(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-heatmapViewer-marginBottom'] === 9,
    null, { timeout: 30000 });
  record('heatmap margin slider reflects agent-set value', (await hmMargin()) === 9);
  const w2 = await runHook(p1, 'widgets', {
    patch: { 'dataspace.expr_heatmap.heatmap_colors': 'Spectral' }
  });
  record('invalid palette rejected per key with suggestions',
    !w2.hook_error && (w2.rejected || []).length === 1 &&
    /BrBG|RdGy|RdBu/.test(w2.rejected[0].reason),
    JSON.stringify(w2.rejected || []));
  record('rejected key leaves current palette unchanged', (await hmColor()) === 'RdGy');
  const w3 = await runHook(p1, 'widgets', {
    patch: { 'dataspace.expr_heatmap.colour_scheme': 'RdGy' }
  });
  record('unknown widget id rejected with a suggestion',
    !w3.hook_error && (w3.rejected || []).length === 1 &&
    /Unknown widget id/.test(w3.rejected[0].reason) &&
    /dataspace\.expr_heatmap\./.test(w3.rejected[0].reason),
    JSON.stringify(w3.rejected || []));

  // ---- 4c. S4 acceptance: heatmap completion ---------------------------
  // Sorting, clustering, and annotation widgets (the multi_select kind's
  // first live consumers) driven through the same store path.
  const hmInput = (id) => p1.evaluate(
    (k) => Shiny.shinyapp.$inputValues[k],
    `app-dataspace-heatmapViewer-${id}`);
  const w4 = await runHook(p1, 'widgets', {
    patch: {
      'dataspace.expr_heatmap.col_sort_by': 'none',
      'dataspace.expr_heatmap.row_sort_by': 'hierarchical cluster',
      'dataspace.expr_heatmap.cluster_row_dist': 'Spearman correlation',
      'dataspace.expr_heatmap.cluster_row_link': 'complete',
      'dataspace.expr_heatmap.annot_col': ['General|All|Cell.line']
    }
  });
  record('S4 sorting/clustering/annotation apply succeeds',
    !w4.hook_error, w4.hook_error || '');
  record('receipt lists all five S4 keys applied',
    ['dataspace.expr_heatmap.col_sort_by', 'dataspace.expr_heatmap.row_sort_by',
     'dataspace.expr_heatmap.cluster_row_dist', 'dataspace.expr_heatmap.cluster_row_link',
     'dataspace.expr_heatmap.annot_col']
      .every(k => (w4.applied || []).includes(k)),
    JSON.stringify(w4.applied || []));
  await p1.waitForFunction(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-heatmapViewer-colSortBy'] === 'none' &&
    Shiny.shinyapp.$inputValues['app-dataspace-heatmapViewer-rowSortBy'] === 'hierarchical cluster',
    null, { timeout: 30000 });
  record('sorting selects reflect agent-set values',
    (await hmInput('colSortBy')) === 'none' &&
    (await hmInput('rowSortBy')) === 'hierarchical cluster');
  record('clustering selects reflect agent-set values',
    (await hmInput('clusterRowDist')) === 'Spearman correlation' &&
    (await hmInput('clusterRowLink')) === 'complete');
  record('annotation multi-select reflects agent-set column',
    JSON.stringify(await hmInput('annotCol')) ===
      JSON.stringify(['General|All|Cell.line']),
    JSON.stringify(await hmInput('annotCol')));
  const w5 = await runHook(p1, 'widgets', {
    patch: { 'dataspace.expr_heatmap.annot_col': ['Cell.line'] }
  });
  record('unknown annotation column rejected with suggestions',
    !w5.hook_error && (w5.rejected || []).length === 1 &&
    /Cell\.line|Unknown value/.test(w5.rejected[0].reason),
    JSON.stringify(w5.rejected || []));
  record('rejected annotation leaves selection unchanged',
    JSON.stringify(await hmInput('annotCol')) ===
      JSON.stringify(['General|All|Cell.line']));
  const w6 = await runHook(p1, 'widgets', {
    patch: { 'dataspace.expr_heatmap.annot_col': ['General|All|MDR'] }
  });
  await p1.waitForFunction(() =>
    JSON.stringify(Shiny.shinyapp.$inputValues['app-dataspace-heatmapViewer-annotCol']) ===
      JSON.stringify(['General|All|MDR']), null, { timeout: 30000 });
  record('multi-select write replaces the selection wholesale',
    JSON.stringify(await hmInput('annotCol')) === JSON.stringify(['General|All|MDR']),
    JSON.stringify(await hmInput('annotCol')));

  // ---- 4d. S4 acceptance: data tables ---------------------------------
  // The Feature table registers its user-editable surface: the
  // multi-row-selection switch and the shown-columns set.
  await runHook(p1, 'state', { data_space_tab: 'Feature table' });
  await waitTab(p1, 'Feature table');
  const tabSwitch = () => p1.evaluate(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-tab_feature-multisel']);
  const tableHeaders = () => p1.evaluate(() =>
    Array.from(document.querySelectorAll('#app-dataspace-tab_feature-table thead th'))
      .filter(el => !el.querySelector('input'))
      .map(el => el.innerText.trim()).filter(t => t));
  const w7 = await runHook(p1, 'widgets', {
    patch: {
      'dataspace.tab_feature.multi_selection': true,
      'dataspace.tab_feature.columns': ['General|All|Gene.name', 'mean|Origin|RE']
    }
  });
  record('table widget apply succeeds', !w7.hook_error, w7.hook_error || '');
  record('receipt lists both table keys applied',
    (w7.applied || []).includes('dataspace.tab_feature.multi_selection') &&
    (w7.applied || []).includes('dataspace.tab_feature.columns'),
    JSON.stringify(w7.applied || []));
  await p1.waitForFunction(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-tab_feature-multisel'] === true,
    null, { timeout: 30000 });
  record('multi-selection switch reflects agent-set value',
    (await tabSwitch()) === true);
  await p1.waitForFunction(() => {
    const hs = Array.from(document.querySelectorAll('#app-dataspace-tab_feature-table thead th'))
      .filter(el => !el.querySelector('input'))
      .map(el => el.innerText.trim());
    return hs.includes('mean|Origin|RE') && !hs.includes('General|All|Protein.ID');
  }, null, { timeout: 30000 });
  record('feature table headers reflect the agent-set columns',
    (await tableHeaders()).includes('mean|Origin|RE'));
  const w8 = await runHook(p1, 'widgets', {
    patch: { 'dataspace.tab_feature.columns': [] }
  });
  record('emptying table columns is rejected with the min rule',
    !w8.hook_error && (w8.rejected || []).length === 1 &&
    /at least 1/.test(w8.rejected[0].reason),
    JSON.stringify(w8.rejected || []));
  record('rejected column set leaves headers unchanged',
    (await tableHeaders()).includes('mean|Origin|RE'));
  const w9 = await runHook(p1, 'widgets', {
    patch: { 'dataspace.tab_pheno.columns': ['General|All|Cell.line'] }
  });
  record('sample table columns are settable through the same tier',
    !w9.hook_error && [w9.applied].flat().includes('dataspace.tab_pheno.columns'),
    JSON.stringify({ applied: w9.applied || [], rejected: w9.rejected || [] }));
  // return to the Feature tab so later isolation assertions keep their
  // baseline
  await runHook(p1, 'state', { data_space_tab: 'Feature' });
  await waitTab(p1, 'Feature');

  // ---- 4e. S4 acceptance: result space --------------------------------
  // The result-space modules join the control plane: the analyst navbar,
  // sample_general's link-variable cascade (survival / contingency
  // switch), feature_general's cascade + plot-type radio, and the shared
  // attr4 panel. Driven through the same widgets hook as every S4 tier.
  const rsInput = (id) => p1.evaluate(
    (k) => Shiny.shinyapp.$inputValues[k], `app-resultspace-${id}`);
  const waitRsInput = (id, v, timeout = 45000) => p1.waitForFunction(
    ({ k, x }) => Shiny.shinyapp.$inputValues[k] === x,
    { k: `app-resultspace-${id}`, x: v }, { timeout });
  const w10 = await runHook(p1, 'widgets', {
    patch: { 'resultspace.analyst_tab': 'Sample' }
  });
  record('result-space tab apply succeeds', !w10.hook_error, w10.hook_error || '');
  record('receipt lists the analyst tab applied',
    [w10.applied].flat().includes('resultspace.analyst_tab'),
    JSON.stringify(w10.applied || []));
  await waitRsInput('analyst', 'Sample');
  record('analyst navbar reflects agent-set tab',
    (await rsInput('analyst')) === 'Sample');

  // sample_general: the Surv category switches to the Kaplan-Meier view
  const w11 = await runHook(p1, 'widgets', {
    patch: {
      'resultspace.sample_general.xax_analysis': 'Surv',
      'resultspace.sample_general.xax_subset': 'all',
      'resultspace.sample_general.xax_variable': 'OS'
    }
  });
  record('sample link-variable cascade applies', !w11.hook_error, w11.hook_error || '');
  record('receipt lists all three cascade keys applied',
    ['resultspace.sample_general.xax_analysis',
     'resultspace.sample_general.xax_subset',
     'resultspace.sample_general.xax_variable']
      .every(k => (w11.applied || []).includes(k)),
    JSON.stringify(w11.applied || []));
  await waitRsInput('sample_general-tris_sample_general-variable', 'OS');
  record('sample cascade select reflects agent-set Surv variable',
    (await rsInput('sample_general-tris_sample_general-variable')) === 'OS');
  await p1.waitForFunction(() =>
    !!document.querySelector('[id*="sample_general_surv"]'), null, { timeout: 45000 });
  record('survival view renders for the Surv link variable', true);

  // a categorical variable switches the same panel to a contingency table
  // (MDR is numeric in the demo data and keeps the beeswarm view)
  const w12 = await runHook(p1, 'widgets', {
    patch: {
      'resultspace.sample_general.xax_analysis': 'General',
      'resultspace.sample_general.xax_subset': 'All',
      'resultspace.sample_general.xax_variable': 'Origin'
    }
  });
  record('categorical sample cascade applies', !w12.hook_error, w12.hook_error || '');
  await waitRsInput('sample_general-tris_sample_general-variable', 'Origin');
  await p1.waitForFunction(() =>
    !!document.querySelector('[id*="sample_general_contab"]'), null, { timeout: 45000 });
  record('contingency view renders for a categorical link variable', true);

  // feature_general: cascade + plot type. Features are selected through
  // the real UI (Feature table rows) because the analysis panel consumes
  // the data-space selection, which only table/scatter interactions feed
  // durably (a pre-existing state-bridge gap: apply_agent_state feature
  // patches are transient and get overwritten by the restore roundtrip).
  const w13 = await runHook(p1, 'widgets', {
    patch: { 'resultspace.analyst_tab': 'Feature' }
  });
  await waitRsInput('analyst', 'Feature');
  await runHook(p1, 'state', { data_space_tab: 'Feature table' });
  await waitTab(p1, 'Feature table');
  const selRows = async () => p1.evaluate(() =>
    (Shiny.shinyapp.$inputValues['app-dataspace-tab_feature-table_rows_selected'] || []).length);
  {
    // multiple-row selection needs the switch on; tab re-renders reset it,
    // so toggle it through the DOM exactly like a user would. Read the
    // rendered bootstrap-switch state (the client input map can be stale
    // right after a tab switch).
    const sw = '.bootstrap-switch-id-app-dataspace-tab_feature-multisel';
    await p1.waitForSelector(sw, { timeout: 30000 });
    await p1.waitForFunction((s) => {
      const el = document.querySelector(s);
      return el && (el.classList.contains('bootstrap-switch-on') ||
                    el.classList.contains('bootstrap-switch-off'));
    }, sw, { timeout: 15000 });
    if (!(await p1.evaluate(s =>
      document.querySelector(s).classList.contains('bootstrap-switch-on'), sw))) {
      await p1.click(sw);
      await p1.waitForFunction(s =>
        document.querySelector(s)?.classList.contains('bootstrap-switch-on'),
        sw, { timeout: 15000 });
    }
    const rows = p1.locator('#app-dataspace-tab_feature-table tbody tr');
    await rows.first().waitFor({ timeout: 30000 });
    // plain clicks accumulate in the DT os-select style (ctrl-click
    // replaces); the DT redraw on tab return can swallow a click, so
    // self-correct against the live selection until the wanted rows hold
    const wanted = [0, 2, 4];
    const getSel = () => p1.evaluate(() =>
      Shiny.shinyapp.$inputValues['app-dataspace-tab_feature-table_rows_selected'] || []);
    const deadline = Date.now() + 45000;
    while (Date.now() < deadline) {
      const sel = await getSel();
      if (sel.length === wanted.length && wanted.every(i => sel.includes(i + 1)))
        break;
      const missing = wanted.find(i => !sel.includes(i + 1));
      const extra = sel.filter(x => !wanted.includes(x - 1));
      const target = missing !== undefined ? missing : extra[0] - 1;
      await rows.nth(target).click().catch(() => {});
      await p1.waitForTimeout(600);
    }
    if ((await selRows()) !== 3)
      throw new Error('feature-table row selection did not converge');
  }
  record('three features selected through the feature table', (await selRows()) === 3);
  await runHook(p1, 'state', { data_space_tab: 'Feature' });
  await waitTab(p1, 'Feature');
  const w14 = await runHook(p1, 'widgets', {
    patch: {
      'resultspace.feature_general.xax_analysis': 'General',
      'resultspace.feature_general.xax_subset': 'All',
      'resultspace.feature_general.xax_variable': 'TP53.Status',
      'resultspace.feature_general.plot_type': 'Curve'
    }
  });
  record('feature link-variable + plot type apply', !w14.hook_error, w14.hook_error || '');
  await waitRsInput('feature_general-tris_feature_general-variable', 'TP53.Status');
  await waitRsInput('feature_general-internal_radio', 'Curve');
  record('feature plot-type radio reflects agent-set value',
    (await rsInput('feature_general-internal_radio')) === 'Curve');
  await p1.waitForFunction(() =>
    !!document.querySelector('[id*="feature_general_roc_pr"]'), null, { timeout: 45000 });
  record('ROC/PR view renders for the Curve plot type', true);

  // attr4: the shared figure-attribute panel (nested namespace)
  const w15 = await runHook(p1, 'widgets', {
    patch: {
      'resultspace.feature_general.attr4.color_analysis': 'General',
      'resultspace.feature_general.attr4.color_subset': 'All',
      'resultspace.feature_general.attr4.color_variable': 'Origin'
    }
  });
  record('attr4 color cascade applies', !w15.hook_error, w15.hook_error || '');
  await waitRsInput('feature_general-a4_gf-selectColorUI-variable', 'Origin');
  record('attr4 color selector reflects agent-set variable',
    (await rsInput('feature_general-a4_gf-selectColorUI-variable')) === 'Origin');

  // negatives: per-key rejections with suggestions, current state intact
  const w16 = await runHook(p1, 'widgets', {
    patch: { 'resultspace.feature_general.plot_type': 'Histogram' }
  });
  record('invalid result-space plot type rejected with allowed values',
    !w16.hook_error && (w16.rejected || []).length === 1 &&
    /Allowed: Bees, Curve/.test(w16.rejected[0].reason),
    JSON.stringify(w16.rejected || []));
  record('rejected plot type leaves the radio unchanged',
    (await rsInput('feature_general-internal_radio')) === 'Curve');
  const w17 = await runHook(p1, 'widgets', {
    patch: { 'resultspace.analyst_tab': 'Response' }
  });
  record('dataset-absent analysis tab rejected with allowed tabs',
    !w17.hook_error && (w17.rejected || []).length === 1 &&
    /StringDB|Sample/.test(w17.rejected[0].reason),
    JSON.stringify(w17.rejected || []));
  record('rejected tab leaves the navbar on Feature',
    (await rsInput('analyst')) === 'Feature');

  // ---- 4f. S4 acceptance: enrichment tables (fGSEA) -------------------
  // The ranking cascade plus the dataTableDownload row selection (the
  // settled S4 decision: row selection registers where it drives a
  // downstream view - here the leading-edge bar plot). The same machinery
  // backs the ORA overlap table (unit-covered in test_agentWidgets.R).
  const w18 = await runHook(p1, 'widgets', {
    patch: { 'resultspace.analyst_tab': 'fGSEA' }
  });
  record('fGSEA tab apply succeeds', !w18.hook_error, w18.hook_error || '');
  await waitRsInput('analyst', 'fGSEA');
  const w19 = await runHook(p1, 'widgets', {
    patch: {
      'resultspace.fgsea.xax_analysis': 'PCA',
      'resultspace.fgsea.xax_subset': 'All',
      'resultspace.fgsea.xax_variable': 'PC1(10.5%)'
    }
  });
  record('fGSEA ranking cascade applies', !w19.hook_error, w19.hook_error || '');
  // re-apply the same patch: every key must now be a no-op (unchanged),
  // proving the whole triple holds the requested values even when a key
  // was already equal before the first apply (diff-only receipts)
  const w19b = await runHook(p1, 'widgets', {
    patch: {
      'resultspace.fgsea.xax_analysis': 'PCA',
      'resultspace.fgsea.xax_subset': 'All',
      'resultspace.fgsea.xax_variable': 'PC1(10.5%)'
    }
  });
  record('all three ranking keys hold the requested values',
    ['resultspace.fgsea.xax_analysis',
     'resultspace.fgsea.xax_subset',
     'resultspace.fgsea.xax_variable']
      .every(k => (w19.applied || []).includes(k) ||
                  (w19b.unchanged || []).includes(k)),
    JSON.stringify({ applied: w19.applied, unchanged: w19b.unchanged }));
  await waitRsInput('fgsea-tris_fgsea-variable', 'PC1(10.5%)');
  record('fGSEA ranking select reflects agent-set variable',
    (await rsInput('fgsea-tris_fgsea-variable')) === 'PC1(10.5%)');
  const fgRows = p1.locator('#app-resultspace-fgsea-stab-table tbody tr');
  await fgRows.nth(1).waitFor({ timeout: 60000 });
  record('fGSEA results table renders for the ranking variable', true);
  const fgPathway = async (i) => {
    await fgRows.nth(i).waitFor({ timeout: 30000 });
    return (await fgRows.nth(i).locator('td').first().innerText()).trim();
  };
  const path1 = await fgPathway(0);
  const path2 = await fgPathway(1);
  // user row click lands in the store (verified via the no-diff receipt)
  await fgRows.nth(1).click();
  await p1.waitForFunction(() =>
    (Shiny.shinyapp.$inputValues['app-resultspace-fgsea-stab-table_rows_selected'] || [])[0] === 2,
    null, { timeout: 30000 });
  const w20 = await runHook(p1, 'widgets', {
    patch: { 'resultspace.fgsea.selected_row': path2 }
  });
  record('user row click reached the store (same-value patch is a no-op)',
    !w20.hook_error && (w20.applied || []).length === 0 &&
      (w20.unchanged || []).includes('resultspace.fgsea.selected_row'),
    JSON.stringify({ applied: w20.applied, unchanged: w20.unchanged }));
  // agent push of another pathway re-renders the table with that row
  // preselected; DT reports it and the write acknowledges
  const w21 = await runHook(p1, 'widgets', {
    patch: { 'resultspace.fgsea.selected_row': path1 }
  });
  record('agent pathway-row push applies',
    !w21.hook_error && (w21.applied || []).includes('resultspace.fgsea.selected_row'),
    JSON.stringify(w21.applied || []));
  await p1.waitForFunction(() =>
    (Shiny.shinyapp.$inputValues['app-resultspace-fgsea-stab-table_rows_selected'] || [])[0] === 1,
    null, { timeout: 30000 });
  record('table re-renders with the agent-selected pathway row', true);
  const w23 = await runHook(p1, 'widgets', {
    patch: { 'resultspace.fgsea.selected_row': 'NOT_A_PATHWAY' }
  });
  record('unknown pathway rejected with allowed values',
    !w23.hook_error && (w23.rejected || []).length === 1 &&
      /Allowed:/.test(w23.rejected[0].reason),
    JSON.stringify(w23.rejected || []));

  // ---- 4g. S4 completion: snapshot round-trip equals store state ------
  // The real .ESS save/restore flow (snapshot modal) must reproduce the
  // saved canonical store state exactly: patch widgets, read the store,
  // save, drift, restore through the modal table, read again, deep-compare.
  const canon = (v) => Array.isArray(v)
    ? v.map(canon).sort()
    : (v !== null && typeof v === 'object')
      ? Object.fromEntries(Object.entries(v).map(([k, x]) => [k, canon(x)])
          .sort((a, b) => (a[0] < b[0] ? -1 : 1)))
      : [v];
  const storeSnapshot = async () => {
    const r = await runHook(p1, 'store', {});
    if (r.hook_error) throw new Error(r.hook_error);
    return r.values;
  };
  const w24 = await runHook(p1, 'widgets', {
    patch: {
      'dataspace.expr_heatmap.heatmap_colors': 'RdGy',
      'dataspace.expr_heatmap.margin_bottom': 7,
      'resultspace.feature_general.plot_type': 'Curve'
    }
  });
  record('round-trip widget writes apply', !w24.hook_error,
    w24.hook_error || '');
  const rtSaved = await storeSnapshot();
  // save through the real modal
  await p1.click('[data-testid="app-snapshot-button"]');
  await p1.fill('#app-snapshot_name', 'rt1');
  await p1.click('#app-snapshot_save');
  await p1.waitForFunction(() =>
    (document.querySelector('#app-tab_saveSS tbody tr') || null) !== null,
    null, { timeout: 30000 });
  record('snapshot modal lists the saved .ESS', true);
  // drift the widgets after saving
  const w25 = await runHook(p1, 'widgets', {
    patch: {
      'dataspace.expr_heatmap.heatmap_colors': 'PiYG',
      'dataspace.expr_heatmap.margin_bottom': 3,
      'resultspace.feature_general.plot_type': 'Bees'
    }
  });
  record('post-save drift applies', !w25.hook_error, w25.hook_error || '');
  const rtMid = await storeSnapshot();
  record('drift changed the live store state',
    rtMid['dataspace.expr_heatmap.heatmap_colors'] === 'PiYG');
  // restore by clicking the saved row in the modal's snapshot table
  await p1.click('[data-testid="app-snapshot-button"]');
  await p1.locator('#app-tab_saveSS tbody tr').first().click();
  await p1.waitForFunction((want) =>
    Shiny.shinyapp.$inputValues['app-dataspace-heatmapViewer-heatmapColors'] === want,
    'RdGy', { timeout: 45000 });
  record('restore reverts the drifted heatmap palette in the widget', true);
  const rtAfter = await storeSnapshot();
  const keys = Object.keys(rtSaved);
  const diffs = keys.filter(k =>
    JSON.stringify(canon(rtSaved[k])) !== JSON.stringify(canon(rtAfter[k])));
  record('restore reproduces the saved store state exactly',
    diffs.length === 0,
    diffs.slice(0, 5).map(k => `${k}: ${JSON.stringify(rtSaved[k])} -> ${JSON.stringify(rtAfter[k])}`).join('; '));
  // clean the test snapshot out of the repo extdata directory
  try {
    for (const f of readdirSync(EXTDATA))
      if (/rt1\.ESS$/i.test(f)) unlinkSync(path.join(EXTDATA, f));
  } catch {}

  // ---- 5. cross-session isolation -------------------------------------
  const s2 = await openSession(browser);
  const p2 = s2.page;
  const res5 = await runHook(p2, 'state', { data_space_tab: 'Heatmap' });
  record('second session applies its own state update', !res5.hook_error, res5.hook_error || '');
  await waitTab(p2, 'Heatmap');
  record('second session visible tab switched to Heatmap', true);
  record('first session tab unchanged by second session',
    (await getTab(p1)) === 'Feature');
  record('session 1 has no uncaught page errors', s1.pageErrors.length === 0, s1.pageErrors[0] || '');
  record('session 2 has no uncaught page errors', s2.pageErrors.length === 0, s2.pageErrors[0] || '');
  await s2.ctx.close();

  await browser.close();
  try { r.kill('SIGTERM'); } catch {}
} catch (e) {
  record('harness completed without fatal errors', false, e.message);
  console.error('--- R process log (tail) ---');
  console.error(rLog.split('\n').slice(-25).join('\n'));
  exitCode = 1;
}

const passed = results.filter(x => x.pass).length;
console.log(`\nTier A summary: ${passed}/${results.length} passed`);
writeFileSync(path.join(here, 'tier_a_results.json'), JSON.stringify({
  finished_at: new Date().toISOString(),
  passed, total: results.length, results
}, null, 2));
process.exit(exitCode === 0 && passed === results.length ? 0 : 1);
