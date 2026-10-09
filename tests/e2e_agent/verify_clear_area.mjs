// Browser-level repro: does "Clear figure selection" erase the drawn
// box/lasso (x/y area) selection on the scatter?
//
// Steps: load volcano -> manual box select of topmost point -> click
// "Clear figure selection" -> inspect bus count, plotly DOM selection
// artifacts, render (DOM mutation) count, plotly event inputs.
//
// Run:  node verify_clear_area.mjs     (from tests/e2e_agent/)

import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import { mkdirSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../..');
const EXTDATA = path.join(REPO, 'inst/extdata');
const ARTIFACTS = path.join(here, 'artifacts');
mkdirSync(ARTIFACTS, { recursive: true });

const PORT = 7781;
const results = [];
const record = (name, pass, detail = '') => {
  results.push({ name, pass, detail: String(detail).slice(0, 400) });
  console.log(`${pass ? 'PASS' : 'FAIL'} | ${name}${detail ? ' | ' + detail : ''}`);
};

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
  return false;
};

const HOOK = {
  op: '#app-agentTestHooks-op', payload: '#app-agentTestHooks-payload',
  run: '#app-agentTestHooks-run', result: '#app-agentTestHooks-result'
};
const setHookSelect = (page, sel, value) => page.evaluate(({ sel, value }) => {
  const el = document.querySelector(sel);
  if (el.selectize) el.selectize.setValue(value);
  else {
    el.value = value;
    el.dispatchEvent(new Event('change', { bubbles: true }));
  }
}, { sel, value });

async function selCount(page) {
  const before = (await page.innerText(HOOK.result)).trim();
  await setHookSelect(page, HOOK.op, 'overview');
  await page.fill(HOOK.payload, JSON.stringify({ sections: [] }));
  const n = (parseInt((await page.evaluate(
    s => document.querySelector(s).dataset.runs || '0', HOOK.run)), 10) + 1);
  await page.evaluate(({ s, n }) => {
    document.querySelector(s).dataset.runs = String(n);
  }, { s: HOOK.run, n });
  await page.click(HOOK.run);
  await page.waitForFunction(({ s, b, n }) => {
    const t = (document.querySelector(s)?.innerText || '').trim();
    return t !== b && t.includes(`"hook_run": ${n}`);
  }, { s: HOOK.result, b: before, n }, { timeout: 30000 });
  const txt = (await page.innerText(HOOK.result)).trim();
  const j = JSON.parse(txt);
  if (j.hook_error) throw new Error('overview hook: ' + j.hook_error);
  const out = j.result || j;
  const sel = (out.state || out).selection || out;
  return {
    n: (sel.features && sel.features.count) != null ? sel.features.count : null,
    raw: sel
  };
}

// census of plotly-side selection artifacts on the VISIBLE feature scatter
const selectionCensus = (page) => page.evaluate(() => {
  const gd = Array.from(document.querySelectorAll('.js-plotly-plot'))
    .filter(e => e.offsetParent !== null && e.data)[0];
  if (!gd) return { error: 'no visible plot' };
  const outlines = gd.querySelectorAll('.select-outline').length;
  const selLayer = gd.querySelectorAll('.selectionlayer *').length;
  const selLayerHtml = selLayer ?
    (gd.querySelector('.selectionlayer').innerHTML || '').slice(0, 120) : '';
  const op = gd.data[0].marker.opacity;
  let solid = null, dim = null;
  if (Array.isArray(op)) {
    solid = 0; dim = 0;
    for (const v of op) { if (v > 0.5) solid++; else dim++; }
  }
  // shiny event inputs currently registered for plotly selections
  const inputs = {};
  for (const [k, v] of Object.entries(Shiny.shinyapp.$inputValues)) {
    if (/plotly_selected/.test(k) && v !== null && v !== undefined &&
        !(Array.isArray(v) && v.length === 0)) {
      const s = JSON.stringify(v);
      if (s && s !== '[]') inputs[k] = s.slice(0, 160);
    }
  }
  return {
    outlines, selLayer, selLayerHtml,
    solid, dim, total: Array.isArray(op) ? op.length : null,
    selInputs: inputs,
    renderId: gd.dataset.renderId || null,
    epoch: gd.__renderEpoch || null
  };
});

const startPlotWatch = (page) => page.evaluate(() => {
  const visiblePlot = () => Array.from(document.querySelectorAll('.js-plotly-plot'))
    .filter(e => e.offsetParent !== null && e._fullLayout)[0];
  window.__plotRenders = 0;
  window.__plotWatchStopped = false;
  let pending = false;
  const settle = () => {
    pending = false;
    if (window.__plotWatchStopped) return;
    const p = visiblePlot();
    if (!p) return;
    window.__plotRenders++;
  };
  const mo = new MutationObserver(() => {
    if (window.__plotWatchStopped || pending) return;
    pending = true;
    requestAnimationFrame(settle);
  });
  const el0 = visiblePlot();
  if (el0) mo.observe(el0, { childList: true, subtree: true, attributes: true });
  window.__plotDisconnect = () => {
    window.__plotWatchStopped = true;
    mo.disconnect();
  };
});
const readPlotRenders = (page) => page.evaluate(() => {
  if (window.__plotDisconnect) window.__plotDisconnect();
  return window.__plotRenders;
});

const waitVolcano = (page, timeout = 90000) => page.waitForFunction(() => {
  const gd = Array.from(document.querySelectorAll('.js-plotly-plot'))
    .filter(e => e.offsetParent !== null && e._fullLayout)[0];
  if (!gd) return false;
  const t = (ax) => (ax && ax.title ? (ax.title.text || '') : '');
  return /mean\.diff/.test(t(gd._fullLayout.xaxis)) &&
         /log\.(fdr|pvalue)/.test(t(gd._fullLayout.yaxis));
}, null, { timeout });

const settle = (ms = 5000) => new Promise(s => setTimeout(s, ms));

let exitCode = 0;
try {
  await waitPort(PORT);
  const browser = await chromium.launch({
    executablePath: '/usr/bin/google-chrome',
    headless: true,
    args: ['--no-sandbox', '--disable-gpu', '--disable-dev-shm-usage',
           '--disable-webgl', '--disable-webgl2']
  });
  const ctx = await browser.newContext({ viewport: { width: 1800, height: 1200 } });
  const page = await ctx.newPage();
  await page.goto(`http://127.0.0.1:${PORT}/`, { timeout: 60000 });
  await page.waitForFunction(() => {
    const tab = Shiny.shinyapp.$inputValues['app-dataspace-eset'];
    const box = document.querySelector('#app-contents');
    return !!tab && box && box.offsetParent !== null;
  }, null, { timeout: 90000 });
  await page.evaluate(s => { document.querySelector(s).style.display = 'block'; },
    '#app-agentTestHooks-container');
  await waitVolcano(page);
  await settle();

  // -- 1. manual box selection of the topmost point ----------------------
  await page.click('.modebar-btn[data-title="Box Select"]');
  const target = await page.evaluate(() => {
    const gd = Array.from(document.querySelectorAll('.js-plotly-plot'))
      .filter(e => e.offsetParent !== null && e.data)[0];
    const x = gd.data[0].x, y = gd.data[0].y;
    let i = 0;
    for (let k = 1; k < y.length; k++) if (y[k] > y[i]) i = k;
    const fl = gd._fullLayout, rect = gd.getBoundingClientRect();
    const cx = rect.left + fl.margin.l + fl.xaxis.l2p(x[i]);
    const cy = rect.top + fl.margin.t + fl.yaxis.l2p(y[i]);
    return { sx: cx - 25, sy: cy - 25, ex: cx + 25, ey: cy + 25 };
  });
  await page.mouse.move(target.sx, target.sy);
  await page.mouse.down();
  await page.mouse.move(target.ex, target.ey, { steps: 10 });
  await page.mouse.up();
  await settle();

  const beforeSel = await selCount(page);
  const beforeCensus = await selectionCensus(page);
  record('setup: manual box selection is live on the bus',
    beforeSel.n === 1, `count=${beforeSel.n}`);
  record('setup: selection artifacts visible in plotly DOM',
    beforeCensus.outlines > 0 || beforeCensus.selLayer > 0,
    JSON.stringify(beforeCensus).slice(0, 300));
  await page.screenshot({ path: path.join(ARTIFACTS, 'cleararea_1_selected.png') });

  // -- 2. click "Clear figure selection" ---------------------------------
  const clearBtn = page.locator('[data-testid$="-clear-selection-button"]').first();
  const btnVisible = await clearBtn.count();
  record('setup: clear button found', btnVisible === 1,
    await clearBtn.getAttribute('data-testid').catch(() => 'n/a'));

  await startPlotWatch(page);
  await clearBtn.click();
  await settle(7000);
  const renders = await readPlotRenders(page);

  const afterSel = await selCount(page);
  const afterCensus = await selectionCensus(page);
  record('clear: bus selection dropped to 0', afterSel.n === 0, `count=${afterSel.n}`);
  record('clear: figure repainted (DOM mutations observed)', renders > 0, `renders=${renders}`);
  record('clear: marker emphasis dropped (server-side render current)',
    afterCensus.solid === 0 || afterCensus.solid === null,
    `solid=${afterCensus.solid} dim=${afterCensus.dim}`);
  record('BUG CHECK: plotly x/y area (box) selection artifact GONE after clear',
    afterCensus.outlines === 0 && afterCensus.selLayer === 0,
    JSON.stringify(afterCensus).slice(0, 300));
  await page.screenshot({ path: path.join(ARTIFACTS, 'cleararea_2_after_clear.png') });

  // -- 3. control: a plain re-selection still works after clear ----------
  await page.click('.modebar-btn[data-title="Box Select"]');
  const t2 = await page.evaluate(() => {
    const gd = Array.from(document.querySelectorAll('.js-plotly-plot'))
      .filter(e => e.offsetParent !== null && e.data)[0];
    const x = gd.data[0].x, y = gd.data[0].y;
    let i = 0;
    for (let k = 1; k < y.length; k++) if (y[k] < y[i]) i = k; // bottommost
    const fl = gd._fullLayout, rect = gd.getBoundingClientRect();
    const cx = rect.left + fl.margin.l + fl.xaxis.l2p(x[i]);
    const cy = rect.top + fl.margin.t + fl.yaxis.l2p(y[i]);
    return { sx: cx - 25, sy: cy - 25, ex: cx + 25, ey: cy + 25 };
  });
  await page.mouse.move(t2.sx, t2.sy);
  await page.mouse.down();
  await page.mouse.move(t2.ex, t2.ey, { steps: 10 });
  await page.mouse.up();
  await settle();
  const reSel = await selCount(page);
  record('control: re-selection after clear still works', reSel.n === 1,
    `count=${reSel.n}`);

  await ctx.close();
  await browser.close();
} catch (e) {
  record('script crashed', false, String(e));
  exitCode = 1;
}

console.log('\n---- summary ----');
const fails = results.filter(r => !r.pass).length;
console.log(`${results.length - fails}/${results.length} checks passed`);
if (rLog) {
  const tail = rLog.split('\n').slice(-25).join('\n');
  console.log('---- R log tail ----\n' + tail);
}
process.exit(exitCode || (fails ? 1 : 0));
