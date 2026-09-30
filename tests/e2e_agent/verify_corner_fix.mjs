// Browser-level verification of the two reported corner-selection defects.
//
// Problem 1: manual square selection, then a volcano view switch - the
//            selection bus must carry the re-derived corner selection (was:
//            empty after the post-render echo cleared it).
// Problem 2: Area volcano -> topleft -> volcano on unchanged axes - the
//            marker emphasis must follow the live rects (was: stale opacity
//            from the previous corner while the bus carried the new ids).
//
// Run:  node verify_corner_fix.mjs     (from tests/e2e_agent/)

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

const PORT = 7779;
const results = [];
const record = (name, pass, detail = '') => {
  results.push({ name, pass, detail: String(detail).slice(0, 300) });
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

// marker opacity census of the visible feature scatter
const opacityCensus = (page) => page.evaluate(() => {
  const gd = Array.from(document.querySelectorAll('.js-plotly-plot'))
    .filter(e => e.offsetParent !== null && e.data)[0];
  if (!gd) return null;
  const op = gd.data[0].marker.opacity;
  if (!Array.isArray(op)) return { solid: null, dim: null, total: gd.data[0].x.length };
  let solid = 0, dim = 0;
  for (const v of op) { if (v > 0.5) solid++; else dim++; }
  return { solid, dim, total: op.length };
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

  // ================= Scenario A: manual selection + view switch ==========
  {
    const ctx = await browser.newContext();
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

    const loadSel = await selCount(page);
    const loadOp = await opacityCensus(page);
    record('A: load-time volcano corner selection is live on the bus',
      loadSel.n > 10, `count=${loadSel.n}`);
    record('A: load-time emphasis matches the corner selection',
      loadOp.solid === loadSel.n, `solid=${loadOp.solid} count=${loadSel.n}`);
    await page.screenshot({ path: path.join(ARTIFACTS, 'cornerfix_A1_load.png'), fullPage: false });

    // manual box selection around one isolated point (topmost): use the
    // real modebar button (the reliable path - a bare Plotly.relayout
    // dragmode flip does not always arm the drag cover)
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
      return { sx: cx - 18, sy: cy - 18, ex: cx + 18, ey: cy + 18 };
    });
    await page.mouse.move(target.sx, target.sy);
    await page.mouse.down();
    await page.mouse.move(target.ex, target.ey, { steps: 8 });
    await page.mouse.up();
    await settle();

    const manualSel = await selCount(page);
    record('A: manual square selection reports exactly one feature',
      manualSel.n === 1, `count=${manualSel.n}`);
    await page.screenshot({ path: path.join(ARTIFACTS, 'cornerfix_A2_manual.png'), fullPage: false });

    // switch to another volcano view
    await page.click('[data-quick-view-id="volcano_RE_vs_LE"]');
    await page.waitForFunction(() => {
      const gd = Array.from(document.querySelectorAll('.js-plotly-plot'))
        .filter(e => e.offsetParent !== null && e._fullLayout)[0];
      if (!gd) return false;
      const t = (ax) => (ax && ax.title ? (ax.title.text || '') : '');
      return /RE_vs_LE/.test(t(gd._fullLayout.xaxis));
    }, null, { timeout: 45000 });
    await settle();
    await settle();

    const switchSel = await selCount(page);
    const switchOp = await opacityCensus(page);
    record('A: view switch re-derives the volcano corner selection (was: empty)',
      switchSel.n > 10, `count=${switchSel.n}`);
    record('A: switched figure emphasizes the re-derived corner',
      switchOp.solid === switchSel.n, `solid=${switchOp.solid} count=${switchSel.n}`);
    await page.screenshot({ path: path.join(ARTIFACTS, 'cornerfix_A3_switch.png'), fullPage: false });
    await ctx.close();
  }

  // ================= Scenario B: Area topleft -> volcano emphasis ========
  {
    const ctx = await browser.newContext();
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

    const bothSel = await selCount(page);
    const bothOp = await opacityCensus(page);
    record('B: load: both volcano corners selected and solid',
      bothOp.solid === bothSel.n && bothSel.n > 10,
      `solid=${bothOp.solid} count=${bothSel.n}`);
    await page.screenshot({ path: path.join(ARTIFACTS, 'cornerfix_B1_volcano.png'), fullPage: false });

    // Area -> topleft (selectize API; the dropdown need not be open)
    await page.evaluate(() => {
      const el = document.querySelector('select[id$="feature_space-a4selector-scorner"]');
      if (!el) throw new Error('scorner select not found');
      if (el.selectize) el.selectize.setValue('topleft');
      else {
        el.value = 'topleft';
        el.dispatchEvent(new Event('change', { bubbles: true }));
      }
    });
    await settle();

    const tlSel = await selCount(page);
    const tlOp = await opacityCensus(page);
    record('B: Area topleft narrows the selection and the emphasis follows',
      tlOp.solid === tlSel.n && tlSel.n < bothSel.n,
      `solid=${tlOp.solid} count=${tlSel.n} (was ${bothSel.n})`);
    await page.screenshot({ path: path.join(ARTIFACTS, 'cornerfix_B2_topleft.png'), fullPage: false });

    // Area -> volcano again (the reported bug: stale topleft emphasis)
    await page.evaluate(() => {
      const el = document.querySelector('select[id$="feature_space-a4selector-scorner"]');
      if (el.selectize) el.selectize.setValue('volcano');
      else {
        el.value = 'volcano';
        el.dispatchEvent(new Event('change', { bubbles: true }));
      }
    });
    await settle();
    await settle();

    const backSel = await selCount(page);
    const backOp = await opacityCensus(page);
    record('B: Area volcano restores BOTH corners in the emphasis (was: stale topleft)',
      backOp.solid === backSel.n && backSel.n === bothSel.n,
      `solid=${backOp.solid} count=${backSel.n} (load: ${bothSel.n})`);
    await page.screenshot({ path: path.join(ARTIFACTS, 'cornerfix_B3_volcano.png'), fullPage: false });
    await ctx.close();
  }

  await browser.close();
} catch (e) {
  console.error('HARNESS ERROR:', e);
  exitCode = 1;
}

const failed = results.filter(x => !x.pass).length;
console.log(`\n${results.length - failed}/${results.length} passed`);
if (failed) { exitCode = 1; console.log('--- app log tail ---\n' + rLog.slice(-2000)); }
process.exit(exitCode);
