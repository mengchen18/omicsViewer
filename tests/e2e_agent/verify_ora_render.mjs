// Browser-level guard against ORA results-table render amplification.
//
// One volcano quick-view switch = one selection-bus report + one ORA
// computation; the ORA results table must re-render essentially once.
// The dataTableDownload R-M2 fallback epoch used to bump on EVERY
// post-render input$table_state NULL window (the browser resets it while
// DT re-initializes), feeding a self-sustaining re-render loop: one
// selection change re-rendered the identical table 12-20 times (visible
// flashing, run length varying with browser timing - user reported
// "renders 2, 3, 4 or more times, not stable").
//
// Counts DOM table (re)initializations under the ORA results panel;
// DataTables.js rebuilds the table node once per render, so the counter
// sees ~2 nodes per actual render.
//
// Run:  node verify_ora_render.mjs     (from tests/e2e_agent/)

import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../..');
const EXTDATA = path.join(REPO, 'inst/extdata');
const PORT = 7791;

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

let exitCode = 0;
try {
  await waitPort(PORT);
  const browser = await chromium.launch({
    executablePath: '/usr/bin/google-chrome',
    headless: true,
    args: ['--no-sandbox', '--disable-gpu', '--disable-dev-shm-usage',
           '--disable-webgl', '--disable-webgl2']
  });
  const ctx = await browser.newContext();
  const page = await ctx.newPage();
  await page.goto(`http://127.0.0.1:${PORT}/`, { timeout: 90000 });
  await page.waitForFunction(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-eset'] &&
    document.querySelector('#app-contents').offsetParent !== null, null, { timeout: 90000 });
  await page.click('.navbar a:has-text("ORA")');
  await page.waitForSelector('[id*="ora-stab-table"] table', { timeout: 90000, state: 'attached' });
  await page.waitForFunction(() =>
    (document.querySelector('[id*="ora-stab-table"] table tbody') || { children: [] }).children.length > 0,
    null, { timeout: 90000 });
  await new Promise(s => setTimeout(s, 6000));

  await page.evaluate(() => {
    window.__renders = 0;
    const mo = new MutationObserver(muts => {
      for (const m of muts) for (const n of m.addedNodes) {
        if (n.nodeType !== 1) continue;
        const tb = n.tagName === 'TABLE' ? n : (n.querySelector ? n.querySelector('table') : null);
        if (tb && tb.closest('[id*="ora-stab"]')) window.__renders++;
      }
    });
    mo.observe(document.body, { childList: true, subtree: true });
  });

  const views = ['volcano_RE_vs_LE', 'volcano_RE_vs_ME', 'volcano_MT_vs_WT', 'volcano_RE_vs_ME'];
  for (const v of views) {
    const before = await page.evaluate(() => window.__renders);
    await page.click(`[data-quick-view-id="${v}"]`);
    await new Promise(s => setTimeout(s, 10000));
    const n = (await page.evaluate(() => window.__renders)) - before;
    const rows = await page.evaluate(() =>
      (document.querySelector('[id*="ora-stab-table"] table tbody') || { children: [] }).children.length);
    // 2 counted nodes == 1 actual render (DataTables rebuilds the node);
    // allow <= 3 for a benign extra init, fail on loops (was 12-20)
    record(`view switch -> ${v}: ORA table renders once (counted ${n}, rows ${rows})`,
      n <= 3, `counted=${n} rows=${rows}`);
  }
  await ctx.close();
  await browser.close();
} catch (e) {
  console.error('HARNESS ERROR:', e);
  exitCode = 1;
}

const failed = results.filter(x => !x.pass).length;
console.log(`\n${results.length - failed}/${results.length} passed`);
if (failed) exitCode = 1;
process.exit(exitCode);
