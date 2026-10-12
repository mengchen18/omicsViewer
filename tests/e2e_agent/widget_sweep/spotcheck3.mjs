// spotcheck3 — WP-B2 (F6): expr-heatmap cases 022/026 live in isolation.
// The plan names this tool but it did not exist; built 2026-10-12 for this WP.
//
// Question: do writes to dataspace.expr_heatmap.row_sort_by (022) and
// .annot_row (026) reach the browser inputs when observed PROPERLY -
// i.e. polling for the SPECIFIC input key (the sweep's waitUiEffect
// returns on the first input diff, which is always the hook's own
// controls: waited_ms 21-22 in the committed evidence, so its real
// observation window was ~20 ms, not 30 s)?
//
// Per case, measures:
//   * receipt (applied/rejected)          * time-to-effect on the target input
//   * store read-back (correct result.values path, post-WP-B harness lesson)
//   * main-thread stall profile (heartbeat timestamp gaps around the write)
//   * the expr-heatmap re-render duration (the alleged 2702x60 blocker)
// Run: node spotcheck3.mjs            (from tests/e2e_agent/widget_sweep/)
//      node spotcheck3.mjs --tabbed   (also activate the Heatmap tab first)
import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import path from 'node:path';
import fs from 'node:fs';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../../..');
const EXTDATA = path.join(REPO, 'inst/extdata');
const PORT = 7793;
const TABBED = process.argv.includes('--tabbed');
const r = spawn('Rscript', ['-e', `
  options(shiny.port = ${PORT}, shiny.host = '127.0.0.1')
  eset <- readRDS(file.path('${EXTDATA}', 'demo.RDS'))
  omicsViewer::omicsViewer(dir = '${EXTDATA}', ESVObj = eset)
`], { cwd: REPO, env: { ...process.env, OMICSVIEWER_TEST_HOOKS: 'true' }, stdio: ['ignore', 'pipe', 'pipe'] });
process.on('exit', () => { try { r.kill('SIGKILL'); } catch {} });
process.on('SIGINT', () => { try { r.kill('SIGKILL'); } catch {}; process.exit(130); });

const waitPort = async (port, ms = 120000) => {
  const t0 = Date.now();
  while (Date.now() - t0 < ms) {
    try { const res = await fetch(`http://127.0.0.1:${port}/`); if (res.status > 0) return true; } catch {}
    await new Promise(s => setTimeout(s, 500));
  }
  return false;
};
const HOOK = { op: '#app-agentTestHooks-op', payload: '#app-agentTestHooks-payload',
  run: '#app-agentTestHooks-run', result: '#app-agentTestHooks-result', box: '#app-agentTestHooks-container' };
async function runHook(page, op, payload) {
  const before = (await page.innerText(HOOK.result)).trim();
  await page.evaluate(({ sel, value }) => {
    const el = document.querySelector(sel);
    if (el.selectize) el.selectize.setValue(value);
    else { el.value = value; el.dispatchEvent(new Event('change', { bubbles: true })); }
  }, { sel: HOOK.op, value: op });
  await page.fill(HOOK.payload, JSON.stringify(payload ?? {}));
  const n = await page.evaluate(sel => {
    const el = document.querySelector(sel);
    const n = (parseInt(el.dataset.runs || '0', 10) + 1);
    el.dataset.runs = String(n); return n;
  }, HOOK.run);
  await page.click(HOOK.run);
  await page.waitForFunction(({ sel, n, before }) => {
    const t = (document.querySelector(sel)?.innerText || '').trim();
    return t !== before && t.includes(`"hook_run": ${n}`);
  }, { sel: HOOK.result, n, before }, { timeout: 120000 });
  return JSON.parse((await page.innerText(HOOK.result)).trim());
}
const inputVal = (page, key) => page.evaluate(k => (Shiny.shinyapp.$inputValues || {})[k] ?? null, key);

// main-thread heartbeat: 50 ms interval; gaps > 500 ms = blocked main thread
const startHeartbeat = page => page.evaluate(() => {
  window.__hb = [];
  window.__hbTimer = setInterval(() => window.__hb.push(Date.now()), 50);
});
const stopHeartbeat = page => page.evaluate(() => {
  clearInterval(window.__hbTimer);
  const t = window.__hb || [];
  const gaps = [];
  for (let i = 1; i < t.length; i++) if (t[i] - t[i - 1] > 500) gaps.push(t[i] - t[i - 1]);
  return { ticks: t.length, stalls: gaps, totalBlockedMs: gaps.reduce((a, b) => a + b, 0) };
});

// poll until the target input deep-equals want (or timeout); returns ms or null
async function waitInput(page, key, want, budgetMs) {
  const eq = v => JSON.stringify(v) === JSON.stringify(want);
  const t0 = Date.now();
  let last = await inputVal(page, key);
  if (eq(last)) return { ms: 0, ok: true, before: last };
  while (Date.now() - t0 < budgetMs) {
    await new Promise(s => setTimeout(s, 50));
    last = await inputVal(page, key);
    if (eq(last)) return { ms: Date.now() - t0, ok: true, before: last };
  }
  return { ms: Date.now() - t0, ok: false, before: last };
}

const CASES = [
  { n: '022', id: 'dataspace.expr_heatmap.row_sort_by', input: 'app-dataspace-heatmapViewer-rowSortBy', write: process.env.ROW_WRITE || 'General|All|Gene.name' },
  { n: '026', id: 'dataspace.expr_heatmap.annot_row', input: 'app-dataspace-heatmapViewer-annotRow', write: ['General|All|Gene.name'] }
];

const out = { generated: new Date().toISOString(), tabbed: TABBED, cases: [] };
try {
  if (!(await waitPort(PORT))) throw new Error('app did not start');
  const browser = await chromium.launch({ executablePath: '/usr/bin/google-chrome', headless: true,
    args: ['--no-sandbox', '--disable-gpu', '--disable-dev-shm-usage', '--disable-webgl', '--disable-webgl2'] });
  const page = await (await browser.newContext()).newPage();
  await page.goto(`http://127.0.0.1:${PORT}/`, { timeout: 90000 });
  await page.waitForFunction(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-eset'] &&
    document.querySelector('#app-contents').offsetParent !== null, null, { timeout: 90000 });
  await page.evaluate(sel => { document.querySelector(sel).style.display = 'block'; }, HOOK.box);
  await new Promise(s => setTimeout(s, 1500));

  if (TABBED) {
    const opened = await page.evaluate(() => {
      const a = Array.from(document.querySelectorAll('a[data-toggle="tab"]'))
        .find(x => x.offsetParent !== null && /^\s*Heatmap\s*$/.test(x.innerText));
      if (a) { a.click(); return true; } return false;
    });
    out.heatmap_tab_opened = opened;
    await new Promise(s => setTimeout(s, 8000)); // let the 2702x60 heatmap render once
  }

  for (const c of CASES) {
    await startHeartbeat(page);
    const before = await inputVal(page, c.input);
    const t0 = Date.now();
    const receipt = await runHook(page, 'widgets', { patch: { [c.id]: c.write } });
    const eff = await waitInput(page, c.input, c.write, 60000);
    const storeRb = (await runHook(page, 'store', {})).values?.[c.id] ?? null;
    const hb = await stopHeartbeat(page);
    out.cases.push({
      case: c.n, binding: c.id, input: c.input, written: c.write,
      applied: receipt.applied ?? [], rejected: receipt.rejected ?? [],
      input_before: before, time_to_effect_ms: eff.ms, effect_observed: eff.ok,
      store_readback: storeRb,
      heartbeat: hb,
      total_case_ms: Date.now() - t0
    });
    console.log(`${c.n} ${c.id}: applied=${(receipt.applied || []).length > 0} ` +
      `effect=${eff.ok ? eff.ms + 'ms' : 'NOT OBSERVED in 60s'} ` +
      `store=${JSON.stringify(storeRb)} stalls=${JSON.stringify(hb.stalls)}`);
    await new Promise(s => setTimeout(s, 1500));
  }
  await browser.close();
} catch (e) {
  out.fatal = String(e);
  console.log('FATAL', String(e));
}
fs.writeFileSync(path.join(here, 'spotcheck3_results.json'), JSON.stringify(out, null, 2));
const ok = out.cases.length === CASES.length && out.cases.every(c => c.effect_observed);
console.log(ok ? 'SPOTCHECK3: both cases GREEN in isolation' : 'SPOTCHECK3: at least one case NOT observed');
process.exit(ok ? 0 : 1);
