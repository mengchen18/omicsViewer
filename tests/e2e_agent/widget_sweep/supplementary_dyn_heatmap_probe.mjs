// Supplementary probe: dynamic-heatmap widget operability vs tab state.
// The sweep recorded dataspace.dyn_heatmap.{col_sort_by,row_sort_by,
// annot_col,annot_row,tooltip_info} as no-probe because their capability
// records expose NO allowed values (choices providers evaluate to NULL
// while the Dynamic heatmap tab is not active). This probe answers the
// follow-up question: can set_widgets operate them at all?
//   A) write while the tab is INACTIVE  -> expect rejection (validator
//      uses the same provider)
//   B) activate the tab, re-read capability records, write again
// Run: node supplementary_dyn_heatmap_probe.mjs
import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import path from 'node:path';
import fs from 'node:fs';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../../..');
const EXTDATA = path.join(REPO, 'inst/extdata');
const PORT = 7795;
const r = spawn('Rscript', ['-e', `
  options(shiny.port = ${PORT}, shiny.host = '127.0.0.1')
  eset <- readRDS(file.path('${EXTDATA}', 'demo.RDS'))
  omicsViewer::omicsViewer(dir = '${EXTDATA}', ESVObj = eset)
`], { cwd: REPO, env: { ...process.env, OMICSVIEWER_TEST_HOOKS: 'true' }, stdio: ['ignore', 'pipe', 'pipe'] });
process.on('exit', () => { try { r.kill('SIGKILL'); } catch {} });

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
  await page.waitForFunction(({ sel, n }) =>
    (document.querySelector(sel)?.innerText || '').trim().includes(`"hook_run": ${n}`),
    { sel: HOOK.result, n }, { timeout: 60000 });
  return JSON.parse((await page.innerText(HOOK.result)).trim());
}
const PROBES = {
  'dataspace.dyn_heatmap.annot_col': ['General|All|Origin'],
  'dataspace.dyn_heatmap.annot_row': ['General|All|Gene.name'],
  'dataspace.dyn_heatmap.tooltip_info': ['General|All|Gene.name'],
  'dataspace.dyn_heatmap.col_sort_by': 'none',
  'dataspace.dyn_heatmap.row_sort_by': 'none'
};
const out = { generated: new Date().toISOString(), probes: [] };
try {
  if (!(await waitPort(PORT))) throw new Error('no app');
  const browser = await chromium.launch({ executablePath: '/usr/bin/google-chrome', headless: true,
    args: ['--no-sandbox', '--disable-gpu', '--disable-dev-shm-usage', '--disable-webgl', '--disable-webgl2'] });
  const page = await (await browser.newContext()).newPage();
  await page.goto(`http://127.0.0.1:${PORT}/`, { timeout: 90000 });
  await page.waitForFunction(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-eset'] &&
    document.querySelector('#app-contents').offsetParent !== null, null, { timeout: 90000 });
  await page.evaluate(sel => { document.querySelector(sel).style.display = 'block'; }, HOOK.box);

  for (const [id, value] of Object.entries(PROBES)) {
    const rec = (await runHook(page, 'capabilities', { query: id, max_results: 5 })).capabilities?.find(c => c.id === id);
    const receiptA = await runHook(page, 'widgets', { patch: { [id]: value } });
    out.probes.push({ id, phase: 'A: tab inactive', allowed_at_query: rec?.allowed_values ?? null,
      applied: receiptA.applied ?? [], rejected: receiptA.rejected ?? [], error: receiptA.hook_error ?? null });
  }
  // activate the Dynamic heatmap tab
  await page.evaluate(() => {
    const a = Array.from(document.querySelectorAll('a[data-toggle="tab"]'))
      .find(x => x.offsetParent !== null && /Dynamic heatmap/.test(x.innerText));
    if (a) a.click();
  });
  await page.waitForTimeout(4000);
  for (const [id, value] of Object.entries(PROBES)) {
    const rec = (await runHook(page, 'capabilities', { query: id, max_results: 5 })).capabilities?.find(c => c.id === id);
    const receiptB = await runHook(page, 'widgets', { patch: { [id]: value } });
    const ui = await page.evaluate(() => {
      const iv = Shiny.shinyapp.$inputValues || {};
      const o = {}; for (const k of Object.keys(iv)) if (/dynheatmap/.test(k)) o[k] = iv[k];
      return o;
    });
    out.probes.push({ id, phase: 'B: tab active', allowed_at_query: rec?.allowed_values ?? null,
      applied: receiptB.applied ?? [], rejected: receiptB.rejected ?? [], error: receiptB.hook_error ?? null,
      dyn_inputs_after: ui });
  }
  await browser.close();
} catch (e) { out.fatal = String(e); }
fs.writeFileSync(path.join(here, 'supplementary_dyn_heatmap.json'), JSON.stringify(out, null, 2));
for (const p of out.probes) {
  const ok = (p.applied || []).length > 0;
  console.log(`${p.phase} | ${p.id} | ${ok ? 'APPLIED' : 'REJECTED'} | allowed: ${JSON.stringify(p.allowed_at_query)?.slice(0, 60)}`);
}
process.exit(0);
