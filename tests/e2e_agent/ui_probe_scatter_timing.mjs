// Probe 3: does the set_scatter_view seam fail on a fresh app?
// Sequence: scatter-hook quick view (fresh) -> widgets warm-up -> retry
// scatter hook -> explicit axes variant. Mirrors what the assistant did.
import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import { mkdirSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../..');
const EXTDATA = path.join(REPO, 'inst/extdata');
const OUT = '/home/chen/.omicsviewer-llm-eval/probe4';
mkdirSync(OUT, { recursive: true });
const PORT = 7789;

const r = spawn('Rscript', ['-e', `
  options(shiny.port = ${PORT}, shiny.host = '127.0.0.1')
  eset <- readRDS(file.path('${EXTDATA}', 'demo.RDS'))
  omicsViewer::omicsViewer(dir = '${EXTDATA}', ESVObj = eset)
`], { cwd: REPO,
  env: { ...process.env, OMICSVIEWER_TEST_HOOKS: 'true' },
  stdio: ['ignore', 'pipe', 'pipe'] });
process.on('exit', () => { try { r.kill('SIGKILL'); } catch {} });

await (async () => {
  const t0 = Date.now();
  while (Date.now() - t0 < 90000) {
    try { const res = await fetch(`http://127.0.0.1:${PORT}/`); if (res.status > 0) return; } catch {}
    await new Promise(s => setTimeout(s, 500));
  }
  throw new Error('app did not start');
})();

const browser = await chromium.launch({
  executablePath: '/usr/bin/google-chrome', headless: true,
  args: ['--no-sandbox', '--disable-gpu', '--disable-dev-shm-usage',
         '--disable-webgl', '--disable-webgl2']
});
const page = await (await browser.newContext()).newPage();
await page.goto(`http://127.0.0.1:${PORT}/`, { timeout: 60000 });
const HOOK = {
  box: '#app-agentTestHooks-container',
  op: '#app-agentTestHooks-op', payload: '#app-agentTestHooks-payload',
  run: '#app-agentTestHooks-run', result: '#app-agentTestHooks-result'
};
await page.waitForFunction(() => !!document.querySelector('#app-contents')?.offsetParent,
  null, { timeout: 90000 });
await page.evaluate(sel => { document.querySelector(sel).style.display = 'block'; }, HOOK.box);

const runHook = async (op, payload) => {
  const before = (await page.innerText(HOOK.result)).trim();
  await page.evaluate(({ sel, value }) => {
    const el = document.querySelector(sel);
    if (el.selectize) el.selectize.setValue(value);
    else { el.value = value; el.dispatchEvent(new Event('change', { bubbles: true })); }
  }, { sel: HOOK.op, value: op });
  await page.fill(HOOK.payload, JSON.stringify(payload ?? {}));
  const n = await page.evaluate(sel => {
    const el = document.querySelector(sel);
    const k = (parseInt(el.dataset.runs || '0', 10) + 1); el.dataset.runs = String(k); return k;
  }, HOOK.run);
  await page.click(HOOK.run);
  await page.waitForFunction(({ sel, before, n }) => {
    const t = (document.querySelector(sel)?.innerText || '').trim();
    return t !== before && t.includes(`"hook_run": ${n}`);
  }, { sel: HOOK.result, before, n }, { timeout: 35000 });
  const raw = (await page.innerText(HOOK.result)).trim();
  try { return JSON.parse(raw); } catch { return raw; }
};
const show = (label, res) =>
  console.log(label, '->', JSON.stringify(res).slice(0, 300));

await page.waitForTimeout(4000);
show('S0 capabilities read (mimic turn 1)', await runHook('capabilities',
  { query: 'download', max_results: 5 }).then(r => Object.keys(r).length));
await page.waitForTimeout(90000);   // mimic batch3: ~90s between load and the tool call
show('S1 fresh scatter quick_view  ', await runHook('scatter',
  { space: 'feature', quick_view_id: 'volcano_MT_vs_WT' }));
show('S2 fresh scatter explicit   ', await runHook('scatter',
  { space: 'feature', x_axis: 'ttest|MT_vs_WT|mean.diff', y_axis: 'ttest|MT_vs_WT|log.fdr' }));
// warm the axis providers through the widgets seam (agent_widget_apply)
show('W  widgets warm-up          ', await runHook('widgets', {
  patch: { 'dataspace.feature_space.x_analysis': 'ttest' } }));
await page.waitForTimeout(3000);
show('S3 scatter quick after warm ', await runHook('scatter',
  { space: 'feature', quick_view_id: 'volcano_MT_vs_WT' }));
show('S4 scatter explicit after   ', await runHook('scatter',
  { space: 'feature', x_axis: 'ttest|MT_vs_WT|mean.diff', y_axis: 'ttest|MT_vs_WT|log.fdr' }));
await browser.close();
process.exit(0);
