// Probe 2: reproduce the fresh-app scatter-axis patch failure at the
// store level via the official test-hook bridge (OMICSVIEWER_TEST_HOOKS).
import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import { mkdirSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../..');
const EXTDATA = path.join(REPO, 'inst/extdata');
const OUT = '/home/chen/.omicsviewer-llm-eval/probe2';
mkdirSync(OUT, { recursive: true });
const PORT = 7787;

const r = spawn('Rscript', ['-e', `
  options(shiny.port = ${PORT}, shiny.host = '127.0.0.1')
  eset <- readRDS(file.path('${EXTDATA}', 'demo.RDS'))
  omicsViewer::omicsViewer(dir = '${EXTDATA}', ESVObj = eset)
`], { cwd: REPO,
  env: { ...process.env, OMICSVIEWER_TEST_HOOKS: 'true' },
  stdio: ['ignore', 'pipe', 'pipe'] });
let rLog = '';
r.stderr.on('data', d => { rLog += d; });
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
await page.waitForFunction(sel => !!document.querySelector('#app-contents')?.offsetParent,
  HOOK.box, { timeout: 90000 });
await page.evaluate(sel => { document.querySelector(sel).style.display = 'block'; }, HOOK.box);

const patch = { 'dataspace.feature_space.x_analysis': 'ttest',
                'dataspace.feature_space.x_subset': 'MT_vs_WT',
                'dataspace.feature_space.x_variable': 'mean.diff' };
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

const attempt = async label => {
  const res = await runHook('widgets', { patch });
  const s = JSON.stringify(res);
  console.log(label, '->',
    /error|unavailable|reject/i.test(s) ? 'REJECTED: ' + s.slice(0, 260) : 'applied: ' + s.slice(0, 120));
};

await page.waitForTimeout(3000);
await attempt('A fresh-app immediate');
await page.waitForTimeout(15000);
await attempt('B +15s idle           ');
const clicked = await page.evaluate(() => {
  for (const a of document.querySelectorAll('.navbar-nav a, .nav-tabs a'))
    if (a.innerText.trim() === 'Feature' && a.offsetParent) { a.click(); return true; }
  return false;
});
console.log('clicked Feature tab:', clicked);
await page.waitForTimeout(5000);
await attempt('C after tab click     ');
const drove = await page.evaluate(() => {
  const el = document.querySelector('#app-dataspace-feature_space-tris-x_analysis');
  if (el && el.selectize) { el.selectize.setValue('ttest'); return true; }
  return false;
});
console.log('selectize setValue x_analysis:', drove);
await page.waitForTimeout(5000);
await attempt('D after selectize drive');
await browser.close();
process.exit(0);
