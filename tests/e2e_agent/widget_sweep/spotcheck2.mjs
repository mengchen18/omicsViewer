import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import path from 'node:path';
import fs from 'node:fs';
import { fileURLToPath } from 'node:url';
const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../../..');
const PORT = 7799;
const r = spawn('Rscript', ['-e', `
  options(shiny.port = ${PORT}, shiny.host = '127.0.0.1')
  eset <- readRDS(file.path('${REPO}/inst/extdata', 'demo.RDS'))
  omicsViewer::omicsViewer(dir = '${REPO}/inst/extdata', ESVObj = eset)
`], { cwd: REPO, env: { ...process.env, OMICSVIEWER_TEST_HOOKS: 'true' }, stdio: ['ignore', 'pipe', 'pipe'] });
process.on('exit', () => { try { r.kill('SIGKILL'); } catch {} });
const waitPort = async (p, ms = 120000) => { const t0 = Date.now(); while (Date.now() - t0 < ms) { try { const res = await fetch(`http://127.0.0.1:${p}/`); if (res.status > 0) return true; } catch {} await new Promise(s => setTimeout(s, 500)); } return false; };
const HOOK = { op: '#app-agentTestHooks-op', payload: '#app-agentTestHooks-payload', run: '#app-agentTestHooks-run', result: '#app-agentTestHooks-result', box: '#app-agentTestHooks-container' };
async function runHook(page, op, payload) {
  const before = (await page.innerText(HOOK.result)).trim();
  await page.evaluate(({ sel, value }) => { const el = document.querySelector(sel); if (el.selectize) el.selectize.setValue(value); else { el.value = value; el.dispatchEvent(new Event('change', { bubbles: true })); } }, { sel: HOOK.op, value: op });
  await page.fill(HOOK.payload, JSON.stringify(payload ?? {}));
  const n = await page.evaluate(sel => { const el = document.querySelector(sel); const n = (parseInt(el.dataset.runs || '0', 10) + 1); el.dataset.runs = String(n); return n; }, HOOK.run);
  await page.click(HOOK.run);
  await page.waitForFunction(({ sel, n }) => (document.querySelector(sel)?.innerText || '').trim().includes(`"hook_run": ${n}`), { sel: HOOK.result, n }, { timeout: 90000 });
  return JSON.parse((await page.innerText(HOOK.result)).trim());
}
const readInputs = (page, re) => page.evaluate(re => { const iv = Shiny.shinyapp.$inputValues || {}; const o = {}; for (const k of Object.keys(iv)) if (new RegExp(re).test(k)) o[k] = iv[k]; return o; }, re);
const wait = ms => new Promise(s => setTimeout(s, ms));
if (!(await waitPort(PORT))) { console.log('NO APP'); process.exit(1); }
const browser = await chromium.launch({ executablePath: '/usr/bin/google-chrome', headless: true, args: ['--no-sandbox', '--disable-gpu', '--disable-dev-shm-usage', '--disable-webgl', '--disable-webgl2'] });
const page = await (await browser.newContext()).newPage();
await page.goto(`http://127.0.0.1:${PORT}/`, { timeout: 90000 });
await page.waitForFunction(() => Shiny.shinyapp.$inputValues['app-dataspace-eset'] && document.querySelector('#app-contents').offsetParent !== null, null, { timeout: 90000 });
await page.evaluate(sel => { document.querySelector(sel).style.display = 'block'; }, HOOK.box);
await wait(2000);
const out = {};
// A) REAL single-key cascade write: color_analysis General -> PCA
out.before = await readInputs(page, 'feature_space-a4selector-selectColorUI-');
const rcA = await runHook(page, 'widgets', { patch: { 'dataspace.feature_space.attr4.color_analysis': 'PCA' } });
out.receipt_A = { applied: rcA.applied, skipped: rcA.unchanged };
await wait(10000);
out.after_analysis_PCA = await readInputs(page, 'feature_space-a4selector-selectColorUI-');
// B) open the Feature table tab, THEN generic columns write
await page.evaluate(() => { const a = Array.from(document.querySelectorAll('a[data-toggle="tab"]')).find(x => x.offsetParent !== null && /Feature table/.test(x.innerText)); if (a) a.click(); });
await wait(3000);
out.tab_open_table_state = await readInputs(page, 'dataspace-tab_feature-table_state');
const rcB = await runHook(page, 'widgets', { patch: { 'dataspace.tab_feature.columns': ['featureId', 'General|All|Gene.name'] } });
out.receipt_B = { applied: rcB.applied, skipped: rcB.unchanged };
await wait(8000);
out.after_columns = await readInputs(page, 'dataspace-tab_feature-table_state');
const rcC = await runHook(page, 'widgets', { patch: { 'dataspace.tab_feature.page': 2 } });
await wait(8000);
out.after_page = await readInputs(page, 'dataspace-tab_feature-table_state');
console.log(JSON.stringify(out, (k, v) => typeof v === 'string' ? v.slice(0, 70) + (v.length > 70 ? '…' : '') : v, 1));
fs.writeFileSync(path.join(here, 'spotcheck2_results.json'), JSON.stringify(out, null, 2));
await browser.close(); process.exit(0);
