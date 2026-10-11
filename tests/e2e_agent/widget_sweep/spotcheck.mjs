// Live spot-check of the 'no-ui-effect' classes: single-key cascade writes,
// axis_mode, generic table writes, survival censor. Decides TRUE no-effect
// vs instrumentation artifact.
import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../../..');
const PORT = 7798;
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
const readInputs = (page, re) => page.evaluate(re => {
  const iv = Shiny.shinyapp.$inputValues || {}; const o = {};
  for (const k of Object.keys(iv)) if (new RegExp(re).test(k)) o[k] = iv[k];
  return o;
}, re);
const wait = ms => new Promise(s => setTimeout(s, ms));

if (!(await waitPort(PORT))) { console.log('NO APP'); process.exit(1); }
const browser = await chromium.launch({ executablePath: '/usr/bin/google-chrome', headless: true, args: ['--no-sandbox', '--disable-gpu', '--disable-dev-shm-usage', '--disable-webgl', '--disable-webgl2'] });
const page = await (await browser.newContext()).newPage();
await page.goto(`http://127.0.0.1:${PORT}/`, { timeout: 90000 });
await page.waitForFunction(() => Shiny.shinyapp.$inputValues['app-dataspace-eset'] && document.querySelector('#app-contents').offsetParent !== null, null, { timeout: 90000 });
await page.evaluate(sel => { document.querySelector(sel).style.display = 'block'; }, HOOK.box);
await wait(2000);

const out = {};
// 1) axis_mode single write
out.axis_mode_before = await readInputs(page, 'feature_space-axisMode$');
await runHook(page, 'widgets', { patch: { 'dataspace.feature_space.axis_mode': 'custom' } });
await wait(8000);
out.axis_mode_after = await readInputs(page, 'feature_space-axisMode$');

// 2) attr4 color_analysis single-key write
out.color_before = await readInputs(page, 'feature_space-a4selector-selectColorUI');
await runHook(page, 'widgets', { patch: { 'dataspace.feature_space.attr4.color_analysis': 'General' } });
await wait(8000);
out.color_after_analysis_only = await readInputs(page, 'feature_space-a4selector-selectColorUI');
// full-triple transaction
await runHook(page, 'widgets', { patch: { 'dataspace.feature_space.attr4.color_analysis': 'General', 'dataspace.feature_space.attr4.color_subset': 'All', 'dataspace.feature_space.attr4.color_variable': 'Gene.name' } });
await wait(8000);
out.color_after_full_triple = await readInputs(page, 'feature_space-a4selector-selectColorUI');

// 3) generic table write: columns + page (feature table = tab_feature)
out.table_before = await readInputs(page, 'dataspace-tab_feature-table_state');
const colBefore = Object.keys(out.table_before || {}).length ? out.table_before : null;
await runHook(page, 'widgets', { patch: { 'dataspace.tab_feature.columns': ['featureId', 'General|All|Gene.name'] } });
await wait(8000);
out.table_after_columns = await readInputs(page, 'dataspace-tab_feature-table_state');
await runHook(page, 'widgets', { patch: { 'dataspace.tab_feature.page': 2 } });
await wait(8000);
out.table_after_page = await readInputs(page, 'dataspace-tab_feature-table_state');

// 4) survival censor (resultspace sample_general)
out.surv_before = await readInputs(page, 'sample_general_surv.*censor');
await runHook(page, 'widgets', { patch: { 'resultspace.sample_general.survival_censor': 500 } });
await wait(8000);
out.surv_after = await readInputs(page, 'sample_general_surv.*censor');

// 5) fgsea xax_analysis single-key
out.fg_before = await readInputs(page, 'resultspace-fgsea-xax');
await runHook(page, 'widgets', { patch: { 'resultspace.fgsea.xax_analysis': 'General' } });
await wait(8000);
out.fg_after = await readInputs(page, 'resultspace-fgsea-xax');

console.log(JSON.stringify(out, (k, v) => typeof v === 'string' ? v.slice(0, 90) : v, 1));

import fs from 'node:fs';
fs.writeFileSync(path.join(here, 'spotcheck_results.json'), JSON.stringify(out, null, 2));
await browser.close(); process.exit(0);
