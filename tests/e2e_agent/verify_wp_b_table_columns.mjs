// WP-B verification: table `columns` write semantics (F4), end to end.
// Tests what the sweep harness cannot observe (columns writes act through
// a table RE-RENDER - an output-level effect - not through an input value):
//   1) INITIALIZED table: array write, 1-element write (the sweep's probe
//      shape), read-back equality, DOM header set, pending/ack echo
//      (store must still hold the written set after the echo settles),
//      page push AFTER a columns-driven re-render, invalid value rejects.
//   2) UNINITIALIZED table: write while the tab was never opened applies
//      in the store; opening the tab renders the stored column set (the
//      deferred-effect design that the sweep records as no-ui-effect).
// Run: node verify_wp_b_table_columns.mjs   (from tests/e2e_agent/)
import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import path from 'node:path';
import fs from 'node:fs';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../..');
const EXTDATA = path.join(REPO, 'inst/extdata');
const PORT = 7797;
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
  }, { sel: HOOK.result, n, before }, { timeout: 60000 });
  return JSON.parse((await page.innerText(HOOK.result)).trim());
}
const storeValue = async (page, id) =>
  (await runHook(page, 'store', {})).values?.[id] ?? null;

const openTab = (page, label) => page.evaluate(label => {
  const a = Array.from(document.querySelectorAll('a[data-toggle="tab"]'))
    .find(x => x.offsetParent !== null && new RegExp(label).test(x.innerText));
  if (a) a.click();
  return !!a;
}, label);
// visible column headers of the dataTable nested under an ancestor id matching re
const headers = (page, containerRe) => page.evaluate(re => {
  for (const t of Array.from(document.querySelectorAll('table.dataTable'))) {
    let el = t;
    while (el && !(el.id && new RegExp(re).test(el.id))) el = el.parentElement;
    const inHost = !!el;
    const visible = t.offsetParent !== null;
    const th = Array.from(t.querySelectorAll('thead th')).map(x => x.innerText.trim())
      .filter(x => x.length && x !== '\u200b');
    if (inHost && visible && th.length) return th;
  }
  return null;
}, containerRe);
const tableState = (page, re) => page.evaluate(re => {
  const iv = Shiny.shinyapp.$inputValues || {};
  for (const k of Object.keys(iv)) if (new RegExp(re).test(k)) return iv[k];
  return null;
}, re);
const wait = ms => new Promise(s => setTimeout(s, ms));

const results = [];
const check = (name, pass, detail = '') => {
  results.push({ name, pass, detail });
  console.log(`${pass ? 'ok  ' : 'FAIL'} ${name}${detail ? ' -- ' + detail : ''}`);
};
const eqSet = (a, b) => Array.isArray(a) && Array.isArray(b) &&
  a.length === b.length && [...a].sort().join('\u0001') === [...b].sort().join('\u0001');

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
  await wait(1500);

  // ---------- part 1: initialized table (Sample table = tab_pheno) ----------
  check('open Sample table tab', await openTab(page, 'Sample table'));
  await wait(2500);
  const h0 = await headers(page, 'dataspace-tab_pheno');
  check('table rendered with headers', Array.isArray(h0) && h0.length >= 3, `${h0?.length} headers`);

  // 1a. two-element array write
  const rc1 = await runHook(page, 'widgets', { patch: {
    'dataspace.tab_pheno.columns': ['General|All|Cell.line', 'General|All|MDR'] } });
  await wait(2500);
  const rb1 = await storeValue(page, 'dataspace.tab_pheno.columns');
  const h1 = await headers(page, 'dataspace-tab_pheno');
  check('2-element write applied', (rc1.applied || []).includes('dataspace.tab_pheno.columns'),
    JSON.stringify(rc1.rejected || []).slice(0, 90));
  check('read-back equals written set (after echo settles)', eqSet(rb1, ['General|All|Cell.line', 'General|All|MDR']),
    `store: ${JSON.stringify(rb1)}`);
  check('DOM headers equal written set', eqSet(h1, ['General|All|Cell.line', 'General|All|MDR']),
    `headers: ${JSON.stringify(h1)}`);

  // 1b. one-element write (the sweep probe shape - displays as scalar in receipts)
  const rc2 = await runHook(page, 'widgets', { patch: {
    'dataspace.tab_pheno.columns': ['General|All|Origin'] } });
  await wait(2500);
  const rb2 = await storeValue(page, 'dataspace.tab_pheno.columns');
  const h2 = await headers(page, 'dataspace-tab_pheno');
  check('1-element write applied', (rc2.applied || []).includes('dataspace.tab_pheno.columns'));
  check('read-back equals written value (1 element)',
    eqSet(rb2, ['General|All|Origin']) || rb2 === 'General|All|Origin',
    `store: ${JSON.stringify(rb2)} (JSON scalar display of a length-1 R vector is expected)`);
  check('DOM headers equal written value (1 element)', eqSet(h2, ['General|All|Origin']),
    `headers: ${JSON.stringify(h2)}`);

  // 1c. page push AFTER a columns-driven re-render (deferred proxy interplay)
  const rc3 = await runHook(page, 'widgets', { patch: { 'dataspace.tab_pheno.page': 2 } });
  await wait(2500);
  const st3 = await tableState(page, 'dataspace-tab_pheno-table_state');
  check('page push lands after columns re-render', st3?.start === 25,
    `start: ${st3?.start}`);

  // 1d. invalid value rejects
  const rc4 = await runHook(page, 'widgets', { patch: {
    'dataspace.tab_pheno.columns': ['NOPE|X|Y'] } });
  await wait(1000);
  const rb4 = await storeValue(page, 'dataspace.tab_pheno.columns');
  check('invalid column rejects', (rc4.rejected || []).length > 0 && !(rc4.applied || []).length,
    (rc4.rejected || [])[0]?.reason?.slice(0, 90));
  check('store unchanged after rejection',
    eqSet(rb4, ['General|All|Origin']) || rb4 === 'General|All|Origin',
    `store: ${JSON.stringify(rb4)}`);

  // ---------- part 2: uninitialized table (Expression = tab_expr) ----------
  const preH = await headers(page, 'dataspace-tab_expr');
  check('precondition: Expression table not rendered yet', preH === null);
  const rc5 = await runHook(page, 'widgets', { patch: {
    'dataspace.tab_expr.columns': ['feature', 'X786O_NCI60'] } });
  await wait(1500);
  const rb5 = await storeValue(page, 'dataspace.tab_expr.columns');
  check('write to uninitialized table applied in store',
    (rc5.applied || []).includes('dataspace.tab_expr.columns') && eqSet(rb5, ['feature', 'X786O_NCI60']),
    `store: ${JSON.stringify(rb5)}`);
  check('open Expression tab', await openTab(page, 'Expression'));
  await wait(3500);
  const h5 = await headers(page, 'dataspace-tab_expr');
  check('stored column set renders when the tab opens (deferred effect)',
    eqSet(h5, ['feature', 'X786O_NCI60']), `headers: ${JSON.stringify(h5)}`);
  const rb5b = await storeValue(page, 'dataspace.tab_expr.columns');
  check('store still equals written set after render echo', eqSet(rb5b, ['feature', 'X786O_NCI60']),
    `store: ${JSON.stringify(rb5b)}`);

  await browser.close();
} catch (e) {
  check('fatal', false, String(e));
}
fs.writeFileSync(path.join(here, 'wp_b_results.json'),
  JSON.stringify({ generated: new Date().toISOString(), results }, null, 2));
const failed = results.filter(x => !x.pass).length;
console.log(`\nWP-B verification: ${results.length - failed}/${results.length} ok`);
process.exit(failed ? 1 : 0);
