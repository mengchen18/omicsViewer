// WP-A acceptance verification: dynamic-heatmap validation/discovery hole (F3).
// Checks, with the Dynamic-heatmap tab INACTIVE for phases 1-3:
//   1) DISCOVERY: capability records for the 5 dyn-heatmap bindings expose
//      non-NULL allowed_values (choices providers are activation-independent).
//   2) GARBAGE REJECTS: the exact F3 garbage writes (annot_col
//      ['GARBAGE|X|Y'], col_sort_by 'NOT_A_SORT_MODE', ...) are rejected
//      instead of applied.
//   3) VALID WRITES APPLY while inactive (store + browser input value).
//   4) After activating the tab, the written values are live in the UI.
// Run: node verify_wp_a_dyn_heatmap.mjs   (from tests/e2e_agent/)
import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import path from 'node:path';
import fs from 'node:fs';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../..');
const EXTDATA = path.join(REPO, 'inst/extdata');
const PORT = 7796;
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
const dynInputs = page => page.evaluate(() => {
  const iv = Shiny.shinyapp.$inputValues || {};
  const o = {}; for (const k of Object.keys(iv)) if (/dynheatmap/i.test(k)) o[k] = iv[k];
  return o;
});

const GARBAGE = {
  'dataspace.dyn_heatmap.annot_col': ['GARBAGE|X|Y'],
  'dataspace.dyn_heatmap.annot_row': ['GARBAGE|X|Y'],
  'dataspace.dyn_heatmap.tooltip_info': ['GARBAGE|X|Y'],
  'dataspace.dyn_heatmap.col_sort_by': 'NOT_A_SORT_MODE',
  'dataspace.dyn_heatmap.row_sort_by': 'NOT_A_SORT_MODE'
};
const VALID = {
  'dataspace.dyn_heatmap.annot_col': ['General|All|Origin'],
  'dataspace.dyn_heatmap.annot_row': ['General|All|Gene.name'],
  'dataspace.dyn_heatmap.tooltip_info': ['General|All|Gene.name'],
  'dataspace.dyn_heatmap.col_sort_by': 'none',
  'dataspace.dyn_heatmap.row_sort_by': 'General|All|Gene.name'
};
const INPUT_OF = {
  'dataspace.dyn_heatmap.annot_col': 'app-dataspace-dynheatmapViewer-annotCol',
  'dataspace.dyn_heatmap.annot_row': 'app-dataspace-dynheatmapViewer-annotRow',
  'dataspace.dyn_heatmap.tooltip_info': 'app-dataspace-dynheatmapViewer-tooltipInfo',
  'dataspace.dyn_heatmap.col_sort_by': 'app-dataspace-dynheatmapViewer-colSortBy',
  'dataspace.dyn_heatmap.row_sort_by': 'app-dataspace-dynheatmapViewer-rowSortBy'
};

const results = [];
const check = (name, pass, detail = '') => {
  results.push({ name, pass, detail });
  console.log(`${pass ? 'ok  ' : 'FAIL'} ${name}${detail ? ' -- ' + detail : ''}`);
};

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
  // sanity: Dynamic heatmap tab really is inactive
  const inactive = await page.evaluate(() => {
    const a = Array.from(document.querySelectorAll('a[data-toggle="tab"]'))
      .find(x => /Dynamic heatmap/.test(x.innerText));
    return !(a && a.offsetParent !== null && a.classList.contains('active'));
  });
  check('precondition: Dynamic-heatmap tab inactive at start', inactive === true);

  // 1) discovery: allowed_values non-NULL while inactive
  for (const id of Object.keys(GARBAGE)) {
    const rec = (await runHook(page, 'capabilities', { query: id, max_results: 5 }))
      .capabilities?.find(c => c.id === id);
    const n = Array.isArray(rec?.allowed_values) ? rec.allowed_values.length : 0;
    check(`discovery (inactive): ${id} exposes allowed values`, n > 0,
      n ? `${n} values, e.g. ${JSON.stringify(rec.allowed_values.slice(0, 3))}` : 'NULL/empty');
  }

  // 2) garbage rejects while inactive
  for (const [id, value] of Object.entries(GARBAGE)) {
    const receipt = await runHook(page, 'widgets', { patch: { [id]: value } });
    const rejected = (receipt.rejected || []).some(rx => rx.id === id.split('.').pop() || rx.id === id);
    const applied = (receipt.applied || []).includes(id);
    const reason = (receipt.rejected || [])[0]?.reason || '';
    check(`garbage rejects (inactive): ${id}`, !applied && (rejected || reason),
      applied ? 'WAS APPLIED (F3 hole still open)' : reason.slice(0, 110));
  }

  // 3) valid writes apply while inactive AND move the browser input
  for (const [id, value] of Object.entries(VALID)) {
    const receipt = await runHook(page, 'widgets', { patch: { [id]: value } });
    const applied = (receipt.applied || []).includes(id);
    await page.waitForTimeout(600);
    const iv = await dynInputs(page);
    const inputId = INPUT_OF[id];
    let uiOk = false;
    if (Array.isArray(value)) {
      const cur = iv[inputId];
      uiOk = Array.isArray(cur) && value.every(v => cur.includes(v));
    } else {
      uiOk = iv[inputId] === value;
    }
    check(`valid write applies (inactive): ${id}`, applied,
      applied ? '' : JSON.stringify(receipt).slice(0, 110));
    check(`valid write reaches UI input (inactive): ${inputId}`, uiOk,
      `input now: ${JSON.stringify(iv[inputId])?.slice(0, 60)}`);
  }

  // 4) activate the tab; the values must be live
  await page.evaluate(() => {
    const a = Array.from(document.querySelectorAll('a[data-toggle="tab"]'))
      .find(x => x.offsetParent !== null && /Dynamic heatmap/.test(x.innerText));
    if (a) a.click();
  });
  await page.waitForTimeout(5000);
  const ivActive = await dynInputs(page);
  for (const [id, value] of Object.entries(VALID)) {
    const inputId = INPUT_OF[id];
    let live = false;
    if (Array.isArray(value)) {
      const cur = ivActive[inputId];
      live = Array.isArray(cur) && value.every(v => cur.includes(v));
    } else {
      live = ivActive[inputId] === value;
    }
    check(`value live after activation: ${inputId}`, live,
      `input now: ${JSON.stringify(ivActive[inputId])?.slice(0, 60)}`);
  }

  await browser.close();
} catch (e) {
  check('fatal', false, String(e));
}
fs.writeFileSync(path.join(here, 'wp_a_results.json'), JSON.stringify({ generated: new Date().toISOString(), results }, null, 2));
const failed = results.filter(x => !x.pass).length;
console.log(`\nWP-A verification: ${results.length - failed}/${results.length} ok`);
process.exit(failed ? 1 : 0);
