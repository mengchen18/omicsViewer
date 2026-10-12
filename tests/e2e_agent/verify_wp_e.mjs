// WP-E verification: F5 dead-end chains and conditional widgets, live.
// Evidence per the WP-B/WP-B2/WP-D methodology: store-level assertions and
// DOM-level side effects (table rows computed and selected, slider value),
// never input-diff grades.
//   chain 1 (stringdb): STRING-connected features selected -> run via
//     run_control -> enrichment table POPULATES -> selected_row choices
//     MATERIALIZE -> set_widgets selects the row -> tr.selected in the DOM
//   chain 2 (batch): samples selected -> show toggle on -> panel visible ->
//     phenotype comparison computes -> batch_phenotype_selected_row choices
//     materialize -> set -> row selected
//   chain 3 (censor): Surv column set through the cascade -> survival view
//     renders the slider -> survival_censor write lands on the slider
// Also asserts the PRECONDITION notes are present in the capability
// records (the discovery surface).
// Run: node verify_wp_e.mjs   (from tests/e2e_agent/)
import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import path from 'node:path';
import fs from 'node:fs';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../..');
const EXTDATA = path.join(REPO, 'inst/extdata');
const PORT = 7788;
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
const rec = async (page, id) =>
  (await runHook(page, 'capabilities', { query: id, max_results: 3 }))
    .capabilities?.find(c => c.id === id) || null;
// scoped tab opening: the result-space navbar (id app-resultspace-analyst)
// materializes only after the dataset-level reactives settle, and the two
// navbars share tab labels (Sample!) - always scope + poll for existence
const openTabIn = async (page, navId, label, ms = 30000) => {
  const t0 = Date.now();
  let dropdownsOpened = false;
  while (Date.now() - t0 < ms) {
    const clicked = await page.evaluate(({ navId, label, openMenus }) => {
      const nav = document.getElementById(navId);
      if (!nav) return false;
      // tabs can hide in the navbar's overflow dropdown (additionalTabs):
      // open every dropdown menu once, then the anchor becomes visible
      if (openMenus) {
        nav.querySelectorAll('.dropdown-toggle, a[data-toggle="dropdown"]')
          .forEach(t => { try { t.click(); } catch (e) {} });
      }
      const a = Array.from(nav.querySelectorAll('a[data-toggle="tab"]'))
        .find(x => x.offsetParent !== null && new RegExp(label, 'i').test(x.innerText));
      if (a) { a.click(); return true; } return false;
    }, { navId, label, openMenus: !dropdownsOpened && Date.now() - t0 > 3000 });
    if (clicked) return true;
    dropdownsOpened = true;
    await wait(500);
  }
  return false;
};
const openTab = (page, label) => openTabIn(page, 'app-resultspace-analyst', label);
// dataTable info+selected under a container id regex
const dtState = (page, containerRe) => page.evaluate(re => {
  let best = null;
  for (const t of Array.from(document.querySelectorAll('table.dataTable'))) {
    let el = t;
    while (el && !(el.id && new RegExp(re).test(el.id))) el = el.parentElement;
    if (el && t.offsetParent !== null) {
      const rows = t.querySelectorAll('tbody tr').length;
      if (!best || rows > best.rows) {
        const w = t.closest('.dataTables_wrapper');
        best = { rows,
                 info: w?.querySelector('.dataTables_info')?.innerText || '',
                 selected: t.querySelectorAll('tbody tr.selected').length };
      }
    }
  }
  return best;
}, containerRe);
const wait = ms => new Promise(s => setTimeout(s, ms));

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
  await wait(1500);

  // ---------- precondition notes on the capability records ----------
  for (const id of ['resultspace.stringdb.selected_row',
                    'resultspace.sample_general.batch_phenotype_selected_row',
                    'resultspace.sample_general.batch_feature_selected_row',
                    'resultspace.sample_general.survival_censor']) {
    const rr = await rec(page, id);
    check(`capability record carries the PRECONDITION note: ${id}`,
      /PRECONDITION/.test(rr?.help_text || rr?.description || ''),
      (rr?.help_text || '').slice(0, 90));
  }

  // ---------- chain 1: stringdb run -> selected_row (criterion 1) ----------
  const GENES = ['X1452.EGFR', 'X1477.MYC', 'X1570.TP53', 'X2033.BRAF', 'X2371.CCND1'];
  await runHook(page, 'state', { features: GENES });
  await wait(1500);
  const ov1 = await runHook(page, 'overview', {});
  check('precondition: 5 STRING-connected features selected',
    ov1.selection?.features?.count === 5, JSON.stringify(ov1.selection?.features?.count));
  const pre = await rec(page, 'resultspace.stringdb.selected_row');
  check('precondition: selected_row choices empty before the run',
    !(pre?.allowed_values || []).length, JSON.stringify((pre?.allowed_values || []).length));
  const run1 = await runHook(page, 'runcontrol', { action_id: 'resultspace.stringdb.run' });
  check('run_control drove resultspace.stringdb.run', run1.ran === 'resultspace.stringdb.run',
    `epoch ${run1.epoch}${run1.error ? ' | ' + run1.error : ''}`);
  let enriched = null, terms = [];
  for (let i = 0; i < 45; i++) {
    await wait(2000);
    enriched = await dtState(page, 'resultspace-stringdb-strtab');
    const rr = await rec(page, 'resultspace.stringdb.selected_row');
    terms = rr?.allowed_values || [];
    if (terms.length && enriched && enriched.rows > 0) break;
  }
  check('selected_row choices MATERIALIZED (F5 unblocked)', terms.length > 0,
    `${terms.length} terms, first: ${JSON.stringify(terms[0])}`);
  check('open StringDB tab (result space, after the run)', await openTab(page, 'StringDB'));
  for (let i = 0; i < 20 && !(enriched && enriched.rows > 0); i++) {
    await wait(1000);
    enriched = await dtState(page, 'resultspace-stringdb-strtab');
  }
  check('enrichment table POPULATED after the run',
    !!(enriched && enriched.rows > 0), `rows: ${enriched?.rows} | info: ${enriched?.info}`);
  const set1 = await runHook(page, 'widgets', {
    patch: { 'resultspace.stringdb.selected_row': [terms[0]] } });
  let sel1 = null;
  for (let i = 0; i < 15; i++) {
    await wait(1000);
    sel1 = await dtState(page, 'resultspace-stringdb-strtab');
    if (sel1 && sel1.selected > 0) break;
  }
  const snap1 = await runHook(page, 'store', {});
  const rb1 = snap1.values?.['resultspace.stringdb.selected_row'];
  check('selected_row WRITE applies and the table row is SELECTED (end-to-end)',
    (set1.applied === 'resultspace.stringdb.selected_row' ||
      (set1.applied || []).includes?.('resultspace.stringdb.selected_row')) &&
      sel1?.selected > 0,
    `applied: ${JSON.stringify(set1.applied)} | tr.selected: ${sel1?.selected} | read-back: ${JSON.stringify(rb1)}`);

  // ---------- chain 2: batch phenotype comparison -> selected_row ----------
  check('open Sample tab (result space, scoped)', await openTab(page, /^\s*Sample\s*$/));
  await runHook(page, 'state', { samples: ['X786O_NCI60', 'A498_NCI60', 'A549_NCI60'] });
  await runHook(page, 'widgets', { patch: {
    'resultspace.sample_general.batch_show_phenotype': true } });
  await wait(2000);
  let btab = null, bterms = [];
  for (let i = 0; i < 30; i++) {
    await wait(2000);
    btab = await dtState(page, 'batch_comp-phenotype_table');
    const rr = await rec(page, 'resultspace.sample_general.batch_phenotype_selected_row');
    bterms = rr?.allowed_values || [];
    if (bterms.length && btab && btab.rows > 0) break;
  }
  check('batch phenotype table computes with preconditions met',
    !!(btab && btab.rows > 0), `rows: ${btab?.rows} | info: ${btab?.info}`);
  check('batch_phenotype_selected_row choices materialize', bterms.length > 0,
    `${bterms.length} rows, first: ${JSON.stringify(bterms[0])}`);
  if (bterms.length) {
    const set2 = await runHook(page, 'widgets', {
      patch: { 'resultspace.sample_general.batch_phenotype_selected_row': [bterms[0]] } });
    let sel2 = null;
    for (let i = 0; i < 15; i++) {
      await wait(1000);
      sel2 = await dtState(page, 'batch_comp-phenotype_table');
      if (sel2 && sel2.selected > 0) break;
    }
    check('batch_phenotype_selected_row write selects the row',
      (set2.applied === 'resultspace.sample_general.batch_phenotype_selected_row' ||
        (set2.applied || []).includes?.('resultspace.sample_general.batch_phenotype_selected_row')) &&
        sel2?.selected > 0,
      `applied: ${JSON.stringify(set2.applied)} | tr.selected: ${sel2?.selected}`);
  } else {
    check('batch_phenotype_selected_row write selects the row', false, 'no choices materialized');
  }

  // ---------- chain 3: survival_censor (criterion 2) ----------
  // set the x-axis cascade to the Surv column through set_widgets
  const set3 = await runHook(page, 'widgets', { patch: {
    'resultspace.sample_general.xax_analysis': 'Surv',
    'resultspace.sample_general.xax_subset': 'all',
    'resultspace.sample_general.xax_variable': 'OS' } });
  check('cascade write to the Surv column applies (3 keys)',
    (set3.applied || []).length === 3 || set3.applied === 'resultspace.sample_general.xax_variable',
    `applied: ${JSON.stringify(set3.applied)}`);
  let slider = null;
  for (let i = 0; i < 20; i++) {
    await wait(1000);
    slider = await page.evaluate(() => {
      const iv = Shiny.shinyapp.$inputValues || {};
      for (const k of Object.keys(iv))
        if (!k.startsWith('.clientdata') && /-censor$/.test(k))
          return { key: k, value: iv[k] };
      return null;
    });
    if (slider) break;
  }
  check('survival view renders the censor slider after the Surv column is set',
    !!slider, JSON.stringify(slider));
  const set4 = await runHook(page, 'widgets', { patch: {
    'resultspace.sample_general.survival_censor': 20 } });
  await wait(2000);
  const sliderAfter = await page.evaluate(() => {
    const iv = Shiny.shinyapp.$inputValues || {};
    for (const k of Object.keys(iv))
      if (!k.startsWith('.clientdata') && /-censor$/.test(k)) return iv[k];
    return null;
  });
  const snap4 = await runHook(page, 'store', {});
  const rb4 = snap4.values?.['resultspace.sample_general.survival_censor'];
  check('censor WRITE lands once the precondition is met (slider + store agree)',
    (set4.applied === 'resultspace.sample_general.survival_censor' ||
      (set4.applied || []).includes?.('resultspace.sample_general.survival_censor')) &&
      sliderAfter !== null && sliderAfter === rb4,
    `applied: ${JSON.stringify(set4.applied)} | slider: ${JSON.stringify(slider?.value)} -> ${JSON.stringify(sliderAfter)} | store: ${JSON.stringify(rb4)} (write 20 clamps into the rendered time range, R-M7)`);

  await browser.close();
} catch (e) {
  check('fatal', false, String(e));
}
fs.writeFileSync(path.join(here, 'widget_sweep', 'wp_e_results.json'),
  JSON.stringify({ generated: new Date().toISOString(), results }, null, 2));
const failed = results.filter(x => !x.pass).length;
console.log(`\nWP-E verification: ${results.length - failed}/${results.length} ok`);
process.exit(failed ? 1 : 0);
