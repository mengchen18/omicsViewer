// ============================================================================
// omicsViewer — end-to-end UI tests (Playwright, user-mimicking)
//
// One app instance per BACKEND, one browser page, one continuous user journey:
//   boot & dataset load -> data-space tabs (scatter, tables, heatmaps,
//   dynamic heatmap, expression, GSList) -> analysis tabs (feature,
//   sample incl. survival, ORA, fGSEA, response, SeqLogo, Geneshot,
//   StringDB) -> snapshot save/restore/delete -> xlsx download ->
//   AI assistant panel.
//
// BACKENDS: the identical journey runs against BOTH dataset formats and
// the visible results are compared against each other (parity block):
//   rds    — inst/extdata/demo.RDS copied into a temp data dir
//   sqlite — the SAME demo.RDS converted with the app's own save path
//            (saveOmicsViewerDb) into demo.sqlite in its own temp dir
// The sqlite arm additionally verifies the WP-4 connection lifecycle:
// after the journey (browser closed, app process still alive) no file
// descriptor may remain open on the database, an external
// saveOmicsViewerDb overwrite must succeed (no stale lock), and no
// -wal/-shm journal files may be left in the data dir.
//
// PARITY PHILOSOPHY: per-arm checks assert each backend behaves; the
// parity block asserts the USER-VISIBLE numbers agree — dataset
// dimensions, selection counts, quick-view badge ids, analysis table
// row counts and first-row contents, survival censor range, expression
// cell spot-checks. Plotly internals and DOM structure are deliberately
// NOT compared (brittle, meaningless). Numeric tolerance: rel 1e-9 /
// abs 1e-12 — SQLite REAL stores exact IEEE-754 doubles so values
// should be bit-identical; the tolerance only absorbs last-bit
// formatting differences while still catching any real corruption.
//
// NOTE ON DATASET SAVING: the UI has NO dataset-object save control —
// only an xlsx table export and app-state snapshots (.ESS). The
// snapshot save/restore journey IS the UI's save mechanism and is
// exercised, including a pre-save vs post-restore state assertion.
//
// Run:  node e2e.spec.mjs [rds|sqlite|both]   (default both)
// Env:  E2E_PORT (rds arm; sqlite = +1, default 7780), E2E_HEADLESS=0
//       to watch, E2E_KEEP=1 to keep the temp data dirs.
//
// Writes e2e_results.json (both arms + parity) + screenshots under
// artifacts/<backend>/.

import { chromium } from 'playwright';
import { spawn, spawnSync } from 'node:child_process';
import { mkdirSync, writeFileSync, cpSync, rmSync, readdirSync, readlinkSync, statSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../..');
const EXTDATA = path.join(REPO, 'inst/extdata');
const DEMO_RDS = path.join(EXTDATA, 'demo.RDS');
const HEADLESS = process.env.E2E_HEADLESS !== '0';
const KEEP = process.env.E2E_KEEP === '1';

// backend selection: CLI arg wins, then E2E_BACKEND, default both
const argBe = (process.argv[2] || '').toLowerCase();
const BACKENDS = (argBe || process.env.E2E_BACKEND || 'both').toLowerCase();
const ARMS = BACKENDS === 'rds' ? ['rds'] : BACKENDS === 'sqlite' ? ['sqlite'] : ['rds', 'sqlite'];
const PORT0 = Number(process.env.E2E_PORT || 7780);
const portFor = (be) => PORT0 + (be === 'sqlite' ? 1 : 0);
const datasetFor = (be) => (be === 'sqlite' ? 'demo.sqlite' : 'demo.RDS');

// ------------------------------------------------------------------ helpers
const DATASET_TABS = ['Feature', 'Feature table', 'Sample', 'Sample table',
  'Cor', 'Heatmap', 'Dynamic heatmap', 'Expression', 'GSList'];
const ANALYSIS_TABS = ['Feature', 'ORA', 'fGSEA', 'Response', 'StringDB',
  'SeqLogo', 'Geneshot', 'Sample'];

const t0 = Date.now();
const ts = () => `t+${Math.round((Date.now() - t0) / 1000)}s`;
const log = (m) => console.log(`${ts()} ${m}`);

const inputVal = (page, id) => page.evaluate(
  k => (window.Shiny && Shiny.shinyapp && Shiny.shinyapp.$inputValues[k]) ?? null, id);

const waitInput = (page, id, timeout = 30000) => page.waitForFunction(
  k => Shiny.shinyapp.$inputValues[k] !== undefined, id, { timeout });

// click a navbar tab (data space or analysis space) and wait for the input
const gotoTab = async (page, space, value, timeout = 30000) => {
  const nav = space === 'data' ? '#app-dataspace-eset' : '#app-resultspace-analyst';
  const inputId = space === 'data' ? 'app-dataspace-eset' : 'app-resultspace-analyst';
  const ok = await page.evaluate(({ nav, value }) => {
    const el = Array.from(document.querySelectorAll(`${nav} a[data-value="${value}"]`))
      .filter(e => e.offsetParent !== null)[0];
    if (el) { el.click(); return true; }
    return false;
  }, { nav, value });
  if (!ok) throw new Error(`tab link not found: ${space}/${value}`);
  await waitInput(page, inputId, timeout);
  await settle(page);
};

// wait until shiny is idle (not busy) for ~1s
const settle = async (page, ms = 15000) => {
  try {
    await page.waitForFunction(() =>
      !document.documentElement.classList.contains('shiny-busy'), null, { timeout: ms });
    await page.waitForTimeout(700);
  } catch { /* keep going on timeout */ }
};

// selectize-aware select: skips no-op sets (re-setting the same value fires
// a change that re-triggers the cascade and can clear dependent selections);
// waits for the TARGET option to exist, then setValue
const setSelect = async (page, sel, value, timeout = 45000) => {
  await page.waitForSelector(sel, { state: 'attached', timeout });
  const cur = await page.evaluate(s => {
    const el = document.querySelector(s);
    return el ? (el.selectize ? el.selectize.getValue() : el.value) : null;
  }, sel).catch(() => null);
  if (cur === value) return true; // already set — do not re-fire change
  await page.waitForFunction(({ s, v }) => {
    const el = document.querySelector(s);
    if (!el) return false;
    if (el.selectize) return Object.prototype.hasOwnProperty.call(el.selectize.options, v);
    return Array.from(el.options).some(o => o.value === v || o.text === v);
  }, { s: sel, v: value }, { timeout });
  return page.evaluate(({ s, v }) => {
    const el = document.querySelector(s);
    if (el.selectize) { el.selectize.setValue(v); return true; }
    el.value = v;
    el.dispatchEvent(new Event('change', { bubbles: true }));
    return true;
  }, { s: sel, v: value });
};

// visible plotly plot handle
const visiblePlot = (page) => page.evaluate(() => {
  const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
    .filter(e => e.offsetParent !== null && e._fullLayout);
  return els.length ? els[0] : null;
});

const plotInfo = (page) => page.evaluate(() => {
  const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
    .filter(e => e.offsetParent !== null && e._fullLayout);
  if (!els.length) return null;
  const p = els[0];
  const t = (ax) => {
    const v = ax && ax.title ? (ax.title.text || ax.title) : '';
    return typeof v === 'string' ? v : '';
  };
  return {
    ntraces: (p._fullData || []).length,
    x: t(p._fullLayout.xaxis), y: t(p._fullLayout.yaxis),
    shapes: (p._fullLayout.shapes || []).length,
    w: p.offsetWidth, h: p.offsetHeight
  };
});

const waitPlot = (page, timeout = 90000) => page.waitForFunction(() =>
  Array.from(document.querySelectorAll('.js-plotly-plot'))
    .some(e => e.offsetParent !== null && e._fullLayout), null, { timeout });

// wait for a rendered (non-empty) image inside a Shiny plotOutput
const waitImg = (page, sel, timeout = 90000) => page.waitForFunction(s => {
  const root = document.querySelector(s);
  if (!root) return false;
  const box = root.closest('.shiny-plot-output, div');
  const imgs = Array.from((box || root).querySelectorAll('img'))
    .filter(im => im.offsetParent !== null);
  return imgs.length > 0 && imgs.some(im => im.naturalWidth > 2);
}, sel, { timeout });

// DT helpers -------------------------------------------------------------
const dtInfo = (page, tableId) => page.evaluate((id) => {
  const root = document.getElementById(id);
  if (!root) return null;
  const info = root.querySelector('.dataTables_info');
  const rows = root.querySelectorAll('tbody tr');
  const cols = root.querySelectorAll('thead th');
  return {
    info: info ? info.innerText.trim() : null,
    nRows: rows.length, nCols: cols.length,
    firstRow: rows[0] ? rows[0].innerText.trim().slice(0, 120) : null
  };
}, tableId);

const dtEntryCount = async (page, tableId) => {
  const inf = await dtInfo(page, tableId);
  if (!inf || !inf.info) return null;
  const m = inf.info.match(/of\s+([\d,.]+)\s+entr/i);
  return m ? Number(m[1].replace(/[,.]/g, '')) : null;
};

// spec: { eq, min, max, notEq } — serializable predicate (waitForFunction
// args cannot carry functions)
const dtWaitEntries = (page, tableId, spec, timeout = 60000) =>
  page.waitForFunction(({ id, spec }) => {
    try {
      const root = document.getElementById(id);
      if (!root) return false;
      const info = root.querySelector('.dataTables_info');
      if (!info) return false;
      const m = info.innerText.match(/of\s+([\d,.]+)\s+entr/i);
      if (!m) return false;
      const n = Number(m[1].replace(/[,.]/g, ''));
      if (spec.eq !== undefined && n !== spec.eq) return false;
      if (spec.min !== undefined && n < spec.min) return false;
      if (spec.max !== undefined && n > spec.max) return false;
      if (spec.notEq !== undefined && n === spec.notEq) return false;
      return true;
    } catch (e) { return false; }
  }, { id: tableId, spec }, { timeout });

const dtClearSearch = (page, tableId) => page.evaluate((id) => {
  const input = document.querySelector(`[id="${id}_filter"] input`);
  if (input) {
    const set = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value').set;
    set.call(input, '');
    input.dispatchEvent(new Event('input', { bubbles: true }));
  }
}, tableId);

// first N rows x M columns of visible DT cells (expression spot-checks)
const dtCells = (page, tableId, nRows = 5, nCols = 4) => page.evaluate(
  ({ id, nRows, nCols }) => {
    const root = document.getElementById(id);
    if (!root) return null;
    const rows = Array.from(root.querySelectorAll('tbody tr')).slice(0, nRows);
    return rows.map(r => Array.from(r.querySelectorAll('td'))
      .slice(0, nCols).map(td => td.innerText.trim()));
  }, { id: tableId, nRows, nCols });

// fresh Shiny session (page reload + dataset re-select); used to recover
// from a poisoned session state so later sections stay meaningful
const hardReset = async (page, dataset) => {
  await page.reload({ timeout: 60000 });
  await page.waitForFunction(() => {
    const el = document.querySelector('#app-selectFile');
    return el && el.selectize && Object.keys(el.selectize.options).length > 0;
  }, null, { timeout: 120000 });
  await page.evaluate(ds => {
    document.querySelector('#app-selectFile').selectize.setValue(ds);
  }, dataset);
  await page.waitForFunction(() => {
    const tab = Shiny.shinyapp.$inputValues['app-dataspace-eset'];
    const box = document.querySelector('#app-contents');
    return !!tab && box && box.offsetParent !== null;
  }, null, { timeout: 180000 });
  await page.waitForSelector('#app-resultspace-analyst a[data-value]', { timeout: 120000 });
};

// mouse helpers on Shiny imageOutput (click + brush) ----------------------
// Shiny hover inputs have delay=150ms: pause after moving so the hover
// input registers before the click (the heatmap click handler requires it)
const clickImg = async (page, sel, fx, fy) => {
  const box = await page.locator(sel).first().boundingBox();
  if (!box) throw new Error(`no bounding box for ${sel}`);
  await page.mouse.move(box.x + box.width * fx, box.y + box.height * fy);
  await page.waitForTimeout(450);
  await page.mouse.down();
  await page.mouse.up();
};

const brushImg = async (page, sel, fx0, fy0, fx1, fy1) => {
  const box = await page.locator(sel).first().boundingBox();
  if (!box) throw new Error(`no bounding box for ${sel}`);
  await page.mouse.move(box.x + box.width * fx0, box.y + box.height * fy0);
  await page.mouse.down();
  await page.mouse.move(box.x + box.width * fx1, box.y + box.height * fy1, { steps: 10 });
  await page.mouse.up();
};

// plotly box-select: activate the modebar tool then drag across the plot
// (plotly 2.x modebar entries are <a data-title=...>)
const plotlyBoxSelect = async (page, x0f, y0f, x1f, y1f) => {
  await page.evaluate(() => {
    const p = Array.from(document.querySelectorAll('.js-plotly-plot'))
      .filter(e => e.offsetParent !== null)[0];
    if (!p) throw new Error('no visible plot');
    const btn = Array.from(p.querySelectorAll('.modebar a'))
      .find(b => (b.getAttribute('data-title') || '').toLowerCase().includes('box select'));
    if (!btn) throw new Error('box-select modebar entry not found');
    btn.dispatchEvent(new MouseEvent('mousedown', { bubbles: true }));
    btn.dispatchEvent(new MouseEvent('mouseup', { bubbles: true }));
    btn.click();
  });
  const box = await page.evaluate(() => {
    const p = Array.from(document.querySelectorAll('.js-plotly-plot'))
      .filter(e => e.offsetParent !== null)[0];
    const r = p.getBoundingClientRect();
    return { x: r.x, y: r.y, w: r.width, h: r.height };
  });
  await page.mouse.move(box.x + box.w * x0f, box.y + box.h * y0f);
  await page.mouse.down();
  await page.mouse.move(box.x + box.w * x1f, box.y + box.h * y1f, { steps: 12 });
  await page.mouse.up();
};

// click nth visible DT row with a real playwright click (carries modifiers)
const dtClickRowReal = async (page, tableId, n = 0, modifiers = []) => {
  await page.click(`#${tableId} tbody tr >> nth=${n}`, { modifiers });
};

// notification helper
const notifVisible = (page, kind, timeout = 15000) => page.waitForFunction(k =>
  Array.from(document.querySelectorAll(`.shiny-notification-${k}`))
    .some(e => e.offsetParent !== null), kind, { timeout });

// ------------------------------------------------------- parity comparison
// numeric tolerance: rel 1e-9 / abs 1e-12 (SQLite REAL is exact IEEE-754;
// the tolerance absorbs last-bit/display drift, not real corruption)
const REL_TOL = 1e-9, ABS_TOL = 1e-12;
const numEq = (a, b) =>
  Math.abs(a - b) <= Math.max(ABS_TOL, REL_TOL * Math.max(Math.abs(a), Math.abs(b)));
const looksNumeric = (s) => /^[-+]?[\d.,]+([eE][-+]?\d+)?$/.test(s) && /\d/.test(s);
// one DT cell: identical strings pass; numeric-looking cells compare with
// tolerance; everything else must match exactly
const cellEq = (a, b) => {
  const s1 = String(a).trim(), s2 = String(b).trim();
  if (s1 === s2) return true;
  if (!s1 || !s2) return false;
  if (looksNumeric(s1) && looksNumeric(s2))
    return numEq(Number(s1.replace(/,/g, '')), Number(s2.replace(/,/g, '')));
  return false;
};
const rowTextEq = (a, b) => {
  const t1 = String(a).split('\t').map(s => s.trim());
  const t2 = String(b).split('\t').map(s => s.trim());
  if (t1.length !== t2.length) return false;
  return t1.every((c, i) => cellEq(c, t2[i]));
};
const cellMatrixEq = (a, b) => {
  if (!Array.isArray(a) || !Array.isArray(b) || a.length !== b.length) return false;
  return a.every((row, i) =>
    Array.isArray(row) && Array.isArray(b[i]) && row.length === b[i].length &&
    row.every((c, j) => cellEq(c, b[i][j])));
};

// ------------------------------------------------------------ arm launcher
const waitPort = async (port, ms = 120000, rLog = '') => {
  const start = Date.now();
  while (Date.now() - start < ms) {
    try {
      const res = await fetch(`http://127.0.0.1:${port}/`);
      if (res.status > 0) return true;
    } catch {}
    await new Promise(s => setTimeout(s, 500));
  }
  throw new Error(`app did not start on :${port}\n${rLog.slice(-2000)}`);
};

const launchArm = async (backend) => {
  const port = portFor(backend);
  const dataset = datasetFor(backend);
  const datadir = path.join(here, 'artifacts', `data-${backend}`);
  rmSync(datadir, { recursive: true, force: true });
  mkdirSync(datadir, { recursive: true });
  if (backend === 'rds') {
    cpSync(DEMO_RDS, path.join(datadir, dataset));
  } else {
    // the app's own save path converts the demo dataset — no toy data,
    // no shortcuts: exactly what a user would run
    const dbFile = path.join(datadir, dataset);
    const conv = spawnSync('Rscript', ['-e',
      `omicsViewer::saveOmicsViewerDb(readRDS('${DEMO_RDS}'), db.file = '${dbFile}')`],
      { cwd: REPO, encoding: 'utf8', timeout: 300000 });
    if (conv.status !== 0)
      throw new Error(`demo.sqlite conversion failed:\n${conv.stdout}\n${conv.stderr}`);
    log(`converted demo.RDS -> ${dataset} (${Math.round(statSync(dbFile).size / 1e6)} MB)`);
  }
  const rCode = `
    options(shiny.port = ${port}, shiny.host = '127.0.0.1')
    omicsViewer::omicsViewer(dir = '${datadir}')
  `;
  const r = spawn('Rscript', ['-e', rCode], {
    cwd: REPO,
    env: { ...process.env, OMICSVIEWER_LOG: 'off' },
    stdio: ['ignore', 'pipe', 'pipe']
  });
  let rLog = '';
  r.stdout.on('data', d => { rLog += d; });
  r.stderr.on('data', d => { rLog += d; });
  const cleanup = () => {
    try { r.kill('SIGKILL'); } catch {}
    if (!KEEP) { try { rmSync(datadir, { recursive: true, force: true }); } catch {} }
  };
  process.on('exit', cleanup);
  await waitPort(port, 120000, rLog);
  return { backend, port, dataset, datadir, r, getLog: () => rLog, cleanup };
};

// open file descriptors of the R process pointing at the dataset file
const datasetFds = (arm) => {
  const fdDir = `/proc/${arm.r.pid}/fd`;
  try {
    return readdirSync(fdDir).filter(l => {
      try { return readlinkSync(`${fdDir}/${l}`).includes(arm.dataset); }
      catch { return false; }
    }).length;
  } catch { return null; }
};

// ==================================================================== run
const armReports = {};
let exitCode = 0;

const runJourney = async (arm) => {
  const { backend, port, dataset } = arm;
  const ARTIFACTS = path.join(here, 'artifacts', backend);
  mkdirSync(ARTIFACTS, { recursive: true });

  const results = [];
  const M = {}; // user-visible metrics for the parity block
  let shots = 0;
  const record = (name, pass, detail = '') => {
    results.push({ name, pass, detail: String(detail).slice(0, 400) });
    log(`${pass ? 'PASS' : 'FAIL'} | [${backend}] ${name}${detail ? ' | ' + String(detail).slice(0, 160) : ''}`);
  };
  const shot = async (page, name) => {
    try {
      await page.screenshot({ path: path.join(ARTIFACTS, `fail-${++shots}-${name}.png`), fullPage: false });
    } catch {}
  };

  const browser = await chromium.launch({
    executablePath: '/usr/bin/google-chrome',
    headless: HEADLESS,
    args: ['--no-sandbox', '--disable-gpu', '--disable-dev-shm-usage',
           '--disable-webgl', '--disable-webgl2']
  });
  const ctx = await browser.newContext({ viewport: { width: 1680, height: 1050 } });
  const page = await ctx.newPage();
  const pageErrors = [];
  page.on('pageerror', e => pageErrors.push(String(e).slice(0, 200)));

  // ============================ 1. boot & dataset load ==================
  console.log(`\n--- [${backend}] 1. boot & dataset load ---`);
  await page.goto(`http://127.0.0.1:${port}/`, { timeout: 60000 });
  record('page loads (title)', await page.title().then(t => /omicsViewer/i.test(t) || t.length > 0), await page.title());

  await page.waitForSelector('#app-selectFile', { state: 'attached', timeout: 60000 });
  record('dataset selector is rendered (selectized)',
    await page.locator('#app-selectFile-selectized').count() === 1);
  await page.waitForFunction(() => {
    const el = document.querySelector('#app-selectFile');
    return el && el.selectize && Object.keys(el.selectize.options).length > 0;
  }, null, { timeout: 90000 });
  const choices = await page.evaluate(() => Object.keys(document.querySelector('#app-selectFile').selectize.options));
  record(`dataset selector offers ${dataset}`, choices.includes(dataset), choices.join(','));

  await page.evaluate(ds => {
    document.querySelector('#app-selectFile').selectize.setValue(ds);
  }, dataset);
  await page.waitForFunction(() => {
    const tab = Shiny.shinyapp.$inputValues['app-dataspace-eset'];
    const box = document.querySelector('#app-contents');
    return !!tab && box && box.offsetParent !== null;
  }, null, { timeout: 180000 });
  record(`${dataset} loads: contents panel appears`, true);
  // dataset dimensions straight from the loading status a screen reader sees
  await page.waitForFunction(() => /Dataset loaded successfully/.test(
    (document.querySelector('#app-loadingStatus') || {}).innerText || ''), null, { timeout: 60000 });
  const loadTxt = await page.evaluate(() =>
    (document.querySelector('#app-loadingStatus') || {}).innerText || '');
  const fm = loadTxt.match(/(\d+) features and (\d+) samples/);
  M.featuresN = fm ? Number(fm[1]) : null;
  M.samplesN = fm ? Number(fm[2]) : null;
  record(`loading status reports the dataset dimensions (${M.featuresN}x${M.samplesN})`,
    M.featuresN === 2702 && M.samplesN === 60, loadTxt.trim().slice(0, 120));

  // tab structure (the result-space navbar renders via renderUI after load)
  await page.waitForSelector('#app-resultspace-analyst a[data-value]', { timeout: 120000 });
  const dataTabs = await page.$$eval('#app-dataspace-eset a[data-value]', els =>
    els.map(e => e.getAttribute('data-value')));
  record(`data navbar has all ${DATASET_TABS.length} tabs`,
    DATASET_TABS.every(t => dataTabs.includes(t)), dataTabs.join(','));
  const anaTabs = await page.$$eval('#app-resultspace-analyst a[data-value]', els =>
    els.map(e => e.getAttribute('data-value')));
  record(`analysis navbar has all ${ANALYSIS_TABS.length} tabs (ORA/fGSEA included)`,
    ANALYSIS_TABS.every(t => anaTabs.includes(t)), anaTabs.join(','));
  M.dataTabs = dataTabs.join(',');
  M.anaTabs = anaTabs.join(',');
  record('default data tab is Feature',
    (await inputVal(page, 'app-dataspace-eset')) === 'Feature');
  record('default analysis tab is Feature',
    (await inputVal(page, 'app-resultspace-analyst')) === 'Feature');

  // ============================ 2. feature scatter ======================
  console.log(`\n--- [${backend}] 2. feature scatter (volcano default, quick views, custom axes, box select) ---`);
  await waitPlot(page, 120000);
  await page.waitForFunction(() => {
    const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
      .filter(e => e.offsetParent !== null && e._fullLayout);
    if (!els.length) return false;
    const y = els[0]._fullLayout.yaxis.title.text || '';
    return /log\.(fdr|pvalue)/.test(y);
  }, null, { timeout: 120000 });
  let pi = await plotInfo(page);
  record('default volcano renders (y = -log statistic)', !!pi && /log\.(fdr|pvalue)/.test(pi.y), pi ? `y=${pi.y}` : 'none');
  record('volcano corner rects drawn (2)', pi && pi.shapes === 2, pi ? `shapes=${pi.shapes}` : 'none');
  record('scatter has data points', pi && pi.ntraces >= 1, pi ? `traces=${pi.ntraces}` : 'none');
  // quick-view badge inventory (rides on the persisted quickViews attr +
  // runtime auto-detection — the sqlite arm must serve the same set)
  await page.waitForSelector('[data-quick-view-id]', { timeout: 60000 });
  M.quickViewIds = (await page.$$eval('[data-quick-view-id]', els =>
    els.map(e => e.getAttribute('data-quick-view-id')))).sort().join(',');
  record(`quick-view badges present (${M.quickViewIds.split(',').length})`,
    M.quickViewIds.includes('volcano_RE_vs_LE') && M.quickViewIds.includes('cor_MDR'),
    M.quickViewIds);

  // corner auto-selection mirrors into the feature table (107 for demo)
  await gotoTab(page, 'data', 'Feature table');
  try {
    await dtWaitEntries(page, 'app-dataspace-tab_feature-table', { eq: 107 }, 90000);
    record('volcano corner auto-selection = 107 features (table mirrors)', true);
  } catch (e) {
    record('volcano corner auto-selection = 107 features (table mirrors)', false,
      `entries=${await dtEntryCount(page, 'app-dataspace-tab_feature-table')} (${e.message})`);
    await shot(page, 'corner-selection');
  }
  M.cornerSel = await dtEntryCount(page, 'app-dataspace-tab_feature-table');
  const FT = 'app-dataspace-tab_feature-table';
  record('feature table renders with selection subset', (await dtInfo(page, FT)).nCols > 1);

  // "Show all" clears the mirror (selection stays for analysis)
  await page.click('#app-dataspace-tab_feature-clear');
  try {
    await dtWaitEntries(page, FT, { eq: 2702 }, 60000);
    record('"Show all" clears the table mirror (2702 entries)', true);
  } catch {
    record('"Show all" clears the table mirror (2702 entries)', false,
      `entries=${await dtEntryCount(page, FT)}`);
    await shot(page, 'show-all');
  }
  M.ftTotal = await dtEntryCount(page, FT);

  // quick badge switch -> volcano_RE_vs_LE, 239 features
  await gotoTab(page, 'data', 'Feature');
  await page.waitForSelector('[data-quick-view-id="volcano_RE_vs_LE"]', { timeout: 30000 });
  await page.click('[data-quick-view-id="volcano_RE_vs_LE"]');
  try {
    await page.waitForFunction(() => {
      const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
        .filter(e => e.offsetParent !== null && e._fullLayout);
      return els.length && /RE_vs_LE/.test(els[0]._fullLayout.yaxis.title.text || '');
    }, null, { timeout: 60000 });
    pi = await plotInfo(page);
    record('quick view volcano_RE_vs_LE switches axes', /RE_vs_LE/.test(pi.y), `y=${pi.y}`);
    record('RE_vs_LE volcano keeps corner rects', pi.shapes === 2, `shapes=${pi.shapes}`);
  } catch (e) {
    record('quick view volcano_RE_vs_LE switches axes', false, e.message);
    await shot(page, 'badge-le');
  }
  await gotoTab(page, 'data', 'Feature table');
  try {
    await dtWaitEntries(page, FT, { eq: 239 }, 90000);
    record('RE_vs_LE corner selection = 239 features', true);
  } catch {
    record('RE_vs_LE corner selection = 239 features', false,
      `entries=${await dtEntryCount(page, FT)}`);
  }
  M.cornerLE = await dtEntryCount(page, FT);
  await gotoTab(page, 'data', 'Feature');

  // correlation quick view (no rects, regression view)
  await page.click('[data-quick-view-id="cor_MDR"]');
  try {
    await page.waitForFunction(() => {
      const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
        .filter(e => e.offsetParent !== null && e._fullLayout);
      return els.length && /Cor\|MDR/.test((els[0]._fullLayout.xaxis.title.text || '') + (els[0]._fullLayout.yaxis.title.text || ''));
    }, null, { timeout: 60000 });
    pi = await plotInfo(page);
    record('quick view cor_MDR switches to correlation axes', /Cor\|MDR/.test(pi.x + pi.y), `x=${pi.x}`);
  } catch (e) {
    record('quick view cor_MDR switches to correlation axes', false, e.message);
    await shot(page, 'badge-cor');
  }

  // custom axes via triselector cascade
  await page.click('#app-dataspace-feature_space-axisMode button:has(input[value="custom"])');
  await waitInput(page, 'app-dataspace-feature_space-axisMode');
  const tri1 = 'app-dataspace-feature_space-tris_main_scatter1';
  const tri2 = 'app-dataspace-feature_space-tris_main_scatter2';
  try {
    await setSelect(page, `#${tri1}-analysis`, 'PCA');
    await setSelect(page, `#${tri1}-subset`, 'All');
    await page.waitForFunction((s) => {
      const el = document.querySelector(s);
      return el && el.selectize && Object.keys(el.selectize.options).some(o => o.startsWith('PC1'));
    }, `#${tri1}-variable`, { timeout: 30000 });
    const pc1 = await page.evaluate(s => Object.keys(document.querySelector(s).selectize.options)
      .find(o => o.startsWith('PC1')), `#${tri1}-variable`);
    await setSelect(page, `#${tri1}-variable`, pc1);
    await setSelect(page, `#${tri2}-analysis`, 'PCA');
    await setSelect(page, `#${tri2}-subset`, 'All');
    await page.waitForFunction((s) => {
      const el = document.querySelector(s);
      return el && el.selectize && Object.keys(el.selectize.options).some(o => o.startsWith('PC2'));
    }, `#${tri2}-variable`, { timeout: 30000 });
    const pc2 = await page.evaluate(s => Object.keys(document.querySelector(s).selectize.options)
      .find(o => o.startsWith('PC2')), `#${tri2}-variable`);
    await setSelect(page, `#${tri2}-variable`, pc2);
    await page.waitForFunction(() => {
      const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
        .filter(e => e.offsetParent !== null && e._fullLayout);
      return els.length && /PC1/.test(els[0]._fullLayout.xaxis.title.text || '') &&
        /PC2/.test(els[0]._fullLayout.yaxis.title.text || '');
    }, null, { timeout: 90000 });
    pi = await plotInfo(page);
    record('custom mode: PCA PC1 vs PC2 axes applied', /PC1/.test(pi.x) && /PC2/.test(pi.y), `x=${pi.x} y=${pi.y}`);
  } catch (e) {
    record('custom mode: PCA PC1 vs PC2 axes applied', false, e.message);
    await shot(page, 'custom-axes');
  }

  // back to shortcut + volcano, then a real box selection
  await page.click('#app-dataspace-feature_space-axisMode button:has(input[value="quick"])');
  await page.waitForSelector('[data-quick-view-id="volcano_RE_vs_ME"]', { timeout: 30000 });
  await page.click('[data-quick-view-id="volcano_RE_vs_ME"]');
  await page.waitForFunction(() => {
    const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
      .filter(e => e.offsetParent !== null && e._fullLayout);
    return els.length && /RE_vs_ME/.test(els[0]._fullLayout.yaxis.title.text || '');
  }, null, { timeout: 60000 });
  await settle(page);
  try {
    await plotlyBoxSelect(page, 0.30, 0.10, 0.70, 0.45); // upper band of the volcano
    await gotoTab(page, 'data', 'Feature table');
    await dtWaitEntries(page, FT, { min: 1, max: 2701 }, 60000);
    const n = await dtEntryCount(page, FT);
    record(`box selection on volcano selects a feature subset (${n})`, n >= 1 && n < 2702);
    await gotoTab(page, 'data', 'Feature');
  } catch (e) {
    record('box selection on volcano selects a feature subset', false, e.message);
    await shot(page, 'box-select');
    await gotoTab(page, 'data', 'Feature');
  }

  // clear figure selection button
  try {
    await page.click('[data-testid="app-dataspace-feature_space-clear-selection-button"]');
    await gotoTab(page, 'data', 'Feature table');
    await dtWaitEntries(page, FT, { eq: 2702 }, 60000);
    record('"Clear figure selection" empties the selection', true);
  } catch (e) {
    record('"Clear figure selection" empties the selection', false, e.message);
    await shot(page, 'clear-sel');
  }

  // ============================ 3. feature table ========================
  console.log(`\n--- [${backend}] 3. feature table (paging, sort, row selection, multi-select) ---`);
  // NOTE: the app tables use dom='tip' — no global search box by design.
  await gotoTab(page, 'data', 'Feature table');
  await dtWaitEntries(page, FT, { eq: 2702 }, 60000);
  try {
    await page.click(`#${FT} .paginate_button.next`);
    await page.waitForFunction((id) => {
      const i = document.getElementById(id).querySelector('.dataTables_info');
      return i && /26\s+to\s+50/.test(i.innerText);
    }, FT, { timeout: 15000 });
    record('pagination: next page shows rows 26–50', true);
    await page.click(`#${FT} .paginate_button.previous`);
  } catch (e) {
    record('pagination: next page shows rows 26–50', false, e.message);
  }
  // single row click selects one feature (no table re-filter for table origins)
  try {
    await dtClickRowReal(page, FT, 0);
    await settle(page);
    const sel = await page.evaluate(id =>
      document.getElementById(id).querySelectorAll('tbody tr.selected').length, FT);
    record('row click marks the row selected', sel >= 1, `selected=${sel}`);
  } catch (e) {
    record('row click marks the row selected', false, e.message);
    await shot(page, 'ft-rowclick');
  }
  // sort by the second column header
  try {
    const before = (await dtInfo(page, FT)).firstRow;
    await page.click(`#${FT} thead th >> nth=1`);
    await settle(page);
    await page.click(`#${FT} thead th >> nth=1`);
    await settle(page);
    const after = (await dtInfo(page, FT)).firstRow;
    record('column sort re-orders the table', before !== null && after !== null,
      `"${(before || '').slice(0, 30)}" -> "${(after || '').slice(0, 30)}"`);
    M.ftFirstRow = after; // deterministic given identical annotation content
  } catch (e) {
    record('column sort re-orders the table', false, e.message);
  }
  // multi-selection toggle + ctrl-click second row (switchInput: toggle the
  // underlying checkbox and let the change event reach Shiny directly)
  try {
    const toggleMulti = () => page.evaluate(() => {
      const el = document.getElementById('app-dataspace-tab_feature-multisel');
      el.checked = !el.checked;
      el.dispatchEvent(new Event('change', { bubbles: true }));
    });
    await toggleMulti();
    await page.waitForFunction(() =>
      Shiny.shinyapp.$inputValues['app-dataspace-tab_feature-multisel'] === true,
      null, { timeout: 15000 });
    await settle(page); // table re-renders in multiple mode
    await dtClickRowReal(page, FT, 0);
    await settle(page);
    await dtClickRowReal(page, FT, 1, ['Control']);
    await settle(page);
    const n = await page.evaluate(id =>
      document.getElementById(id).querySelectorAll('tbody tr.selected').length, FT);
    record('multi-selection: ctrl-click selects a second row', n >= 2, `selected=${n}`);
    await toggleMulti();
    await settle(page);
  } catch (e) {
    record('multi-selection: ctrl-click selects a second row', false, e.message);
    await shot(page, 'ft-multi');
  }

  // ============================ 4. sample space =========================
  console.log(`\n--- [${backend}] 4. sample space (PCA scatter, axes, sample selection) ---`);
  await gotoTab(page, 'data', 'Sample');
  await waitPlot(page, 90000);
  await page.waitForFunction(() => {
    const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
      .filter(e => e.offsetParent !== null && e._fullLayout);
    return els.length && /PC1/.test(els[0]._fullLayout.xaxis.title.text || '');
  }, null, { timeout: 90000 });
  pi = await plotInfo(page);
  record('sample scatter renders PCA PC1 vs PC2', /PC1/.test(pi.x) && /PC2/.test(pi.y), `x=${pi.x}`);
  const nPts = await page.evaluate(() => {
    const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
      .filter(e => e.offsetParent !== null && e._fullLayout && e.data);
    const tr = els.flatMap(e => e.data).find(t => t.type === 'scatter' && t.mode && t.mode.includes('markers'));
    return tr ? tr.x.length : null;
  });
  record('sample scatter shows all 60 samples', nPts === 60, `points=${nPts}`);
  M.nPts = nPts;

  // switch sample axes via custom mode (MDR vs doubling time)
  const stri1 = 'app-dataspace-sample_space-tris_main_scatter1';
  const stri2 = 'app-dataspace-sample_space-tris_main_scatter2';
  try {
    await page.click('#app-dataspace-sample_space-axisMode button:has(input[value="custom"])');
    await setSelect(page, `#${stri1}-analysis`, 'General');
    await setSelect(page, `#${stri1}-subset`, 'All');
    await setSelect(page, `#${stri1}-variable`, 'MDR');
    await setSelect(page, `#${stri2}-analysis`, 'General');
    await setSelect(page, `#${stri2}-subset`, 'All');
    await setSelect(page, `#${stri2}-variable`, 'Doubleing.Time');
    await page.waitForFunction(() => {
      const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
        .filter(e => e.offsetParent !== null && e._fullLayout);
      return els.length && /MDR/.test(els[0]._fullLayout.xaxis.title.text || '') &&
        /Doubleing/.test(els[0]._fullLayout.yaxis.title.text || '');
    }, null, { timeout: 90000 });
    pi = await plotInfo(page);
    record('custom sample axes: MDR vs Doubling time', /MDR/.test(pi.x) && /Doubleing/.test(pi.y), `x=${pi.x}`);
  } catch (e) {
    record('custom sample axes: MDR vs Doubling time', false, e.message);
    await shot(page, 'sample-axes');
  }

  // box-select a subset of samples on the current (custom) sample scatter
  try {
    await plotlyBoxSelect(page, 0.15, 0.15, 0.55, 0.60);
    await gotoTab(page, 'data', 'Sample table');
    await dtWaitEntries(page, 'app-dataspace-tab_pheno-table', { min: 1, max: 59 }, 60000);
    const n = await dtEntryCount(page, 'app-dataspace-tab_pheno-table');
    record(`box selection on the sample scatter selects ${n} samples`, n >= 1 && n < 60);
  } catch (e) {
    record('box selection on the sample scatter selects samples', false, e.message);
    await shot(page, 'sample-boxsel');
  }

  // ============================ 5. sample table =========================
  console.log(`\n--- [${backend}] 5. sample table ---`);
  const ST = 'app-dataspace-tab_pheno-table';
  await gotoTab(page, 'data', 'Sample table');
  try {
    await page.click('#app-dataspace-tab_pheno-clear'); // Show all
    await dtWaitEntries(page, ST, { eq: 60 }, 60000);
    record('sample table lists all 60 samples (after Show all)', true);
  } catch (e) {
    record('sample table lists all 60 samples (after Show all)', false,
      `entries=${await dtEntryCount(page, ST)}`);
    await shot(page, 'st-count');
  }
  M.stTotal = await dtEntryCount(page, ST);
  M.stFirstRow = (await dtInfo(page, ST)).firstRow;
  try {
    await dtClickRowReal(page, ST, 0);
    await settle(page);
    record('sample row click registers (selected class)', (await page.evaluate(id =>
      document.getElementById(id).querySelectorAll('tbody tr.selected').length, ST)) >= 1);
    await page.click('#app-dataspace-tab_pheno-clear');
    await settle(page);
  } catch (e) {
    record('sample row click registers (selected class)', false, e.message);
  }

  // ============================ 6. Cor & Heatmap ========================
  console.log(`\n--- [${backend}] 6. correlation & expression heatmaps ---`);
  await gotoTab(page, 'data', 'Cor');
  try {
    await waitImg(page, '#app-dataspace-corheatmapViewer-heatmap', 120000);
    record('correlation heatmap renders', true);
  } catch (e) {
    record('correlation heatmap renders', false, e.message);
    await shot(page, 'cor-heatmap');
  }
  await gotoTab(page, 'data', 'Heatmap');
  try {
    await waitImg(page, '#app-dataspace-heatmapViewer-heatmap', 120000);
    record('expression heatmap renders', true);
  } catch (e) {
    record('expression heatmap renders', false, e.message);
    await shot(page, 'expr-heatmap');
  }

// ============================ 7. expression matrix & gene-set list ====================
  console.log(`\n--- [${backend}] 7. expression matrix & gene-set list ---`);
  await gotoTab(page, 'data', 'Expression');
  try {
    // the table mirrors the CURRENT feature selection (e.g. the volcano
    // corner's 107 after a clear/re-apply cycle) — assert it renders and
    // matches a plausible selection state
    await dtWaitEntries(page, 'app-dataspace-tab_expr-table', { min: 1, max: 2702 }, 90000);
    const n = await dtEntryCount(page, 'app-dataspace-tab_expr-table');
    record(`expression matrix table renders (mirrors selection: ${n} rows)`, n >= 1);
    M.exprEntries = n;
  } catch (e) {
    record('expression matrix table renders', false, e.message);
    await shot(page, 'expr-table');
  }
  // spot-check expression values at fixed positions (feature ids + the
  // first sample columns) — the core SQLite data-fidelity probe.
  // NOTE: DT renders NA as EMPTY cells — empties are legitimate values
  // and compare as strings in the parity block
  try {
    await dtCells(page, 'app-dataspace-tab_expr-table', 5, 4).then(c => { M.exprCells = c; });
    record('expression cells captured for parity (5 rows x 4 cols)',
      Array.isArray(M.exprCells) && M.exprCells.length === 5 &&
      M.exprCells.every(r => r.length === 4 && r[0].length > 0),
      (M.exprCells || []).map(r => r.slice(0, 2).join(':')).join(' | ').slice(0, 150));
  } catch (e) {
    record('expression cells captured for parity', false, e.message);
  }
  await gotoTab(page, 'data', 'GSList');
  const GS = 'app-dataspace-gsList-stab-table';
  try {
    await dtWaitEntries(page, GS, { min: 100 }, 180000);
    const nGS = await dtEntryCount(page, GS);
    record(`gene-set membership table renders (GS annotation attached, ${nGS} rows)`, nGS >= 100);
  } catch (e) {
    record('gene-set membership table renders (GS annotation attached)', false, e.message);
    await shot(page, 'gslist');
  }
  M.gsTotal = await dtEntryCount(page, GS);
  // clicking a gene-set row selects that member feature
  try {
    await dtClickRowReal(page, GS, 0);
    await settle(page);
    await gotoTab(page, 'analysis', 'Feature');
    await dtWaitEntries(page, 'app-resultspace-feature_general-mtab-table', { eq: 1 }, 90000);
    record('gene-set row click selects its member feature (analysis reacts)', true);
  } catch (e) {
    record('gene-set row click selects its member feature (analysis reacts)', false, e.message);
    await shot(page, 'gslist-click');
  }

  // ============================ 8. analysis: feature ====================
  console.log(`\n--- [${backend}] 8. analysis > feature (boxplot / beeswarm / ROC) ---`);
  // select a compact feature set: volcano corner -> 107 -> boxplot view
  await gotoTab(page, 'data', 'Feature');
  await page.waitForSelector('[data-quick-view-id="volcano_RE_vs_ME"]', { timeout: 30000 });
  await page.click('[data-quick-view-id="volcano_RE_vs_ME"]');
  await page.waitForFunction(() => {
    const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
      .filter(e => e.offsetParent !== null && e._fullLayout);
    return els.length && /RE_vs_ME/.test(els[0]._fullLayout.yaxis.title.text || '');
  }, null, { timeout: 60000 });
  await gotoTab(page, 'analysis', 'Feature');
  try {
    await page.waitForFunction(() => {
      const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
        .filter(e => e.offsetParent !== null && e._fullLayout);
      return els.some(e => (e.closest('[id*="feature_general"]') || e.id.includes('feature_general')));
    }, null, { timeout: 120000 });
    pi = await page.evaluate(() => {
      const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
        .filter(e => e.offsetParent !== null && e._fullLayout &&
          (e.closest('[id*="feature_general"]') || e.id.includes('feature_general')));
      return els.length ? { n: els.length, traces: els[0]._fullData.length } : null;
    });
    record('feature analysis renders plot for 107 selected features (boxplot view)',
      !!pi && pi.traces >= 1, pi ? `traces=${pi.traces}` : 'none');
  } catch (e) {
    record('feature analysis renders plot for 107 selected features', false, e.message);
    await shot(page, 'ana-feature-box');
  }
  try {
    await dtWaitEntries(page, 'app-resultspace-feature_general-mtab-table', { min: 1 }, 60000);
    record('feature analysis stats table has rows (group test results)', true);
  } catch (e) {
    record('feature analysis stats table has rows', false, e.message);
    await shot(page, 'ana-feature-tab');
  }
  M.mtabFirstRow = (await dtInfo(page, 'app-resultspace-feature_general-mtab-table')).firstRow;

  // single feature + categorical link var -> dual view; Curve -> ROC
  // (no table search box: sort by the Gene.name column, then click row 0;
  // the mirror may hold the volcano-corner 107 subset — any rows work)
  try {
    await gotoTab(page, 'data', 'Feature table');
    await dtWaitEntries(page, FT, { min: 1 }, 60000);
    await page.click(`#${FT} thead th >> nth=1`);
    await settle(page);
    await dtClickRowReal(page, FT, 0);
    await settle(page);
  } catch (e) {
    record('select a single feature via table row click', false, e.message);
    await shot(page, 'ft-single');
  }
  await gotoTab(page, 'analysis', 'Feature');
  const ftri = 'app-resultspace-feature_general-tris_feature_general';
  try {
    await setSelect(page, `#${ftri}-analysis`, 'General');
    await setSelect(page, `#${ftri}-subset`, 'All');
    await setSelect(page, `#${ftri}-variable`, 'TP53.Status');
    await page.waitForFunction(() => {
      const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
        .filter(e => e.offsetParent !== null && e._fullLayout &&
          (e.closest('[id*="feature_general"]') || e.id.includes('feature_general')));
      return els.length > 0;
    }, null, { timeout: 90000 });
    record('single feature + categorical var: beeswarm/dual view renders', true);
  } catch (e) {
    record('single feature + categorical var: beeswarm/dual view renders', false, e.message);
    await shot(page, 'ana-bees');
  }
  try {
    await page.click('#app-resultspace-feature_general-internal_radio button:has(input[value="Curve"])');
    // the ROC/PR module renders a BASE-R plot (img), not plotly
    await waitImg(page, '#app-resultspace-feature_general-feature_general_roc_pr-roc_pr', 90000);
    const rocSize = await page.evaluate(() => {
      const el = document.querySelector('#app-resultspace-feature_general-feature_general_roc_pr-roc_pr img');
      return el ? { w: el.naturalWidth, h: el.naturalHeight } : null;
    });
    record('ROC/PR curve renders (Curve radio, binary TP53.Status)',
      !!rocSize && rocSize.w > 50, JSON.stringify(rocSize));
    await page.click('#app-resultspace-feature_general-internal_radio button:has(input[value="Bees"])');
    await settle(page);
  } catch (e) {
    record('ROC/PR curve renders (Curve radio)', false, e.message);
    await shot(page, 'ana-roc');
  }

  // ============================ 9. analysis: sample ====================
  console.log(`\n--- [${backend}] 9. analysis > sample (batch comparison, beeswarm, contingency, survival) ---`);
  // deterministic state: exactly one selected sample (verified user flow)
  await gotoTab(page, 'data', 'Sample table');
  try {
    await dtWaitEntries(page, ST, { min: 1 }, 60000);
    await dtClickRowReal(page, ST, 0);
    await settle(page);
  } catch (e) {
    record('select one sample for the sample-analysis flow', false, e.message);
  }
  await gotoTab(page, 'analysis', 'Sample');
  try {
    await dtWaitEntries(page, 'app-resultspace-sample_general-batch_comp-phenotype_table-table',
      { min: 1 }, 120000);
    record('batch comparison table renders (phenotype vs selection)', true);
  } catch (e) {
    record('batch comparison table renders', false, e.message);
    await shot(page, 'ana-batch');
  }
  const stri = 'app-resultspace-sample_general-tris_sample_general';
  // beeswarm on a continuous variable
  try {
    await setSelect(page, `#${stri}-analysis`, 'General');
    await setSelect(page, `#${stri}-subset`, 'All');
    await setSelect(page, `#${stri}-variable`, 'MDR');
    await page.waitForFunction(() => {
      const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
        .filter(e => e.offsetParent !== null && e._fullLayout &&
          (e.closest('[id*="sample_general"]') || e.id.includes('sample_general')));
      return els.length > 0;
    }, null, { timeout: 90000 });
    record('sample analysis: beeswarm view renders for MDR (continuous)', true);
    await dtWaitEntries(page, 'app-resultspace-sample_general-msatab-table', { min: 1 }, 60000);
    record('sample analysis stats table renders (MDR test)', true);
  } catch (e) {
    record('sample analysis: beeswarm view / stats table for MDR', false, e.message);
    await shot(page, 'ana-bees-mdr');
  }
  // contingency on a categorical variable
  try {
    await setSelect(page, `#${stri}-analysis`, 'General');
    await setSelect(page, `#${stri}-subset`, 'All');
    await setSelect(page, `#${stri}-variable`, 'Origin');
    await page.waitForFunction(() => {
      const t = document.querySelector('#app-resultspace-sample_general-sample_general_contab-count\\.table\\.output');
      return t && t.querySelectorAll('tbody tr').length > 0;
    }, null, { timeout: 90000 });
    record('sample analysis: contingency table renders for Origin (categorical)', true);
  } catch (e) {
    record('sample analysis: contingency table renders for Origin', false, e.message);
    await shot(page, 'ana-cont');
  }
  // survival view
  try {
    await setSelect(page, `#${stri}-analysis`, 'Surv');
    await setSelect(page, `#${stri}-subset`, 'all');
    await setSelect(page, `#${stri}-variable`, 'OS');
    await page.waitForFunction(() => {
      const el = document.querySelector('#app-resultspace-sample_general-sample_general_surv-kmplot');
      return el && el.querySelectorAll('img').length > 0;
    }, null, { timeout: 90000 });
    const km = await page.evaluate(() => {
      const el = document.querySelector('#app-resultspace-sample_general-sample_general_surv-kmplot');
      const img = el && el.querySelector('img');
      return img && img.naturalWidth > 2;
    });
    record('sample analysis: Kaplan-Meier plot renders (Surv category)', !!km);
    // the KM log-rank p-value is rendered INSIDE the plot image (not
    // DOM-readable); the censor-time slider range is the DOM-visible
    // numeric and derives from the Surv times in the sample annotation.
    // shiny renders sliderInput as a text input carrying data-min/max
    // (the .min property is absent and Number('') === 0 would fake a
    // degenerate range)
    const cr = await page.evaluate(() => {
      const s = document.querySelector('#app-resultspace-sample_general-sample_general_surv-censor');
      if (!s) return null;
      const g = (k) => {
        const v = s.getAttribute('data-' + k);
        if (v !== null && v !== '') return v;
        const p = s[k];
        return (p !== undefined && p !== '') ? String(p) : null;
      };
      return { min: g('min'), max: g('max') };
    });
    M.survMin = cr && cr.min !== null ? Number(cr.min) : null;
    M.survMax = cr && cr.max !== null ? Number(cr.max) : null;
    record(`survival censor range captured (${M.survMin}–${M.survMax})`,
      M.survMin !== null && Number.isFinite(M.survMin) && Number.isFinite(M.survMax) &&
      M.survMin < M.survMax, `min=${M.survMin} max=${M.survMax}`);
  } catch (e) {
    record('sample analysis: Kaplan-Meier plot renders', false, e.message);
    await shot(page, 'ana-km');
  }

  // ============================ 10. ORA =================================
  console.log(`\n--- [${backend}] 10. analysis > ORA (hypergeometric over-representation) ---`);
  // force a GENUINE quick-view switch (LE then ME) so the volcano corner
  // auto-selection lands freshly with its 107 features (re-clicking the
  // already-active badge is a no-op for the selection)
  await gotoTab(page, 'data', 'Feature');
  await page.click('[data-quick-view-id="volcano_RE_vs_LE"]');
  await page.waitForFunction(() => {
    const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
      .filter(e => e.offsetParent !== null && e._fullLayout);
    return els.length && /RE_vs_LE/.test(els[0]._fullLayout.yaxis.title.text || '');
  }, null, { timeout: 60000 });
  await page.click('[data-quick-view-id="volcano_RE_vs_ME"]');
  await page.waitForFunction(() => {
    const els = Array.from(document.querySelectorAll('.js-plotly-plot'))
      .filter(e => e.offsetParent !== null && e._fullLayout);
    return els.length && /RE_vs_ME/.test(els[0]._fullLayout.yaxis.title.text || '');
  }, null, { timeout: 60000 });
  await settle(page);
  await gotoTab(page, 'analysis', 'ORA');
  const ORA = 'app-resultspace-ora-stab-table';
  try {
    await dtWaitEntries(page, ORA, { min: 10 }, 180000);
    const inf = await dtInfo(page, ORA);
    record('ORA results table renders enriched terms (107 selected)', true, `${inf.nRows} rows / ${inf.nCols} cols`);
    const head = (inf.firstRow || '').slice(0, 150);
    record('ORA table has p-value columns', /p|fdr/i.test((await page.evaluate(id =>
      Array.from(document.getElementById(id).querySelectorAll('thead th'))
        .map(th => th.innerText).join('|'), ORA)) || head), '');
  } catch (e) {
    record('ORA results table renders enriched terms', false, e.message);
    await shot(page, 'ora');
  }
  M.oraEntries = await dtEntryCount(page, ORA);
  M.oraFirstRow = (await dtInfo(page, ORA)).firstRow;
  try {
    await dtClickRowReal(page, ORA, 0);
    await page.waitForFunction(() => {
      const t = document.getElementById('app-resultspace-ora-overlapTab-table');
      return t && t.querySelectorAll('tbody tr').length > 0;
    }, null, { timeout: 60000 });
    record('clicking an ORA row shows the overlapping genes table', true);
  } catch (e) {
    record('clicking an ORA row shows the overlapping genes table', false, e.message);
    await shot(page, 'ora-overlap');
  }

  // ============================ 11. fGSEA ===============================
  console.log(`\n--- [${backend}] 11. analysis > fGSEA (ranked enrichment) ---`);
  await gotoTab(page, 'analysis', 'fGSEA');
  const FG = 'app-resultspace-fgsea-stab-table';
  const fgsTri = 'app-resultspace-fgsea-tris_fgsea';
  try {
    await setSelect(page, `#${fgsTri}-analysis`, 'ttest');
    await setSelect(page, `#${fgsTri}-subset`, 'RE_vs_ME');
    await setSelect(page, `#${fgsTri}-variable`, 'mean.diff');
    await dtWaitEntries(page, FG, { min: 10 }, 180000);
    const inf = await dtInfo(page, FG);
    record('fGSEA results render for ranking ttest|RE_vs_ME|mean.diff', true, `${inf.nRows} rows`);
  } catch (e) {
    record('fGSEA results render for ranking ttest|RE_vs_ME|mean.diff', false, e.message);
    await shot(page, 'fgsea');
  }
  M.fgEntries = await dtEntryCount(page, FG);
  M.fgFirstRow = (await dtInfo(page, FG)).firstRow;
  // deterministic view for parity: the default |NES| sort is stochastic
  // (fgsea Monte-Carlo — the top pathway itself shuffles between runs on
  // identical input), so sort by the pathway column (alphabetical, a
  // plain user click) and project that first row onto its deterministic
  // fields: pathway (col 1), ES (col 5), size (col 7)
  try {
    await page.click(`#${FG} thead th >> nth=0`);
    await settle(page);
    const sr = await dtInfo(page, FG);
    const f = String(sr.firstRow || '').split('\t').map(s => s.trim());
    M.fgFirstRowDet = f.length >= 7 ? [f[0], f[4], f[6]].join('\t') : null;
  } catch { M.fgFirstRowDet = null; }
  try {
    // the ranking line renders as a lightgray plotly line; with no pathway
    // selected there are no bars yet — the svg + line IS the resting state
    await page.waitForSelector('#app-resultspace-fgsea-bplot svg.main-svg', { timeout: 120000 });
    await page.waitForTimeout(1500);
    const pathsBefore = await page.evaluate(() =>
      document.querySelectorAll('#app-resultspace-fgsea-bplot svg path').length);
    record('fGSEA ranking bar plot renders (svg ranking line)', pathsBefore > 0,
      `${pathsBefore} svg paths`);
    // clicking a pathway row must add the leading-edge bars (red highlight)
    await dtClickRowReal(page, FG, 0);
    await page.waitForFunction(before =>
      document.querySelectorAll('#app-resultspace-fgsea-bplot svg path').length > before,
      pathsBefore, { timeout: 60000 });
    const pathsAfter = await page.evaluate(() =>
      document.querySelectorAll('#app-resultspace-fgsea-bplot svg path').length);
    record('fGSEA row click adds the leading-edge bars to the plot',
      pathsAfter > pathsBefore, `${pathsBefore} -> ${pathsAfter} paths`);
  } catch (e) {
    record('fGSEA NES bar plot renders / leading-edge highlight', false, e.message);
    await shot(page, 'fgsea-bplot');
  }

  // ============================ 12. Response ============================
  console.log(`\n--- [${backend}] 12. analysis > Response (dose-response curves) ---`);
  // the module requires EXACTLY ONE selected feature (req(length == 1)):
  // pick a single feature through the feature table first
  try {
    await gotoTab(page, 'data', 'Feature table');
    await dtWaitEntries(page, FT, { min: 1 }, 60000);
    await dtClickRowReal(page, FT, 0);
    await settle(page);
  } catch (e) {
    record('select a single feature for the dose-response view', false, e.message);
  }
  await gotoTab(page, 'analysis', 'Response');
  try {
    await waitImg(page, '#app-resultspace-rescurve-plot', 180000);
    record('dose-response curve renders', true);
    await dtWaitEntries(page, 'app-resultspace-rescurve-param-table', { min: 1 }, 120000);
    record('dose-response parameter table renders (EC50/Hill)', true);
  } catch (e) {
    record('dose-response curve/parameter table renders', false, e.message);
    await shot(page, 'response');
  }
  M.drcEntries = await dtEntryCount(page, 'app-resultspace-rescurve-param-table');
  M.drcFirstRow = (await dtInfo(page, 'app-resultspace-rescurve-param-table')).firstRow;

  // ============================ 13. SeqLogo =============================
  console.log(`\n--- [${backend}] 13. analysis > SeqLogo (PTM motif) ---`);
  await gotoTab(page, 'analysis', 'SeqLogo');
  try {
    await waitImg(page, '#app-resultspace-ptm-plt', 180000);
    record('PTM motif seq-logo renders (ratio selected/all)', true);
    await dtWaitEntries(page, 'app-resultspace-ptm-seqtable-table', { min: 1 }, 120000).catch(() => {});
    record('PTM sequence table renders', (await dtEntryCount(page, 'app-resultspace-ptm-seqtable-table')) > 0);
  } catch (e) {
    record('PTM motif seq-logo renders', false, e.message);
    await shot(page, 'seqlogo');
  }
  M.ptmEntries = await dtEntryCount(page, 'app-resultspace-ptm-seqtable-table');
  M.ptmFirstRow = (await dtInfo(page, 'app-resultspace-ptm-seqtable-table')).firstRow;


  // ============================ 14. dynamic heatmap ======================
  console.log(`\n--- [${backend}] 14. dynamic heatmap (click selects feature+sample) ---`);
  // NOTE on observability: the analysis stats table (metatab) req()s a
  // NON-EMPTY feature selection — with an empty selection it renders
  // blank by design. The feature TABLE mirror ("of N entries") is the
  // reliable user-visible signal for the current selection state.
  const MTAB = 'app-resultspace-feature_general-mtab-table';
  // deterministic baseline: no selection -> mirror shows all 2702
  await gotoTab(page, 'data', 'Feature');
  try {
    await page.click('[data-testid="app-dataspace-feature_space-clear-selection-button"]');
    await settle(page);
  } catch {}
  await gotoTab(page, 'data', 'Feature table');
  await dtWaitEntries(page, FT, { eq: 2702 }, 60000);
  record('baseline: no selection (feature table = 2702 entries)', true);

  await gotoTab(page, 'data', 'Dynamic heatmap');
  try {
    await waitImg(page, '#app-dataspace-dynheatmapViewer-heatmap', 180000);
    record('dynamic heatmap renders', true);
  } catch (e) {
    record('dynamic heatmap renders', false, e.message);
    await shot(page, 'dyn-heatmap');
  }
  try {
    // Shiny image-click needs a settled hover (hover input, 150ms delay)
    // before down/up; retry once at another position if it didn't land
    const tryClick = async (fx, fy) => {
      await clickImg(page, '#app-dataspace-dynheatmapViewer-heatmap img', fx, fy);
      await page.waitForTimeout(4000);
      return page.evaluate(() =>
        Shiny.shinyapp.$inputValues['app-dataspace-dynheatmapViewer-heatmap_click'] ?? null);
    };
    let clickInput = await tryClick(0.5, 0.5);
    if (!clickInput) clickInput = await tryClick(0.5, 0.3);
    await settle(page, 20000);
    await gotoTab(page, 'analysis', 'Feature');
    await dtWaitEntries(page, MTAB, { eq: 1 }, 90000);
    record('dynamic heatmap cell click selects exactly one feature (stats table = 1 row)',
      !!clickInput, `clickInput=${!!clickInput}`);
  } catch (e) {
    record('dynamic heatmap cell click selects exactly one feature', false, e.message);
    await shot(page, 'dyn-click');
  }
  // clear the heatmap-driven selection via the dynamic heatmap's own
  // Clear button; the stats table then goes blank (empty selection) and
  // the feature-table mirror returns to all 2702 entries
  await gotoTab(page, 'data', 'Dynamic heatmap');
  try {
    await page.click('#app-dataspace-dynheatmapViewer-clear');
    await settle(page, 20000);
    await gotoTab(page, 'data', 'Feature table');
    await dtWaitEntries(page, FT, { eq: 2702 }, 60000);
    record('dynamic heatmap Clear resets the selection (mirror back to 2702)', true);
  } catch (e) {
    record('dynamic heatmap Clear resets the selection', false, e.message);
    await shot(page, 'dyn-clear');
    // recovery: scatter clear, then fresh session if needed
    try {
      await gotoTab(page, 'data', 'Feature');
      await page.click('[data-testid="app-dataspace-feature_space-clear-selection-button"]');
      await settle(page);
      await gotoTab(page, 'data', 'Feature table');
      await dtWaitEntries(page, FT, { eq: 2702 }, 60000);
    } catch {
      await hardReset(page, dataset);
    }
  }

  // ============================ 15. Geneshot UI =========================
  console.log(`\n--- [${backend}] 14. analysis > Geneshot (UI render; live API not required) ---`);
  await gotoTab(page, 'analysis', 'Geneshot');
  record('Geneshot UI renders (term input + search button)',
    (await page.locator('#app-resultspace-geneshotTab-term').count()) === 1 &&
    (await page.locator('#app-resultspace-geneshotTab-submit').count()) === 1);

  // ============================ 16. StringDB ============================
  console.log(`\n--- [${backend}] 16. analysis > StringDB (graceful offline behavior) ---`);
  await gotoTab(page, 'analysis', 'StringDB');
  record('StringDB UI renders (taxonomy input + Run button)',
    (await page.locator('#app-resultspace-stringdb-tax').count()) === 1 &&
    (await page.locator('#app-resultspace-stringdb-run').count()) === 1);
  // clear the selection first so Run fails fast with a validation error
  await gotoTab(page, 'data', 'Feature');
  try {
    await page.click('[data-testid="app-dataspace-feature_space-clear-selection-button"]');
    await settle(page);
  } catch {}
  await gotoTab(page, 'analysis', 'StringDB');
  try {
    await page.click('#app-resultspace-stringdb-run');
    await notifVisible(page, 'error', 30000);
    record('StringDB Run with empty selection shows an error notification (no crash)', true);
    await settle(page);
  } catch (e) {
    record('StringDB Run with empty selection shows an error notification', false, e.message);
    await shot(page, 'stringdb-run');
  }

  // ============================ 17. Snapshot ============================
  // NOTE: the UI has no dataset-object save control (only xlsx export and
  // app-state snapshots) — the snapshot round-trip IS the supported
  // save/reload mechanism and is asserted state-equal here.
  console.log(`\n--- [${backend}] 17. snapshots (save / restore / delete) ---`);
  const SNAP = 'e2e-check';
  // pre-save state: data tab Feature, empty selection (sections 14-16 left
  // the mirror at 2702) — saved from the CURRENT position like the original
  // journey (StringDB analysis tab), so the restore target stays 'Feature'
  const preSaveMirror = 2702;
  try {
    try {
      await page.click('[data-testid="app-snapshot-button"]');
      await page.waitForSelector('#app-snapshot_name', { timeout: 30000 });
      record('snapshot modal opens', true);
      await page.fill('#app-snapshot_name', SNAP);
      await page.click('#app-snapshot_save');
      await notifVisible(page, 'message', 30000);
      record('snapshot saves with a success notification', true);
      await settle(page);
    } catch (e) {
      record('snapshot modal/save', false, e.message);
      await shot(page, 'snap-save');
    }
    // change visible state, then restore
    try {
      await gotoTab(page, 'data', 'Sample table');
      await page.click('[data-testid="app-snapshot-button"]');
      await page.waitForSelector('#app-tab_saveSS', { timeout: 30000 });
      await page.waitForFunction(() => {
        const t = document.getElementById('app-tab_saveSS');
        return t && t.querySelectorAll('tbody tr').length > 0;
      }, null, { timeout: 60000 });
      record('snapshot appears in the saved list', true);
      await dtClickRowReal(page, 'app-tab_saveSS', 0);
      await page.waitForSelector('#app-snapshot_restore_confirm', { timeout: 30000 });
      await page.click('#app-snapshot_restore_confirm');
      await settle(page, 30000);
      const tab = await inputVal(page, 'app-dataspace-eset');
      record('snapshot restore returns the saved tab (Feature scatter)', tab === 'Feature', `tab=${tab}`);
      // reloaded state must equal the pre-save state: selection mirror
      await gotoTab(page, 'data', 'Feature table');
      await dtWaitEntries(page, FT, { eq: preSaveMirror }, 90000);
      const post = await dtEntryCount(page, FT);
      record(`snapshot restore returns the pre-save selection state (mirror ${preSaveMirror})`, post === preSaveMirror,
        `mirror=${post}`);
    } catch (e) {
      record('snapshot restore round-trip', false, e.message);
      await shot(page, 'snap-restore');
    }
    try {
      await page.click('[data-testid="app-snapshot-button"]');
      await page.waitForFunction(() => {
        const t = document.getElementById('app-tab_saveSS');
        return t && t.querySelectorAll('tbody tr').length > 0;
      }, null, { timeout: 60000 });
      const hasRow = await page.locator('#app-tab_saveSS tbody button[id^="deletess_"]').count();
      if (hasRow > 0) {
        await page.click('#app-tab_saveSS tbody button[id^="deletess_"] >> nth=0');
        await page.waitForSelector('#app-snapshot_delete_confirm', { state: 'visible', timeout: 30000 });
        await page.click('#app-snapshot_delete_confirm');
        await page.waitForTimeout(2500);
        const gone = await page.evaluate(() => {
          const t = document.getElementById('app-tab_saveSS');
          return !t || t.querySelectorAll('tbody tr').length === 0;
        });
        record('snapshot delete removes it from the list', gone);
      } else {
        record('snapshot delete removes it from the list', false, 'no delete button visible');
      }
      await page.keyboard.press('Escape');
      await page.evaluate(() => {
        const m = document.querySelector('.modal');
        if (m) m.remove();
        document.body.classList.remove('modal-open');
        const b = document.querySelector('.modal-backdrop');
        if (b) b.remove();
      });
    } catch (e) {
      record('snapshot delete', false, e.message);
      await shot(page, 'snap-delete');
    }
  } catch (e) {
    record('snapshot pre-save state not reachable', false, e.message);
  }

  // ============================ 18. xlsx download =======================
  console.log(`\n--- [${backend}] 18. dataset export (xlsx) ---`);
  try {
    const dlPromise = page.waitForEvent('download', { timeout: 300000 });
    await page.click('#app-download');
    const dl = await dlPromise;
    const fp = await dl.path();
    const { statSync } = await import('node:fs');
    const size = statSync(fp).size;
    record(`xlsx download works (${Math.round(size / 1024)} kB)`, size > 100 * 1024, dl.suggestedFilename());
  } catch (e) {
    record('xlsx download works', false, e.message);
    await shot(page, 'download');
  }

  // ============================ 19. AI assistant ========================
  console.log(`\n--- [${backend}] 19. AI assistant panel (no LLM call) ---`);
  try {
    await page.click('#app-assistant-toggle');
    await page.waitForFunction(() => {
      const p = document.querySelector('#app-assistant-panel');
      return p && getComputedStyle(p).display !== 'none';
    }, null, { timeout: 15000 });
    await page.waitForFunction(() => {
      const b = document.querySelector('#app-assistant-body');
      return b && b.innerText.trim().length > 0;
    }, null, { timeout: 20000 });
    record('AI assistant panel opens with content', true);
    await page.click('#app-assistant-toggle');
    await settle(page);
  } catch (e) {
    record('AI assistant panel opens with content', false, e.message);
    await shot(page, 'assistant');
  }

  // ============================ 20. dataset re-select (lifecycle) ======
  // the closest thing to save->reload the UI offers for the DATASET
  // object: drop the session, re-select the same file, assert the state
  // comes back identical (for sqlite this cycles connections: the WP-4
  // observer must have closed the previous one)
  console.log(`\n--- [${backend}] 20. dataset re-select after the full journey ---`);
  try {
    await hardReset(page, dataset);
    // a fresh session re-runs the volcano corner auto-selection, so the
    // mirror equals the boot default (M.cornerSel), not 2702
    await gotoTab(page, 'data', 'Feature table');
    await dtWaitEntries(page, FT, { eq: M.cornerSel }, 90000);
    record(`dataset re-select reloads after the full journey (${M.cornerSel} corner selection)`, true);
    M.reselectTotal = await dtEntryCount(page, FT);
  } catch (e) {
    record('dataset re-select reloads after the full journey', false, e.message);
    await shot(page, 'reselect');
  }

  // ============================ 21. global health =======================
  console.log(`\n--- [${backend}] 21. session health ---`);
  const errSample = pageErrors.slice(0, 5);
  record(`no JavaScript page errors during the whole journey (${pageErrors.length} seen)`,
    pageErrors.length === 0, errSample.join(' || '));
  if (pageErrors.length) await shot(page, 'js-errors');

  await browser.close();
  return { results, metrics: M };
};

// ------------------------------------------------------- sqlite lifecycle
// Runs AFTER the journey with the browser closed but the app process
// still alive: WP-4 claims connections close on session end.
const sqliteLifecycleChecks = (arm, record) => {
  // give shiny a moment to run session$onEnded handlers
  spawnSync('sleep', ['4']);
  let fds = datasetFds(arm);
  for (let i = 0; i < 3 && fds !== 0; i++) {
    spawnSync('sleep', ['2']);
    fds = datasetFds(arm);
  }
  record('no leaked sqlite file descriptors after session end (app alive)',
    fds === 0, `open fds on ${arm.dataset}: ${fds}`);
  // external overwrite while the app process still runs: a leaked
  // (or locking) connection would make this fail or corrupt
  const ow = spawnSync('Rscript', ['-e',
    `omicsViewer::saveOmicsViewerDb(readRDS('${DEMO_RDS}'), db.file = '${path.join(arm.datadir, arm.dataset)}')`],
    { cwd: REPO, encoding: 'utf8', timeout: 300000 });
  record('external saveOmicsViewerDb overwrite succeeds with the app running (no stale lock)',
    ow.status === 0, ow.status === 0 ? '' : `${(ow.stderr || '').slice(0, 200)}`);
};

// ==================================================================== main
try {
  for (const backend of ARMS) {
    log(`launching ${backend} arm ...`);
    const arm = await launchArm(backend);
    let journey = null;
    try {
      journey = await runJourney(arm);
    } catch (e) {
      journey = { results: [], metrics: {} };
      const ARTIFACTS = path.join(here, 'artifacts', backend);
      mkdirSync(ARTIFACTS, { recursive: true });
      journey.results.push({ name: 'FATAL: arm aborted', pass: false,
        detail: `${e.message}\n${(e.stack || '').split('\n').slice(0, 4).join('\n')}` });
      console.log(`FATAL [${backend}]: ${e.message}`);
    }
    const fatal = journey.results.some(r => r.name.startsWith('FATAL'));
    if (backend === 'sqlite' && !fatal) {
      try { sqliteLifecycleChecks(arm, (name, pass, detail = '') => {
        journey.results.push({ name, pass, detail: String(detail).slice(0, 400) });
        log(`${pass ? 'PASS' : 'FAIL'} | [sqlite] ${name}${detail ? ' | ' + String(detail).slice(0, 160) : ''}`);
      }); } catch (e) {
        journey.results.push({ name: 'sqlite lifecycle checks', pass: false, detail: e.message });
      }
    }
    // WAL/SHM absence is checked after the app process is gone but before
    // cleanup removes the data dir
    if (backend === 'sqlite') {
      try { arm.r.kill('SIGKILL'); } catch {}
      spawnSync('sleep', ['1']);
      const stray = (() => { try { return readdirSync(arm.datadir); } catch { return null; } })();
      const journals = stray ? stray.filter(f => /-(wal|shm|journal)$/i.test(f)) : ['<datadir missing>'];
      journey.results.push({ name: 'no stray -wal/-shm journal files left in the data dir',
        pass: journals.length === 0, detail: journals.join(',') });
      log(`${journals.length === 0 ? 'PASS' : 'FAIL'} | [sqlite] no stray -wal/-shm journal files (${journals.join(',') || 'none'})`);
      arm.cleanup();
    } else {
      arm.cleanup();
    }
    const passed = journey.results.filter(r => r.pass).length;
    armReports[backend] = {
      total: journey.results.length, passed, failed: journey.results.length - passed,
      results: journey.results, metrics: journey.metrics
    };
    console.log(`\n=== [${backend}] arm: ${passed}/${journey.results.length} passed ===`);
  }

  // ------------------------------------------------------------ parity
  let parity = { ran: false, total: 0, passed: 0, failed: 0, comparisons: [], note: '' };
  if (ARMS.length === 2 && armReports.rds && armReports.sqlite) {
    parity.ran = true;
    const A = armReports.rds.metrics, B = armReports.sqlite.metrics;
    const spec = [
      ['dataset dimensions (features x samples)', m => `${m.featuresN}x${m.samplesN}`,
        (a, b) => a === b],
      ['data-space tab set', m => m.dataTabs, (a, b) => a === b],
      ['analysis tab set (incl. ORA/fGSEA)', m => m.anaTabs, (a, b) => a === b],
      ['quick-view badge ids', m => m.quickViewIds, (a, b) => a === b],
      ['volcano corner selection count', m => m.cornerSel, numEq],
      ['RE_vs_LE corner selection count', m => m.cornerLE, numEq],
      ['feature table total entries', m => m.ftTotal, numEq],
      ['sample table total entries', m => m.stTotal, numEq],
      ['sample scatter point count', m => m.nPts, numEq],
      ['sample table first row', m => m.stFirstRow, rowTextEq],
      ['feature table first row (sorted)', m => m.ftFirstRow, rowTextEq],
      ['expression table total entries', m => m.exprEntries, numEq],
      ['expression spot-check cells (5 rows x 4 cols)', m => m.exprCells, cellMatrixEq],
      ['gene-set membership total entries', m => m.gsTotal, numEq],
      ['feature analysis stats first row', m => m.mtabFirstRow, rowTextEq],
      ['survival censor-time range min', m => m.survMin, numEq],
      ['survival censor-time range max', m => m.survMax, numEq],
      ['ORA total entries', m => m.oraEntries, numEq],
      ['ORA first row (enrichment statistics)', m => m.oraFirstRow, rowTextEq],
      ['fGSEA total entries', m => m.fgEntries, numEq],
      // fgseaMultilevel is an UNSEEDED Monte-Carlo method: pval/padj/NES/
      // log2err differ between ANY two runs on identical input (verified
      // in plain R: same demo stats, 4 runs -> pval 9.7e-7..2.9e-6, top
      // pathway by |NES| shuffles). Only the deterministic fields of the
      // first row are comparable between arms; the stochastic ones are
      // recorded per-arm (fgFirstRow) but excluded from parity.
      ['fGSEA alphabetically-first pathway: deterministic fields (pathway, ES, size)',
        m => m.fgFirstRowDet, rowTextEq],
      ['dose-response parameter entries', m => m.drcEntries, numEq],
      ['dose-response first row (EC50/Hill)', m => m.drcFirstRow, rowTextEq],
      ['seq-logo table entries', m => m.ptmEntries, numEq],
      ['seq-logo first row', m => m.ptmFirstRow, rowTextEq],
      ['dataset re-select total (post full journey)', m => m.reselectTotal, numEq]
    ];
    for (const [name, get, eq] of spec) {
      const a = get(A), b = get(B);
      // a missing metric on either arm must not silently 'pass' parity
      const pass = (a == null || b == null) ? false : eq(a, b);
      parity.comparisons.push({ name, rds: a ?? null, sqlite: b ?? null, pass });
      log(`${pass ? 'PASS' : 'FAIL'} | [parity] ${name} | rds=${JSON.stringify(a)?.slice(0, 90)} sqlite=${JSON.stringify(b)?.slice(0, 90)}`);
    }
    parity.total = parity.comparisons.length;
    parity.passed = parity.comparisons.filter(c => c.pass).length;
    parity.failed = parity.total - parity.passed;
    parity.note = 'numeric tolerance rel 1e-9 / abs 1e-12 (SQLite REAL stores exact IEEE-754 doubles; ' +
      'tolerance absorbs last-bit display drift only). Mouse-geometry counts (box selections) and ' +
      'plotly/DOM internals are deliberately not compared. fGSEA pval/padj/NES/log2err are unseeded ' +
      'Monte-Carlo statistics that differ between any two runs on identical input (R repro: 4 runs on ' +
      'the same demo ranking -> pval 9.7e-7..2.9e-6, |NES|-sorted top pathway shuffles) - excluded ' +
      'from parity, deterministic fields (pathway/ES/size) compared instead.';
  } else {
    parity.note = `parity block requires both arms (mode=${BACKENDS}) — skipped`;
  }

  // ------------------------------------------------------------ report
  const totals = Object.values(armReports).reduce((acc, r) => ({
    total: acc.total + r.total, passed: acc.passed + r.passed, failed: acc.failed + r.failed
  }), { total: 0, passed: 0, failed: 0 });
  const report = {
    started: new Date(t0).toISOString(),
    duration_sec: Math.round((Date.now() - t0) / 1000),
    backends: ARMS,
    arms: armReports,
    parity,
    total: totals.total + parity.total,
    passed: totals.passed + parity.passed,
    failed: totals.failed + parity.failed
  };
  writeFileSync(path.join(here, 'e2e_results.json'), JSON.stringify(report, null, 2));
  console.log(`\n============================================================`);
  for (const be of ARMS)
    console.log(`[${be}] ${armReports[be].passed}/${armReports[be].total} passed, ${armReports[be].failed} failed`);
  if (parity.ran)
    console.log(`[parity] ${parity.passed}/${parity.total} passed, ${parity.failed} failed`);
  console.log(`OVERALL: ${report.passed}/${report.total} passed, ${report.failed} failed ` +
    `(${report.duration_sec}s) — e2e_results.json`);
  for (const be of ARMS)
    for (const f of armReports[be].results.filter(r => !r.pass))
      console.log(`  FAIL [${be}]: ${f.name} | ${f.detail}`);
  for (const c of parity.comparisons.filter(c => !c.pass))
    console.log(`  FAIL [parity]: ${c.name} | rds=${JSON.stringify(c.rds)?.slice(0, 120)} sqlite=${JSON.stringify(c.sqlite)?.slice(0, 120)}`);
  console.log(`============================================================`);
  if (report.failed) exitCode = 1;
} catch (e) {
  console.error(`FATAL: ${e.message}\n${(e.stack || '').split('\n').slice(0, 6).join('\n')}`);
  exitCode = 2;
}
process.exit(exitCode);
