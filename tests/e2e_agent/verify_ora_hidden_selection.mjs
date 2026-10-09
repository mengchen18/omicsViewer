// Browser-level guard for the R-H3 stale-ORA bug (module_ora.R).
//
// Repro: load demo.RDS (default volcano selection) -> stay on the Feature
// tab -> change the scatter selection while the ORA tab is HIDDEN -> click
// the ORA tab. The ORA results must reflect the NEW selection. The old
// collapse observer was an observeEvent whose handler ran inside isolate()
// (bindEvent.Observer), so the output_visible() clientData read could never
// re-trigger it: hidden selection changes were dropped and the tab kept
// showing the previous (default) selection's results.
//
// The hidden selection change is driven through a volcano quick-view
// switch - the established harness idiom - which flows through the exact
// same selection bus -> ri -> reactive_i() path a mouse rect selection
// uses; the bug was mechanism-agnostic.
//
// Run:  node verify_ora_hidden_selection.mjs     (from tests/e2e_agent/)

import { chromium } from 'playwright';
import { spawn, execSync } from 'node:child_process';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../..');
const EXTDATA = path.join(REPO, 'inst/extdata');
const PORT = 7792;

// orphaned R processes poison subsequent runs (AGENTS.md): sweep first
try {
  const out = execSync('ss -tlnp', { encoding: 'utf8', stdio: ['pipe', 'pipe', 'pipe'] });
  for (const pid of out.split('\n').filter(l => l.includes(`:${PORT}`))
       .flatMap(l => [...l.matchAll(/pid=\K[0-9]+/g)].map(m => m[0]))) {
    try { process.kill(parseInt(pid, 10), 'SIGKILL'); } catch {}
  }
} catch {}

const results = [];
const record = (name, pass, detail = '') => {
  results.push({ name, pass, detail: String(detail).slice(0, 300) });
  console.log(`${pass ? 'PASS' : 'FAIL'} | ${name}${detail ? ' | ' + detail : ''}`);
};

const r = spawn('Rscript', ['-e', `
  options(shiny.port = ${PORT}, shiny.host = '127.0.0.1')
  eset <- readRDS(file.path('${EXTDATA}', 'demo.RDS'))
  omicsViewer::omicsViewer(dir = '${EXTDATA}', ESVObj = eset)
`], { cwd: REPO, env: { ...process.env, OMICSVIEWER_TEST_HOOKS: 'true' }, stdio: ['ignore', 'pipe', 'pipe'] });
process.on('exit', () => { try { r.kill('SIGKILL'); } catch {} });
process.on('SIGINT', () => { try { r.kill('SIGKILL'); } catch {}; process.exit(130); });

const waitPort = async (port, ms = 90000) => {
  const t0 = Date.now();
  while (Date.now() - t0 < ms) {
    try {
      const res = await fetch(`http://127.0.0.1:${port}/`);
      if (res.status > 0) return true;
    } catch {}
    await new Promise(s => setTimeout(s, 500));
  }
  return false;
};

const HOOK = {
  box: '#app-agentTestHooks-container',
  op: '#app-agentTestHooks-op',
  payload: '#app-agentTestHooks-payload',
  run: '#app-agentTestHooks-run',
  result: '#app-agentTestHooks-result'
};

async function runHook(page, op, payload) {
  const before = (await page.innerText(HOOK.result)).trim();
  await page.evaluate(({ sel, value }) => {
    const el = document.querySelector(sel);
    if (el.selectize) el.selectize.setValue(value);
    else {
      el.value = value;
      el.dispatchEvent(new Event('change', { bubbles: true }));
    }
  }, { sel: HOOK.op, value: op });
  await page.fill(HOOK.payload, JSON.stringify(payload ?? {}));
  const runCount = await page.evaluate(sel => {
    const el = document.querySelector(sel);
    const n = (parseInt(el.dataset.runs || '0', 10) + 1);
    el.dataset.runs = String(n);
    return n;
  }, HOOK.run);
  await page.click(HOOK.run);
  await page.waitForFunction(({ sel, before, n }) => {
    const t = (document.querySelector(sel)?.innerText || '').trim();
    return t !== before && t.includes(`"hook_run": ${n}`);
  }, { sel: HOOK.result, before, n: runCount }, { timeout: 30000 });
  return JSON.parse((await page.innerText(HOOK.result)).trim());
}

const selectionCount = (page) => runHook(page, 'overview', {})
  .then(ov => (ov && ov.selection && ov.selection.features) ? ov.selection.features.count : null);

const oraSig = (page) => page.evaluate(() => {
  const tb = document.querySelector('[id*="ora-stab-table"] table tbody');
  if (!tb || !tb.children.length) return null;
  const rows = Array.from(tb.children).slice(0, 8).map(tr =>
    Array.from(tr.cells).map(td => td.innerText.trim()).join('|'));
  return { n: tb.children.length, sig: rows.join('\n') };
});

let exitCode = 0;
try {
  if (!await waitPort(PORT)) throw new Error('app did not start');
  const browser = await chromium.launch({
    executablePath: '/usr/bin/google-chrome',
    headless: true,
    args: ['--no-sandbox', '--disable-gpu', '--disable-dev-shm-usage',
           '--disable-webgl', '--disable-webgl2']
  });
  const ctx = await browser.newContext();
  const page = await ctx.newPage();
  await page.goto(`http://127.0.0.1:${PORT}/`, { timeout: 90000 });
  await page.waitForFunction(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-eset'] &&
    document.querySelector('#app-contents').offsetParent !== null, null, { timeout: 90000 });
  await page.evaluate(sel => {
    document.querySelector(sel).style.display = 'block';
  }, HOOK.box);

  // 1. baseline: ORA results for the DEFAULT volcano selection
  // (analyst navbar links carry leading whitespace: " Feature" - use
  // has-text scoped to the analyst navbar ul, never :text-is)
  const analystTab = (t) => `#app-resultspace-analyst a:has-text("${t}")`;
  await page.waitForSelector('#app-resultspace-analyst', { timeout: 60000 });
  await page.click(analystTab('ORA'));
  await page.waitForSelector('[id*="ora-stab-table"] table', { timeout: 90000, state: 'attached' });
  await page.waitForFunction(() =>
    (document.querySelector('[id*="ora-stab-table"] table tbody') || { children: [] }).children.length > 0,
    null, { timeout: 90000 });
  await new Promise(s => setTimeout(s, 3000));
  const sigA = await oraSig(page);
  const selA = await selectionCount(page);
  record('baseline: ORA table renders for the default selection',
    !!sigA && sigA.n > 0, sigA ? `rows=${sigA.n} selected=${selA}` : 'no table');

  // 2. leave the ORA tab (it becomes hidden)
  await page.click(analystTab('Feature'));
  await new Promise(s => setTimeout(s, 2000));

  // 3. change the selection while the ORA tab is hidden (rect-selection
  //    equivalent): switch to another volcano quick view. NOTE: the view
  //    must settle a NON-EMPTY selection - an empty one suspends the
  //    collapse observer at req(reactive_i()) (pre-existing, orthogonal
  //    edge), which would also leave the old table on screen.
  await page.click('[data-quick-view-id="volcano_RE_vs_LE"]');
  let selB = null;
  {
    const counts = [];
    const t0 = Date.now();
    while (Date.now() - t0 < 30000) {
      const c = await selectionCount(page);
      if (typeof c === 'number') counts.push(c);
      const n = counts.length;
      if (n >= 2 && typeof c === 'number' && c !== selA && c > 0 &&
          counts[n - 1] === counts[n - 2]) {
        selB = c;
        break;
      }
      await new Promise(s => setTimeout(s, 1000));
    }
  }
  record('hidden selection change lands in the app state',
    typeof selB === 'number' && selB !== selA,
    `selected ${selA} -> ${selB}`);
  if (!(typeof selB === 'number' && selB !== selA))
    throw new Error('selection never changed; cannot verify ORA catch-up');

  // 4. switch to the ORA tab: the results must catch up with the hidden
  //    change (the buggy build keeps the default selection's table)
  await page.click(analystTab('ORA'));
  let sigB = null;
  {
    const t0 = Date.now();
    while (Date.now() - t0 < 45000) {
      sigB = await oraSig(page);
      if (sigB && sigB.n > 0 && sigB.sig !== sigA.sig) break;
      await new Promise(s => setTimeout(s, 1000));
    }
  }
  record('ORA table reflects the hidden selection change after the tab switch',
    !!(sigB && sigB.n > 0 && sigB.sig !== sigA.sig),
    sigB ? `rows ${sigA.n} -> ${sigB.n}; head "${(sigB.sig || '').split('\n')[0].slice(0, 120)}"`
         : 'no table after flip');

  // 5. stability: no re-render churn after the catch-up
  await page.evaluate(() => {
    window.__ora_renders = 0;
    const mo = new MutationObserver(muts => {
      for (const m of muts) for (const n of m.addedNodes) {
        if (n.nodeType !== 1) continue;
        const tb = n.tagName === 'TABLE' ? n : (n.querySelector ? n.querySelector('table') : null);
        if (tb && tb.closest('[id*="ora-stab"]')) window.__ora_renders++;
      }
    });
    mo.observe(document.body, { childList: true, subtree: true });
  });
  await new Promise(s => setTimeout(s, 6000));
  const churn = await page.evaluate(() => window.__ora_renders);
  record('ORA table stable after the catch-up (no re-render loop)',
    churn <= 4, `counted=${churn}`);

  await ctx.close();
  await browser.close();
} catch (e) {
  console.error('HARNESS ERROR:', e);
  exitCode = 1;
}

const failed = results.filter(x => !x.pass).length;
console.log(`\n${results.length - failed}/${results.length} passed`);
if (failed) exitCode = 1;
process.exit(exitCode);
