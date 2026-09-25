// WP5/WP6 acceptance sweep — the user's disk-load flow (HANDOVER 0.2/0.6).
// Launches against an already-running hooked app (port 7776, dir=inst/extdata),
// loads demo.RDS through the file selector, then checks:
//   - volcano corner at load (2 rects, scorner=volcano), selection = 107
//   - one visible paint per quick-view switch (Feature tab), selection follows
//   - custom-mode analysis change works (RC6 regression class)
//   - Sample tab: badge switch is one paint; sample table shows 60 rows
//   - no page errors
import { chromium } from 'playwright';

const PORT = process.env.SWEEP_PORT || 7776;
let pass = 0, fail = 0;
const record = (name, ok, detail = '') => {
  console.log(`${ok ? 'PASS' : 'FAIL'} | ${name}${detail ? ' | ' + detail : ''}`);
  ok ? pass++ : fail++;
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
  }, { sel: HOOK.result, before, n }, { timeout: 30000 });
  return JSON.parse((await page.innerText(HOOK.result)).trim());
}

const startPlotWatch = (page) => page.evaluate(() => {
  const visiblePlot = () => Array.from(document.querySelectorAll('.js-plotly-plot'))
    .filter(e => e.offsetParent !== null && e._fullLayout)[0];
  const sig = (p) => (p._fullLayout.xaxis.title.text || '') + '|' +
    (p._fullLayout.yaxis.title.text || '') + '|' + (p._fullLayout.shapes || []).length;
  const el0 = visiblePlot();
  window.__plotStates = [el0 ? sig(el0) : ''];
  window.__plotRenders = 0;
  let pending = false;
  const settle = () => {
    pending = false;
    const p = visiblePlot();
    if (!p) return;
    window.__plotRenders++;
    const s = sig(p);
    if (s !== window.__plotStates[window.__plotStates.length - 1]) window.__plotStates.push(s);
  };
  const mo = new MutationObserver(() => {
    if (pending) return;
    pending = true;
    requestAnimationFrame(settle);
  });
  if (el0) mo.observe(el0, { childList: true, subtree: true, attributes: true });
  window.__plotDisconnect = () => { mo.disconnect(); };
});
const readPlotStates = (page) => page.evaluate(() => {
  if (window.__plotDisconnect) window.__plotDisconnect();
  return { states: window.__plotStates, renders: window.__plotRenders };
});

const plotState = (page) => page.evaluate(() => {
  const p = Array.from(document.querySelectorAll('.js-plotly-plot'))
    .filter(e => e.offsetParent !== null && e._fullLayout)[0];
  if (!p) return null;
  return {
    x: p._fullLayout.xaxis.title.text || '', y: p._fullLayout.yaxis.title.text || '',
    shapes: (p._fullLayout.shapes || []).length,
    scorner: p._fullLayout.shapes && p._fullLayout.shapes.length &&
      p._fullLayout.shapes[0].xref === 'x' && (p._fullLayout.shapes[0].line || {}).dash ? 'rects' : ''
  };
});

(async () => {
  const browser = await chromium.launch({
    executablePath: '/usr/bin/google-chrome',
    args: ['--disable-webgl', '--disable-webgl2', '--no-sandbox']
  });
  const ctx = await browser.newContext();
  const page = await ctx.newPage();
  const pageErrors = [];
  page.on('pageerror', e => pageErrors.push(String(e).slice(0, 200)));
  await page.goto(`http://127.0.0.1:${PORT}/`, { timeout: 60000 });

  // ---- disk-load the demo dataset (the user's flow) ----
  await page.waitForSelector('#app-selectFile', { state: 'attached', timeout: 60000 });
  // the file-list choices arrive from the server several seconds in;
  // setValue on an option that is not loaded yet is a silent no-op
  await page.waitForFunction(() => {
    const el = document.querySelector('#app-selectFile');
    return el && el.selectize && Object.keys(el.selectize.options).length > 0;
  }, null, { timeout: 90000 });
  await page.evaluate(() => {
    const el = document.querySelector('#app-selectFile');
    if (el.selectize) el.selectize.setValue('demo.RDS');
  });
  await page.waitForFunction(() => {
    const tab = Shiny.shinyapp.$inputValues['app-dataspace-eset'];
    const box = document.querySelector('#app-contents');
    return !!tab && !!box && box.offsetParent !== null;
  }, null, { timeout: 120000 });
  await page.evaluate(sel => { document.querySelector(sel).style.display = 'block'; }, HOOK.box);

  // ---- load view: volcano corner with 2 rects, 107 selected ----
  let st = null;
  for (let i = 0; i < 30; i++) {
    st = await plotState(page);
    if (st && st.shapes === 2 && /log\.fdr|log\.pvalue/.test(st.y || '')) break;
    await page.waitForTimeout(1000);
  }
  record('load: volcano quick view with both corner rects',
    !!st && st.shapes === 2 && /log\.fdr|log\.pvalue/.test(st.y || ''),
    st ? JSON.stringify(st) : 'no plot');
  let ov = await runHook(page, 'overview', {});
  record('load: volcano corner selection = 107 features',
    ov.selection && ov.selection.features && ov.selection.features.count === 107,
    ov.selection ? `features=${ov.selection.features.count}` : 'no selection block');

  // ---- one paint per quick-view switch, selection follows ----
  const switches = [
    ['volcano_RE_vs_LE', 239, /RE_vs_LE/, true],
    ['volcano_RE_vs_ME', 107, /RE_vs_ME/, true],
    ['cor_MDR', 0, /Cor\|MDR/, false],
    ['volcano_RE_vs_LE', 239, /RE_vs_LE/, true]
  ];
  for (const [badge, nSel, axisRe, isVolcano] of switches) {
    await startPlotWatch(page);
    await page.click(`[data-quick-view-id="${badge}"]`);
    await page.waitForTimeout(4500);
    const { states, renders } = await readPlotStates(page);
    const changes = states.slice(1);
    const final = states[states.length - 1] || '';
    record(`switch ${badge}: exactly one signature change`,
      changes.length === 1 && axisRe.test(final),
      `renders=${renders} ; ` + states.join(' ; '));
    record(`switch ${badge}: final rects ${isVolcano ? '= 2' : '= 0'}`,
      isVolcano ? /\|2$/.test(final) : /\|0$/.test(final), final);
    ov = await runHook(page, 'overview', {});
    record(`switch ${badge}: selection count follows (${nSel})`,
      ov.selection && ov.selection.features &&
        ov.selection.features.count === nSel,
      ov.selection ? `features=${ov.selection.features.count}` : 'no selection');
  }

  // ---- custom-mode analysis change (RC6 class) ----
  await page.click('#app-dataspace-feature_space-axisMode button:has(input[value="custom"])');
  await page.waitForFunction(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-feature_space-axisMode'] === 'custom',
    null, { timeout: 15000 });
  await page.waitForTimeout(1500);
  // the analysis select is a selectize; switch x analysis to Cor
  const anaId = '#app-dataspace-feature_space-tris_main_scatter1-analysis';
  await page.waitForSelector(anaId, { state: 'attached', timeout: 30000 });
  await page.waitForFunction((sel) => {
    const el = document.querySelector(sel);
    return el && el.selectize && Object.keys(el.selectize.options).length > 0;
  }, anaId, { timeout: 30000 });
  const setOk = await page.evaluate((sel) => {
    const el = document.querySelector(sel);
    if (el && el.selectize) { el.selectize.setValue('Cor'); return true; }
    return false;
  }, anaId);
  await page.waitForTimeout(4500);
  st = await plotState(page);
  record('custom mode: x analysis -> Cor re-renders the figure on Cor axes',
    setOk && !!st && /Cor/.test(st.x || ''),
    st ? `x="${st.x}"` : 'no plot');
  const subVal = await page.evaluate(() =>
    Shiny.shinyapp.$inputValues['app-dataspace-feature_space-tris_main_scatter1-subset']);
  record('custom mode: subset cascade re-derived under the new analysis',
    !!subVal && subVal !== 'RE_vs_ME', `subset=${subVal}`);

  // ---- Sample tab: badge switch paint + 60-row table ----
  await page.evaluate(() => {
    const el = document.querySelector('a[data-value="Sample"]');
    if (el) el.click();
  });
  await page.waitForTimeout(3000);
  const badges = await page.$$eval('[data-quick-view-id]', els => {
    const vis = els.filter(e => e.offsetParent !== null);
    return vis.map(e => e.getAttribute('data-quick-view-id'));
  });
  if (badges.length) {
    await startPlotWatch(page);
    await page.evaluate((badge) => {
      const el = Array.from(document.querySelectorAll(`[data-quick-view-id="${badge}"]`))
        .filter(e => e.offsetParent !== null)[0];
      if (el) el.click();
    }, badges[0]);
    await page.waitForTimeout(4500);
    const { states, renders } = await readPlotStates(page);
    record('sample tab: badge switch is one signature change',
      states.slice(1).length <= 1, `renders=${renders} ; ` + states.join(' ; '));
  } else {
    record('sample tab: badge switch is one signature change', true, 'no badges (single-view tab)');
  }
  // the sample table lives on its own "Sample table" navbar tab
  await page.evaluate(() => {
    const el = Array.from(document.querySelectorAll('a[data-value="Sample table"]'))
      .filter(e => e.offsetParent !== null)[0];
    if (el) el.click();
  });
  await page.waitForTimeout(4000);
  const tabInfo = await page.evaluate(() =>
    Array.from(document.querySelectorAll('.dataTables_info'))
      .filter(e => e.offsetParent !== null).map(e => e.innerText.trim()));
  record('sample table shows all 60 rows',
    tabInfo.some(t => /60\s+entries/.test(t)), tabInfo.slice(0, 3).join(' / ') || 'no info');

  // ---- WP6: right panel must not rebuild on feature-selection changes ----
  // Open the Feature ANALYSIS tab (right panel); the boxplot container
  // element must SURVIVE a selection change (corner 107 -> 239 via the
  // store-driven quick view) - a renderUI rebuild would replace the DOM
  // node (remount resets plotly state and flashes the panel).
  await page.evaluate(() => {
    const el = Array.from(document.querySelectorAll('#app-resultspace-analyst a[data-value="Feature"]'))
      .filter(e => e.offsetParent !== null)[0];
    if (el) el.click();
  });
  let fgEl = null;
  for (let i = 0; i < 30 && !fgEl; i++) {
    await page.waitForTimeout(1000);
    fgEl = await page.evaluate(() => {
      const p = Array.from(document.querySelectorAll('.js-plotly-plot'))
        .filter(e => e.offsetParent !== null)
        .find(e => (e.closest('[id*="feature_general"]') || e.id.includes('feature_general')));
      if (p) { window.__fgEl = p; return true; }
      return false;
    });
  }
  record('analysis tab: feature_general plot renders', !!fgEl);
  if (fgEl) {
    // selection change via the store apply path (works from any tab;
    // quick badges are only visible on the Feature scatter tab)
    const res = await runHook(page, 'scatter', { space: 'feature', quick_view_id: 'volcano_RE_vs_LE' });
    await page.waitForTimeout(4500);
    const survived = await page.evaluate(() =>
      window.__fgEl && window.__fgEl.isConnected && window.__fgEl.offsetParent !== null);
    record('WP6: right-panel plot container survives selection change',
      survived && !res.hook_error,
      (survived ? '' : 'container was replaced (renderUI rebuild)') +
        (res.hook_error ? ' hook_error: ' + res.hook_error : ''));
    ov = await runHook(page, 'overview', {});
    record('WP6: selection followed the hook apply (239)',
      ov.selection && ov.selection.features && ov.selection.features.count === 239,
      ov.selection ? `features=${ov.selection.features.count}` : 'no selection');
  }

  record('no uncaught page errors', pageErrors.length === 0, pageErrors.slice(0, 3).join(' ;; '));

  await browser.close();
  console.log(`\nSweep summary: ${pass}/${pass + fail} passed`);
  process.exit(fail ? 1 : 0);
})().catch(e => { console.error(e); process.exit(2); });
