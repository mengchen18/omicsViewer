// UI mimic driver — plays a realistic user session against the LLM
// assistant through the real chat UI (one app process, one browser
// context, many turns), with durable LLM diagnostic logs.
//
// Run:  node ui_mimic.mjs <turns.json> <outDirName>
//       turns.json = [{ id, prompt, toggle_logging_before? }, ...]

import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import { mkdirSync, copyFileSync, existsSync, readFileSync, writeFileSync,
         readdirSync, statSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../..');
const EXTDATA = path.join(REPO, 'inst/extdata');
const TURNS_FILE = process.argv[2];
const SESSION = process.argv[3] || `session-${Date.now()}`;
const OUT = `/home/chen/.omicsviewer-llm-eval/${SESSION}`;
const LOGDIR = path.join(OUT, 'llm-logs');
const ART = path.join(OUT, 'artifacts');
mkdirSync(LOGDIR, { recursive: true });
mkdirSync(ART, { recursive: true });

const TURNS = JSON.parse(readFileSync(TURNS_FILE, 'utf8'));
console.log(`session ${SESSION}: ${TURNS.length} turns -> ${OUT}`);

// ---- provider credentials from provider.env --------------------------------
const env = { ...process.env };
for (const line of readFileSync(path.join(here, 'provider.env'), 'utf8').split('\n')) {
  const m = line.match(/^\s*([A-Z_][A-Z0-9_]*)\s*=\s*(.*?)\s*$/);
  if (m && !line.trim().startsWith('#')) env[m[1]] = m[2];
}
if (!env.OMICSVIEWER_LLM_API_KEY) { console.error('FATAL: no API key'); process.exit(2); }
console.log(`provider: ${env.OMICSVIEWER_LLM_PROVIDER} model: ${env.OMICSVIEWER_LLM_MODEL}`);

// ---- spawn app --------------------------------------------------------------
const PORT = 7785;
const r = spawn('Rscript', ['-e', `
  options(shiny.port = ${PORT}, shiny.host = '127.0.0.1')
  eset <- readRDS(file.path('${EXTDATA}', 'demo.RDS'))
  omicsViewer::omicsViewer(dir = '${EXTDATA}', ESVObj = eset)
`], { cwd: REPO, env: { ...env, OMICSVIEWER_LLM_LOG: 'true', OMICSVIEWER_LLM_LOG_DIR: LOGDIR },
     stdio: ['ignore', 'pipe', 'pipe'] });
let rLog = '';
r.stdout.on('data', d => { rLog += d; });
r.stderr.on('data', d => { rLog += d; });
const cleanup = () => { try { r.kill('SIGKILL'); } catch {} };
process.on('exit', cleanup);
process.on('SIGINT', () => { cleanup(); process.exit(130); });

const waitPort = async (port, ms = 90000) => {
  const t0 = Date.now();
  while (Date.now() - t0 < ms) {
    try { const res = await fetch(`http://127.0.0.1:${port}/`); if (res.status > 0) return true; } catch {}
    await new Promise(s => setTimeout(s, 500));
  }
  throw new Error(`app did not start:\n${rLog}`);
};

const logRequestInFlight = () => {
  try {
    const logs = readFileSync(path.join(LOGDIR, latestLogName()), 'utf8')
      .split('\n').filter(l => l.trim());
    if (!logs.length) return false;
    return JSON.parse(logs[logs.length - 1]).event === 'provider_request_start';
  } catch { return false; }
};
const latestLogName = () => readdirSync(LOGDIR).filter(f => f.endsWith('.jsonl'))
  .map(f => ({ f, m: statSync(path.join(LOGDIR, f)).mtimeMs }))
  .sort((a, b) => b.m - a.m)[0].f;

const chatBusy = (pg) => pg.evaluate(() => {
  const b = document.querySelector('[class*="cancel"]');
  return !!(b && b.offsetParent !== null && !b.disabled);
});

// wait until the assistant is truly idle (no visible stream AND no
// provider request in flight) before delivering the next prompt -
// the composer silently swallows submissions while a chain is active
const waitIdle = async (pg, ms = 420000) => {
  const dl = Date.now() + ms;
  while (Date.now() < dl) {
    if (!(await chatBusy(pg)) && !logRequestInFlight()) return true;
    await new Promise(s => setTimeout(s, 2000));
  }
  return false;
};

const t0 = Date.now();
const ts = () => `t+${Math.round((Date.now() - t0) / 1000)}s`;

try {
  await waitPort(PORT);
  const browser = await chromium.launch({
    executablePath: '/usr/bin/google-chrome', headless: true,
    args: ['--no-sandbox', '--disable-gpu', '--disable-dev-shm-usage',
           '--disable-webgl', '--disable-webgl2']
  });
  const page = await (await browser.newContext({ viewport: { width: 1600, height: 1000 } })).newPage();
  const pageErrors = [];
  page.on('pageerror', e => pageErrors.push(String(e).slice(0, 200)));

  await page.goto(`http://127.0.0.1:${PORT}/`, { timeout: 60000 });
  await page.waitForFunction(() => {
    const box = document.querySelector('#app-contents');
    return !!box && box.offsetParent !== null;
  }, null, { timeout: 90000 });
  console.log(ts(), 'app loaded');

  await page.click('#app-assistant-toggle');
  await page.waitForFunction(() => {
    const p = document.querySelector('#app-assistant-panel');
    return p && getComputedStyle(p).display !== 'none';
  }, null, { timeout: 15000 });
  await page.waitForFunction(() => {
    const body = document.querySelector('#app-assistant-body');
    return body && body.innerText.length > 0;
  }, null, { timeout: 15000 });
  const bt = await page.innerText('#app-assistant-body');
  if (/not installed|Configure a model|Select and load a dataset/i.test(bt))
    throw new Error('assistant not ready: ' + bt.slice(0, 200));
  console.log(ts(), 'assistant drawer open');

  const uiSnapshot = () => page.evaluate(() => {
    const vals = Shiny.shinyapp.$inputValues;
    const out = {};
    for (const k of Object.keys(vals)) {
      if (/\.clientdata/.test(k)) continue;
      const v = vals[k];
      if (typeof v === 'string' || typeof v === 'number' || typeof v === 'boolean' || v === null)
        out[k] = v;
    }
    const tabs = {};
    for (const li of document.querySelectorAll('.nav-tabs > .active > a, .navbar-nav > .active > a'))
      tabs[li.getAttribute('data-value') || li.innerText.trim()] = true;
    return { inputs: out, active_tabs: Object.keys(tabs) };
  });

  const sendTurn = async (turn) => {
    const rec = { id: turn.id, prompt: turn.prompt, started: new Date().toISOString() };
    const tTurn = Date.now();

    // optional mid-session logging checkbox exercise (off -> on, no chat
    // in between so no events are lost)
    if (turn.toggle_logging_before) {
      try {
        await page.click('#app-assistant-settings');
        await page.waitForTimeout(1500);
        const cb = page.locator('#app-assistant-enable_logging, input[id*="enable_logging"]');
        if (await cb.count()) { await cb.first().click(); await page.waitForTimeout(800); await cb.first().click(); }
        const cancel = page.locator('button:has-text("Cancel"), .modal button[data-dismiss]');
        if (await cancel.count()) await cancel.first().click().catch(() => {});
        await page.waitForTimeout(800);
        rec.logging_toggled = true;
      } catch (e) { rec.logging_toggled = 'failed: ' + String(e).slice(0, 80); }
    }

    const before = await uiSnapshot();
    const textBefore = await page.innerText('#app-assistant-body').catch(() => '');
    const composer = page.locator(
      '#app-assistant-panel textarea, #app-assistant-panel [contenteditable="true"]');
    await composer.first().waitFor({ state: 'visible', timeout: 20000 });

    // deliver the prompt; if no user bubble appears within 25s the
    // submission was swallowed (busy composer) - wait for idle and retry once
    const deliver = async () => {
      await composer.first().click();
      await composer.first().fill(turn.prompt);
      await page.keyboard.press('Enter');
    };
    await deliver();
    console.log(ts(), `>> [${turn.id}] "${turn.prompt.slice(0, 70)}"`);
    const delivered = await page.waitForFunction(prevLen =>
      (document.querySelector('#app-assistant-body')?.innerText || '').length > prevLen + 5,
      textBefore.length, { timeout: 25000 }).then(() => true).catch(() => false);
    if (!delivered) {
      rec.warn_submission = 'no user bubble in 25s; waiting idle + retrying';
      console.log(ts(), `   [${turn.id}] submission swallowed - waiting for idle...`);
      const ok = await waitIdle(page, 420000);
      rec.waited_idle_before_retry = ok;
      await deliver();
      const ok2 = await page.waitForFunction(prevLen =>
        (document.querySelector('#app-assistant-body')?.innerText || '').length > prevLen + 5,
        textBefore.length, { timeout: 30000 }).then(() => true).catch(() => false);
      if (!ok2) { rec.warn_submission += '; retry also failed'; rec.submission_failed = true; return rec; }
    }

    // settle: streaming / in-flight request / stable text
    let prev = '', stableRounds = 0, figureReady = false;
    const deadline = Date.now() + (turn.budget_ms || 300000);
    while (Date.now() < deadline) {
      await new Promise(s => setTimeout(s, 3000));
      const st = await page.evaluate(() => ({
        text: (document.querySelector('#app-assistant-body')?.innerText || ''),
        figures: document.querySelectorAll('.omicsviewer-ai-figure img').length,
        figLoaded: Array.from(document.querySelectorAll('.omicsviewer-ai-figure img'))
          .every(i => i.complete && i.naturalWidth > 0),
        streaming: (() => {
          const b = document.querySelector('[class*="cancel"]');
          return !!(b && b.offsetParent !== null && !b.disabled);
        })()
      }));
      if (st.figures > 0 && st.figLoaded) figureReady = true;
      if (st.streaming || logRequestInFlight()) { stableRounds = 0; prev = st.text; continue; }
      if (st.text === prev && st.text.length > 200) {
        stableRounds++;
        if (stableRounds >= (figureReady ? 6 : 8)) {
          rec.final_text = st.text.slice(textBefore.length).slice(0, 4000);
          rec.figures = st.figures;
          break;
        }
      } else stableRounds = 0;
      prev = st.text;
    }
    if (!rec.final_text) { rec.final_text = '(TIMEOUT)'; rec.timed_out = true; }

    rec.duration_ms = Date.now() - tTurn;
    const after = await uiSnapshot();
    const changed = {};
    for (const k of Object.keys(after.inputs))
      if (JSON.stringify(after.inputs[k]) !== JSON.stringify(before.inputs[k]))
        changed[k] = { from: before.inputs[k], to: after.inputs[k] };
    rec.ui_changes = changed;
    rec.active_tabs_after = after.active_tabs;
    await page.screenshot({ path: path.join(ART, `${String(turn.idx).padStart(2, '0')}-${turn.id}.png`) });
    console.log(ts(), `<< [${turn.id}] settled in ${(rec.duration_ms / 1000).toFixed(0)}s, ` +
      `${Object.keys(changed).length} ui changes, figures=${rec.figures || 0}` +
      (rec.timed_out ? ' (TIMEOUT)' : ''));
    return rec;
  };

  const records = [];
  for (let i = 0; i < TURNS.length; i++) {
    const turn = { ...TURNS[i], idx: i + 1 };
    try {
      // never deliver into a busy chat: wait for a previous chain to close
      const idle = await waitIdle(page, 420000);
      if (!idle) records.push({ id: turn.id, error: 'chat never went idle' });
      records.push(await sendTurn(turn));
    }
    catch (e) { records.push({ id: turn.id, error: String(e).slice(0, 400) }); }
  }

  records.push({ page_errors: pageErrors });
  writeFileSync(path.join(OUT, 'session_report.json'), JSON.stringify(records, null, 2));
  // copy the session's diagnostic log next to the report
  try { copyFileSync(path.join(LOGDIR, latestLogName()), path.join(OUT, 'session_llm_log.jsonl')); }
  catch {}
  writeFileSync(path.join(OUT, 'app_stderr.log'), rLog.slice(-20000));
  console.log(ts(), `DONE. ${records.length - 1} turns recorded -> ${OUT}`);
  await browser.close();
} catch (e) {
  console.error('FATAL:', e);
  writeFileSync(path.join(OUT, 'fatal.txt'), String(e) + '\n---\n' + rLog.slice(-5000));
  process.exitCode = 1;
} finally { cleanup(); }
