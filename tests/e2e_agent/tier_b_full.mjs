// WP7 gate — full Tier B benchmark driver (12 core tasks, multi-prompt aware).
//
// Extends tier_b.mjs: drives the whole core task list (plan §3 WP5) against
// demo.RDS with a FRESH app process + fresh chat per task, supports multi-
// prompt tasks (setup + probe in one conversation: tasks 8/10/11), records
// DOM observations + the read-only state overview after every prompt, and
// archives per-task diagnostic logs + screenshots for the gate scorer
// (score_tier_b_gate.R).
//
// Run:   node tier_b_full.mjs --run 1 [--tasks 10,11] [--label gate]
// Creds: tests/e2e_agent/provider.env (gitignored), sourced into each app.
//
// Task 11 note: the demo dataset has no "category" feature column (only
// Gene.name / Protein.ID / StringDB IDs are discrete), so the probe colors
// by a real column (Intensity) — same revision skill (aesthetic remap via
// update_figure), dataset-faithful wording. Documented in tier_b_tasks.md.

import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import { mkdirSync, copyFileSync, existsSync, readFileSync, writeFileSync, readdirSync, statSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../..');
const EXTDATA = path.join(REPO, 'inst/extdata');

// ---- CLI -------------------------------------------------------------------
const argv = process.argv.slice(2);
const arg = (name, dflt) => {
  const i = argv.indexOf(`--${name}`);
  if (i < 0) return dflt;
  const v = argv[i + 1];
  return (v && !v.startsWith('--')) ? v : true;
};
const RUN = parseInt(arg('run', '1'), 10);
const LABEL = arg('label', 'gate');
const TASKS = arg('tasks', null); // e.g. "10,11"

// ---- benchmark tasks (plan §3 WP5 core set; dataset-adapted wording) ------
const ALL_TASKS = [
  { id: 1,  prompts: ['What dataset is loaded?'] },
  { id: 2,  prompts: ['How many genes are currently selected?'] },
  { id: 3,  prompts: ['Switch to the Sample tab'] },
  { id: 4,  prompts: ['Make a volcano plot of the current results'] },
  { id: 5,  prompts: ['Plot log FDR versus mean difference for the RE vs ME comparison with custom axes'] },
  { id: 6,  prompts: ["Find genes matching 'kinase' and select the first five"] },
  { id: 7,  prompts: ['Summarize the sample group annotation'] },
  { id: 8,  prompts: ["Find genes matching 'kinase' and select the first five",
                      'Make a boxplot of the expression of the selected genes grouped by sample group'] },
  { id: 9,  prompts: ['Create a histogram of the t-test log FDR'] },
  { id: 10, prompts: ['Create a histogram of the t-test log FDR',
                      "Change the last figure's title to 'Expression-wide t-test significance'"] },
  { id: 11, prompts: ['Create a volcano plot of the t-test RE vs ME results',
                      'Color the volcano points by intensity'] },
  { id: 12, prompts: ['What can you change in this app?'] },
  // dataset-faithful variants of tasks 6/8 (run 1 exposed that 'kinase'
  // matches nothing in demo.RDS ids/annotations — search returns 0 hits
  // and the honest model asks for clarification instead of selecting;
  // 'MAPK' matches 10 real ids, so these measure the actual find+select
  // skill). Not part of the default set: request with --tasks 106,108.
  { id: 106, prompts: ["Find genes matching 'MAPK' and select the first five"] },
  { id: 108, prompts: ["Find genes matching 'MAPK' and select the first five",
                       'Make a boxplot of the expression of the selected genes grouped by sample group'] }
];
const tasks = TASKS
  ? ALL_TASKS.filter(t => TASKS.split(',').map(s => parseInt(s, 10)).includes(t.id))
  : ALL_TASKS;

// ---- provider credentials ---------------------------------------------------
const env = { ...process.env };
if (existsSync(path.join(here, 'provider.env'))) {
  for (const line of readFileSync(path.join(here, 'provider.env'), 'utf8').split('\n')) {
    const m = line.match(/^\s*([A-Z_][A-Z0-9_]*)\s*=\s*(.*?)\s*$/);
    if (m && !line.trim().startsWith('#')) env[m[1]] = m[2];
  }
}
if (!env.OMICSVIEWER_LLM_API_KEY) {
  console.error('FATAL: no OMICSVIEWER_LLM_API_KEY (fill tests/e2e_agent/provider.env)');
  process.exit(2);
}
console.log(`provider: ${env.OMICSVIEWER_LLM_PROVIDER} model: ${env.OMICSVIEWER_LLM_MODEL} ` +
  `base: ${env.OMICSVIEWER_LLM_BASE_URL || 'default'}`);

const OUTDIR = path.join(here, 'artifacts', LABEL, `run${RUN}`);
mkdirSync(OUTDIR, { recursive: true });

const HOOK = {
  box: '#app-agentTestHooks-container',
  op: '#app-agentTestHooks-op',
  payload: '#app-agentTestHooks-payload',
  run: '#app-agentTestHooks-run',
  result: '#app-agentTestHooks-result'
};

const newestLog = (dir) => {
  try {
    const logs = readdirSync(dir).filter(f => f.endsWith('.jsonl'))
      .map(f => ({ f, mtime: statSync(path.join(dir, f)).mtimeMs }))
      .sort((a, b) => b.mtime - a.mtime);
    return logs.length ? path.join(dir, logs[0].f) : null;
  } catch { return null; }
};

const logRequestInFlight = (dir) => {
  const f = newestLog(dir);
  if (!f) return false;
  try {
    const lines = readFileSync(f, 'utf8').split('\n').filter(l => l.trim());
    if (!lines.length) return false;
    // request still active: the newest event is a request start with no
    // assistant_response yet, OR the stream is still marked streaming
    const last = JSON.parse(lines[lines.length - 1]);
    if (last.event === 'provider_request_start') return true;
    for (let i = lines.length - 1; i >= 0; i--) {
      const ev = JSON.parse(lines[i]);
      if (ev.event === 'stream_status') return ev.details && ev.details.status === 'streaming';
      if (ev.event === 'assistant_response') return false;
      if (ev.event === 'user_message') return false;
    }
    return false;
  } catch { return false; }
};

// A turn is COMPLETE when the newest request has an assistant_response and
// the stream went idle — the log is the authoritative completion signal, so
// the fixed quiet window shrinks to a short double-check (2 rounds / 6 s).
const logTurnComplete = (dir) => {
  const f = newestLog(dir);
  if (!f) return false;
  try {
    const lines = readFileSync(f, 'utf8').split('\n').filter(l => l.trim());
    let sawResponse = false;
    for (let i = lines.length - 1; i >= 0; i--) {
      const ev = JSON.parse(lines[i]);
      if (ev.event === 'user_message') return sawResponse;
      if (ev.event === 'assistant_response') sawResponse = true;
      if (ev.event === 'provider_request_start') return false;
      if (ev.event === 'stream_status' && sawResponse && ev.details && ev.details.status === 'idle') return true;
    }
    return false;
  } catch { return false; }
};

// ---- one task = one fresh app + fresh chat ----------------------------------
async function runTask(browser, task, taskIndex) {
  const PORT = 7780 + taskIndex;
  const TASKDIR = path.join(OUTDIR, `task${String(task.id).padStart(2, '0')}`);
  const LOGDIR = path.join(TASKDIR, 'logs');
  mkdirSync(LOGDIR, { recursive: true });

  const r = spawn('Rscript', ['-e', `
    options(shiny.port = ${PORT}, shiny.host = '127.0.0.1')
    eset <- readRDS(file.path('${EXTDATA}', 'demo.RDS'))
    omicsViewer::omicsViewer(dir = '${EXTDATA}', ESVObj = eset)
  `], {
    cwd: REPO,
    env: { ...env, OMICSVIEWER_LLM_LOG: 'true', OMICSVIEWER_LLM_LOG_DIR: LOGDIR,
           OMICSVIEWER_TEST_HOOKS: 'true', OMICSVIEWER_LLM_MAX_REQUESTS: '8' },
    stdio: ['ignore', 'pipe', 'pipe']
  });
  let rLog = '';
  r.stdout.on('data', d => { rLog += d; });
  r.stderr.on('data', d => { rLog += d; });
  const cleanup = () => { try { r.kill('SIGKILL'); } catch {} };
  process.on('exit', cleanup);

  const record = { task_id: task.id, prompts: task.prompts, observations: [], status: 'ok', error: null };
  const t0 = Date.now();
  const ts = () => `t+${String(Math.round((Date.now() - t0) / 100) / 10)}s`;

  let page = null, ctx = null;
  try {
    // wait for the app port
    const t1 = Date.now();
    while (Date.now() - t1 < 90000) {
      try { const res = await fetch(`http://127.0.0.1:${PORT}/`); if (res.status > 0) break; } catch {}
      await new Promise(s => setTimeout(s, 500));
    }
    if (Date.now() - t1 >= 90000) throw new Error('app did not start');

    ctx = await browser.newContext();
    page = await ctx.newPage();
    await page.goto(`http://127.0.0.1:${PORT}/`, { timeout: 60000 });
    await page.waitForFunction(() => {
      const tab = Shiny.shinyapp.$inputValues['app-dataspace-eset'];
      const box = document.querySelector('#app-contents');
      return !!tab && !!box && box.offsetParent !== null;
    }, null, { timeout: 90000 });
    console.log(ts(), `task ${task.id}: app loaded`);

    // expose the (hidden) read-only test-hook panel for observation calls
    await page.evaluate(sel => { const el = document.querySelector(sel); if (el) el.style.display = 'block'; }, HOOK.box);

    const observe = async () => page.evaluate(async () => {
      const tab = (window.Shiny && Shiny.shinyapp && Shiny.shinyapp.$inputValues['app-dataspace-eset']) || null;
      const figs = Array.from(document.querySelectorAll('.omicsviewer-ai-figure img'));
      const pt = Array.from(document.querySelectorAll('.js-plotly-plot'))
        .filter(e => e.offsetParent !== null && e._fullLayout);
      const ax = (a) => { const v = a && a.title ? (a.title.text || a.title) : ''; return typeof v === 'string' ? v : ''; };
      // read-only state overview through the test hook (selection counts)
      const el = document.querySelector('#app-agentTestHooks-op');
      if (el && el.selectize) el.selectize.setValue('overview');
      const payload = document.querySelector('#app-agentTestHooks-payload');
      const result = document.querySelector('#app-agentTestHooks-result');
      let overview = null;
      if (el && payload && result) {
        const before = result.innerText;
        payload.value = '{}';
        document.querySelector('#app-agentTestHooks-run').click();
        for (let i = 0; i < 40; i++) {
          await new Promise(s => setTimeout(s, 250));
          const t = result.innerText;
          if (t !== before && t.includes('hook_run')) {
            try { overview = JSON.parse(t.trim()); } catch { overview = { parse_error: true }; }
            break;
          }
        }
      }
      return {
        tab,
        figures: figs.length,
        figures_loaded: figs.every(i => i.complete && i.naturalWidth > 0),
        axis_titles: pt.length && pt[0]._fullLayout ? [ax(pt[0]._fullLayout.xaxis), ax(pt[0]._fullLayout.yaxis)] : null,
        feature_selected: overview && overview.selection && overview.selection.features
          ? overview.selection.features.count : null,
        sample_selected: overview && overview.selection && overview.selection.samples
          ? overview.selection.samples.count : null
      };
    });

    // open the assistant drawer
    await page.click('#app-assistant-toggle');
    await page.waitForFunction(() => {
      const p = document.querySelector('#app-assistant-panel');
      return p && getComputedStyle(p).display !== 'none';
    }, null, { timeout: 15000 });
    await page.waitForFunction(() => {
      const body = document.querySelector('#app-assistant-body');
      return body && body.innerText.length > 0;
    }, null, { timeout: 15000 });
    const bodyText = await page.innerText('#app-assistant-body');
    if (/not installed|Configure a model|Select and load a dataset/i.test(bodyText)) {
      throw new Error('assistant not ready: ' + bodyText.slice(0, 200));
    }
    console.log(ts(), `task ${task.id}: drawer open`);

    record.baseline = await observe();

    // send prompts sequentially, settling between them
    for (let pi = 0; pi < task.prompts.length; pi++) {
      const prompt = task.prompts[pi];
      const composer = page.locator(
        '#app-assistant-panel textarea, #app-assistant-panel [contenteditable="true"]');
      await composer.first().waitFor({ state: 'visible', timeout: 20000 });
      await composer.first().click();
      await composer.first().fill(prompt);
      await page.keyboard.press('Enter');
      console.log(ts(), `task ${task.id}: prompt ${pi + 1}/${task.prompts.length} sent: "${prompt.slice(0, 60)}"`);

      // settle: the diagnostic log is the authoritative turn-completion
      // signal (assistant_response + stream idle); once it says complete,
      // confirm with 2 short stable-text polls (6 s) so late figure PNGs
      // and trailing tokens are captured. Streams/requests in flight reset.
      let prev = '', stableRounds = 0, figureReady = false;
      const needed = 2;
      const deadline = Date.now() + 300000;
      while (Date.now() < deadline) {
        await new Promise(s => setTimeout(s, 3000));
        const st = await page.evaluate(() => {
          const figs = Array.from(document.querySelectorAll('.omicsviewer-ai-figure img'));
          const cancel = document.querySelector('[class*="cancel"]');
          return {
            text: (document.querySelector('#app-assistant-body')?.innerText || ''),
            figures: figs.length,
            figLoaded: figs.every(i => i.complete && i.naturalWidth > 0),
            streaming: !!(cancel && cancel.offsetParent !== null && !cancel.disabled)
          };
        });
        if (st.figures > 0 && st.figLoaded) figureReady = true;
        const busy = st.streaming || logRequestInFlight(LOGDIR) || !logTurnComplete(LOGDIR);
        if (busy) { stableRounds = 0; prev = st.text; continue; }
        if (st.text === prev && st.text.length > 200) {
          stableRounds++;
          if (stableRounds >= needed) break;
        } else stableRounds = 0;
        prev = st.text;
      }

      const obs = await observe();
      obs.prompt = prompt;
      obs.chat_tail = prev.slice(-1500);
      obs.figure_ready = figureReady;
      record.observations.push(obs);
      await page.screenshot({ path: path.join(TASKDIR, `prompt${pi + 1}.png`) });
      console.log(ts(), `task ${task.id}: prompt ${pi + 1} settled — tab=${obs.tab} figures=${obs.figures} featSel=${obs.feature_selected}`);
    }

    await page.screenshot({ path: path.join(TASKDIR, 'final.png'), fullPage: false });
  } catch (e) {
    record.status = 'error';
    record.error = String(e.message || e);
    console.error(ts(), `task ${task.id} FAILED:`, record.error);
    console.error('--- R log tail ---\n' + rLog.split('\n').slice(-15).join('\n'));
    try { if (page) await page.screenshot({ path: path.join(TASKDIR, 'error.png'), fullPage: false }); } catch {}
  } finally {
    try { if (ctx) await ctx.close(); } catch {}
    cleanup();
    process.removeListener('exit', cleanup);
    // archive the diagnostic log for this task's app session
    const f = newestLog(LOGDIR);
    if (f) { try { copyFileSync(f, path.join(TASKDIR, 'log.jsonl')); } catch {} }
    writeFileSync(path.join(TASKDIR, 'record.json'), JSON.stringify(record, null, 2));
    // give the OS a beat to release the port
    await new Promise(s => setTimeout(s, 1500));
  }
  return record;
}

// ---- main -------------------------------------------------------------------
const results = [];
const browser = await chromium.launch({
  executablePath: '/usr/bin/google-chrome',
  headless: true,
  args: ['--no-sandbox', '--disable-gpu', '--disable-dev-shm-usage',
         '--disable-webgl', '--disable-webgl2']
});
try {
  for (let i = 0; i < tasks.length; i++) {
    console.log(`\n=== run ${RUN} task ${tasks[i].id} (${i + 1}/${tasks.length}) ===`);
    const rec = await runTask(browser, tasks[i], i);
    results.push(rec);
    writeFileSync(path.join(OUTDIR, 'records.json'), JSON.stringify(results, null, 2));
  }
} finally {
  await browser.close();
}
console.log(`\nrun ${RUN} complete: ${results.length} tasks, ${results.filter(r => r.status === 'ok').length} ok`);
