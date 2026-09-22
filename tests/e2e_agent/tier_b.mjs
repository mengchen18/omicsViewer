// WP0 Tier B — end-to-end agent validation (real provider, env-gated).
//
// Drives the actual chat UI: opens the assistant drawer, sends a prompt,
// waits for the tool-driven UI effect / figure, and archives screenshots +
// the session's diagnostic log for inspection.
//
// Run:  node tier_b.mjs "your prompt here"
// Creds: tests/e2e_agent/provider.env (gitignored) is sourced into the
//        spawned app process; nothing is hardcoded here.

import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import { mkdirSync, copyFileSync, existsSync, readFileSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../..');
const EXTDATA = path.join(REPO, 'inst/extdata');
const ARTIFACTS = path.join(here, 'artifacts');
mkdirSync(ARTIFACTS, { recursive: true });

const PROMPT = process.argv[2] ||
  'for the volcano plot on left panel, generate a high quality png';
const PORT = 7780;
const LOGDIR = '/tmp/tierb-llm-logs';
mkdirSync(LOGDIR, { recursive: true });

// ---- provider credentials from provider.env --------------------------------
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

// ---- spawn app --------------------------------------------------------------
const r = spawn('Rscript', ['-e', `
  options(shiny.port = ${PORT}, shiny.host = '127.0.0.1')
  eset <- readRDS(file.path('${EXTDATA}', 'demo.RDS'))
  omicsViewer::omicsViewer(dir = '${EXTDATA}', ESVObj = eset)
`], {
  cwd: REPO,
  env: {
    ...env,
    OMICSVIEWER_LLM_LOG: 'true',
    OMICSVIEWER_LLM_LOG_DIR: LOGDIR
  },
  stdio: ['ignore', 'pipe', 'pipe']
});
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

const t0 = Date.now();
const ts = () => `t+${String(Math.round((Date.now() - t0) / 100) / 10)}s`;

try {
  await waitPort(PORT);
  const browser = await chromium.launch({
    executablePath: '/usr/bin/google-chrome',
    headless: true,
    args: ['--no-sandbox', '--disable-gpu', '--disable-dev-shm-usage',
           '--disable-webgl', '--disable-webgl2']
  });
  const page = await (await browser.newContext()).newPage();
  page.on('pageerror', e => console.log(ts(), 'PAGEERROR:', String(e).slice(0, 150)));

  await page.goto(`http://127.0.0.1:${PORT}/`, { timeout: 60000 });
  await page.waitForFunction(() => {
    const box = document.querySelector('#app-contents');
    return !!box && box.offsetParent !== null;
  }, null, { timeout: 90000 });
  console.log(ts(), 'app loaded, dataset preloaded');

  // open the assistant drawer
  await page.click('#app-assistant-toggle');
  // the drawer panel is position:fixed, so offsetParent is always null;
  // shinyjs::show flips display:none -> block
  await page.waitForFunction(() => {
    const p = document.querySelector('#app-assistant-panel');
    return p && getComputedStyle(p).display !== 'none';
  }, null, { timeout: 15000 });
  console.log(ts(), 'assistant drawer open');

  // the drawer shows either the chat or a setup notice
  await page.waitForFunction(() => {
    const body = document.querySelector('#app-assistant-body');
    return body && body.innerText.length > 0;
  }, null, { timeout: 15000 });
  const bodyText = await page.innerText('#app-assistant-body');
  if (/not installed|Configure a model|Select and load a dataset/i.test(bodyText)) {
    throw new Error('assistant not ready: ' + bodyText.slice(0, 200));
  }
  const status = await page.innerText('#app-assistant-status').catch(() => '');
  console.log(ts(), 'assistant status:', status);

  // find the composer (shinychat web component) and type the prompt
  const composer = page.locator(
    '#app-assistant-panel textarea, #app-assistant-panel [contenteditable="true"]');
  await composer.first().waitFor({ state: 'visible', timeout: 20000 });
  await composer.first().click();
  await composer.first().fill(PROMPT);
  await page.keyboard.press('Enter');
  console.log(ts(), `prompt sent: "${PROMPT}"`);

  await page.screenshot({ path: path.join(ARTIFACTS, 'tier_b_prompt_sent.png') });

  // confirm submission: a user bubble must appear beyond the greeting
  await page.waitForFunction(() => {
    const body = document.querySelector('#app-assistant-body');
    if (!body) return false;
    const users = body.querySelectorAll('.chat-message-user, [class*="user"]');
    return users.length > 0 || /volcano plot on left panel/.test(body.innerText);
  }, null, { timeout: 20000 }).catch(() => console.log(ts(), 'WARN: user bubble not detected'));

  // settle: poll chat text; done when a figure rendered or text is stable
  // across consecutive polls (greeting excluded) with a generous LLM budget.
  // The shinychat cancel-control probe misses some active streams (observed
  // with glm flash: long silent continuation after tool results), so also
  // treat "newest diagnostic log ends with provider_request_start" as
  // active - that event has no matching assistant_response yet.
  const { readdirSync: lsSync, statSync: stSync } = await import('node:fs');
  const logRequestInFlight = () => {
    try {
      const logs = lsSync(LOGDIR).filter(f => f.endsWith('.jsonl'))
        .map(f => ({ f, mtime: stSync(path.join(LOGDIR, f)).mtimeMs }))
        .sort((a, b) => b.mtime - a.mtime);
      if (!logs.length) return false;
      const lines = readFileSync(path.join(LOGDIR, logs[0].f), 'utf8')
        .split('\n').filter(l => l.trim());
      if (!lines.length) return false;
      const last = JSON.parse(lines[lines.length - 1]);
      return last.event === 'provider_request_start';
    } catch { return false; }
  };
  let prev = '', stableRounds = 0, figureCount = 0, figureReady = false;
  const deadline = Date.now() + 300000;
  while (Date.now() < deadline) {
    await new Promise(s => setTimeout(s, 3000));
    const st = await page.evaluate(() => ({
      text: (document.querySelector('#app-assistant-body')?.innerText || ''),
      figures: document.querySelectorAll('.omicsviewer-ai-figure img').length,
      figLoaded: Array.from(document.querySelectorAll('.omicsviewer-ai-figure img'))
        .every(i => i.complete && i.naturalWidth > 0),
      // shinychat shows a cancel control while a stream is active
      streaming: (() => {
        const b = document.querySelector('[class*="cancel"]');
        return !!(b && b.offsetParent !== null && !b.disabled);
      })()
    }));
    figureCount = st.figures;
    if (st.figures > 0 && st.figLoaded) figureReady = true;
    if (st.streaming || logRequestInFlight()) { stableRounds = 0; prev = st.text; continue; }
    if (st.text === prev && st.text.length > 200) {
      stableRounds++;
      // once a figure is rendered, settle after 18s of quiet so a pending
      // revision (second figure) is still captured; text-only tasks use
      // the original 24s stability rule
      const needed = figureReady ? 6 : 8;
      if (stableRounds >= needed) {
        console.log(ts(), figureReady ? 'chat settled (figure ready, 18s quiet)'
                                      : 'chat settled (24s stable)');
        break;
      }
    } else stableRounds = 0;
    prev = st.text;
  }
  const finalState = { figureCount, chatText: prev.slice(-600) };
  console.log(ts(), 'figure count in drawer:', finalState.figureCount);
  console.log(ts(), 'chat tail:', JSON.stringify(finalState.chatText.slice(-300)));

  await page.screenshot({ path: path.join(ARTIFACTS, 'tier_b_final.png'), fullPage: false });
  await browser.close();
  cleanup();

  // archive the diagnostic log of this run
  const { readdirSync, statSync } = await import('node:fs');
  const logs = readdirSync(LOGDIR).filter(f => f.endsWith('.jsonl'))
    .map(f => ({ f, mtime: statSync(path.join(LOGDIR, f)).mtimeMs }))
    .sort((a, b) => b.mtime - a.mtime);
  if (logs.length) {
    copyFileSync(path.join(LOGDIR, logs[0].f), path.join(ARTIFACTS, 'tier_b_last_log.jsonl'));
    console.log('diagnostic log archived to artifacts/tier_b_last_log.jsonl');
  }
  process.exit(finalState.figureCount > 0 ? 0 : 1);
} catch (e) {
  console.error(ts(), 'FAILED:', e.message);
  console.error('--- R log tail ---\n' + rLog.split('\n').slice(-20).join('\n'));
  cleanup();
  process.exit(1);
}
