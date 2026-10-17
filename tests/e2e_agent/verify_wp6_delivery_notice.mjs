// WP-6 (F2/F8) verification: composer delivery watchdog + Retry.
//
// Scenario A (notice): arm the watchdog via the env-gated JS hook (no
// provider needed) - the inline notice must appear after ~10 s.
// Scenario B (retry clears): click Retry -> the server re-drives the
// submission through update_user_input(submit = TRUE) -> the delivery
// ack clears the notice and the message reaches the server. Uses a
// black-hole provider endpoint: the ack + user_message still fire
// (delivery succeeded); the provider request fails visibly instead.
// Scenario C (no false positive): a normal submit acks and never shows
// the notice.
//
// Run: NOT_CRAN=true node verify_wp6_delivery_notice.mjs   (exit 0 = pass)
import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import { mkdirSync, readFileSync, readdirSync, statSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../..');
const EXTDATA = path.join(REPO, 'inst/extdata');
const OUT = '/home/chen/.omicsviewer-llm-eval/wp6_verify';
mkdirSync(OUT, { recursive: true });

const env = { ...process.env };
const PORT = 7789;
const r = spawn('Rscript', ['-e', `
  Sys.setenv(
    OMICSVIEWER_LLM_PROVIDER = 'openai',
    OMICSVIEWER_LLM_MODEL = 'test-model',
    OMICSVIEWER_LLM_API_KEY = 'sk-test-only',
    OMICSVIEWER_LLM_BASE_URL = 'http://127.0.0.1:9/v1',
    OMICSVIEWER_LLM_ALLOW_USER_ENDPOINT = 'true')
  options(shiny.port = ${PORT}, shiny.host = '127.0.0.1')
  eset <- readRDS(file.path('${EXTDATA}', 'demo.RDS'))
  omicsViewer::omicsViewer(dir = '${EXTDATA}', ESVObj = eset)
`], { cwd: REPO, env: { ...env, OMICSVIEWER_TEST_HOOKS: 'true',
     OMICSVIEWER_LLM_LOG: 'true', OMICSVIEWER_LLM_LOG_DIR: OUT },
     stdio: ['ignore', 'pipe', 'pipe'] });
let rLog = '';
r.stdout.on('data', d => { rLog += d; });
r.stderr.on('data', d => { rLog += d; });
process.on('exit', () => { try { r.kill('SIGKILL'); } catch {} });
const waitPort = async (ms = 90000) => {
  const t0 = Date.now();
  while (Date.now() - t0 < ms) {
    try { const res = await fetch(`http://127.0.0.1:${PORT}/`); if (res.status > 0) return; } catch {}
    await new Promise(s => setTimeout(s, 500));
  }
  throw new Error(`app did not start:\n${rLog}`);
};
// NOTE: import fs functions at the top (ESM require trap).
const liveLog = () => {
  const fs = readdirSync(OUT).filter(f => f.endsWith('.jsonl'))
    .map(f => ({ f, m: statSync(path.join(OUT, f)).mtimeMs }))
    .sort((a, b) => b.m - a.m);
  return fs.length ? path.join(OUT, fs[0].f) : null;
};
const events = () => {
  const f = liveLog();
  if (!f) return [];
  try { return readFileSync(f, 'utf8').split('\n').filter(Boolean).map(l => JSON.parse(l)); }
  catch { return []; }
};

let failures = 0;
const check = (name, cond) => {
  console.log(`${cond ? 'PASS' : 'FAIL'} - ${name}`);
  if (!cond) failures++;
};

await waitPort();
const browser = await chromium.launch({
  executablePath: '/usr/bin/google-chrome', headless: true,
  args: ['--no-sandbox', '--disable-gpu', '--disable-dev-shm-usage',
         '--disable-webgl', '--disable-webgl2']
});
const page = await (await browser.newContext({ viewport: { width: 1600, height: 1000 } })).newPage();
await page.goto(`http://127.0.0.1:${PORT}/`, { timeout: 60000 });
await page.waitForFunction(() => !!document.querySelector('#app-contents')?.offsetParent, null, { timeout: 90000 });
await page.click('#app-assistant-toggle');
await page.waitForFunction(() => getComputedStyle(document.querySelector('#app-assistant-panel')).display !== 'none', null, { timeout: 15000 });
await page.waitForFunction(() => {
  const b = document.querySelector('#app-assistant-body');
  return b && b.innerText.length > 0;
}, null, { timeout: 15000 });
console.log('app + assistant ready');

// hook must exist (TEST_HOOKS on)
check('delivery watchdog hook exposed',
  await page.evaluate(() => typeof window.__omicsviewerAiDelivery === 'object'));

// --- Scenario A: arm -> notice within ~12 s --------------------------------
await page.evaluate(() => window.__omicsviewerAiDelivery.arm('probe word zeta'));
const noticed = await page.waitForFunction(
  () => !!document.querySelector('.omicsviewer-ai-delivery-notice'),
  null, { timeout: 14000 }).then(() => true).catch(() => false);
check('notice appears after 10 s without ack', noticed);

// --- Scenario B: Retry re-drives through the server ------------------------
const msgBefore = events().filter(e => e.event === 'user_message').length;
if (noticed) await page.click('.omicsviewer-ai-delivery-notice button');
// update_user_input -> composer submit -> ack clears the notice; the
// provider call then fails fast against the dead endpoint
await page.waitForTimeout(9000);
const noticeGone = await page.evaluate(
  () => !document.querySelector('.omicsviewer-ai-delivery-notice'));
check('retry clears the notice (ack round-trip)', noticeGone);
const msgAfter = events().filter(e => e.event === 'user_message').length;
check('retried message reached the server (user_message logged)',
  msgAfter > msgBefore);
const retried = events().some(e =>
  e.event === 'user_message' &&
  JSON.stringify(e.details || {}).includes('probe word zeta'));
check('retried message text matches', retried);

// --- Scenario C: normal submit acks, no false notice ------------------------
// wait for scenario B's request to settle against the dead endpoint first
await page.waitForFunction(() => {
  const b = document.querySelector('[class*="cancel"]');
  return !(b && b.offsetParent !== null && !b.disabled);
}, null, { timeout: 60000 });
await page.waitForTimeout(2000);
const composer = page.locator('#app-assistant-panel textarea, #app-assistant-panel [contenteditable="true"]');
await composer.first().click();
await composer.first().fill('normal delivery probe');
await page.keyboard.press('Enter');
await page.waitForTimeout(4000);
// the ack arrives within ~1 s; give the 10 s timer room to prove it stays hidden
await page.waitForTimeout(9000);
const noNotice = await page.evaluate(
  () => !document.querySelector('.omicsviewer-ai-delivery-notice'));
check('normal submit never shows the notice', noNotice);
const delivered = events().some(e =>
  e.event === 'user_message' &&
  JSON.stringify(e.details || {}).includes('normal delivery probe'));
check('normal submit delivered', delivered);

console.log(failures === 0 ? '\nWP-6 verification: ALL PASS' : `\nWP-6 verification: ${failures} FAILURES`);
await browser.close();
process.exit(failures === 0 ? 0 : 1);
