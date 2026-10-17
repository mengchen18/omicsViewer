// WP-8 (F7) verification: progress line under the composer.
//
// Scenario A (server pipeline, no real provider): with a black-hole
// provider endpoint, a submitted message still fires the server-side
// submit push (Working - waiting for the model...), the request-start
// push (model phase), and the settle push (line hides) - proving the
// whole custom-message plumbing without a live model.
// Scenario B (step rendering): env-gated hook drives the exact tool-step
// render path (Working - step N: tool) and the idle clear.
//
// Run: NOT_CRAN=true node verify_wp8_progress.mjs   (exit 0 = pass)
import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import { mkdirSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../..');
const EXTDATA = path.join(REPO, 'inst/extdata');

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
`], { cwd: REPO, env: { ...env, OMICSVIEWER_TEST_HOOKS: 'true' },
     stdio: ['ignore', 'pipe', 'pipe'] });
let rLog = '';
r.stdout.on('data', d => { rLog += d; });
r.stderr.on('data', d => { rLog += d; });
process.on('exit', () => { try { r.kill('SIGKILL'); } catch {} });
const waitPort = async (ms = 90000) => {
  const t0 = Date.now();
  while (Date.now() - t0 < ms) {
    try { const res = await fetch(`http://127.0.0.1:${PORT}/`); if (res.status > 0) return true; } catch {}
    await new Promise(s => setTimeout(s, 500));
  }
  throw new Error(`app did not start:\n${rLog}`);
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

const lineState = () => page.evaluate(() => {
  const l = document.getElementById('app-assistant-progress');
  if (!l) return { exists: false };
  return { exists: true, shown: getComputedStyle(l).display !== 'none',
           txt: l.textContent };
});

// --- initial state: line exists, hidden -----------------------------------
const init = await lineState();
check('progress line exists in the panel', init.exists);
check('progress line hidden before any turn', init.exists && !init.shown);

// --- Scenario A: submit against the dead endpoint -------------------------
// The server pushes {working, submit} on receipt and {working, model} on
// request start; the failed request settles the stream -> idle -> hidden.
const composer = page.locator('#app-assistant-panel textarea, #app-assistant-panel [contenteditable="true"]');
// record message-driven renders WITHOUT breaking the real one: Shiny
// keeps ONE handler per type, so re-register with a wrapper that records
// and then calls the module's original render (captured from the hook)
await page.evaluate(() => {
  const orig = window.__omicsviewerAiProgress.update;
  window.__wp8Msgs = [];
  const rec = m => { window.__wp8Msgs.push(m); return orig(m); };
  window.__omicsviewerAiProgress.update = rec;
  Shiny.addCustomMessageHandler('omicsviewerAssistantProgress', rec);
});
await composer.first().click();
await composer.first().fill('progress pipeline probe');
await page.keyboard.press('Enter');
let sawWorking = false, sawModel = false;
const deadline = Date.now() + 30000;
while (Date.now() < deadline) {
  const msgs = await page.evaluate(() => window.__wp8Msgs);
  sawWorking = sawWorking || msgs.some(m => m.state === 'working' && m.phase === 'submit');
  sawModel = sawModel || msgs.some(m => m.state === 'working' && m.phase === 'model');
  if (sawWorking && sawModel) break;
  await page.waitForTimeout(250);
}
check('submit pushes the working/submit message', sawWorking);
check('request start pushes the working/model message', sawModel);
check('submit message precedes the model message',
  await page.evaluate(() => {
    const m = window.__wp8Msgs;
    const i1 = m.findIndex(x => x.state === 'working' && x.phase === 'submit');
    const i2 = m.findIndex(x => x.state === 'working' && x.phase === 'model');
    return i1 >= 0 && i2 >= 0 && i1 < i2;
  }));
// the request fails against the dead endpoint; the settle must hide it
let hidden = false;
const settleDeadline = Date.now() + 45000;
while (Date.now() < settleDeadline) {
  const s = await lineState();
  if (!s.shown) { hidden = true; break; }
  await page.waitForTimeout(500);
}
check('stream settle hides the line (idle)', hidden);
check('settle pushes the idle message',
  await page.evaluate(() =>
    window.__wp8Msgs.some(m => m.state === 'idle')));
check('message sequence ends working...working...idle',
  await page.evaluate(() => {
    const m = window.__wp8Msgs.map(x => x.state);
    return m.length >= 3 && m[0] === 'working' && m[1] === 'working' &&
      m[m.length - 1] === 'idle';
  }));

// --- Scenario B: step rendering through the hook --------------------------
const stepTxt = await page.evaluate(() => {
  window.__omicsviewerAiProgress.update({ state: 'working', phase: 'tool', step: 3, tool: 'set_enrichment_parameters' });
  const l = document.getElementById('app-assistant-progress');
  return l.textContent;
});
check('tool step renders as Working - step N: tool',
  /Working .*step 3: set_enrichment_parameters/.test(stepTxt), stepTxt);
await page.evaluate(() => window.__omicsviewerAiProgress.update({ state: 'idle' }));
check('idle clears through the hook', !(await lineState()).shown);

console.log(failures === 0 ? '\nWP-8 verification: ALL PASS' : `\nWP-8 verification: ${failures} FAILURES`);
await browser.close();
process.exit(failures === 0 ? 0 : 1);
