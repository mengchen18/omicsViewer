// Probe: does the chat transport die, or does the composer swallow
// messages? Sends short prompts at intervals, watches the websocket
// lifecycle + shiny connection state + server receipt (log file growth).
import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import { mkdirSync, existsSync, readFileSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../..');
const EXTDATA = path.join(REPO, 'inst/extdata');
const OUT = '/home/chen/.omicsviewer-llm-eval/probe1';
mkdirSync(OUT, { recursive: true });

const env = { ...process.env };
for (const line of readFileSync(path.join(here, 'provider.env'), 'utf8').split('\n')) {
  const m = line.match(/^\s*([A-Z_][A-Z0-9_]*)\s*=\s*(.*?)\s*$/);
  if (m && !line.trim().startsWith('#')) env[m[1]] = m[2];
}
const PORT = 7786;
const r = spawn('Rscript', ['-e', `
  options(shiny.port = ${PORT}, shiny.host = '127.0.0.1')
  eset <- readRDS(file.path('${EXTDATA}', 'demo.RDS'))
  omicsViewer::omicsViewer(dir = '${EXTDATA}', ESVObj = eset)
`], { cwd: REPO, env: { ...env, OMICSVIEWER_LLM_LOG: 'true', OMICSVIEWER_LLM_LOG_DIR: OUT },
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
  throw new Error('app did not start');
};

const logLines = () => {
  try {
    const f = require('node:fs').readdirSync(OUT).filter(x => x.endsWith('.jsonl'))[0];
    return readFileSync(path.join(OUT, f), 'utf8').split('\n').filter(Boolean).length;
  } catch { return 0; }
};
const userMsgCount = () => {
  try {
    const fs2 = require('node:fs');
    const f = fs2.readdirSync(OUT).filter(x => x.endsWith('.jsonl'))[0];
    return fs2.readFileSync(path.join(OUT, f), 'utf8').split('\n').filter(Boolean)
      .map(l => JSON.parse(l)).filter(l => l.event === 'user_message').length;
  } catch { return 0; }
};

await waitPort();
const browser = await chromium.launch({
  executablePath: '/usr/bin/google-chrome', headless: true,
  args: ['--no-sandbox', '--disable-gpu', '--disable-dev-shm-usage',
         '--disable-webgl', '--disable-webgl2']
});
const page = await (await browser.newContext()).newPage();
const wsEvents = [];
page.on('websocket', ws => {
  wsEvents.push({ t: Date.now(), ev: 'open', url: ws.url().slice(0, 60) });
  ws.on('close', () => wsEvents.push({ t: Date.now(), ev: 'CLOSE', url: ws.url().slice(0, 60) }));
  ws.on('socketerror', () => wsEvents.push({ t: Date.now(), ev: 'SOCKETERROR', url: ws.url().slice(0, 60) }));
});
await page.goto(`http://127.0.0.1:${PORT}/`, { timeout: 60000 });
await page.waitForFunction(() => !!document.querySelector('#app-contents')?.offsetParent, null, { timeout: 90000 });
await page.click('#app-assistant-toggle');
await page.waitForFunction(() => getComputedStyle(document.querySelector('#app-assistant-panel')).display !== 'none');
await page.waitForTimeout(2000);

const prompts = [
  'Reply with just the word apple.',
  'Reply with just the word banana.',
  'Now reply with just the word cherry.',
  'Reply with just the word date.',
  'Reply with just the word elderberry.',
  'Reply with just the word fig.'
];
const t0 = Date.now();
const ts = () => `t+${Math.round((Date.now() - t0) / 1000)}s`;
for (let i = 0; i < prompts.length; i++) {
  const composer = page.locator('#app-assistant-panel textarea, #app-assistant-panel [contenteditable="true"]');
  await composer.first().click();
  await composer.first().fill(prompts[i]);
  await page.keyboard.press('Enter');
  const before = userMsgCount();
  await page.waitForTimeout(40000);   // generous: tiny prompt, tiny reply
  const after = userMsgCount();
  const shiny = await page.evaluate(() => ({
    overlay: !!document.getElementById('shiny-disconnected-overlay'),
    busy: !!document.querySelector('.shiny-busy'),
    inputs: Object.keys(Shiny.shinyapp.$inputValues).length
  }));
  console.log(ts(), `prompt ${i + 1}: serverReceipts=${after}(${after - before ? '+' + (after - before) : 'NO-DELTA'})`,
    `overlay=${shiny.overlay} inputs=${shiny.inputs}`,
    `lastText=`, (await page.innerText('#app-assistant-body')).slice(-40).replace(/\n/g, '|'));
}
console.log('--- ws events ---');
for (const e of wsEvents) console.log(new Date(e.t).toISOString().slice(11, 19), e.ev, e.url);
await browser.close();
process.exit(0);
