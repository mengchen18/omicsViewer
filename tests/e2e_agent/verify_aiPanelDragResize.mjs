// Branch verification: AI assistant panel drag + native resize + persistence.
// Spawns the app with the demo dataset preloaded (mirrors tier_a.mjs).
import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../..');
const EXTDATA = path.join(REPO, 'inst/extdata');
const PORT = 7779;

const results = [];
const record = (name, pass, detail = '') => {
  results.push({ name, pass });
  console.log(`${pass ? 'ok' : 'FAIL'} - ${name}${detail ? ' [' + detail + ']' : ''}`);
};

const rCode = `
  pkgload::load_all('${REPO}', quiet = TRUE)
  options(shiny.port = ${PORT}, shiny.host = '127.0.0.1')
  eset <- readRDS(file.path('${EXTDATA}', 'demo.RDS'))
  omicsViewer::omicsViewer(dir = '${EXTDATA}', ESVObj = eset)
`;
const r = spawn('Rscript', ['-e', rCode], { cwd: REPO, stdio: ['ignore', 'pipe', 'pipe'] });
let rLog = '';
r.stdout.on('data', d => { rLog += d; });
r.stderr.on('data', d => { rLog += d; });
process.on('exit', () => { try { r.kill('SIGKILL'); } catch {} });

const waitPort = async (ms = 120000) => {
  const t0 = Date.now();
  while (Date.now() - t0 < ms) {
    try {
      const res = await fetch(`http://127.0.0.1:${PORT}/`);
      if (res.status > 0) return true;
    } catch {}
    await new Promise(s => setTimeout(s, 500));
  }
  throw new Error(`app did not start on :${PORT}\n${rLog}`);
};

const browser = await chromium.launch({
  executablePath: '/usr/bin/google-chrome', headless: true,
  args: ['--no-sandbox', '--disable-gpu', '--disable-dev-shm-usage', '--disable-webgl', '--disable-webgl2']
});
try {
  await waitPort();
  const page = await browser.newPage({ viewport: { width: 1400, height: 900 } });
  await page.goto(`http://127.0.0.1:${PORT}/`, { timeout: 60000 });
  // dataset is preloaded: wait for the data-space tab + contents panel, then
  // give the server a beat so the assistant-module observers are registered
  // (an action click delivered before registration is silently dropped).
  await page.waitForFunction(() => {
    const tab = window.Shiny && Shiny.shinyapp && Shiny.shinyapp.$inputValues['app-dataspace-eset'];
    const box = document.querySelector('#app-contents');
    return !!tab && !!box && box.offsetParent !== null;
  }, null, { timeout: 90000 });
  await page.waitForTimeout(3000);

  const panel = page.locator('.omicsviewer-ai-panel');
  const hiddenAtStart = await panel.evaluate(el => getComputedStyle(el).display === 'none');
  record('panel hidden at start', hiddenAtStart);

  await page.locator('.omicsviewer-ai-launcher').click();
  let chatOk = false;
  for (let i = 0; i < 30; i++) {
    await page.waitForTimeout(1000);
    const st = await page.evaluate(() => ({
      cls: document.querySelector('.omicsviewer-ai-panel')?.className,
      bodyLen: (document.querySelector('.omicsviewer-ai-body')?.innerHTML || '').length,
      chat: !!document.querySelector('.omicsviewer-ai-body shiny-chat-container')
    }));
    if (i % 5 === 0 || st.chat) console.log('poll', i, JSON.stringify(st));
    if (st.chat) { chatOk = true; break; }
  }
  record('panel opens and chat renders with preloaded dataset', chatOk);

  // header buttons must click, not drag: gear opens the settings modal
  await page.locator('.omicsviewer-ai-header button[title="Configure the AI model"]').click();
  await page.waitForSelector('.modal-dialog', { timeout: 10000 });
  record('header gear button still opens settings modal', true);
  await page.locator('#app-assistant-settings_cancel').click();
  await page.waitForSelector('.modal-dialog', { state: 'detached', timeout: 10000 });

  // drag via the header title
  const before = await panel.boundingBox();
  const title = await page.locator('.omicsviewer-ai-title').boundingBox();
  await page.mouse.move(title.x + title.width / 2, title.y + title.height / 2);
  await page.mouse.down();
  for (let i = 1; i <= 12; i++) {
    await page.mouse.move(
      title.x + (140 - title.x) * i / 12,
      title.y + (240 - title.y) * i / 12);
  }
  await page.mouse.up();
  await page.waitForTimeout(300);
  const afterDrag = await panel.boundingBox();
  record('panel moved by header drag',
    Math.abs(afterDrag.x - before.x) > 80 && Math.abs(afterDrag.y - before.y) > 80,
    `${Math.round(before.x)},${Math.round(before.y)} -> ${Math.round(afterDrag.x)},${Math.round(afterDrag.y)}`);
  record('anchored left/top after drag', await panel.evaluate(el =>
    el.style.right === 'auto' && el.style.bottom === 'auto'));

  // native CSS resize grip (bottom-right corner)
  const b4 = await panel.boundingBox();
  await page.mouse.move(b4.x + b4.width - 4, b4.y + b4.height - 4);
  await page.mouse.down();
  await page.mouse.move(b4.x + b4.width + 120, b4.y + b4.height + 100, { steps: 8 });
  await page.mouse.up();
  await page.waitForTimeout(400);
  const b5 = await panel.boundingBox();
  record('native resize grip grows panel',
    b5.width > b4.width + 60 && b5.height > b4.height + 60,
    `${Math.round(b4.width)}x${Math.round(b4.height)} -> ${Math.round(b5.width)}x${Math.round(b5.height)}`);

  // chat flex-fills the resized panel body
  const chatFills = await page.evaluate(() => {
    const c = document.querySelector('.omicsviewer-ai-body shiny-chat-container');
    const b = document.querySelector('.omicsviewer-ai-body');
    if (!c || !b) return false;
    const cr = c.getBoundingClientRect(), br = b.getBoundingClientRect();
    return cr.height > 300 && cr.height <= br.height + 1;
  });
  record('chat container flex-fills panel body', chatFills);

  // persistence
  const stored = await page.evaluate(() =>
    window.localStorage.getItem('omicsviewerAiPanelGeometry'));
  record('geometry persisted to localStorage', !!stored, stored || 'null');
  const saved = JSON.parse(stored || 'null');

  // reload: geometry restored while panel hidden, applied on reopen
  await page.reload({ waitUntil: 'networkidle' });
  await page.waitForTimeout(2000);
  const restoredInline = await page.evaluate(() => {
    const el = document.querySelector('.omicsviewer-ai-panel');
    return { left: el.style.left, top: el.style.top, width: el.style.width };
  });
  record('geometry applied on reload (inline styles present)',
    !!restoredInline.left && !!restoredInline.top,
    JSON.stringify(restoredInline));

  await page.locator('.omicsviewer-ai-launcher').click();
  await page.waitForSelector('.omicsviewer-ai-body shiny-chat-container', { timeout: 30000 });
  const reopenedBox = await panel.boundingBox();
  record('reopened panel at restored spot',
    saved && Math.abs(reopenedBox.x - saved.left) < 40 && Math.abs(reopenedBox.width - saved.width) < 40,
    `${Math.round(reopenedBox.x)},${Math.round(reopenedBox.width)} vs saved ${saved && saved.left},${saved && saved.width}`);

  await page.keyboard.press('Escape');
  let closed = false;
  for (let i = 0; i < 10; i++) {
    await page.waitForTimeout(500);
    if (await panel.evaluate(el => getComputedStyle(el).display === 'none')) { closed = true; break; }
  }
  record('Escape still closes panel', closed);
} finally {
  await browser.close();
  try { r.kill('SIGKILL'); } catch {}
}

const failed = results.filter(x => !x.pass);
console.log(`\n${results.length - failed.length}/${results.length} passed`);
process.exit(failed.length ? 1 : 0);
