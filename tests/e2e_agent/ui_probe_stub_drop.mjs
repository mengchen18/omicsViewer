// WP-4 (F2/F8) probe: force context stubbing with heavy turns, then
// deliver 5 short probe messages and record, per message, whether the
// server received it (chat_submit_seen / user_message in the jsonl),
// what the browser actually sent on the websocket, the cancel-control
// state, and the jsonl tail. A DROP = user bubble renders client-side
// but no server receipt within 30 s while the websocket stayed open.
//
// Run: node ui_probe_stub_drop.mjs [outName]
// Env: provider.env (real provider - the bug may need real timings).
// Exit: 0 always; the last line is the verdict (drops N of M).

import { chromium } from 'playwright';
import { spawn } from 'node:child_process';
import { mkdirSync, readFileSync, readdirSync, statSync, writeFileSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(here, '../..');
const EXTDATA = path.join(REPO, 'inst/extdata');
const SESSION = process.argv[2] || `stub_probe_${new Date().toISOString().replace(/[:.]/g, '')}`;
const OUT = `/home/chen/.omicsviewer-llm-eval/${SESSION}`;
mkdirSync(OUT, { recursive: true });

const env = { ...process.env };
for (const line of readFileSync(path.join(here, 'provider.env'), 'utf8').split('\n')) {
  const m = line.match(/^\s*([A-Z_][A-Z0-9_]*)\s*=\s*(.*?)\s*$/);
  if (m && !line.trim().startsWith('#')) env[m[1]] = m[2];
}
if (!env.OMICSVIEWER_LLM_API_KEY) { console.error('FATAL: no API key'); process.exit(2); }

const PORT = 7787;
const r = spawn('Rscript', ['-e', `
  options(shiny.port = ${PORT}, shiny.host = '127.0.0.1')
  eset <- readRDS(file.path('${EXTDATA}', 'demo.RDS'))
  omicsViewer::omicsViewer(dir = '${EXTDATA}', ESVObj = eset)
`], { cwd: REPO, env: { ...env, OMICSVIEWER_LLM_LOG: 'true', OMICSVIEWER_LLM_LOG_DIR: OUT },
     stdio: ['ignore', 'pipe', 'pipe'] });
let rLog = '';
r.stdout.on('data', d => { rLog += d; });
r.stderr.on('data', d => { rLog += d; });
const cleanup = () => { try { r.kill('SIGKILL'); } catch {} };
process.on('exit', cleanup);
process.on('SIGINT', () => { cleanup(); process.exit(130); });

const waitPort = async (ms = 90000) => {
  const t0 = Date.now();
  while (Date.now() - t0 < ms) {
    try { const res = await fetch(`http://127.0.0.1:${PORT}/`); if (res.status > 0) return; } catch {}
    await new Promise(s => setTimeout(s, 500));
  }
  throw new Error(`app did not start:\n${rLog}`);
};

const logFile = () => {
  // multiple shiny sessions may write logs here (ws reconnects spawn a
  // fresh session + logger); follow the LIVE one by mtime, not [0]
  const fs = readdirSync(OUT).filter(f => f.endsWith('.jsonl'))
    .map(f => ({ f, m: statSync(path.join(OUT, f)).mtimeMs }))
    .sort((a, b) => b.m - a.m);
  return fs.length ? path.join(OUT, fs[0].f) : null;
};
const events = () => {
  const f = logFile();
  if (!f) return [];
  try {
    return readFileSync(f, 'utf8').split('\n').filter(Boolean).map(l => JSON.parse(l));
  } catch { return []; }
};
const countEvent = (name) => events().filter(e => e.event === name).length;
const requestInFlight = () => {
  const ev = events();
  return ev.length > 0 && ev[ev.length - 1].event === 'provider_request_start';
};
const chatBusy = async (pg) => pg.evaluate(() => {
  const b = document.querySelector('[class*="cancel"]');
  return !!(b && b.offsetParent !== null && !b.disabled);
});
const waitIdle = async (pg, ms = 420000) => {
  const dl = Date.now() + ms;
  while (Date.now() < dl) {
    if (!(await chatBusy(pg)) && !requestInFlight()) return true;
    await new Promise(s => setTimeout(s, 2000));
  }
  return false;
};

const t0 = Date.now();
const ts = () => `t+${Math.round((Date.now() - t0) / 1000)}s`;

await waitPort();
const browser = await chromium.launch({
  executablePath: '/usr/bin/google-chrome', headless: true,
  args: ['--no-sandbox', '--disable-gpu', '--disable-dev-shm-usage',
         '--disable-webgl', '--disable-webgl2']
});
const page = await (await browser.newContext({ viewport: { width: 1600, height: 1000 } })).newPage();

// ws frame recorder: every client->server frame, timestamped
const wsFrames = [];
let wsOpen = 0, wsClosed = 0;
page.on('websocket', ws => {
  wsOpen++;
  ws.on('close', () => wsClosed++);
  ws.on('socketerror', () => wsFrames.push({ t: Date.now(), dir: 'SOCKETERROR' }));
  ws.on('framesent', f => {
    const p = typeof f.payload === 'string' ? f.payload : '(binary)';
    wsFrames.push({ t: Date.now(), dir: 'c>s', payload: p.slice(0, 300) });
  });
});

await page.goto(`http://127.0.0.1:${PORT}/`, { timeout: 60000 });
await page.waitForFunction(() => !!document.querySelector('#app-contents')?.offsetParent,
  null, { timeout: 90000 });
await page.click('#app-assistant-toggle');
await page.waitForFunction(() => getComputedStyle(document.querySelector('#app-assistant-panel')).display !== 'none',
  null, { timeout: 15000 });
await page.waitForFunction(() => {
  const b = document.querySelector('#app-assistant-body');
  return b && b.innerText.length > 0;
}, null, { timeout: 15000 });
console.log(ts(), 'app + assistant ready');

const submit = async (prompt) => {
  const composer = page.locator(
    '#app-assistant-panel textarea, #app-assistant-panel [contenteditable="true"]');
  await composer.first().click();
  await composer.first().fill(prompt);
  await page.keyboard.press('Enter');
};

const deliverTurn = async (prompt, label, budget_ms = 300000) => {
  const before = events().length;
  const textBefore = (await page.innerText('#app-assistant-body')).length;
  await submit(prompt);
  console.log(ts(), `>> [${label}] "${prompt.slice(0, 60)}"`);
  const tTurn = Date.now();
  let prev = '', stable = 0, timedOut = true;
  while (Date.now() - tTurn < budget_ms) {
    await new Promise(s => setTimeout(s, 3000));
    const st = await page.evaluate(() => ({
      text: (document.querySelector('#app-assistant-body')?.innerText || ''),
      streaming: (() => {
        const b = document.querySelector('[class*="cancel"]');
        return !!(b && b.offsetParent !== null && !b.disabled);
      })()
    }));
    if (st.streaming || requestInFlight()) { stable = 0; prev = st.text; continue; }
    if (st.text === prev && st.text.length > textBefore + 20) {
      if (++stable >= 6) { timedOut = false; break; }
    } else stable = 0;
    prev = st.text;
  }
  return { timedOut, newEvents: events().length - before };
};

// ---- Phase A: heavy turns to force stub(s) --------------------------------
const RACE = (process.env.PROBE_MODE || 'idle') === 'race';
const MIDSTREAM = (process.env.PROBE_MODE || 'idle') === 'midstream';
const heavy = [
  ['h1-volcano', 'Create a volcano plot for the RE vs ME comparison.'],
  ['h2-pca-fig', 'Create a second figure: a static PCA scatter of the samples colored by Origin, with axis titles.'],
  ['h3-ora', 'Run over-representation analysis on the RE vs ME result (FDR < 0.05, |logFC| > 1) and summarize the top enriched gene sets.'],
  ['h4-update', 'Update the volcano figure: color points by whether their log FDR is below 0.05 and label the top 10 genes.'],
];
if (!MIDSTREAM || /^(1|true|yes|on)$/i.test(process.env.PROBE_PHASE_A || '')) {
for (const [id, prompt] of heavy) {
  const rec = await deliverTurn(prompt, id, 300000);
  console.log(ts(), `<< [${id}] ${rec.timedOut ? 'TIMEOUT' : 'settled'}; stubs so far: ${countEvent('history_stubbed')}`);
}
}
console.log(ts(), `phase A done: history_stubbed=${countEvent('history_stubbed')} ` +
  `compactions=${countEvent('history_compacted')} requests=${countEvent('provider_request_start')}`);

// ---- Phase B: probe messages ----------------------------------------------
const words = ['alpha', 'bravo', 'charlie', 'delta', 'echo'];
const probes = [];
if (MIDSTREAM) {
  // Reproduces the batch2 death shape: prompts submitted WHILE a long
  // tool chain streams (queued by shinychat), the chain then aborts at
  // the session request limit, and we watch whether the queued prompts
  // (and fresh ones) ever reach the server.
  console.log(ts(), 'MIDSTREAM MODE: queuing probes during a long chain');
  const queued = ['alpha', 'bravo'];
  const after = ['charlie', 'delta', 'echo'];
  const queuedRecords = [];
  const seenBefore = countEvent('chat_submit_seen');
  const msgBefore = countEvent('user_message');
  // long-chaining turn (ORA flails = many requests); with a low
  // OMICSVIEWER_LLM_MAX_REQUESTS the chain aborts mid-turn
  const chain = page.locator(
    '#app-assistant-panel textarea, #app-assistant-panel [contenteditable="true"]');
  await chain.first().click();
  await chain.first().fill('Run over-representation analysis on the RE vs ME result (FDR < 0.05, |logFC| > 1), then describe the top enriched gene sets in detail.');
  await page.keyboard.press('Enter');
  console.log(ts(), '>> long chain turn submitted');
  const chainT = Date.now();
  for (let i = 0; i < queued.length; i++) {
    await new Promise(s => setTimeout(s, 15000 + i * 30000));
    const busy = await chatBusy(page);
    const textBeforeQ = (await page.innerText('#app-assistant-body')).length;
    const qSubmitT = Date.now();
    const composer = page.locator(
      '#app-assistant-panel textarea, #app-assistant-panel [contenteditable="true"]');
    await composer.first().click();
    await composer.first().fill(`Reply with just the word ${queued[i]}.`);
    await page.keyboard.press('Enter');
    // bubble render check at +3 s and +12 s (optimistic render?)
    await new Promise(s => setTimeout(s, 3000));
    const text3 = (await page.innerText('#app-assistant-body')).length;
    await new Promise(s => setTimeout(s, 9000));
    const text12 = (await page.innerText('#app-assistant-body')).length;
    const qFrames = wsFrames.filter(f => f.t >= qSubmitT - 1000 && f.t <= qSubmitT + 12000 &&
      f.dir === 'c>s' && f.payload.includes(queued[i]));
    console.log(ts(), `   queued [${queued[i]}] during chain (busy=${busy}) ` +
      `bubbleDelta3s=${text3 - textBeforeQ} bubbleDelta12s=${text12 - textBeforeQ} ` +
      `inputFrames=${qFrames.length}`);
    queuedRecords.push({ word: queued[i], busy, qSubmitT,
      bubbleDelta3s: text3 - textBeforeQ, bubbleDelta12s: text12 - textBeforeQ,
      inputFrames: qFrames.length });
  }
  // let the chain finish or abort (~2.5 min covers a full ORA-style chain)
  await new Promise(s => setTimeout(s, 150000));
  const limitHits0 = countEvent('request_limit_reached');
  // WP-6: a queued submission that was discarded must have surfaced the
  // delivery notice; clicking Retry must deliver it through the server
  let noticeSeen = false, recovered = null;
  const noticeEl = page.locator('.omicsviewer-ai-delivery-notice');
  noticeSeen = await noticeEl.count() > 0;
  if (noticeSeen) {
    const noticeText = (await page.innerText('.omicsviewer-ai-delivery-notice').catch(() => ''));
    console.log(ts(), `   delivery notice surfaced: "${noticeText.replace(/\n/g, ' ').slice(0, 80)}"`);
    const word = queued.find(w => noticeText.includes(w) || true) || null;
    // the notice does not carry the text; retry whatever is pending
    const msgBefore2 = countEvent('user_message');
    await noticeEl.locator('button').click();
    let ok = false;
    const dl2 = Date.now() + 45000;
    while (Date.now() < dl2) {
      await new Promise(s => setTimeout(s, 1500));
      if (countEvent('user_message') > msgBefore2) { ok = true; break; }
    }
    recovered = ok;
    console.log(ts(), `   retry clicked -> ${ok ? 'DELIVERED (user_message logged)' : 'not delivered'}`);
  } else {
    console.log(ts(), '   no delivery notice visible after chain');
  }
  const limitHits = limitHits0;
  const qSeen = countEvent('chat_submit_seen') - seenBefore;
  const qMsg = countEvent('user_message') - msgBefore;
  console.log(ts(), `chain over: limitHits=${limitHits} receipts so far: submit_seen=${qSeen} user_msg=${qMsg}`);
  const allWords = [...queued, ...after];
  for (let i = 0; i < after.length; i++) {
    await waitIdle(page, 420000);
    const sB = countEvent('chat_submit_seen'); const mB = countEvent('user_message');
    const submitT = Date.now();
    const busyAtSubmit = await chatBusy(page);
    await submit(`Reply with just the word ${after[i]}.`);
    let receipt = false;
    const dl = Date.now() + 35000;
    while (Date.now() < dl) {
      await new Promise(s => setTimeout(s, 1500));
      if (countEvent('chat_submit_seen') > sB || countEvent('user_message') > mB) { receipt = true; break; }
    }
    const framesNear = wsFrames.filter(f => f.t >= submitT - 2000 && f.t <= submitT + 10000 &&
      f.dir === 'c>s' && /chat|Reply|word/i.test(f.payload));
    probes.push({ word: after[i], receipt, busyAtSubmit, frames: framesNear.length,
      wsOpen: wsOpen - wsClosed, queuedDuringChain: false });
    console.log(ts(), `post-probe ${i + 1}/${after.length} [${receipt ? 'OK ' : 'DROP'}] word=${after[i]} ` +
      `busyAtSubmit=${busyAtSubmit} frames=${framesNear.length} wsOpen=${wsOpen - wsClosed}`);
    if (!receipt) {
      const ev = events();
      const around = ev.filter(e => Date.parse(e.timestamp) >= chainT - 5000 && Date.parse(e.timestamp) <= Date.now());
      console.log('    --- timeline since chain start (last 30) ---');
      for (const e of around.slice(-30))
        console.log(`    ${e.timestamp.slice(11, 19)} seq${e.sequence} ${e.event} ${JSON.stringify(e.details || {}).slice(0, 100)}`);
    }
    await new Promise(s => setTimeout(s, 15000));
  }
  // queued-word receipts: search the log for each queued word
  const evAll = events();
  for (const w of queued) {
    const found = evAll.some(e => (e.event === 'user_message' || e.event === 'chat_submit_seen') &&
      JSON.stringify(e.details || {}).includes(w));
    probes.push({ word: w, receipt: found, queuedDuringChain: true });
    console.log(ts(), `queued [${w}] receipt=${found}`);
  }
  writeFileSync(path.join(OUT, 'probe_report.json'), JSON.stringify({
    session: SESSION, mode: 'midstream', limitHits, noticeSeen, recovered,
    probes }, null, 2));
  const drops = probes.filter(p => !p.receipt).length;
  console.log(`\nVERDICT (midstream): drops ${drops} of ${probes.length}; limitHits ${limitHits}; ` +
    `notice=${noticeSeen} retryRecovered=${recovered}; ` +
    `report: ${path.join(OUT, 'probe_report.json')}`);
  await browser.close();
  process.exit(0);
}
if (RACE) {
  // H-b discriminator: submit within ~1 s of the idle flip (and of any
  // history_stubbed event), racing the janitor flush instead of waiting
  // for verified idle. Heavy turn -> settle-watch (1 s stability, NOT
  // log-verified idle) -> immediate probe; repeat.
  console.log(ts(), 'RACE MODE: probes submit immediately after apparent settle');
  for (let i = 0; i < words.length; i++) {
    const heavyPrompt = i % 2 === 0
      ? 'Update the volcano figure: label the top ' + (5 + i) + ' genes and adjust the title.'
      : 'Create another static figure: a sample PCA scatter colored by Origin.';
    await deliverTurn(heavyPrompt, `race-h${i + 1}`, 300000);
    // settle-watch: 1 s of stable text, no log verification
    let prev = '', stable = 0;
    const dl = Date.now() + 240000;
    while (Date.now() < dl && stable < 2) {
      await new Promise(s => setTimeout(s, 500));
      const txt = await page.innerText('#app-assistant-body').catch(() => '');
      if (txt === prev && txt.length > 20) stable++; else stable = 0;
      prev = txt;
    }
    const seenBefore = countEvent('chat_submit_seen');
    const msgBefore = countEvent('user_message');
    const stubsBefore = countEvent('history_stubbed');
    const submitT = Date.now();
    const busyAtSubmit = await chatBusy(page);
    await submit(`Reply with just the word ${words[i]}.`);
    await new Promise(s => setTimeout(s, 35000));
    const receipt = countEvent('chat_submit_seen') > seenBefore ||
      countEvent('user_message') > msgBefore;
    const stubbedAfter = countEvent('history_stubbed') - stubsBefore;
    const ev = events();
    const tail3 = ev.slice(-3).map(e => `${e.timestamp.slice(11, 19)}/${e.event}`).join(' ');
    const framesNear = wsFrames.filter(f => f.t >= submitT - 2000 && f.t <= submitT + 10000 &&
      f.dir === 'c>s' && /chat|Reply|word/i.test(f.payload));
    probes.push({ word: words[i], receipt, busyAtSubmit, race: true, stubbedAfter,
      tail: tail3, frames: framesNear.length, wsOpen: wsOpen - wsClosed });
    console.log(ts(), `race-probe ${i + 1}/${words.length} [${receipt ? 'OK ' : 'DROP'}] word=${words[i]} ` +
      `busyAtSubmit=${busyAtSubmit} stubsAfter=${stubbedAfter} frames=${framesNear.length} wsOpen=${wsOpen - wsClosed}`);
    if (!receipt) {
      const around = ev.filter(e => {
        const et = Date.parse(e.timestamp);
        return et >= submitT - 30000 && et <= Date.now();
      });
      console.log('    --- timeline around drop (last 25) ---');
      for (const e of around.slice(-25))
        console.log(`    ${e.timestamp.slice(11, 19)} seq${e.sequence} ${e.event} ${JSON.stringify(e.details || {}).slice(0, 110)}`);
    }
  }
} else {
for (let i = 0; i < words.length; i++) {
  await waitIdle(page, 420000);
  const seenBefore = countEvent('chat_submit_seen');
  const msgBefore = countEvent('user_message');
  const submitT = Date.now();
  const busyAtSubmit = await chatBusy(page);
  await submit(`Reply with just the word ${words[i]}.`);
  // wait up to 35 s for the server receipt marker
  let receipt = false;
  const dl = Date.now() + 35000;
  while (Date.now() < dl) {
    await new Promise(s => setTimeout(s, 1500));
    if (countEvent('chat_submit_seen') > seenBefore || countEvent('user_message') > msgBefore) {
      receipt = true; break;
    }
  }
  const ev = events();
  const tail3 = ev.slice(-3).map(e => `${e.timestamp.slice(11, 19)}/${e.event}`).join(' ');
  const framesNear = wsFrames.filter(f => f.t >= submitT - 2000 && f.t <= submitT + 10000 &&
    f.dir === 'c>s' && /chat|Reply|word/i.test(f.payload));
  const rec = {
    word: words[i], receipt, busyAtSubmit,
    tail: tail3,
    frames: framesNear.length,
    frameSamples: framesNear.slice(0, 2).map(f => f.payload.slice(0, 120)),
    wsOpen: wsOpen - wsClosed,
  };
  probes.push(rec);
  const mark = receipt ? 'OK ' : 'DROP';
  console.log(ts(), `probe ${i + 1}/${words.length} [${mark}] word=${words[i]} ` +
    `busyAtSubmit=${rec.busyAtSubmit} wsOpen=${rec.wsOpen} frames=${rec.frames}`);
  console.log(`    tail: ${rec.tail}`);
  if (rec.frames) console.log(`    frame: ${rec.frameSamples.join(' | ').slice(0, 220)}`);
  if (!receipt) {
    // dump the event timeline around the drop + janitor correlation
    const around = ev.filter(e => {
      const et = Date.parse(e.timestamp);
      return et >= submitT - 30000 && et <= Date.now();
    });
    console.log('    --- timeline around drop (last 25) ---');
    for (const e of around.slice(-25))
      console.log(`    ${e.timestamp.slice(11, 19)} seq${e.sequence} ${e.event} ${JSON.stringify(e.details || {}).slice(0, 110)}`);
  }
  // spacing: next probe at >= 30 s after this submit
  await new Promise(s => setTimeout(s, Math.max(0, 30000 - (Date.now() - submitT))));
}
}

const drops = probes.filter(p => !p.receipt).length;
const report = {
  session: SESSION, started: new Date(t0).toISOString(),
  stubs: countEvent('history_stubbed'), compactions: countEvent('history_compacted'),
  requests: countEvent('provider_request_start'), wsOpenTotal: wsOpen, wsClosed,
  probes,
};
writeFileSync(path.join(OUT, 'probe_report.json'), JSON.stringify(report, null, 2));
console.log(`\nVERDICT: drops ${drops} of ${probes.length}; stubs ${report.stubs}; ` +
  `report: ${path.join(OUT, 'probe_report.json')}`);
await browser.close();
process.exit(0);
