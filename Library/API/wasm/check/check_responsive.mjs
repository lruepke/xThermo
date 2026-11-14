// Responsive regression check of the Explorer in a real browser (manual; needs Chrome, not run in CI).
//
//   "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" --headless=new --remote-debugging-port=9333 \
//       --user-data-dir=/tmp/xt-chrome about:blank &
//   node check_responsive.mjs <site base URL ending in />     e.g. http://127.0.0.1:8000/ (mkdocs build + make_wasm.sh dist
//                                                             in site/app) or https://lruepke.github.io/xThermo/
//
// For phone widths 320-430 px, light and dark scheme, the standalone page (app/) and the homepage embed, and the
// states liquid (default), V+L, V+H, L+H and V+L+H: evaluates the state and fails if the page or the embedded frame
// is wider than the viewport, if an element of the result passes the right edge outside a scroll wrapper, or if a
// cell stays out of reach after scrolling its table to the end. Tables wider than their card scroll inside their own
// wrapper (.xt-scroll); those cases are listed, not failed. Exit status 1 on any failure.
const BASE = process.argv[2];
if (!BASE) { console.error('usage: node check_responsive.mjs <site base URL ending in />'); process.exit(2); }
const t = (await (await fetch('http://127.0.0.1:9333/json/list')).json()).find(t => t.type === 'page');
const ws = new WebSocket(t.webSocketDebuggerUrl); await new Promise(r => ws.addEventListener('open', r));
let id = 0; const pend = new Map();
ws.addEventListener('message', e => { const m = JSON.parse(e.data); if (pend.has(m.id)) { pend.get(m.id)(m); pend.delete(m.id); } });
const send = (method, params = {}) => new Promise(r => { pend.set(++id, r); ws.send(JSON.stringify({ id, method, params })); });
const sleep = ms => new Promise(r => setTimeout(r, ms));
const ev = async x => { const r = await send('Runtime.evaluate', { expression: x, returnByValue: true, awaitPromise: true }); return r.result.result.value; };
const STATES = [
  ['liquid (default)', null, 'Liquid'],
  ['V+L', { T: '400', p: '200', X: '3.2' }, 'V+L'],
  ['V+H', { T: '600', p: '100', X: '3.2' }, 'V+H'],
  ['L+H', { T: '300', p: '100', X: '40' }, 'L+H'],
  ['V+L+H', { mode: 'hpx', H: '2200', p: '100', X: '20' }, 'V+L+H'],
];
await send('Network.enable'); await send('Network.setCacheDisabled', { cacheDisabled: true });   // never a stale stylesheet
const DOC = embed => embed ? "document.getElementById('xt-explorer-frame').contentDocument" : 'document';
let fails = 0, scrolled = 0, n = 0;
for (const w of [320, 360, 375, 390, 430]) {
  await send('Emulation.setDeviceMetricsOverride', { width: w, height: 1400, deviceScaleFactor: 2, mobile: true });
  for (const scheme of ['light', 'dark']) {
    await send('Emulation.setEmulatedMedia', { features: [{ name: 'prefers-color-scheme', value: scheme }] });
    for (const embed of [false, true]) {
      await send('Page.navigate', { url: BASE + (embed ? '' : 'app/index.html') });
      await sleep(2500);
      for (const [name, s, expect] of STATES) {
        const m = await ev(`(async () => { const d = ${DOC(embed)}; const win = d.defaultView;
          const set = (i, v) => { d.getElementById(i).value = v; };
          const s = ${JSON.stringify(s)};
          const mode = d.querySelector('input[name=mode][value=' + ((s && s.mode) || 'tpx') + ']');
          mode.checked = true; mode.dispatchEvent(new win.Event('change', { bubbles: true }));
          if (s) { if (s.T) set('xt-in-T', s.T); if (s.H) set('xt-in-H', s.H); set('xt-in-p', s.p); set('xt-in-X', s.X); }
          d.getElementById('xt-evaluate').click();
          for (let i = 0; i < 100; i++) { await new Promise(r => setTimeout(r, 100)); if (!d.getElementById('xt-result').hidden) break; }
          await new Promise(r => setTimeout(r, 200));
          const W = document.documentElement.clientWidth, fw = d.documentElement.clientWidth;
          const past = [...d.querySelectorAll('#xt-result *')].filter(e => e.getBoundingClientRect().right > fw + 0.5 && !e.closest('.xt-scroll'))
            .map(e => e.id || e.tagName.toLowerCase());
          const wraps = [...d.querySelectorAll('#xt-result .xt-scroll')];
          const inner = wraps.filter(x => x.scrollWidth > x.clientWidth + 0.5).length;
          const unreachable = wraps.map(x => { x.scrollLeft = x.scrollWidth; const r = x.getBoundingClientRect();
            const k = [...x.querySelectorAll('td, th')].filter(c => c.offsetParent && c.getBoundingClientRect().right > r.right + 0.5).length;
            x.scrollLeft = 0; return k; }).reduce((a, b) => a + b, 0);
          return { phase: d.getElementById('xt-phase').textContent, page: document.documentElement.scrollWidth - W,
                   frame: d.documentElement.scrollWidth - fw, past, inner, unreachable }; })()`);
        n++;
        const bad = m.phase !== expect || m.page > 0 || m.frame > 0 || m.past.length || m.unreachable;
        if (bad) fails++;
        if (m.inner) scrolled++;
        console.log(`${bad ? 'FAIL' : 'PASS'} ${w}px ${scheme} ${embed ? 'embed' : 'app'} ${name}: phase ${m.phase}, page +${m.page}, frame +${m.frame}` +
                    `${m.past.length ? ', past edge ' + m.past.slice(0, 4) : ''}${m.unreachable ? ', unreachable cells ' + m.unreachable : ''}` +
                    `${m.inner ? ', table scrolls inside its wrapper' : ''}`);
      }
    }
  }
}
ws.close();
console.log(`${n - fails}/${n} pass; ${scrolled} with a table scrolling inside its wrapper`);
process.exit(fails ? 1 : 0);
