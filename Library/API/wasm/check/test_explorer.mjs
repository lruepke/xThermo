// W2b: interface tests of the Explorer page (docs/app; review/W2B_EXPLORER_PAGE.md). Node only, no packages.
//
//   node test_explorer.mjs <dist dir of make_wasm.sh: xthermo.js, xthermo_wasm.js/.wasm>
//
// 1. page structure: every element the script uses exists; no external scripts, styles or fonts; the About links
//    name exactly the licence files make_wasm.sh writes; the links into the documentation resolve to pages
//    (and the anchor) of docs/
// 2. unit conversion and formatting (explorer.js): input units to SI as documented, empty/non-numeric -> NaN
// 3. worker evaluation (worker.js evaluate) against direct calls of xthermo.js: every value of successful states
//    bitwise identical (both backends, T-P-X and H-P-X, every phase region type), typed failures with the same
//    class, status and library message (InvalidInput, OutOfRange, NotConverged)
// 4. the worker's message protocol in a simulated worker scope: init with a compiled module -> ready with the
//    build information, eval -> result; a missing module -> unavailable (the documentation stays usable)
// Prints PASS/FAIL per test; exit status 1 if any fails.
import fs from 'fs';
import path from 'path';
import vm from 'vm';
import { createRequire } from 'module';
import { fileURLToPath } from 'url';

const require = createRequire(import.meta.url);
const HERE = path.dirname(fileURLToPath(import.meta.url));
const ROOT = path.resolve(HERE, '../../../..');
const APP = path.join(ROOT, 'docs/app');
const dist = path.resolve(process.argv[2] || '');
let fails = 0;
const results = [];
async function test(name, fn) {
  try { await fn(); results.push(`PASS ${name}`); }
  catch (e) { fails++; results.push(`FAIL ${name}: ${e.message}`); }
  console.log(results[results.length - 1]);
}
function assert(c, msg) { if (!c) throw new Error(msg); }
const same = (a, b) => Object.is(a, b) || (Number.isNaN(a) && Number.isNaN(b));

const html = fs.readFileSync(path.join(APP, 'index.html'), 'utf8');
const js = fs.readFileSync(path.join(APP, 'explorer.js'), 'utf8');
const page = require(path.join(APP, 'explorer.js'));
const { evaluate } = require(path.join(APP, 'worker.js'));

// --- 1. page structure
await test('page: every element id used by explorer.js exists in index.html', () => {
  const used = [...js.matchAll(/\$\('([\w-]+)'\)/g)].map(m => m[1]);
  const ids = new Set([...html.matchAll(/\sid="([\w-]+)"/g)].map(m => m[1]));
  const missing = [...new Set(used)].filter(id => !ids.has(id));
  assert(used.length > 10 && !missing.length, `missing ids: ${missing}`);
});
await test('page: no external scripts, stylesheets, fonts or images; links only to the project on GitHub', () => {
  const refs = [...html.matchAll(/\s(src|href)="([^"]+)"/g)].map(m => m[2]);
  const ext = refs.filter(r => /^[a-z]+:/i.test(r) && !r.startsWith('https://github.com/lruepke/xThermo'));
  assert(!ext.length, `external references: ${ext}`);
  for (const f of ['explorer.css', 'explorer.js', 'worker.js']) assert(fs.existsSync(path.join(APP, f)), `${f} missing`);
  assert(!/@import|url\(\s*['"]?https?:/.test(fs.readFileSync(path.join(APP, 'explorer.css'), 'utf8')), 'external CSS resource');
});
await test('layout: result grid tracks can shrink below the tables (phones: no page overflow in multiphase results)', () => {
  // A grid track without minmax(0, ...) grows to the tables' unwrapped width (min-width: auto) and pushes V+L, V+H and
  // L+H results past the right edge on 320-390 px screens; the .xt-scroll wrappers then never scroll. Every
  // grid-template-columns of .xt-res-grid, at every width, must start each track at 0. Browser check:
  // Library/API/wasm/check/check_responsive.mjs.
  const css = fs.readFileSync(path.join(APP, 'explorer.css'), 'utf8').replace(/\/\*[\s\S]*?\*\//g, '');
  const rules = [...css.matchAll(/\.xt-res-grid\s*\{([^}]*)\}/g)].map(m => m[1]);
  const base = rules.find(r => /display\s*:\s*grid/.test(r));
  assert(base, 'no base .xt-res-grid grid rule');
  const cols = rules.map(r => (r.match(/grid-template-columns\s*:\s*([^;]+)/) || [])[1]).filter(Boolean);
  assert(/grid-template-columns/.test(base), 'base .xt-res-grid rule sets no grid-template-columns (track min-width: auto)');
  for (const c of cols) {
    const tracks = c.trim().split(/\s+(?![^(]*\))/);
    assert(tracks.every(t => /^minmax\(\s*0\s*,/.test(t)), `.xt-res-grid track not minmax(0, ...): ${c}`);
  }
  assert(/\.xt-scroll\s*\{[^}]*overflow-x\s*:\s*auto/.test(css), '.xt-scroll lost overflow-x: auto');
});
await test('page: About links name exactly the licence files and notice of make_wasm.sh', () => {
  const sh = fs.readFileSync(path.join(ROOT, 'Library/API/wasm/make_wasm.sh'), 'utf8');
  const written = new Set([...sh.matchAll(/\$OUT\/dist\/(licenses\/[\w.-]+\.txt|NOTICE\.txt)/g)].map(m => m[1]));
  const linked = new Set([...html.matchAll(/href="((?:licenses\/)?[\w.-]+\.txt)"/g)].map(m => m[1]));
  const notWritten = [...linked].filter(l => !written.has(l)), notLinked = [...written].filter(w => !linked.has(w));
  assert(written.size >= 8 && !notWritten.length && !notLinked.length, `linked not written: ${notWritten}; written not linked: ${notLinked}`);
  if (dist) for (const w of written) assert(fs.existsSync(path.join(dist, w)), `${w} not in ${dist}`);
});
await test('page: links into the documentation resolve (pages, anchor, favicon)', () => {
  const docs = path.join(ROOT, 'docs');
  const rel = [...html.matchAll(/href="(\.\.\/[^"]*)"/g)].map(m => m[1]);
  assert(rel.length >= 4, 'expected documentation links');
  for (const r of rel) {
    const [p, anchor] = r.slice(3).split('#');
    const src = p === '' ? 'index.md' : p.endsWith('/') ? p.slice(0, -1) + '.md' : p;
    assert(fs.existsSync(path.join(docs, src)), `${r}: docs/${src} missing`);
    if (anchor) {
      const heads = fs.readFileSync(path.join(docs, src), 'utf8').match(/^#+ .*$/gm) || [];
      const slugs = heads.map(h => h.replace(/^#+ /, '').toLowerCase().replace(/[^\w\s-]/g, '').trim().replace(/\s+/g, '-'));
      assert(slugs.includes(anchor), `${r}: no heading for #${anchor}`);
    }
  }
});
await test('page: module outputs are not committed in docs/app (ignored by git; absent from the public snapshot)', () => {
  const outputs = ['xthermo_wasm.js', 'xthermo_wasm.wasm', 'xthermo.js', 'NOTICE.txt', 'licenses/', 'source/'];
  const gi = path.join(APP, '.gitignore');
  if (fs.existsSync(gi)) {                       // development repository
    const lines = fs.readFileSync(gi, 'utf8').split('\n');
    for (const f of outputs) assert(lines.includes(f), `${f} not in docs/app/.gitignore`);
  } else {                                       // public snapshot: the publish job drops every .gitignore (DROP_GLOBS)
    for (const f of outputs) assert(!fs.existsSync(path.join(APP, f)), `${f} present in docs/app before deployment`);
  }
});

await test('site: Explorer either hidden (excluded, not in nav) or released (in nav, not excluded), never half', () => {
  const y = fs.readFileSync(path.join(ROOT, 'mkdocs.yml'), 'utf8');
  const inNav = /^  - Explorer: explorer\.md$/m.test(y);
  const ex = (y.match(/^exclude_docs: \|\n((?:  .*\n)+)/m) || ['', ''])[1].split('\n').map(l => l.trim());
  const hidden = ex.includes('/explorer.md') && ex.includes('/app/');
  assert(inNav !== hidden && (hidden || (!ex.includes('/explorer.md') && !ex.includes('/app/'))), `nav ${inNav}, excluded ${ex}`);
  console.log(`  Explorer in the site: ${inNav ? 'released' : 'hidden until release'}`);
});

// --- 2. conversions and formatting
await test('units: inputs to SI (°C + 273.15, bar x 1e5, kJ/kg x 1e3, wt % / 100 or kg/kg)', () => {
  const t = page.toSI({ backend: 'IAPS84', mode: 'tpx', T: ' 350 ', p: '200', X: '3.2', Xunit: 'wt' });
  assert(t.backend === 'IAPS84' && t.mode === 'tpx' && t.a === 350 + 273.15 && t.p === 200 * 1e5 && t.X === 3.2 / 100, JSON.stringify(t));
  const h = page.toSI({ backend: 'IAPWS95', mode: 'hpx', H: '1500', p: '1e3', X: '0.032', Xunit: 'frac' });
  assert(h.backend === 'IAPWS95' && h.mode === 'hpx' && h.a === 1.5e6 && h.p === 1e8 && h.X === 0.032, JSON.stringify(h));
  const d = page.toSI({ T: '25', p: '1', X: '0' });
  assert(d.backend === 'IAPS84' && d.mode === 'tpx' && d.X === 0, 'defaults');
  for (const bad of ['', '  ', 'abc', '1,5']) assert(Number.isNaN(page.toSI({ mode: 'tpx', T: bad, p: '1', X: '0' }).a), `'${bad}' -> NaN`);
});
await test('format: 7 significant digits, dash for NaN', () => {
  const cases = [[NaN, '—'], [0, '0'], [1234.5678912, '1234.568'], [623.15, '623.15'], [1e-5, '1e-5'],
    [2.5e-10, '2.5e-10'], [-4.0566e-4, '-4.0566e-4'], [12345678, '1.234568e+7'], [0.0012345678, '0.001234568']];
  for (const [v, s] of cases) assert(page.fmt(v) === s, `fmt(${v}) = ${page.fmt(v)}, expected ${s}`);
});

// --- 3./4. with the module
if (!dist || !fs.existsSync(path.join(dist, 'xthermo_wasm.wasm'))) {
  console.log(`FAIL module tests: no module in '${dist}'`); fails++;
} else {
  const XThermo = require(path.join(dist, 'xthermo.js'));
  const xt = await XThermo.load(require(path.join(dist, 'xthermo_wasm.js')));
  const K = 273.15;
  const OK = [
    ['IAPS84', 'tpx', 350 + K, 200e5, 0.032], ['IAPS84', 'tpx', 25 + K, 1e5, 0], ['IAPS84', 'tpx', 400 + K, 250e5, 0.1],
    ['IAPS84', 'tpx', 500 + K, 300e5, 0.26], ['IAPS84', 'tpx', 600 + K, 100e5, 0], ['IAPS84', 'tpx', 300 + K, 200e5, 0.5],
    ['IAPS84', 'tpx', 380 + K, 250e5, 0], ['IAPS84', 'tpx', 900 + K, 1000e5, 1],
    ['IAPS84', 'hpx', 1.5e6, 200e5, 0.032], ['IAPS84', 'hpx', 2.6e6, 50e5, 0.1],
    ['IAPWS95', 'tpx', 350 + K, 200e5, 0.032], ['IAPWS95', 'tpx', 600 + K, 1000e5, 0.2], ['IAPWS95', 'hpx', 1.5e6, 200e5, 0.032],
  ];
  const FAIL = [
    ['IAPS84', 'tpx', 1300 + K, 100e5, 0.1, 'OutOfRangeError'], ['IAPS84', 'tpx', 350 + K, 200e5, -0.001, 'InvalidInputError'],
    ['IAPS84', 'tpx', NaN, 200e5, 0.032, 'InvalidInputError'], ['IAPWS95', 'hpx', 1e8, 100e5, 0.1, 'OutOfRangeError'],
    ['IAPS84', 'hpx', 2141025.641025641, 22.0e6, 1e-6, 'NotConvergedError'],
  ];
  await test('worker: successful states bitwise identical to direct xthermo.js calls (both backends, T-P-X and H-P-X)', () => {
    const models = {}, phases = new Set();
    for (const [backend, mode, a, p, X] of OK) {
      const r = evaluate(xt, models, { backend, mode, a, p, X });
      assert(r.ok, `${backend} ${mode} ${a} ${p} ${X}: ${JSON.stringify(r.error)}`);
      const d = xt.model(backend)[mode](a, p, X, { errors: 'raise', absent: 'nan' });
      assert(r.phase === d.phase[0] && r.phaseName === d.phaseName[0], 'phase');
      for (const f of xt.fields.concat(['massFraction_l', 'massFraction_v', 'massFraction_h'])) {
        assert(same(r.values[f], d[f][0]), `${backend} ${mode} ${f}: ${r.values[f]} vs ${d[f][0]}`);
      }
      phases.add(r.phaseName);
    }
    assert(['Liquid', 'Vapor', 'Supercritical', 'V+L', 'L+H', 'V+H', 'V+L+H'].every(ph => phases.has(ph)), `phase regions covered: ${[...phases]}`);
  });
  await test('worker: typed failures with the class, status and message of the library', () => {
    const models = {};
    for (const [backend, mode, a, p, X, cls] of FAIL) {
      const r = evaluate(xt, models, { backend, mode, a, p, X });
      let e = null;
      try { xt.model(backend)[mode](a, p, X); } catch (err) { e = err; }
      assert(!r.ok && e && r.error.name === cls && e.name === cls, `${backend} ${mode} ${a}: ${JSON.stringify(r.error)}, direct ${e && e.name}`);
      assert(r.error.status === e.status && r.error.statusName === e.statusName && r.error.message === e.libraryMessage && r.error.message.length > 10,
        `${cls}: ${JSON.stringify(r.error)}`);
    }
  });
  await test('page round trip: form -> SI -> worker -> displayed rows; multiphase placeholders shown as text', () => {
    const req = page.toSI({ backend: 'IAPS84', mode: 'tpx', T: '350', p: '200', X: '3.2', Xunit: 'wt' });
    const st = (q, req) => ({ phase: q.phase, mode: req.mode, X: req.X, S_l: q.values.S_l, S_h: q.values.S_h });
    const r = evaluate(xt, {}, req), rows = page.rows(r.values, st(r, req));
    const get = (label, unit) => rows.bulk.find(x => x[0] === label && x[1] === unit)[2];
    assert(get('Temperature', '°C') === '350' && get('Pressure', 'bar') === '200' && get('Salinity', 'wt % NaCl') === '3.2', JSON.stringify(rows.bulk));
    assert(rows.phases.every(row => row.length === 5), 'phase rows: label, unit, liquid, vapour, halite');
    assert(rows.bulk.every(x => /^-?[\d.]+(e[-+]\d+)?$/.test(x[2])), 'single-phase liquid: every bulk value numeric');
    // multiphase placeholders of the library (docs/limitations.md) are shown as text, not as numbers
    const bulkOf = (b, mode, a, p, X) => { const q = evaluate(xt, {}, { backend: b, mode, a, p, X }); return [q.phaseName, Object.fromEntries(page.rows(q.values, st(q, { mode, X })).bulk.map(x => [x[0] + ' ' + x[1], x[2]]))]; };
    let [ph, b] = bulkOf('IAPS84', 'tpx', 400 + K, 250e5, 0.1);
    assert(ph === 'V+L' && b['Dynamic viscosity Pa·s'] === 'not defined' && b['Isobaric expansivity 1/K'] === 'not computed'
      && /^\d/.test(b['Isobaric heat capacity J/(kg·K)']) && /^\d/.test(b['Isothermal compressibility 1/Pa']), `V+L: ${JSON.stringify(b)}`);
    [ph, b] = bulkOf('IAPS84', 'hpx', 1.5e6, 50e5, 0);
    assert(ph === 'V+L' && b['Isothermal compressibility 1/Pa'] === 'not computed' && b['Dynamic viscosity Pa·s'] === 'not defined', `pure-water V+L (H-P-X): ${JSON.stringify(b)}`);
    for (const f of ['Mu', 'IsobaricExpansivity', 'Cp', 'IsothermalCompressibility']) {
      for (const phase of [0, 1, 2, 3, 4, 5, 6]) assert(page.placeholder(f, { phase, mode: 'hpx', X: 0, S_l: 1, S_h: 0 }) === null, `single phase ${phase} ${f}`);
    }
    const empty = evaluate(xt, {}, page.toSI({ mode: 'tpx', T: '', p: '200', X: '3.2' }));
    assert(!empty.ok && empty.error.name === 'InvalidInputError', 'empty temperature -> InvalidInput from the library');
  });
  await test('placeholder rule (phase region, mode, composition) matches where the library returns placeholders', () => {
    // over T-P-X grids and H-P-X enthalpy scans, both backends: for every multiphase state the rule marks
    // IsobaricExpansivity (library 0), Mu (library 1 except V+L+H and pure-water V+L), Cp and
    // IsothermalCompressibility exactly where the library returns 0; no single-phase field is marked
    let n = 0; const seen = new Set();
    for (const backend of ['IAPS84', 'IAPWS95']) {
      const m = xt.model(backend);
      const cases = [];
      const T = [], P = [], Xs = [];
      for (let t = 1; t <= 1000; t += 13) for (const p of [1, 10, 50, 100, 200, 300, 500, 1000]) for (const x of [0, 1e-6, 0.01, 0.1, 0.3, 0.6, 1]) { T.push(t + K); P.push(p * 1e5); Xs.push(x); }
      cases.push(['tpx', T, P, Xs]);
      const H = [], P2 = [], X2 = [];
      for (const [p, x] of [[50e5, 0], [200e5, 0], [100e5, 0.3], [200e5, 0.1], [300e5, 0.5], [1e5, 1], [100e5, 1], [100e5, 1e-6]]) for (let h = 1e4; h < 4e6; h += 4e3) { H.push(h); P2.push(p); X2.push(x); }
      cases.push(['hpx', H, P2, X2]);
      for (const [mode, a, p, X] of cases) {
        const s = m[mode](a, p, X, { errors: 'nan' });
        for (let i = 0; i < s.n; i++) {
          if (s.status[i]) continue;
          const stt = { phase: s.phase[i], mode, X: X[i], S_l: s.S_l[i], S_h: s.S_h[i] };
          const mp = stt.phase >= 7 && stt.phase <= 11;
          for (const f of ['Cp', 'IsothermalCompressibility', 'IsobaricExpansivity']) {
            const marked = page.placeholder(f, stt) !== null;
            const libZero = s[f][i] === 0 && mp;
            assert(marked === libZero, `${backend} ${mode} #${i} phase ${s.phaseName[i]} X=${X[i]} ${f}=${s[f][i]}: marked ${marked}`);
          }
          const muMarked = page.placeholder('Mu', stt) !== null;
          assert(muMarked === mp, `${backend} ${mode} #${i} Mu`);
          if (mp && stt.phase !== 11 && !(mode === 'hpx' && X[i] === 0)) assert(s.Mu[i] === 1, `${backend} ${mode} #${i} Mu=${s.Mu[i]} in ${s.phaseName[i]}`);
          if (mp) { n++; seen.add(s.phaseName[i] + (mode === 'hpx' && (X[i] === 0 || X[i] === 1) ? ` X=${X[i]} hpx` : '')); }
        }
      }
    }
    for (const need of ['V+L', 'L+H', 'V+H', 'V+L+H', 'V+L X=0 hpx', 'L+H X=1 hpx']) assert(seen.has(need), `scan did not reach ${need}: ${[...seen]}`);
    assert(n > 5000, `multiphase states checked: ${n}`);
  });
  await test('worker protocol: init -> ready (build information), eval -> result; missing module -> unavailable', async () => {
    const compiled = await WebAssembly.compile(fs.readFileSync(path.join(dist, 'xthermo_wasm.wasm')));
    function scope(dir) {
      const out = [];
      const ctx = { console, WebAssembly, performance, setTimeout, clearTimeout, TextDecoder, URL,
        location: { href: 'file://' + path.join(dir, 'worker.js') }, postMessage: m => out.push(m) };
      ctx.self = ctx;
      ctx.importScripts = (...files) => { for (const f of files) vm.runInContext(fs.readFileSync(path.join(dir, f), 'utf8'), ctx, { filename: f }); };
      vm.createContext(ctx);
      vm.runInContext(fs.readFileSync(path.join(APP, 'worker.js'), 'utf8'), ctx, { filename: 'worker.js' });
      const wait = async n => { for (let i = 0; i < 400 && out.length < n; i++) await new Promise(r => setTimeout(r, 10)); return out; };
      return { send: m => ctx.self.onmessage({ data: m }), wait };
    }
    const w = scope(dist);
    w.send({ type: 'init', module: compiled });
    let out = await w.wait(1);
    assert(out[0] && out[0].type === 'ready' && out[0].buildInfo && out[0].buildInfo.abi === XThermo.ABI_VERSION, JSON.stringify(out[0]));
    w.send({ type: 'eval', id: 7, backend: 'IAPWS95', mode: 'tpx', a: 350 + K, p: 200e5, X: 0.032 });
    out = await w.wait(2);
    const d = xt.model('IAPWS95').tpx(350 + K, 200e5, 0.032, { absent: 'nan' });
    assert(out[1].type === 'result' && out[1].id === 7 && out[1].ok && same(out[1].values.Rho, d.Rho[0]) && out[1].ms >= 0, JSON.stringify(out[1]).slice(0, 200));
    w.send({ type: 'eval', id: 8, backend: 'IAPS84', mode: 'tpx', a: 1300 + K, p: 100e5, X: 0.1 });
    out = await w.wait(3);
    assert(out[2].id === 8 && !out[2].ok && out[2].error.statusName === 'OutOfRange', JSON.stringify(out[2]));
    const empty = fs.mkdtempSync(path.join(fs.realpathSync(process.env.TMPDIR || '/tmp'), 'xt-nomodule-'));
    const u = scope(empty);
    u.send({ type: 'init', module: compiled });
    out = await u.wait(1);
    fs.rmdirSync(empty);
    assert(out[0] && out[0].type === 'unavailable', JSON.stringify(out[0]));
  });
}
console.log(`explorer interface tests: ${results.filter(r => r.startsWith('PASS')).length} passed, ${fails} failed`);
process.exit(fails ? 1 : 0);
