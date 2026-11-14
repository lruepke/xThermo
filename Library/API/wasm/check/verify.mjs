// W2a: native-versus-WebAssembly verification (review/W2A_WASM_MODULE.md). Evaluates the inputs of the native
// reference (native_reference.py, same commit and container) with the WebAssembly module through its JavaScript
// interface xthermo.js, and compares point by point. Node only, no packages.
//
//   node verify.mjs <dist dir with xthermo.js, xthermo_wasm.js/.wasm> <reference.json.gz> [report.txt]
//
// Pass criteria (every point of every set, both backends):
//   status, phase code, phase name and the Python/JS error class identical;
//   failure messages identical after replacing every number by a placeholder, and the numbers equal within 1e-9
//     relative (the messages print values with 17 digits, whose last digits are platform round-off);
//   every field bitwise equal, both NaN, or within 1e-9 relative; the saturations S_l, S_v, S_h (volume fractions,
//     from the lever rule: absolute round-off on a fraction of unity, large relative differences when a phase
//     fraction is tiny) in addition within 1e-10 absolute; IsothermalCompressibility in addition within
//     1e-11 1/Pa absolute (a 1 Pa central difference of the bulk density in L+H and V+H: absolute round-off noise of
//     about 1e-12 1/Pa, a large relative difference where the compressibility is close to zero; review/W1);
//   absent="nan": the fields of absent phases NaN exactly where Python's has() says the phase is absent, the
//     present ones unchanged; mass fractions of an evaluated point sum to 1 within 1e-12;
//   build information, field list, status names and the name of every PhaseRegion code as in the native library.
// Exit status 0 only if everything passes. The report lists every violation (up to 30 per category) and, per field,
// the number of non-bitwise points and the largest relative difference.
import fs from 'fs';
import zlib from 'zlib';
import path from 'path';
import { createRequire } from 'module';

const [dist, refPath, reportPath] = process.argv.slice(2);
const require = createRequire(import.meta.url);
const XThermo = require(path.resolve(dist, 'xthermo.js'));
const xt = await XThermo.load(require(path.resolve(dist, 'xthermo_wasm.js')));
const ref = JSON.parse(zlib.gunzipSync(fs.readFileSync(refPath)));

const RTOL = 1e-9, S_ATOL = 1e-10, KAPPA_ATOL = 1e-11, MSG_RTOL = 1e-9, MASS_TOL = 1e-12;
const lines = [];
const out = s => { lines.push(s); console.log(s); };
const violations = {};
const violate = (cat, text) => { (violations[cat] = violations[cat] || []).push(text); };
const hex = s => {
  if (s === 'nan') return NaN; if (s === 'inf') return Infinity; if (s === '-inf') return -Infinity;
  const m = /^(-?)0x([01])\.?([0-9a-f]*)p([+-]\d+)$/.exec(s);
  const v = parseInt(m[2]) + (m[3] ? parseInt(m[3], 16) / 16 ** m[3].length : 0);
  return (m[1] ? -1 : 1) * v * 2 ** parseInt(m[4]);
};
const NUM = /[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?|nan|inf/g;
function sameMessage(a, b) {
  if (a.replace(NUM, '#') !== b.replace(NUM, '#')) return false;
  const x = a.match(NUM) || [], y = b.match(NUM) || [];
  return x.every((s, i) => { const u = Number(s), v = Number(y[i]);
    return s === y[i] || Math.abs(u - v) <= MSG_RTOL * Math.abs(u); });
}

// --- static checks
if (JSON.stringify(xt.fields) !== JSON.stringify(ref.fields)) violate('interface', `fields ${xt.fields} != native ${ref.fields}`);
for (const [code, name] of Object.entries(ref.phase_names)) {
  if (xt.phaseName(Number(code)) !== name) violate('interface', `phase name of ${code}: wasm '${xt.phaseName(Number(code))}', native '${name}'`);
}
const STATUS = ['Ok', 'InvalidInput', 'OutOfRange', 'NotConverged', 'NotAvailable', 'Failed'];
if (JSON.stringify(xt.statusNames) !== JSON.stringify(STATUS)) violate('interface', `status names ${xt.statusNames}`);
out(`module: ${JSON.stringify(xt.buildInfo)}`);

// --- per backend and set
const fields = ref.fields;
const totals = { points: 0, bitwisePoints: 0 };
for (const [backend, recs] of Object.entries(ref.records)) {
  const m = xt.model(backend);
  const sets = [...new Set(recs.map(r => r.set))];
  for (const set of sets) {
    const rs = recs.filter(r => r.set === set), mode = rs[0].mode;
    const [a, p, X] = [0, 1, 2].map(j => Float64Array.from(rs, r => hex(r.in[j])));
    const t0 = performance.now();
    const s = mode === 0 ? m.tpx(a, p, X, { errors: 'nan' }) : m.hpx(a, p, X, { errors: 'nan' });
    const dt = (performance.now() - t0) / 1000;
    const stat = Object.fromEntries(fields.map(f => [f, { n: 0, max: 0, maxAbs: 0, where: '' }]));
    let bitwise = 0;
    const msgs = Object.fromEntries(s.failures.map(f => [f.index, f.message]));
    rs.forEach((r, i) => {
      const id = `${backend} ${set} #${i} (${mode ? 'H' : 'T'}=${a[i]}, p=${p[i]}, X=${X[i]})`;
      if (s.status[i] !== r.status) violate('status', `${id}: wasm ${STATUS[s.status[i]]}, native ${STATUS[r.status]}`);
      if (s.phase[i] !== r.phase) violate('phase', `${id}: wasm ${s.phaseName[i]}, native ${ref.phase_names[r.phase]}`);
      if (r.status && s.status[i] === r.status && !sameMessage(msgs[i] || '', r.message))
        violate('message', `${id}:\n      wasm   ${msgs[i]}\n      native ${r.message}`);
      let same = true;
      fields.forEach((f, k) => {
        const w = s[f][i], n = hex(r.f[k]);
        if (Object.is(w, n) || (Number.isNaN(w) && Number.isNaN(n))) return;
        same = false;
        const d = Math.abs(w - n), rel = d / Math.abs(n);
        const st = stat[f]; st.n++; if (d > st.maxAbs) st.maxAbs = d;
        if (rel > st.max || Number.isNaN(rel)) { st.max = rel; st.where = id; }
        const ok = d <= RTOL * Math.abs(n) || (f.startsWith('S_') && d <= S_ATOL)
          || (f === 'IsothermalCompressibility' && d <= RTOL * Math.abs(n) + KAPPA_ATOL);
        if (!ok) violate('field', `${id} ${f}: wasm ${w}, native ${n} (rel ${rel.toExponential(2)})`);
      });
      if (same) bitwise++;
      if (s.status[i] === 0) {
        const sum = s.massFraction_l[i] + s.massFraction_v[i] + s.massFraction_h[i];
        if (!(Math.abs(sum - 1) <= MASS_TOL)) violate('mass fraction', `${id}: sum ${sum}`);
      }
    });
    // absent="nan" rule (T-P-X sets; the mask depends only on the phase code)
    if (mode === 0) {
      const sn = m.tpx(a, p, X, { errors: 'nan', absent: 'nan' });
      rs.forEach((r, i) => {
        fields.forEach(f => {
          if (!/_[lvh]$/.test(f) || f.startsWith('S_')) return;
          const present = r.present.includes(f[f.length - 1]);
          const okv = present ? Object.is(sn[f][i], s[f][i]) || (Number.isNaN(sn[f][i]) && Number.isNaN(s[f][i])) : Number.isNaN(sn[f][i]);
          if (!okv) violate('absent', `${backend} ${set} #${i} ${f}: present=${present}, value ${sn[f][i]}`);
        });
      });
    }
    // errors="raise": the error class of every failed point
    rs.forEach((r, i) => {
      if (!r.status) return;
      let cls = null;
      try { mode === 0 ? m.tpx(a[i], p[i], X[i]) : m.hpx(a[i], p[i], X[i]); } catch (e) { cls = e.name; }
      if (cls !== r.raise) violate('error class', `${backend} ${set} #${i}: wasm ${cls}, native ${r.raise}`);
    });
    totals.points += rs.length; totals.bitwisePoints += bitwise;
    const nat = ref.timing[`${backend} ${set}`];
    out(`${backend.padEnd(8)} ${set.padEnd(13)} n=${String(rs.length).padStart(5)}  failed=${String(s.failures.length).padStart(4)}` +
        `  all fields bitwise ${String(bitwise).padStart(5)} (${(100 * bitwise / rs.length).toFixed(1)} %)` +
        `  wasm ${dt.toFixed(2)} s, native ${nat[1].toFixed(2)} s`);
    for (const f of fields) {
      const st = stat[f];
      if (st.n) out(`           ${f.padEnd(26)} not bitwise ${String(st.n).padStart(5)}  max rel ${st.max.toExponential(2)}` +
                    `  max abs ${st.maxAbs.toExponential(2)}  at ${st.where}`);
    }
  }
  m.free();
}

out(`\npoints ${totals.points}, all fields bitwise ${totals.bitwisePoints} (${(100 * totals.bitwisePoints / totals.points).toFixed(1)} %)`);
const cats = Object.keys(violations);
if (!cats.length) {
  out('RESULT: PASS (status, phase, messages, error classes, fields, absent-phase rule, mass fractions, interface)');
} else {
  out('RESULT: FAIL');
  for (const c of cats) {
    out(`FAIL ${c}: ${violations[c].length}`);
    violations[c].slice(0, 30).forEach(v => out('   ' + v));
  }
}
if (reportPath) fs.writeFileSync(reportPath, lines.join('\n') + '\n');
process.exit(cats.length ? 1 : 0);
