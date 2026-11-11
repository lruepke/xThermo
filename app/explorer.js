/*
 * xThermo Explorer: page logic (W2b; review/W2B_EXPLORER_PAGE.md). Plain JavaScript, no dependencies.
 * Reads the form, converts the inputs to SI, sends one state to the evaluation worker (worker.js) and shows the
 * library's result or its typed error. The WebAssembly module is compiled once; Cancel terminates the worker and
 * starts a new one from the compiled module.
 * The unit conversions and formatting (UNITS, toSI, fmt, rows) are pure and exported for the Node tests.
 */
(function (root, factory) {
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.addEventListener('DOMContentLoaded', () => api.start(root.document));
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';
  const K = 273.15;

  // inputs: display unit -> SI
  const UNITS = {
    T: { unit: '°C', toSI: v => v + K },
    H: { unit: 'kJ/kg', toSI: v => v * 1e3 },
    p: { unit: 'bar', toSI: v => v * 1e5 },
    X: { wt: { unit: 'wt % NaCl', toSI: v => v / 100 }, frac: { unit: 'kg/kg (mass fraction)', toSI: v => v } },
  };

  // an empty or non-numeric field becomes NaN; the library then reports InvalidInput with its own message
  function num(s) { const t = String(s).trim(); return t === '' ? NaN : Number(t); }

  function toSI(form) {
    const mode = form.mode === 'hpx' ? 'hpx' : 'tpx';
    return {
      backend: form.backend === 'IAPWS95' ? 'IAPWS95' : 'IAPS84',
      mode,
      a: mode === 'tpx' ? UNITS.T.toSI(num(form.T)) : UNITS.H.toSI(num(form.H)),
      p: UNITS.p.toSI(num(form.p)),
      X: (form.Xunit === 'frac' ? UNITS.X.frac : UNITS.X.wt).toSI(num(form.X)),
    };
  }

  // seven significant digits; NaN (absent phase or no such property) as a dash
  function fmt(v) {
    if (typeof v !== 'number' || Number.isNaN(v)) return '—';
    if (v === 0) return '0';
    const a = Math.abs(v);
    if (a >= 1e-3 && a < 1e7) return String(Number(v.toPrecision(7)));
    return v.toExponential(6).replace(/\.?0+e/, 'e');
  }

  // output rows: label, unit, SI -> display, field (bulk) or fields per phase (liquid, vapour, halite; null = none)
  const BULK = [
    ['Temperature', '°C', v => v - K, 'T'],
    ['Temperature', 'K', v => v, 'T'],
    ['Pressure', 'bar', v => v / 1e5, 'p'],
    ['Salinity', 'wt % NaCl', v => v * 100, 'X'],
    ['Specific enthalpy', 'kJ/kg', v => v / 1e3, 'H'],
    ['Density', 'kg/m³', v => v, 'Rho'],
    ['Isobaric heat capacity', 'J/(kg·K)', v => v, 'Cp'],
    ['Dynamic viscosity', 'Pa·s', v => v, 'Mu'],
    ['Isothermal compressibility', '1/Pa', v => v, 'IsothermalCompressibility'],
    ['Isobaric expansivity', '1/K', v => v, 'IsobaricExpansivity'],
  ];
  const PHASES = [
    ['Volume fraction', '–', v => v, ['S_l', 'S_v', 'S_h']],
    ['Mass fraction', '–', v => v, ['massFraction_l', 'massFraction_v', 'massFraction_h']],
    ['Salinity', 'wt % NaCl', v => v * 100, ['X_l', 'X_v', null]],
    ['Density', 'kg/m³', v => v, ['Rho_l', 'Rho_v', 'Rho_h']],
    ['Specific enthalpy', 'kJ/kg', v => v / 1e3, ['H_l', 'H_v', 'H_h']],
    ['Isobaric heat capacity', 'J/(kg·K)', v => v, ['Cp_l', 'Cp_v', 'Cp_h']],
    ['Dynamic viscosity', 'Pa·s', v => v, ['Mu_l', 'Mu_v', null]],
  ];
  // Bulk fields that the library fills with a placeholder (docs/limitations.md, "Multiphase properties"), decided by
  // phase region, input mode and composition (never by the returned number); the values themselves are not altered:
  //   every multiphase state (codes 7-11): Mu has no effective multiphase meaning (1, or a volume-weighted average
  //     in V+L+H and pure-water V+L), IsobaricExpansivity is not computed;
  //   V+L+H (11): Cp and IsothermalCompressibility are not computed;
  //   pure water (X = 0) in vapour-liquid states from H-P-X: IsothermalCompressibility is not computed;
  //   pure NaCl (X = 1) with coexisting liquid and halite from H-P-X (the melting interval): IsothermalCompressibility
  //     is not computed.
  const MULTIPHASE = new Set([7, 8, 9, 10, 11]);
  const VL = new Set([7, 10]);
  function placeholder(f, st) {
    const ph = st.phase;
    if (!MULTIPHASE.has(ph)) return null;
    if (f === 'Mu') return 'not defined';
    if (f === 'IsobaricExpansivity') return 'not computed';
    if (f === 'Cp' && ph === 11) return 'not computed';
    if (f === 'IsothermalCompressibility') {
      if (ph === 11) return 'not computed';
      if (st.mode === 'hpx' && st.X === 0 && VL.has(ph)) return 'not computed';
      if (st.mode === 'hpx' && st.X === 1 && ph === 8 && st.S_l > 0 && st.S_h > 0) return 'not computed';
    }
    return null;
  }
  // st: {phase, mode, X, S_l, S_h} of the evaluated state
  function rows(values, st) {
    const bulk = BULK.map(([label, unit, conv, f]) => [label, unit, (st && placeholder(f, st)) || fmt(conv(values[f]))]);
    const phases = PHASES.map(([label, unit, conv, fs]) => [label, unit, ...fs.map(f => (f ? fmt(conv(values[f])) : '—'))]);
    return { bulk, phases };
  }

  // --- page
  function start(doc) {
    const $ = id => doc.getElementById(id);
    const form = $('xt-form'), statusLine = $('xt-status'), evalBtn = $('xt-evaluate'), cancelBtn = $('xt-cancel');
    let compiled = null, worker = null, busy = null, nextId = 1;

    function setStatus(text, kind) { statusLine.textContent = text; statusLine.dataset.kind = kind || ''; }
    function setBusy(b) { evalBtn.disabled = !!b || !worker; cancelBtn.disabled = !b; form.setAttribute('aria-busy', b ? 'true' : 'false'); }

    function readForm() {
      const fd = new FormData(form), o = {};
      for (const [k, v] of fd.entries()) o[k] = v;
      return o;
    }

    function syncMode() {
      const hpx = form.elements.mode.value === 'hpx';
      $('xt-row-T').hidden = hpx; $('xt-row-H').hidden = !hpx;
    }

    function unavailable(reason) {
      worker = null; setBusy(false); evalBtn.disabled = true;
      setStatus('Explorer not available in this build: ' + reason + '. The documentation is unaffected.', 'error');
    }

    function startWorker() {
      worker = new Worker('worker.js');
      worker.onmessage = ev => onMessage(ev.data);
      worker.onerror = ev => { ev.preventDefault(); unavailable(ev.message || 'the worker failed to start'); };
      worker.postMessage({ type: 'init', module: compiled });
    }

    function onMessage(m) {
      if (m.type === 'ready') {
        showBuild(m.buildInfo); setBusy(false); setStatus('Ready.', 'ok');
      } else if (m.type === 'unavailable') {
        unavailable(m.reason);
      } else if (m.type === 'result' && busy && m.id === busy.id) {
        const req = busy; busy = null; setBusy(false);
        if (m.ok) { showResult(req, m); setStatus(`Evaluated in ${m.ms < 1 ? '< 1' : Math.round(m.ms)} ms.`, 'ok'); }
        else { showError(m.error); setStatus('The library reported an error for this state.', 'error'); }
      }
    }

    function clear(el) { while (el.firstChild) el.removeChild(el.firstChild); }
    function cell(tr, tag, text) { const c = doc.createElement(tag); c.textContent = text; tr.appendChild(c); return c; }

    function showResult(req, m) {
      $('xt-error').hidden = true; $('xt-result').hidden = false;
      $('xt-phase').textContent = m.phaseName;
      $('xt-backend-used').textContent = req.backend === 'IAPWS95' ? 'IAPWS-95' : 'IAPS-84';
      $('xt-mode-used').textContent = req.mode === 'hpx' ? 'H-P-X' : 'T-P-X';
      const r = rows(m.values, { phase: m.phase, mode: req.mode, X: req.X, S_l: m.values.S_l, S_h: m.values.S_h });
      const tb = $('xt-bulk'); clear(tb);
      for (const [label, unit, v] of r.bulk) { const tr = doc.createElement('tr'); cell(tr, 'th', label); cell(tr, 'td', v).className = 'num'; cell(tr, 'td', unit); tb.appendChild(tr); }
      const tp = $('xt-phases'); clear(tp);
      for (const [label, unit, ...vs] of r.phases) {
        const tr = doc.createElement('tr'); cell(tr, 'th', label); for (const v of vs) cell(tr, 'td', v).className = 'num'; cell(tr, 'td', unit); tp.appendChild(tr);
      }
    }

    function showError(e) {
      $('xt-result').hidden = true; $('xt-error').hidden = false;
      $('xt-error-kind').textContent = `${e.statusName}${e.name && e.name !== 'Error' ? ' (' + e.name + ')' : ''}`;
      $('xt-error-message').textContent = e.message;
    }

    function showBuild(b) {
      $('xt-build').textContent = `xThermo ${b.xthermo}, commit ${String(b.commit).slice(0, 12)}, GSL ${b.gsl}, Emscripten ${b.emscripten}, interface ${b.abi}`;
      const src = $('xt-source-commit');
      src.href = `https://github.com/lruepke/xThermo/tree/${encodeURIComponent(b.commit)}`;
      src.textContent = `xThermo source at commit ${String(b.commit).slice(0, 12)}`;
      $('xt-build-script').href = `https://github.com/lruepke/xThermo/blob/${encodeURIComponent(b.commit)}/Library/API/wasm/make_wasm.sh`;
      const gsl = $('xt-source-gsl');
      gsl.href = `source/gsl-${encodeURIComponent(b.gsl)}.tar.gz`; gsl.textContent = `GSL ${b.gsl} source (gsl-${b.gsl}.tar.gz)`;
    }

    form.addEventListener('submit', ev => {
      ev.preventDefault();
      if (!worker || busy) return;
      const req = toSI(readForm());
      busy = Object.assign({ id: nextId++ }, req);
      setBusy(true); setStatus('Evaluating…', 'busy');
      worker.postMessage(Object.assign({ type: 'eval' }, busy));
    });
    cancelBtn.addEventListener('click', () => {
      if (!busy) return;
      busy = null; worker.terminate(); worker = null;
      setBusy(false); evalBtn.disabled = true; setStatus('Cancelled. Restarting the module…', 'busy');
      startWorker();
    });
    for (const r of form.elements.mode) r.addEventListener('change', syncMode);
    syncMode();

    // compile the module once (streaming where the server sends application/wasm); a missing module only disables
    // the Explorer
    setBusy(false); evalBtn.disabled = true; setStatus('Loading the xThermo module…', 'busy');
    if (typeof WebAssembly !== 'object' || typeof Worker !== 'function') { unavailable('this browser does not support WebAssembly workers'); return; }
    fetch('xthermo_wasm.wasm').then(resp => {
      if (!resp.ok) throw new Error('the module files are not part of this build');
      return (WebAssembly.compileStreaming && resp.headers.get('content-type') === 'application/wasm')
        ? WebAssembly.compileStreaming(resp) : resp.arrayBuffer().then(b => WebAssembly.compile(b));
    }).then(mod => { compiled = mod; startWorker(); }, err => unavailable(err.message || String(err)));
  }

  return { UNITS, toSI, fmt, rows, placeholder, BULK, PHASES, start };
});
