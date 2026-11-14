/*
 * xThermo Explorer: evaluation worker (W2b; review/W2B_EXPLORER_PAGE.md).
 * Loads the WebAssembly module (xthermo_wasm.js + .wasm, compiled once by the page and passed in) and its
 * JavaScript interface xthermo.js, and evaluates one state per request with the library's T-P-X or H-P-X
 * function. Inputs and results are in SI units; the page converts. No approximation, no fallback: every value
 * and every error comes from the library.
 *
 * Messages in:  {type: 'init', module: WebAssembly.Module}
 *               {type: 'eval', id, backend: 'IAPS84'|'IAPWS95', mode: 'tpx'|'hpx', a: T [K] or H [J/kg], p [Pa], X}
 * Messages out: {type: 'ready', buildInfo} | {type: 'unavailable', reason}
 *               {type: 'result', id, ms, ok: true, phase, phaseName, values: {field: number}}
 *               {type: 'result', id, ms, ok: false, error: {name, status, statusName, message}}
 */
(function () {
  'use strict';
  const EXTRA = ['massFraction_l', 'massFraction_v', 'massFraction_h'];

  // One state. errors 'raise': a failed state throws the typed error of the library status (InvalidInputError,
  // OutOfRangeError, NotConvergedError, ...) carrying the library's own message. absent 'nan': the fields of a phase
  // that is not present are NaN.
  function evaluate(xt, models, req) {
    try {
      const m = models[req.backend] || (models[req.backend] = xt.model(req.backend));
      if (req.mode !== 'tpx' && req.mode !== 'hpx') throw new Error("mode must be 'tpx' or 'hpx'");
      const s = m[req.mode](req.a, req.p, req.X, { errors: 'raise', absent: 'nan' });
      const values = {};
      for (const f of xt.fields.concat(EXTRA)) values[f] = s[f][0];
      return { ok: true, phase: s.phase[0], phaseName: s.phaseName[0], values };
    } catch (e) {
      if (e && typeof e.status === 'number') {
        return { ok: false, error: { name: e.name, status: e.status, statusName: e.statusName, message: e.libraryMessage } };
      }
      return { ok: false, error: { name: (e && e.name) || 'Error', status: null, statusName: 'Error', message: String(e && e.message || e) } };
    }
  }

  if (typeof module === 'object' && module.exports) { module.exports = { evaluate }; return; }

  let xt = null;
  const models = {};
  self.onmessage = async function (ev) {
    const msg = ev.data;
    if (msg.type === 'init') {
      try {
        importScripts('xthermo_wasm.js', 'xthermo.js');
      } catch (e) {
        self.postMessage({ type: 'unavailable', reason: 'the module files are not part of this build' });
        return;
      }
      try {
        xt = await self.XThermo.load(self.createXThermoModule, {
          instantiateWasm(imports, receive) {
            WebAssembly.instantiate(msg.module, imports).then(inst => receive(inst, msg.module),
              err => self.postMessage({ type: 'unavailable', reason: String(err) }));
            return {};
          },
        });
        self.postMessage({ type: 'ready', buildInfo: xt.buildInfo });
      } catch (e) {
        self.postMessage({ type: 'unavailable', reason: String(e && e.message || e) });
      }
    } else if (msg.type === 'eval') {
      const t0 = performance.now();
      const r = xt ? evaluate(xt, models, msg)
        : { ok: false, error: { name: 'Error', status: null, statusName: 'Error', message: 'module not loaded' } };
      self.postMessage(Object.assign({ type: 'result', id: msg.id, ms: performance.now() - t0 }, r));
    }
  };
})();
