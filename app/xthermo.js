/*
 * xthermo.js: JavaScript interface of the xThermo WebAssembly module (xthermo_wasm.js + xthermo_wasm.wasm).
 * No dependencies. Browser (window or Worker, after importScripts('xthermo_wasm.js')) and Node (require).
 *
 *   const xt = await XThermo.load(createXThermoModule);        // Node: XThermo.load(require('./xthermo_wasm.js'))
 *   const m = xt.model('IAPS84');                               // or 'IAPWS95'
 *   const s = m.tpx(573.15, 200e5, 0.032);                      // T [K], p [Pa], X [kg/kg]; numbers or arrays
 *   s.phaseName[0], s.Rho[0], s.status[0], s.failures            // Float64Array per field, SI units
 *   m.hpx(1.2e6, 200e5, 0.032, {errors: 'nan', absent: 'nan'});
 *   m.free();
 *
 * Semantics as the Python functions xthermo.tpx / xthermo.hpx (states.py):
 *   errors 'raise' (default): throws the EvaluationError subclass of the first failed point (InvalidInputError,
 *     OutOfRangeError, NotConvergedError, NotAvailableError; EvaluationError for Failed) with .status, .index,
 *     .inputs, .message; 'nan': failed points have NaN fields, phase Unknown, their status and message in failures.
 *   absent 'zero' (default): the library's values; 'nan': fields of a phase that is absent at a point are NaN
 *     (saturations and bulk fields unchanged).
 * Saturations S_l, S_v, S_h are volume fractions (library convention). massFraction_l/_v/_h are derived from them as
 * S_i Rho_i / Rho (0 for an absent phase, NaN for a failed point).
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.XThermo = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';
  const ABI_VERSION = 1;
  const FIELDS = ['T', 'p', 'X', 'H', 'S_l', 'S_v', 'S_h', 'X_l', 'X_v', 'Rho', 'Rho_l', 'Rho_v', 'Rho_h', 'H_l', 'H_v',
    'H_h', 'Cp', 'Cp_l', 'Cp_v', 'Cp_h', 'Mu', 'Mu_l', 'Mu_v', 'IsothermalCompressibility', 'IsobaricExpansivity'];
  const BACKENDS = { IAPS84: 0, IAPWS95: 1 };
  const STATUS = ['Ok', 'InvalidInput', 'OutOfRange', 'NotConverged', 'NotAvailable', 'Failed'];
  // phases present per PhaseRegion code, as states.py _PRESENT (single-phase fluid states in the liquid fields)
  const PRESENT = { 0: 'l', 1: 'v', 2: 'h', 3: 'l', 4: 'l', 5: 'l', 6: 'l', 7: 'lv', 8: 'lh', 9: 'vh', 10: 'lv', 11: 'lvh' };
  const PHASE_FIELDS = FIELDS.filter(f => /_[lvh]$/.test(f) && !f.startsWith('S_'));

  class EvaluationError extends Error {
    constructor(text, status, index, inputs, message) {
      super(text);
      this.name = this.constructor.name;
      this.status = status; this.statusName = STATUS[status]; this.index = index; this.inputs = inputs;
      this.libraryMessage = message;
    }
  }
  class InvalidInputError extends EvaluationError {}
  class OutOfRangeError extends EvaluationError {}
  class NotConvergedError extends EvaluationError {}
  class NotAvailableError extends EvaluationError {}
  const ERRORS = [null, InvalidInputError, OutOfRangeError, NotConvergedError, NotAvailableError, EvaluationError];

  function toArray(v, n) {
    if (typeof v === 'number') return new Float64Array(n).fill(v);
    const a = Float64Array.from(v);
    if (a.length !== n) throw new RangeError('the inputs must be numbers or arrays of the same length');
    return a;
  }

  async function load(createModule, options) {
    const M = await createModule(options || {});
    const str = ptr => M.UTF8ToString(ptr);
    if (M._xt_abi_version() !== ABI_VERSION) throw new Error(`xthermo_wasm ABI ${M._xt_abi_version()}, xthermo.js expects ${ABI_VERSION}`);
    const nf = M._xt_nfields();
    const names = Array.from({ length: nf }, (_, k) => str(M._xt_field_name(k)));
    if (names.join() !== FIELDS.join()) throw new Error('xthermo_wasm field list does not match xthermo.js');
    const buildInfo = JSON.parse(str(M._xt_build_info()));
    const phaseNames = {};
    const phaseName = code => (code in phaseNames ? phaseNames[code] : (phaseNames[code] = str(M._xt_phase_name(code))));

    function evaluate(handle, mode, a, p, X, opts) {
      const errors = (opts && opts.errors) || 'raise', absent = (opts && opts.absent) || 'zero';
      if (errors !== 'raise' && errors !== 'nan') throw new Error("errors must be 'raise' or 'nan'");
      if (absent !== 'zero' && absent !== 'nan') throw new Error("absent must be 'zero' or 'nan'");
      const n = typeof a === 'number' ? (typeof p === 'number' ? (typeof X === 'number' ? 1 : X.length) : p.length) : a.length;
      const va = toArray(a, n), vp = toArray(p, n), vX = toArray(X, n);
      const bytes = 8 * Math.max(n, 1);
      const ptrs = [M._malloc(bytes), M._malloc(bytes), M._malloc(bytes), M._malloc(bytes * nf), M._malloc(4 * Math.max(n, 1)),
        M._malloc(4 * Math.max(n, 1))];
      try {
        const [pa, pp, pX, po, ps, ph] = ptrs;
        M.HEAPF64.set(va, pa / 8); M.HEAPF64.set(vp, pp / 8); M.HEAPF64.set(vX, pX / 8);
        const nfail = M._xt_eval(handle, mode, n, pa, pp, pX, po, ps, ph);
        if (nfail < 0) throw new Error('xt_eval: ' + str(M._xt_last_error()));
        const s = { n, status: M.HEAP32.slice(ps / 4, ps / 4 + n), phase: M.HEAP32.slice(ph / 4, ph / 4 + n), failures: [] };
        FIELDS.forEach((f, k) => { s[f] = M.HEAPF64.slice(po / 8 + k * n, po / 8 + (k + 1) * n); });
        for (let i = 0; i < n && s.failures.length < nfail; i++) {
          if (s.status[i] !== 0) s.failures.push({ index: i, status: s.status[i], statusName: STATUS[s.status[i]], message: str(M._xt_message(i)) });
        }
        s.phaseName = Array.from(s.phase, phaseName);
        s.statusName = Array.from(s.status, c => STATUS[c]);
        if (absent === 'nan') {
          for (const f of PHASE_FIELDS) {
            const slot = f[f.length - 1], v = s[f];
            for (let i = 0; i < n; i++) if (!(PRESENT[s.phase[i]] || '').includes(slot)) v[i] = NaN;
          }
        }
        for (const slot of 'lvh') {
          const w = new Float64Array(n), S = s['S_' + slot], R = s['Rho_' + slot];
          for (let i = 0; i < n; i++) {
            if (s.status[i] !== 0) w[i] = NaN;
            else w[i] = (PRESENT[s.phase[i]] || '').includes(slot) ? S[i] * R[i] / s.Rho[i] : 0;
          }
          s['massFraction_' + slot] = w;
        }
        if (errors === 'raise' && s.failures.length) {
          const f = s.failures[0], name1 = mode === 0 ? 'T' : 'H';
          const inputs = { [name1]: va[f.index], p: vp[f.index], X: vX[f.index] };
          const more = s.failures.length - 1;
          const where = (n > 1 ? `point ${f.index}: ` : '') + `${name1}=${inputs[name1]}, p=${inputs.p}, X=${inputs.X}`;
          throw new ERRORS[f.status](`${where}: ${f.message}` + (more ? ` [${more} further failure${more > 1 ? 's' : ''}]` : ''),
            f.status, f.index, inputs, f.message);
        }
        return s;
      } finally {
        ptrs.forEach(q => M._free(q));
      }
    }

    function model(backend) {
      const name = backend || 'IAPS84';
      if (!(name in BACKENDS)) throw new Error("backend must be 'IAPS84' or 'IAPWS95'");
      let handle = M._xt_new(BACKENDS[name]);
      if (!handle) throw new Error('xt_new: ' + str(M._xt_last_error()));
      const live = () => { if (!handle) throw new Error('model was freed'); return handle; };
      return {
        backend: name,
        tpx: (T, p, X, opts) => evaluate(live(), 0, T, p, X === undefined ? 0 : X, opts),
        hpx: (H, p, X, opts) => evaluate(live(), 1, H, p, X === undefined ? 0 : X, opts),
        free: () => { if (handle) { M._xt_free(handle); handle = 0; } },
      };
    }

    return { buildInfo, fields: FIELDS.slice(), statusNames: STATUS.slice(), phaseName, model };
  }

  return { load, ABI_VERSION, FIELDS, STATUS, EvaluationError, InvalidInputError, OutOfRangeError, NotConvergedError,
    NotAvailableError };
});
