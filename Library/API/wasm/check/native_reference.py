"""W2a: native reference for the WebAssembly verification (Library/API/wasm/check/verify.mjs; review/W2A_WASM_MODULE.md).

Builds the input sets with the native Python module (the same commit, built in the same container) and records the
native results of every input: status, phase code and name, the library message of failed points, the phases present,
every field of the module (xthermo.js FIELDS) as float.hex, and for failed points the Python exception class of
errors="raise". Both water backends.

Sets (per backend):
  grid       T-P-X: 19 temperatures 1-1000 degC x 16 pressures 1-5000 bar x 15 salinities 0-1 (pure water and
             pure NaCl included): every phase region, and out-of-range corners
  boundary   T-P-X on both sides of every phase-region change along 63 isobars/isopleths (T scanned in 300 steps,
             the change bracketed by bisection to 1e-7 K), at 1e-6 K and 1e-3 K outside the bracket. Points closer
             than 1e-6 K would test the last bits of the platform, not the module.
  hpx        H-P-X: regular enthalpy grid at 10 pressures x 8 salinities, plus the enthalpies of the boundary states
             +- dH and +- 1 kJ/kg, dH = max(1 J/kg, 3e-5 / X J/kg) (100 x the dilute V+L enthalpy round-off floor
             of XT-13, docs/limitations.md, so that a point is not within platform noise of a boundary)
  dilute     the three dilute V+L cells of the XT-13 diagnosis (review/XT13_DILUTE_VL.md), T-P-X and H-P-X
  failures   typed failures: InvalidInput (NaN inputs, X < 0, X > 1, p <= 0), OutOfRange (T, p, H beyond the model
             range) and NotConverged (H-P-X enthalpy-check failures near the critical point and at X = 1e-8, 1 bar,
             failing on both macOS and Linux in the O-3 fixtures)
Also: model.phase_name for every PhaseRegion code.

Run: XTHERMAL_PYTHONPATH=<native python dir> python native_reference.py <out.json.gz>
"""
import gzip
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.environ["XTHERMAL_PYTHONPATH"])
import xthermo                       # noqa: E402
from xthermo import H2ONaCl          # noqa: E402

FIELDS = ["T", "p", "X", "H", "S_l", "S_v", "S_h", "X_l", "X_v", "Rho", "Rho_l", "Rho_v", "Rho_h", "H_l", "H_v", "H_h",
          "Cp", "Cp_l", "Cp_v", "Cp_h", "Mu", "Mu_l", "Mu_v", "IsothermalCompressibility", "IsobaricExpansivity"]
BACKENDS = ("IAPS84", "IAPWS95")
K = 273.15
NAN = float("nan")

GRID_T = [1, 25, 50, 100, 150, 200, 250, 300, 350, 374, 380, 400, 450, 500, 600, 700, 800, 900, 1000]
GRID_P = [1, 5, 10, 25, 50, 100, 150, 200, 220, 221, 250, 300, 500, 1000, 2000, 5000]
GRID_X = [0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 0.01, 0.032, 0.1, 0.2, 0.26, 0.4, 0.6, 0.8, 1.0]
BND_P = [5, 25, 100, 200, 220.5, 250, 300, 500, 1000]
BND_X = [1e-6, 1e-3, 0.01, 0.032, 0.1, 0.26, 0.5]
HPX_P = [1, 10, 50, 100, 200, 250, 300, 500, 1000, 2000]
HPX_X = [0, 1e-6, 1e-4, 0.01, 0.032, 0.1, 0.26, 0.5]
# XT-13 diagnosis cells (review/XT13_DILUTE_VL.md): p [Pa], X, T [K]
DILUTE = [(3863796.9168229131, 0.030552286206211313, 522.51375577066767),
          (3627821.540650486, 4.4125655421018963e-4, 517.82979115326157),
          (2825688.4637794737, 5.7759820365636842e-6, 503.74485869465201)]
# H-P-X NotConverged on macOS and Linux (O-3 fixtures; HPX-CHK): backend, H [J/kg], p [Pa], X
NOT_CONVERGED = [("IAPS84", 2141025.641025641, 22.0e6, 1e-6), ("IAPS84", 2141025.641025641, 22.0e6, 1e-4),
                 ("IAPS84", 2141025.641025641, 22.04e6, 1e-6), ("IAPS84", 489743.58974358975, 1e5, 1e-8),
                 ("IAPWS95", 2046666.666666667, 22.0e6, 1e-3), ("IAPWS95", 2046666.666666667, 22.04e6, 1e-6),
                 ("IAPWS95", 2046666.666666667, 22.04e6, 1e-4)]
INVALID = [(NAN, 1e7, 0.1), (600.0, NAN, 0.1), (600.0, 1e7, NAN), (600.0, 1e7, -0.1), (600.0, 1e7, 1.5),
           (600.0, 0.0, 0.1), (600.0, -1e5, 0.1)]
OUT_OF_RANGE_TPX = [(1300.0, 1e7, 0.1), (250.0, 1e7, 0.1), (600.0, 6e8, 0.1), (600.0, 10.0, 0.1)]
OUT_OF_RANGE_HPX = [(1e8, 1e7, 0.1), (-1e7, 1e7, 0.1), (2e6, 6e8, 0.1)]


def tpx(sw, T, p, X):
    return xthermo.tpx(sw, np.asarray(T, float), np.asarray(p, float), np.asarray(X, float), errors="nan")


def boundary_temperatures(sw):
    """T on both sides of every phase change along each (p, X), 1e-6 K and 1e-3 K outside a 1e-7 K bracket."""
    out = []
    T = np.linspace(K + 1, K + 1000, 300)
    for pb in BND_P:
        for x in BND_X:
            p = pb * 1e5
            s = tpx(sw, T, np.full(T.size, p), np.full(T.size, x))
            ph, st = np.asarray(s.phase), np.asarray(s.status)
            for i in range(T.size - 1):
                if st[i] or st[i + 1] or ph[i] == ph[i + 1]:
                    continue
                lo, hi, plo = T[i], T[i + 1], ph[i]
                while hi - lo > 1e-7:
                    mid = 0.5 * (lo + hi)
                    m = tpx(sw, [mid], [p], [x])
                    if int(m.status[0]) == 0 and int(m.phase[0]) == plo:
                        lo = mid
                    else:
                        hi = mid
                out += [(lo - 1e-3, p, x), (lo - 1e-6, p, x), (hi + 1e-6, p, x), (hi + 1e-3, p, x)]
    return out


def hpx_inputs(sw, backend, bnd):
    out = []
    for pb in HPX_P:
        for x in HPX_X:
            out += [(h, pb * 1e5, x) for h in np.linspace(1e5, 4e6, 20 if backend == "IAPS84" else 10)]
    # enthalpies of the boundary states, moved away from the boundary by dH and 1 kJ/kg
    T, p, X = (np.array(v) for v in zip(*bnd))
    s = tpx(sw, T, p, X)
    for k in range(0, T.size, 4):   # (lo - 1e-3, lo - 1e-6, hi + 1e-6, hi + 1e-3) per boundary
        if int(s.status[k + 1]) or int(s.status[k + 2]):
            continue
        dh = max(1.0, 3e-5 / X[k]) if X[k] > 0 else 1.0
        hl, hh = float(s.H[k + 1]), float(s.H[k + 2])
        lo, hi = min(hl, hh), max(hl, hh)
        out += [(lo - 1e3, p[k], X[k]), (lo - dh, p[k], X[k]), (hi + dh, p[k], X[k]), (hi + 1e3, p[k], X[k])]
    return out


def record(sw, mode, name, inputs):
    a, p, X = (np.array(v, dtype=float) for v in zip(*inputs))
    f = xthermo.tpx if mode == 0 else xthermo.hpx
    t0 = time.perf_counter()
    s = f(sw, a, p, X, errors="nan")
    dt = time.perf_counter() - t0
    msgs = {i[0]: m for i, st, m in s.failures}
    present = ["".join(c for c, ph in (("l", "liquid"), ("v", "vapour"), ("h", "halite")) if s.has(ph)[k])
               for k in range(a.size)]
    recs = []
    for k in range(a.size):
        r = {"set": name, "mode": mode, "in": [float(a[k]).hex(), float(p[k]).hex(), float(X[k]).hex()],
             "status": int(s.status[k]), "phase": int(s.phase[k]), "message": msgs.get(k, ""), "present": present[k],
             "f": [float(np.asarray(s[fl])[k]).hex() for fl in FIELDS]}
        if r["status"]:
            try:
                f(sw, float(a[k]), float(p[k]), float(X[k]))
                r["raise"] = None
            except xthermo.EvaluationError as e:
                r["raise"] = type(e).__name__
        recs.append(r)
    return recs, dt


def main(path):
    out = {"fields": FIELDS, "records": {}, "timing": {}, "phase_names": {}}
    sw0 = H2ONaCl.cH2ONaCl("IAPS84")
    out["phase_names"] = {str(c): sw0.phase_name(c) for c in range(-1, 15)}
    for b in BACKENDS:
        sw = H2ONaCl.cH2ONaCl(b)
        bnd = boundary_temperatures(sw)
        sets = [
            ("grid", 0, [(t + K, pb * 1e5, x) for t in GRID_T for pb in GRID_P for x in GRID_X]),
            ("boundary", 0, bnd),
            ("hpx", 1, hpx_inputs(sw, b, bnd)),
            ("dilute_tpx", 0, [(T, p, x) for p, x, T in DILUTE]),
            ("dilute_hpx", 1, [(float(tpx(sw, [T], [p], [x]).H[0]), p, x) for p, x, T in DILUTE]),
            ("failures_tpx", 0, INVALID + OUT_OF_RANGE_TPX),
            ("failures_hpx", 1, INVALID + OUT_OF_RANGE_HPX + [(h, p, x) for bb, h, p, x in NOT_CONVERGED if bb == b]),
        ]
        recs = []
        for name, mode, inputs in sets:
            r, dt = record(sw, mode, name, inputs)
            recs += r
            out["timing"][f"{b} {name}"] = [len(r), dt]
            print(f"{b:8s} {name:13s} n={len(r):6d} failed={sum(x['status'] != 0 for x in r):5d}  {dt:7.2f} s", flush=True)
        out["records"][b] = recs
    with gzip.open(path, "wt") as f:
        json.dump(out, f)


if __name__ == "__main__":
    main(sys.argv[1])
