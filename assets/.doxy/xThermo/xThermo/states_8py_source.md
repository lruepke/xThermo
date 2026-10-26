

# File states.py

[**File List**](files.md) **>** [**API**](dir_4b630c9df1aec2cd2fe9acaf75644364.md) **>** [**python**](dir_a5d5f7c295cec99c8e29471acd4df0e8.md) **>** [**xthermo**](dir_c7c9957b02953ad15ede12b0595b5fff.md) **>** [**states.py**](states_8py.md)

[Go to the documentation of this file](states_8py.md)


```Python
"""NumPy interface for state evaluations: xthermo.tpx and xthermo.hpx.

    import numpy as np, xthermo
    from xthermo import H2ONaCl
    sw = H2ONaCl.cH2ONaCl("IAPS84")
    s = xthermo.tpx(sw, T=np.linspace(300, 900, 200)[:, None] + 273.15, p=300e5, X=np.logspace(-4, 0, 150))
    s.Rho.shape          # (200, 150): the broadcast shape of the inputs
    s.phase_name, s.has("vapour"), s.ok, s.status, s.failures

The inputs broadcast like NumPy arrays (scalars give NumPy scalars). Each point is evaluated by the library's batch
evaluation (cxThermal::evaluate_TPX / evaluate_HPX, i.e. tryUpdateState_* per point), so the values are those of the
scalar and vector API. SI units throughout.
"""
import enum

import numpy as np

__all__ = ["tpx", "hpx", "States", "Status", "EvaluationError", "InvalidInputError", "OutOfRangeError",
           "NotConvergedError", "NotAvailableError"]


class Status(enum.IntEnum):
    """Outcome of the evaluation of a point (xThermal::EvaluationStatus, Exception.h)."""
    Ok = 0
    InvalidInput = 1     # non-finite input, X < 0 or p <= 0
    OutOfRange = 2       # outside the range in which the model evaluates
    NotConverged = 3     # a numerical solution failed
    NotAvailable = 4     # not implemented or not available for this model or state
    Failed = 5           # any other error of the library


class EvaluationError(Exception):
    """A point could not be evaluated (errors="raise"). Attributes: index (tuple, () for scalar inputs), status
    (Status), inputs (dict of the point's inputs), message (the library's message)."""

    def __init__(self, text, index, status, inputs, message):
        super().__init__(text)
        self.index, self.status, self.inputs, self.message = index, status, inputs, message


class InvalidInputError(EvaluationError, ValueError):
    pass


class OutOfRangeError(EvaluationError, ValueError):
    pass


class NotConvergedError(EvaluationError, RuntimeError):
    pass


class NotAvailableError(EvaluationError, ValueError):
    pass


_ERRORS = {Status.InvalidInput: InvalidInputError, Status.OutOfRange: OutOfRangeError,
           Status.NotConverged: NotConvergedError, Status.NotAvailable: NotAvailableError, Status.Failed: EvaluationError}

FIELDS = ("T", "p", "X", "H", "S_l", "S_v", "S_h", "X_l", "X_v", "Rho_l", "Rho_v", "Rho_h", "H_l", "H_v", "H_h", "Cp_l",
          "Cp_v", "Cp_h", "Mu_l", "Mu_v", "Rho", "Cp", "Mu", "dRhodP", "dRhodP_l", "dRhodP_v", "dRhodT", "dRhodT_l",
          "dRhodT_v", "IsothermalCompressibility", "IsothermalCompressibility_l", "IsothermalCompressibility_v",
          "IsobaricExpansivity", "IsobaricExpansivity_l", "IsobaricExpansivity_v")

# phases present per PhaseRegion code (DataStructures.h): l = liquid, v = vapour, h = halite. Single-phase fluid states
# of water above the critical point are reported in the liquid fields (library convention). Unknown, MixPhaseRegion
# and NotImposed: none.
_PRESENT = {0: "l", 1: "v", 2: "h", 3: "l", 4: "l", 5: "l", 6: "l", 7: "lv", 8: "lh", 9: "vh", 10: "lv", 11: "lvh"}
_NCODES = 14
_PRESENT_LUT = {slot: np.array([slot in _PRESENT.get(code, "") for code in range(_NCODES)]) for slot in "lvh"}
_SLOT = {"liquid": "l", "vapour": "v", "vapor": "v", "halite": "h"}
# fields of a phase (set to NaN for an absent phase with absent="nan"); saturations stay (0 is their value)
_PHASE_FIELDS = {f: f[-1] for f in FIELDS if f[-2:] in ("_l", "_v", "_h") and not f.startswith("S_")}


class States:
    """States with the shape of the broadcast inputs (attribute access s.Rho or mapping access s["Rho"]).

    Fields: the members of ThermodynamicProperties (SI units), as NumPy arrays converted on first access, with the
    library's values and conventions (DataStructures.h). Also:
      phase       PhaseRegion code per point; phase_name its name (model.phase_name)
      status      Status per point (an int array; a Status for scalar inputs); ok = status == Status.Ok
      failures    [(index, Status, message)] for the failed points
      has(phase)  "liquid", "vapour" or "halite" present at an evaluated point
    A failed point has NaN numbers (its inputs are kept) and phase Unknown. With absent="nan", the fields of a phase
    that is absent at a point (Rho_v, Mu_v, X_v, ... without vapour) are NaN instead of the library's 0; bulk fields
    and saturations are never changed. The pure NaCl model reports its properties in the bulk fields only.
    """

    def __init__(self, raw, model, shape, absent):
        self._raw, self._model, self.shape, self._absent, self._cache = raw, model, shape, absent, {}
        self.fluid = raw.fluidName

    def _array(self, name, dtype):   # the library's values, read-only, cached
        key = ("raw", name)
        if key not in self._cache:
            a = np.frombuffer(self._raw._field_bytes(name), dtype=dtype).reshape(self.shape)
            a.flags.writeable = False
            self._cache[key] = a
        return self._cache[key]

    def _get(self, name):
        cache = self._cache
        if name not in cache:
            if name in ("phase", "status"):
                a = self._array(name, np.intc)
            elif name in FIELDS:
                a = self._array(name, np.float64)
                if self._absent == "nan" and name in _PHASE_FIELDS:
                    a = np.where(self._present(_PHASE_FIELDS[name]), a, np.nan)
                    a.flags.writeable = False
            else:
                raise KeyError(name)
            cache[name] = a
        a = cache[name]
        if self.shape:
            return a
        return Status(int(a)) if name == "status" else a[()]

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        try:
            return self._get(name)
        except KeyError:
            raise AttributeError(f"'States' object has no attribute {name!r}") from None

    def __getitem__(self, name):
        return self._get(name)

    def __dir__(self):
        return sorted(set(super().__dir__()) | set(self.keys()))

    def __repr__(self):
        return f"<xthermo.States {self.fluid} shape={self.shape} failed={len(self._raw.failed_message)}>"

    def keys(self):
        return list(FIELDS) + ["phase", "status"]

    @property
    def ok(self):
        return self._array("status", np.intc) == 0 if self.shape else self.status == Status.Ok

    def _present(self, slot):
        phase = np.clip(self._array("phase", np.intc), 0, _NCODES - 1)
        return _PRESENT_LUT[slot][phase] & (self._array("status", np.intc) == 0)

    def has(self, phase):
        """Mask: the phase ("liquid", "vapour"/"vapor", "halite") is present at an evaluated point."""
        if phase not in _SLOT:
            raise ValueError(f"phase must be one of {sorted(_SLOT)}")
        m = self._present(_SLOT[phase])
        return m if self.shape else bool(m)

    @property
    def phase_name(self):
        phase = self._array("phase", np.intc)
        codes = np.unique(phase)
        names = np.empty(phase.shape, dtype=object)
        for c in codes:
            names[phase == c] = self._model.phase_name(int(c))
        return names if self.shape else names[()]

    @property
    def failures(self):
        status = self._array("status", np.intc).ravel()
        flat = np.flatnonzero(status)
        index = [tuple(int(j) for j in np.unravel_index(k, self.shape)) for k in flat]
        return [(i, Status(int(status[k])), m) for i, k, m in zip(index, flat, self._raw.failed_message)]


def _evaluate(model, name1, batch, a, p, X, errors, absent):
    if errors not in ("raise", "nan"):
        raise ValueError("errors must be 'raise' or 'nan'")
    if absent not in ("zero", "nan"):
        raise ValueError("absent must be 'zero' or 'nan'")
    a, p, X = np.broadcast_arrays(np.asarray(a, dtype=np.float64), np.asarray(p, dtype=np.float64),
                                  np.asarray(X, dtype=np.float64))
    flat = [v.ravel().tolist() for v in (a, p, X)]
    raw = getattr(model, batch)(*flat)
    s = States(raw, model, a.shape, absent)
    if errors == "raise" and len(raw.failed_message):
        status = s._array("status", np.intc).ravel()
        k = int(np.flatnonzero(status)[0])
        code = Status(int(status[k]))
        index = tuple(int(j) for j in np.unravel_index(k, a.shape))
        inputs = {name1: flat[0][k], "p": flat[1][k], "X": flat[2][k]}
        more = len(raw.failed_message) - 1
        where = (f"point {index}: " if a.shape else "") + ", ".join(f"{n}={v!r}" for n, v in inputs.items())
        text = f"{where}: {raw.failed_message[0]}" + (f" [{more} further failure{'s' if more > 1 else ''}]" if more else "")
        raise _ERRORS[code](text, index, code, inputs, raw.failed_message[0])
    return s


def tpx(model, T, p, X=0.0, *, errors="raise", absent="zero"):
    """States at temperature T [K], pressure p [Pa] and bulk salinity X [kg/kg] (ignored by the water models).

    T, p, X: scalars or arrays that broadcast together. errors="raise" (default) raises an EvaluationError subclass for
    the first point that cannot be evaluated; errors="nan" returns failed points as NaN with their status.
    absent="nan": fields of absent phases are NaN in the result (default "zero": the library's values).
    """
    return _evaluate(model, "T", "evaluate_TPX", T, p, X, errors, absent)


def hpx(model, H, p, X=0.0, *, errors="raise", absent="zero"):
    """States at specific enthalpy H [J/kg], pressure p [Pa] and bulk salinity X [kg/kg]; as tpx."""
    return _evaluate(model, "H", "evaluate_HPX", H, p, X, errors, absent)
```


