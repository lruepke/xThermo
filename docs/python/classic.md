# Classic SWIG interface

The model classes also expose the library's C++ methods directly, through SWIG. `xthermo.tpx` and `xthermo.hpx` are
built on them and are usually more convenient. The classic methods remain useful for single points in a loop (they
avoid the small overhead of creating a `States` object) and for existing scripts.

## Single points

`UpdateState_TPX(T, p, X)` and `UpdateState_HPX(H, p, X)` return one `ThermodynamicProperties` object. Its fields
are the same as those of `States` (see the [field table](evaluate.md#properties)), as plain Python numbers. The water
and NaCl models take `UpdateState_TPX(T, p)`.

```python
from xthermo import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")

props = sw.UpdateState_TPX(300 + 273.15, 300e5, 0.1)       # 300 °C, 300 bar, 10 wt% NaCl
print(sw.phase_name(props.phase), f"{props.Rho:.1f} kg/m³", f"{props.H/1e6:.4f} MJ/kg")

props = sw.UpdateState_HPX(1.5e6, 300e5, 0.05)            # H in J/kg
print(f"{props.T - 273.15:.1f} °C", sw.phase_name(props.phase))
```

Errors raise Python exceptions: `ValueError` for invalid or out-of-range input, `RuntimeError` for a failed solution,
`NotImplementedError` for an unknown backend name. These are the base classes of the more specific exceptions that
`tpx` and `hpx` raise ([Errors and status](errors.md)).

## 1-D vectors

The same methods take three 1-D lists or arrays of **equal length** and return an object whose fields are SWIG
vectors; convert them with `np.array(...)`. There is no broadcasting: pass every input at full length. If any point
fails, the whole call raises and returns nothing.

```python
import numpy as np
from xthermo import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")

T = np.linspace(100 + 273.15, 600 + 273.15, 200)
p = np.full_like(T, 300e5)
X = np.full_like(T, 0.1)

state = sw.UpdateState_TPX(T, p, X)
rho = np.array(state.Rho)
phase = np.array(state.phase)
```

The equivalent with the NumPy interface is `xthermo.tpx(sw, T, 300e5, 0.1)`, which broadcasts the scalars, returns
NumPy arrays and offers `errors="nan"`. Its values are identical to the classic calls.

The vector form of `UpdateState_TPX` is not available on the pure-water classes in Python; use `xthermo.tpx` there.

## Other methods

Every model object also has, among others:

| Method | Returns |
|---|---|
| `phase_name(code)` | the name of a phase-region code |
| `Tmin()`, `Tmax()`, `pmin()`, `pmax()` | the range the model declares (K, Pa); see [Model and limitations](../limitations.md#domain-input-checks-and-errors) for what is checked |
| `T_critical()`, `p_critical()` | the critical point of water (water models) |
| `name()` | the model's name |

The H₂O–NaCl model adds the phase-boundary functions ([Phase boundaries](phase_boundaries.md)) and the sections
`Slice_constP` / `Slice_constT` ([Phase diagrams](phase_diagrams.md)).

## The former name xThermal

The package used to be imported as `xThermal`. That name still works and gives the same modules:

```python
import xthermo
from xThermal import H2ONaCl          # the same module as xthermo.H2ONaCl

print(H2ONaCl is xthermo.H2ONaCl)    # True
```

New code should use `import xthermo`.
