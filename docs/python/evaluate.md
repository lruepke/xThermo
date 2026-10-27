# Evaluating states

Two functions evaluate the equation of state on scalars or NumPy arrays:

```text
xthermo.tpx(model, T, p, X=0.0, *, errors="raise", absent="zero")   # temperature input
xthermo.hpx(model, H, p, X=0.0, *, errors="raise", absent="zero")   # enthalpy input
```

- `model` is any model object: `H2ONaCl.cH2ONaCl(...)`, `H2O.cIAPS84()`, `H2O.cIAPWS95()` or `NaCl.cNaCl(...)`.
- `T` in K, `p` in Pa, `H` in J/kg, `X` in kg/kg. The water models ignore `X`.
- The result is a `xthermo.States` object (below). `hpx` also returns the temperature, in `T`.
- `errors` and `absent` choose how failed points and absent phases are reported: see
  [Errors and status](errors.md) and [Absent phases](#absent-phases) below.

## Scalars, arrays and broadcasting

The inputs can be numbers, lists or NumPy arrays of any shape that
[broadcast](https://numpy.org/doc/stable/user/basics.broadcasting.html) together. Every result array has the
broadcast shape. A typical grid puts one input along each axis:

```python
import numpy as np
import xthermo
from xthermo import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")
T = np.linspace(100, 600, 6) + 273.15           # 6 temperatures, K
X = np.array([0.001, 0.032, 0.1, 0.3])          # 4 salinities, kg/kg

s = xthermo.tpx(sw, T[:, None], 300e5, X[None, :])
print(s.shape)                                  # (6, 4): temperature along axis 0, salinity along 1
print(s.Rho.round(1))                           # bulk density, kg/m³
print(s.phase_name[:, 1])                       # phases of 3.2 wt% NaCl from 100 to 600 °C
```

With scalar inputs every field is a single NumPy number:

```python
import xthermo
from xthermo import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")
s = xthermo.tpx(sw, 573.15, 300e5, 0.1)          # 300 °C, 300 bar, 10 wt% NaCl
print(s.phase_name, round(float(s.Rho), 1), "kg/m³", round(float(s.H) / 1e6, 4), "MJ/kg")
```

All points are evaluated in one call into the C++ library. With the IAPS-84 backend a brine state takes a few
microseconds; IAPWS-95 brines are slower (tenths of a millisecond), because the water equation of state is solved
iteratively.

## The States result

Fields are read with attribute or key access, `s.Rho` or `s["Rho"]`; `s.keys()` lists them. Each field is a
read-only NumPy array; use `s.Rho.copy()` for a writable one. Arrays stay valid after the model is deleted.

| Attribute | Content |
|---|---|
| `shape` | shape of the broadcast inputs |
| `fluid` | name of the model's fluid, e.g. `"H2O-NaCl"` |
| `phase` | phase-region code of each point (integer, table below) |
| `phase_name` | the name of each code, e.g. `"V+L"` (array of strings, or a string for scalar input) |
| `has(phase)` | mask: `"liquid"`, `"vapour"` (or `"vapor"`) or `"halite"` is present at an evaluated point |
| `status` | `xthermo.Status` of each point; `Status.Ok` for evaluated points ([Errors and status](errors.md)) |
| `ok` | mask of the evaluated points (`status == Status.Ok`) |
| `failures` | list of `(index, Status, message)` for the points that were not evaluated |

### Phase identification

`phase` holds the library's phase-region code. Its name is `s.phase_name`, or `model.phase_name(code)`.

| Code | Name | Phases present (`has`) |
|---|---|---|
| 0 | Liquid | liquid |
| 1 | Vapor | vapour |
| 2 | Solid | halite (pure NaCl) |
| 3, 4, 5 | Supercritical, Sup.crit. vapor, Sup.crit. liquid | one fluid phase, reported in the liquid fields (`_l`) |
| 6 | Critical point | as 3 |
| 7 | V+L(water) | liquid, vapour (pure water at its boiling curve) |
| 8 | L+H | liquid, halite |
| 9 | V+H | vapour, halite |
| 10 | V+L | liquid, vapour |
| 11 | V+L+H | liquid, vapour, halite |
| 12 | Unknown | none: the point was not evaluated |

The labels "Liquid" and "Vapor" for a single-phase H₂O–NaCl fluid follow the library's classification. Above the
critical curve of the mixture the fluid changes continuously between the two; there is no phase boundary there.

```python
import numpy as np
import xthermo
from xthermo import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")
T = np.linspace(300, 600, 7) + 273.15
s = xthermo.tpx(sw, T, 250e5, 0.032)             # seawater salinity, heated at 250 bar

two_phase = s.has("vapour") & s.has("liquid")    # where seawater has split into vapour and brine
print(T[two_phase] - 273.15)
print(s.X_l[two_phase], s.X_v[two_phase])        # brine and vapour salinity, kg/kg
```

### Properties

Each property comes as a bulk value and, where it applies, as phase values with the suffixes `_l` (liquid), `_v`
(vapour) and `_h` (halite). In a single-phase state the bulk value is the value of that phase.

| Field | Unit | Meaning |
|---|---|---|
| `T`, `p`, `X` | K, Pa, kg/kg | temperature, pressure, bulk salinity (`T` evaluated by `hpx`) |
| `H`, `H_l`, `H_v`, `H_h` | J/kg | specific enthalpy |
| `Rho`, `Rho_l`, `Rho_v`, `Rho_h` | kg/m³ | density |
| `Cp`, `Cp_l`, `Cp_v`, `Cp_h` | J/(kg K) | isobaric heat capacity (∂h/∂T)<sub>p,X</sub>; phases at fixed composition |
| `Mu`, `Mu_l`, `Mu_v` | Pa s | dynamic viscosity; NaN where the viscosity correlation does not evaluate |
| `S_l`, `S_v`, `S_h` | – | phase saturations (below) |
| `X_l`, `X_v` | kg/kg | salinity of the liquid and the vapour |
| `dRhodT`, `dRhodT_l`, `dRhodT_v` | kg/(m³ K) | ∂ρ/∂T at fixed p |
| `dRhodP`, `dRhodP_l`, `dRhodP_v` | kg/(m³ Pa) | ∂ρ/∂p at fixed T |
| `IsothermalCompressibility`, `IsothermalCompressibility_l`, `IsothermalCompressibility_v` | 1/Pa | (1/ρ) ∂ρ/∂p at fixed T |
| `IsobaricExpansivity`, `IsobaricExpansivity_l`, `IsobaricExpansivity_v` | 1/K | −(1/ρ) ∂ρ/∂T at fixed p |

**Multiphase states** (V+L, V+H, L+H, V+L+H): `Rho` and `H` are the bulk values of the phase assemblage, and the
saturations are **volume fractions** (the water models report the vapour **mass** fraction in their two-phase
states). The bulk `Cp` contains no latent heat, the bulk `Mu` is not a mixture viscosity (use `Mu_l` and `Mu_v`), and
some bulk derivatives are not computed. The definitions, and what the fields hold in each case, are on
[Model and limitations](../limitations.md#multiphase-properties). Convert saturations to mass fractions with
w<sub>i</sub> = S<sub>i</sub> ρ<sub>i</sub> / ρ. The pure NaCl model reports its properties in the bulk fields only.

```python
import xthermo
from xthermo import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")
s = xthermo.tpx(sw, 400 + 273.15, 250e5, 0.1)    # 400 °C, 250 bar, 10 wt% NaCl: vapour + liquid

w_v = s.S_v * s.Rho_v / s.Rho                     # vapour mass fraction from the volume fraction
print(s.phase_name, f"vapour: {float(s.S_v):.3f} of the volume, {float(w_v):.3f} of the mass")
print(f"liquid {float(s.X_l):.4f} kg/kg, {float(s.Rho_l):.1f} kg/m³")
print(f"vapour {float(s.X_v):.2e} kg/kg, {float(s.Rho_v):.1f} kg/m³")
```

## Absent phases

In the library's results the fields of a phase that is not present hold 0: `Rho_v` of a pure liquid, for example.
That is the default, `absent="zero"`. With `absent="nan"` these phase fields are NaN instead, which keeps absent
phases out of plots and statistics:

```python
import numpy as np
import xthermo
from xthermo import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")
T = np.array([300.0, 400.0, 500.0]) + 273.15
s = xthermo.tpx(sw, T, 250e5, 0.032, absent="nan")
print(s.phase_name)                               # Liquid, V+L, V+H
print(s.Rho_v)                                    # NaN where there is no vapour
print(s.Rho_l)                                    # NaN where there is no liquid
print(s.S_v)                                      # saturations are never changed: 0 where absent
```

`absent="nan"` affects only the phase fields (`_l`, `_v`, `_h`, except the saturations) of the `States` result. Bulk
fields and the library's own results are unchanged. Whether a phase is present follows the phase code; `has()` gives
the same mask.

## Enthalpy-based input

`hpx` inverts the enthalpy: it finds the state with the given specific enthalpy, pressure and salinity, and returns
it with its temperature.

```python
import numpy as np
import xthermo
from xthermo import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")
H = np.array([1.0e6, 2.0e6, 3.0e6])               # J/kg
s = xthermo.hpx(sw, H, 300e5, 0.05)               # 300 bar, 5 wt% NaCl
print(s.T - 273.15)                               # °C
print(s.phase_name)
```

The field `H` of the result is the enthalpy of the returned state. It agrees with the input to within the solver
tolerance, except at X = 1 (pure NaCl in the H₂O–NaCl model) inside the melting range, where the returned state can
have a different enthalpy ([Under investigation](../limitations.md#under-investigation)). For plots on an enthalpy
grid, use the input values as coordinates (see [Property maps](property_maps.md#enthalpy-grids)). The pure NaCl model has no enthalpy-based input (`hpx` raises
`NotAvailableError`).
