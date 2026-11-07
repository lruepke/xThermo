# Python guide

The Python package `xthermo` gives you the H₂O–NaCl equation of state of Driesner & Heinrich (2007) and Driesner
(2007), pure water (IAPS-84, IAPWS-95) and pure NaCl. You evaluate states on NumPy arrays of any shape and draw phase
diagrams and property maps with a few Matplotlib helpers.

```python
import numpy as np
import xthermo
from xthermo import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")                         # H2O-NaCl with the IAPS-84 water backend

T = np.array([300.0, 400.0, 500.0]) + 273.15            # K
s = xthermo.tpx(sw, T, 250e5, 0.032)                    # p = 250 bar in Pa, X = 3.2 wt% NaCl
print(s.phase_name)                                     # ['Liquid' 'V+L' 'V+H']
print(s.Rho)                                            # bulk density, kg/m³
print(s.X_v)                                            # vapour salinity, kg/kg (0 where absent)
```

## The models

| Module | Class | System |
|---|---|---|
| `xthermo.H2ONaCl` | `cH2ONaCl("IAPS84")`, `cH2ONaCl("IAPWS95")` | H₂O–NaCl; the argument selects the water equation of state |
| `xthermo.H2O` | `cIAPS84()`, `cIAPWS95()` | pure water |
| `xthermo.NaCl` | `cNaCl("IAPS84")` | pure NaCl (halite and liquid) |

The Driesner correlations were fitted with IAPS-84, so **`"IAPS84"` is the reference backend** for H₂O–NaCl.
`"IAPWS95"` gives slightly different brine properties; see
[Model and limitations](../limitations.md#water-backends).

Create a model once and reuse it for all your evaluations. A model object is not meant to be shared between threads;
use one per thread.

## Where to go next

| Page | Content |
|---|---|
| [Evaluating states](evaluate.md) | `xthermo.tpx` and `xthermo.hpx`, broadcasting, the `States` result, phases and the meaning of each property in multiphase states |
| [Errors and status](errors.md) | `errors="raise"` and `errors="nan"`, `xthermo.Status`, the exception classes |
| [Phase diagrams](phase_diagrams.md) | isobaric and isothermal sections from the library's phase-region geometry |
| [Property maps](property_maps.md) | maps and contours of any property, without interpolation across phase boundaries |
| [Phase boundaries](phase_boundaries.md) | the V+L+H surface, the critical curve and the compositions of coexisting phases |
| [Pure water and NaCl](water_nacl.md) | saturation of water, the water models, pure NaCl |
| [Classic SWIG interface](classic.md) | `UpdateState_TPX` / `UpdateState_HPX` on single points and 1-D vectors; the former name `xThermal` |

## Units

Inputs and outputs are in SI units everywhere:

| Quantity | Unit |
|---|---|
| temperature `T` | K |
| pressure `p` | Pa (1 bar = 10⁵ Pa) |
| salinity `X` | kg NaCl per kg of fluid, a mass fraction (0.032 = 3.2 wt% NaCl) |
| density | kg/m³ |
| specific enthalpy `H` | J/kg |
| heat capacity | J/(kg K) |
| viscosity | Pa s |

The plotting helpers label their axes in °C, bar, wt% NaCl and MJ/kg; their inputs stay in SI.

## Installation

`pip install .` in a clone of the repository installs `xthermo` with NumPy; `pip install ".[plotting]"` adds
Matplotlib. See [Installation](../installation.md). `import xthermo` loads only the compiled modules; NumPy and
Matplotlib are imported when you first use `tpx` / `hpx` or a plotting helper.
