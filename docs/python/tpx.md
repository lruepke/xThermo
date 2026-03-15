# T-p-X calculation

`UpdateState_TPX` computes all thermodynamic properties for a given temperature, pressure, and salinity.

## Single point

```python
from xThermal import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")

T = 400 + 273.15   # K
p = 250e5          # Pa  (250 bar)
X = 0.1            # kg/kg  (10 wt% NaCl)

props = sw.UpdateState_TPX(T, p, X)

print(f"Phase:    {sw.phase_name(props.phase)}")
print(f"Density:  {props.Rho:.1f} kg/m³")
print(f"Enthalpy: {props.H/1e6:.4f} MJ/kg")
print(f"Cp:       {props.Cp:.1f} J/kg/K")
print(f"Viscosity (liquid):{props.Mu_l:.3e} Pa·s")
print(f"Viscosity (vapor):{props.Mu_v:.3e} Pa·s")
```

## Vectorised (NumPy arrays)

`UpdateState_TPX` accepts 1-D NumPy arrays and evaluates each point `i` independently (T[i], p[i], X[i]).
All three arrays **must have the same length**. Output fields are SWIG vectors — wrap with `np.array()` to use them as NumPy arrays.

```python
import numpy as np
from xThermal import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")

T = np.linspace(100 + 273.15, 600 + 273.15, 200)
p = np.full_like(T, 300e5)   # 300 bar
X = np.full_like(T, 0.1)     # 10 wt% NaCl

state = sw.UpdateState_TPX(T, p, X)

rho   = np.array(state.Rho)
phase = np.array(state.phase)
```

## 2-D grids

For a 2-D parameter grid, flatten the arrays with `.ravel()` before passing them in, then reshape the outputs:

```python
import numpy as np
from xThermal import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")

T = np.linspace(100 + 273.15, 600 + 273.15, 100)
X = np.linspace(0.0, 0.6, 80)
TT, XX = np.meshgrid(T, X)
PP = np.full_like(TT, 300e5)   # 300 bar

state = sw.UpdateState_TPX(TT.ravel(), PP.ravel(), XX.ravel())

Phase = np.array(state.phase).reshape(TT.shape)
Rho   = np.array(state.Rho).reshape(TT.shape)
```

!!! tip
    Phase names can be retrieved with `sw.phase_name(int(phase_index))`, e.g. `"Liquid"`, `"V+L"`, `"Vapor"`, `"L+H"`.
