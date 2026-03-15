# H-p-X flash

`UpdateState_HPX` performs an enthalpy flash: given specific enthalpy, pressure, and salinity it returns all thermodynamic properties including temperature.

## Single point

```python
from xThermal import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")

H = 1.5e6    # J/kg
p = 300e5    # Pa  (300 bar)
X = 0.05     # kg/kg  (5 wt% NaCl)

props = sw.UpdateState_HPX(H, p, X)

print(f"Temperature: {props.T - 273.15:.1f} °C")
print(f"Phase:       {sw.phase_name(props.phase)}")
print(f"Density:     {props.Rho:.1f} kg/m³")
```

## Vectorised (NumPy arrays)

All three arrays must have the **same length** and be **1-D**. Wrap output fields with `np.array()`.

```python
import numpy as np
from xThermal import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")

H = np.linspace(0.1e6, 5.0e6, 200)   # J/kg
p = np.full_like(H, 300e5)            # 300 bar
X = np.full_like(H, 0.1)              # 10 wt% NaCl

state = sw.UpdateState_HPX(H, p, X)

T_C = np.array(state.T) - 273.15
rho = np.array(state.Rho)
```

For 2-D grids, flatten with `.ravel()` and reshape outputs as shown in the [T-p-X page](tpx.md#2-d-grids).

!!! note
    The HPX flash is more computationally intensive than TPX because it involves an iterative solve. Use `sw.showProgressBar(True)` to display a progress indicator for large arrays.
