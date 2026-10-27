# Evaluating states

A state is given by temperature, pressure and bulk salinity (T–p–X), or by specific enthalpy, pressure and salinity
(H–p–X). Every model offers the same functions:

| Function | Input | On failure |
|---|---|---|
| `UpdateState_TPX(props, T, p, X)`, `UpdateState_HPX(props, H, p, X)` | one state; fills `props` | throws (see [Errors and threads](errors.md)) |
| `UpdateState_TPX(T, p, X)`, `UpdateState_HPX(H, p, X)` | one state; returns the properties | throws |
| `tryUpdateState_TPX(props, T, p, X, &message)`, `tryUpdateState_HPX(...)` | one state | returns an `EvaluationStatus` |
| `evaluate_TPX(T, p, X)`, `evaluate_HPX(H, p, X)` | vectors of points | records a status per point |
| `UpdateState_TPX(T, p, X)`, `UpdateState_HPX(H, p, X)` with vectors | vectors of points | throws at the first failing point |

For pure water and pure NaCl, X is ignored (it may be omitted).

## One state

```cpp
#include "H2ONaCl.h"
#include <iostream>

using namespace xThermal;

int main()
{
    H2ONaCl::cH2ONaCl brine("IAPS84");

    ThermodynamicProperties props;
    // T [K], p [Pa], X [kg/kg]: seawater at 400 °C, 250 bar
    brine.UpdateState_TPX(props, 400 + 273.15, 250e5, 0.032);

    std::cout << brine.phase_name(props.phase) << "\n"          // V+L
              << "bulk density " << props.Rho << " kg/m3\n"     // 194.068
              << "liquid: " << props.Rho_l << " kg/m3, X_l = " << props.X_l << "\n"
              << "vapour: " << props.Rho_v << " kg/m3, X_v = " << props.X_v << "\n"
              << "vapour volume fraction " << props.S_v << "\n";

    // the same state from its enthalpy: T is a result
    ThermodynamicProperties same = brine.UpdateState_HPX(props.H, 250e5, 0.032);
    std::cout << "T = " << same.T - 273.15 << " °C\n";          // 400
}
```

Each call resets the result object, so no value of a previous call survives in `props`.

## Without exceptions

`tryUpdateState_TPX` and `tryUpdateState_HPX` evaluate the same state but return its outcome instead of throwing. On
failure the inputs are kept, every other number is NaN, the phase is `Unknown`, and the optional message receives the
library's explanation:

```cpp
#include "H2ONaCl.h"
#include <iostream>
#include <string>

using namespace xThermal;

int main()
{
    H2ONaCl::cH2ONaCl brine("IAPS84");
    ThermodynamicProperties props;
    std::string message;

    EvaluationStatus status = brine.tryUpdateState_TPX(props, 1100 + 273.15, 300e5, 0.1, &message);
    if (status != EvaluationStatus::Ok)
        std::cout << evaluation_status_name(status) << "\n";   // OutOfRange: above 1000 °C
    std::cout << props.Rho << "\n";                              // nan
}
```

The statuses are listed under [Errors and threads](errors.md#status-codes).

## Many states

`evaluate_TPX` and `evaluate_HPX` take vectors of equal length (point i is `T[i], p[i], X[i]`) and never stop at a
failing point. The result holds one vector per field, `status[i]` (0 = ok) and the messages of the failed points in
`failed_message`, in index order:

```cpp
#include "H2ONaCl.h"
#include <iostream>
#include <vector>

using namespace xThermal;

int main()
{
    H2ONaCl::cH2ONaCl brine("IAPS84");
    // K; the second point is above the model's range
    std::vector<double> T = {300 + 273.15, 2000, 400 + 273.15};
    std::vector<double> p = {300e5, 300e5, 250e5};                 // Pa
    std::vector<double> X = {0.1, 0.1, 0.032};                     // kg/kg

    ThermodynamicPropertiesVector r = brine.evaluate_TPX(T, p, X);
    for (size_t i = 0; i < T.size(); i++)
        std::cout << i << ": status " << r.status[i] << ", Rho " << r.Rho[i] << "\n";
    // statuses 0, 2 (OutOfRange; Rho NaN), 0
    std::cout << r.failed_message.size() << " failed point(s)\n";
}
```

The vector forms of `UpdateState_TPX`/`UpdateState_HPX` evaluate the same points but throw at the first failure. With
`isMeshGrid = true` they evaluate every combination of the three vectors.

## Result fields

`ThermodynamicProperties` (and the vectors of `ThermodynamicPropertiesVector`):

| Field | Meaning |
|---|---|
| `T`, `p`, `X`, `H` | temperature [K], pressure [Pa], bulk salinity [kg/kg], bulk specific enthalpy [J/kg]; T is a result of H–p–X, H of T–p–X |
| `phase` | phase region, a `PhaseRegion`; `phase_name(phase)` gives its name ([Phase region](phase_region.md)) |
| `Rho`, `H`, `Cp`, `Mu` | bulk density [kg/m³], enthalpy, heat capacity [J/(kg K)], viscosity [Pa s] |
| `S_l`, `S_v`, `S_h` | phase fractions: **volume** fractions for H₂O–NaCl; the vapour **mass** fraction (quality) for the water models |
| `X_l`, `X_v` | salinity of the liquid and the vapour [kg/kg] |
| `Rho_l`, `Rho_v`, `Rho_h`; `H_l`, `H_v`, `H_h`; `Cp_l`, `Cp_v`, `Cp_h`; `Mu_l`, `Mu_v` | properties of the liquid, vapour and halite phases |
| `dRhodT`, `dRhodP` | ∂ρ/∂T [kg/(m³ K)], ∂ρ/∂p [kg/(m³ Pa)]; with `_l`, `_v` for the phases |
| `IsothermalCompressibility`, `IsobaricExpansivity` | (1/ρ) ∂ρ/∂p [1/Pa], −(1/ρ) ∂ρ/∂T [1/K]; with `_l`, `_v` for the phases |

In a single-phase state, the bulk fields are those of the phase; fields of an absent phase are 0.

**Multiphase states**: `Rho` and `H` are the bulk values of the phase assemblage. The bulk `Cp` contains no latent
heat, the bulk `Mu` is not a mixture viscosity (use `Mu_l` and `Mu_v`), and some bulk derivatives are not computed.
The definitions, and what each field holds in each case, are on
[Model and limitations](../limitations.md#multiphase-properties).

**Unavailable properties.** A property that cannot be evaluated at an otherwise valid state is NaN while the others
are valid: for example the IAPS-84 viscosity outside its validity range (above 3500 bar at more than 150 °C, among
others), and brine viscosities that depend on it.
