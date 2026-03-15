# T-p-X calculation

`UpdateState_TPX` computes all thermodynamic properties for a given temperature, pressure, and salinity.

```cpp
#include "H2ONaCl.h"
#include <iostream>

using namespace xThermal;

int main()
{
    H2ONaCl::cH2ONaCl eos("IAPS84");   // backend: "IAPS84" or "IAPWS95"

    double T = 400 + 273.15;   // K
    double p = 250e5;           // Pa  (250 bar)
    double X = 0.1;             // kg/kg  (10 wt% NaCl)

    ThermodynamicProperties props;
    eos.UpdateState_TPX(props, T, p, X);

    std::cout << "Phase:    " << eos.phase_name(props.phase) << "\n"
              << "Density:  " << props.Rho                   << " kg/m³\n"
              << "Enthalpy: " << props.H / 1e6               << " MJ/kg\n"
              << "Cp:       " << props.Cp                    << " J/kg/K\n"
              << "Viscosity (liquid): " << props.Mu_l        << " Pa·s\n"
              << "Viscosity (vapor):  " << props.Mu_v        << " Pa·s\n";
}
```

## Key fields of ThermodynamicProperties

| Field | Description |
|---|---|
| `Rho` | Bulk density (kg/m³) |
| `H` | Specific enthalpy (J/kg) |
| `Cp` | Isobaric heat capacity (J/kg/K) |
| `Mu` | Bulk dynamic viscosity (Pa·s) |
| `Mu_l` / `Mu_v` | Liquid / vapour phase viscosity (Pa·s) |
| `Rho_l` / `Rho_v` | Liquid / vapour phase density (kg/m³) |
| `T` | Temperature (K) — especially useful from HPX flash |
| `phase` | Phase region index |

## Enthalpy flash

```cpp
using namespace xThermal;

H2ONaCl::cH2ONaCl eos("IAPS84");
ThermodynamicProperties props;
eos.UpdateState_HPX(props, H, p, X);   // H [J/kg], p [Pa], X [kg/kg]
std::cout << "T = " << props.T - 273.15 << " °C\n";
```
