# Phase region

`findPhaseRegion_TPX` returns the phase region index without computing all thermodynamic properties.

```cpp
#include "H2ONaCl.h"
#include <iostream>

using namespace xThermal;

int main()
{
    H2ONaCl::cH2ONaCl eos("IAPS84");

    double T = 400 + 273.15;   // K
    double p = 200e5;           // Pa  (200 bar)
    double X = 0.032;           // kg/kg  (3.2 wt% NaCl)

    H2ONaCl::PhaseRegion region = eos.findPhaseRegion_TPX(T, p, X);
    std::string name = eos.getPhaseRegionName(region);

    std::cout << "Phase region: " << name << " (index " << region << ")\n";
}
```

## Phase region names

| Name | Description |
|---|---|
| `Liquid` | Single liquid phase |
| `Vapor` | Single vapour phase |
| `V+L` | Vapour–liquid two-phase |
| `L+H` | Liquid–halite two-phase |
| `V+H` | Vapour–halite two-phase |
| `V+L+H` | Three-phase coexistence surface |

## Phase name from an UpdateState result

```cpp
using namespace xThermal;

H2ONaCl::cH2ONaCl eos("IAPS84");
ThermodynamicProperties props;
eos.UpdateState_TPX(props, T, p, X);
std::string name = eos.phase_name(props.phase);
```
