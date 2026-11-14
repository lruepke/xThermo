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

    PhaseRegion region = eos.findPhaseRegion_TPX(T, p, X);
    std::string name = eos.phase_name(region);

    std::cout << "Phase region: " << name << " (index " << region << ")\n";
}
```

## Phase region names

Names returned by `phase_name` for H₂O–NaCl states:

| Name | Enum value | Description |
|---|---|---|
| `Liquid` | `SinglePhase_L` | Single liquid phase. Single-phase "liquid-like" and "vapour-like" fluid above the critical curve are also reported as `Liquid`. |
| `Vapor` | `SinglePhase_V` | Single vapour phase |
| `V+L` | `TwoPhase_VL` | Vapour–liquid coexistence |
| `L+H` | `TwoPhase_LH` | Liquid–halite coexistence |
| `V+H` | `TwoPhase_VH` | Vapour–halite coexistence |
| `V+L+H` | `ThreePhase_VLH` | Three-phase coexistence (reached through the H–p–X flash) |

For pure water, the names of the water backend are used instead, for example `Supercritical` and
`Sup.crit. liquid`, and `V+L(water)` in the two-phase region. This applies to the pure-water classes and also to
`cH2ONaCl` at X = 0, which evaluates pure water directly.

## Phase name from an UpdateState result

```cpp
#include "H2ONaCl.h"
#include <iostream>

using namespace xThermal;

int main()
{
    H2ONaCl::cH2ONaCl eos("IAPS84");
    ThermodynamicProperties props;
    eos.UpdateState_TPX(props, 400 + 273.15, 250e5, 0.1);
    std::cout << eos.phase_name(props.phase) << "\n";   // V+L
}
```
