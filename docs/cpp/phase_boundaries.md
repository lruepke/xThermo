# Phase boundaries

## Pressure of the V+L+H surface

`P_VLH` returns the pressure on the three-phase V+L+H surface for a given temperature.

```cpp
#include "H2ONaCl.h"
#include <iostream>

using namespace xThermal;

int main()
{
    H2ONaCl::cH2ONaCl eos("IAPS84");

    int N = 100;
    double T_min = H2ONaCl::T_MIN_VLH;
    double T_max = H2ONaCl::T_MAX_VLH;
    double dT = (T_max - T_min) / (N - 1);

    std::cout << "T(°C)    P_VLH(bar)\n";
    for (int i = 0; i < N; ++i) {
        double T = T_min + i * dT;   // K
        double P = eos.P_VLH(T);     // Pa
        std::cout << T - 273.15 << "    " << P / 1e5 << "\n";
    }
}
```

## Critical curve

```cpp
#include "H2ONaCl.h"
#include <iostream>

using namespace xThermal;

int main()
{
    H2ONaCl::cH2ONaCl eos("IAPS84");

    double T = 500 + 273.15;   // K
    double P_crit, X_crit;
    eos.P_X_Critical(T, P_crit, X_crit);   // P_crit [Pa], X_crit [kg/kg]
    std::cout << P_crit / 1e5 << " bar, " << X_crit << " kg/kg\n";
}
```

## Liquid and vapour compositions on the V+L+H surface

```cpp
#include "H2ONaCl.h"
#include <iostream>

using namespace xThermal;

int main()
{
    H2ONaCl::cH2ONaCl eos("IAPS84");

    double T = 600 + 273.15;   // K, within [T_MIN_VLH, T_MAX_VLH]
    double P = eos.P_VLH(T);
    double X_L, X_V;
    eos.X_VLH(T, P, X_L, X_V);   // X_L: liquid (halite liquidus), X_V: vapour [kg/kg]
    std::cout << P / 1e5 << " bar: X_L = " << X_L << ", X_V = " << X_V << "\n";
}
```
