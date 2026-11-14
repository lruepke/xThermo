# C++ API overview

xThermo is a C++ library; the Python package and the command-line tool are built on the same classes. This guide
covers the stable public interface: constructing a model, evaluating states, handling errors, and using the library
from several threads. The model itself, its verification and its limitations are described on
[Model and limitations](../limitations.md).

## Installing the library

Build and install xThermo with CMake (details and options on [Installation](../installation.md)):

```bash
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=install
cmake --build build -j
cmake --install build
```

The installation (here `/path/to/xThermo/install`) contains the headers in `include/xThermo/`, the library
`libxThermal` in `lib/` (`lib64` on some Linux distributions), the command-line tool `xThermal` in `bin/`, and a CMake
package in `lib/cmake/xThermo/`.

## Using xThermo in a CMake project

`find_package(xThermo)` provides the target `xThermo::xThermal`, which carries the include directory and the GSL
dependency. A complete project has two files:

```cmake title="CMakeLists.txt"
cmake_minimum_required(VERSION 3.18)
project(my_app LANGUAGES CXX)

find_package(xThermo REQUIRED)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE xThermo::xThermal)
```

```cpp title="main.cpp"
#include "H2ONaCl.h"
#include <iostream>

int main()
{
    // H2O-NaCl with the IAPS-84 water model
    xThermal::H2ONaCl::cH2ONaCl brine("IAPS84");
    // T [K], p [Pa], X [kg/kg]
    xThermal::ThermodynamicProperties s = brine.UpdateState_TPX(300 + 273.15, 300e5, 0.1);
    std::cout << brine.phase_name(s.phase) << ", " << s.Rho << " kg/m3\n";   // Liquid, 849.034 kg/m3
}
```

Configure with the installation prefix in `CMAKE_PREFIX_PATH`, then build and run:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/xThermo/install
cmake --build build
./build/my_app
```

The test suite builds this project against a fresh installation.

??? note "Compiling without CMake"
    The library needs the GSL. With a system GSL (`gsl-config` in the PATH):

    ```bash
    g++ -std=c++11 -I/path/to/xThermo/install/include/xThermo $(gsl-config --cflags) \
        -o my_app main.cpp -L/path/to/xThermo/install/lib -lxThermal $(gsl-config --libs)
    ```

    If xThermo was built with `-DXTHERMO_BUNDLED_GSL=ON`, the GSL is installed with it: use
    `-I/path/to/xThermo/install/lib/xThermo/gsl/include` and
    `-L/path/to/xThermo/install/lib/xThermo/gsl/lib -lgsl -lgslcblas` instead of
    `gsl-config`. CMake does this for you.

## Models

| Class | Model | Construction |
|---|---|---|
| `xThermal::H2ONaCl::cH2ONaCl` | H₂O–NaCl (Driesner & Heinrich 2007; Driesner 2007) | `cH2ONaCl brine("IAPS84")` or `("IAPWS95")`: the water model used inside |
| `xThermal::PROST::cIAPS84` | pure water, IAPS-84 | `cIAPS84 water;` |
| `xThermal::IAPWS95::cIAPWS95` | pure water, IAPWS-95 | `cIAPWS95 water;` |
| `xThermal::NaCl::cNaCl` | pure NaCl | `cNaCl salt("IAPS84");` |

Including `H2ONaCl.h` makes all four available. All derive from `xThermal::cxThermal`, so a function can take any
model as `cxThermal&`. An unknown water-model name throws `NotImplementedError` at construction. `"IAPS84"` is the
reference: the H₂O–NaCl correlations were fitted with it (see [Water backends](../limitations.md#water-backends)).

## Units

All functions use SI units, with temperature in kelvin:

| Quantity | Unit |
|---|---|
| temperature T | K |
| pressure p | Pa |
| salinity X | kg NaCl / kg fluid (mass fraction, 0–1) |
| specific enthalpy H | J/kg |
| density | kg/m³ |
| heat capacity | J/(kg K) |
| viscosity | Pa s |

`Tmin()`, `Tmax()`, `pmin()` and `pmax()` give the range a model declares: for H₂O–NaCl 273.16–1273.15 K and
1e5–5e8 Pa. Which inputs outside it are rejected is listed under
[Errors and threads](errors.md#input-checks).

## In this guide

- [Evaluating states](tpx.md): T–p–X and H–p–X evaluation, single points and batches, the result fields
- [Errors and threads](errors.md): exceptions, status codes, thread safety, OpenMP
- [Phase region](phase_region.md): the phase of a state without evaluating its properties
- [Phase boundaries](phase_boundaries.md): the critical curve, the three-phase surface, coexisting compositions
