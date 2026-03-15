# C++ API overview

xThermo is a C++ library. The Python bindings are a thin SWIG wrapper around the same C++ classes.

## Header

```cpp
#include "H2ONaCl.h"
```

After `make install`, all public headers are installed to `<install>/include/`.

## Linking with CMake

```cmake
cmake_minimum_required(VERSION 3.18)
project(my_app LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)

set(xThermo_DIR "$ENV{xThermo_DIR}" CACHE PATH "xThermo install prefix")

include_directories(${xThermo_DIR}/include)
link_directories(${xThermo_DIR}/STATIC)

add_executable(my_app main.cpp)
target_link_libraries(my_app xThermal)
```

Set the environment variable before configuring:

```bash
export xThermo_DIR=/path/to/xThermo/install
cmake ..
```

## Linking with GCC directly

`libxThermal.a` is a static archive but it does **not** bundle GSL internally — it only references GSL symbols. You must therefore pass GSL explicitly on the linker command line:

```bash
export xThermo_DIR=/path/to/xThermo/install
GSL_LIB=${xThermo_DIR}/../ThirdParties/install/gsl_$(uname -s)_$(uname -m)/lib

g++ -std=c++11 \
    -I${xThermo_DIR}/include \
    -L${xThermo_DIR}/STATIC \
    -L${GSL_LIB} \
    -o my_app main.cpp \
    -lxThermal -lgsl -lgslcblas -lm
```

!!! note "Why GSL?"
    xThermo uses GSL for multidimensional root-finding (saturation pressure, phase boundaries).
    When you build with CMake, `target_link_libraries` resolves the GSL dependency automatically.
    When linking manually with `g++`, you must add `-lgsl -lgslcblas` yourself and point
    `-L` at the GSL library built in `ThirdParties/`.

The `$(uname -s)_$(uname -m)` part expands to the platform-specific install directory, e.g.:

| Platform | Path suffix |
|---|---|
| macOS Apple Silicon | `gsl_Darwin_arm64` |
| macOS Intel | `gsl_Darwin_x86_64` |
| Linux x86_64 | `gsl_Linux_x86_64` |

## Key classes

| Class | Header | Description |
|---|---|---|
| `xThermal::H2ONaCl::cH2ONaCl` | `H2ONaCl.h` | Full H₂O–NaCl EOS |
| `xThermal::ThermodynamicProperties` | `thermo/DataStructures.h` | Property struct returned by UpdateState |
| `xThermal::H2O::cIAPS84` | `H2O/IAPS84.h` | Pure water (IAPS-84) |
| `xThermal::H2O::cIAPWS95` | `H2O/IAPWS95.h` | Pure water (IAPWS-95) |
| `xThermal::NaCl::cNaCl` | `NaCl/NaCl.h` | Pure NaCl EOS |

## Units

All C++ functions use SI units unless noted otherwise in the Doxygen comments:

| Quantity | Unit |
|---|---|
| Temperature | K |
| Pressure | Pa |
| Composition X | kg/kg (mass fraction) |
| Density | kg/m³ |
| Enthalpy | J/kg |
