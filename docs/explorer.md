# Explorer

The [xThermo Explorer](app/index.html) evaluates single states of H₂O–NaCl fluids in the web browser. It runs the
xThermo C++ library, compiled to WebAssembly, on your own computer; no data is sent to a server.

## What it does

- **Input:** temperature, pressure and salinity (T-P-X), or specific enthalpy, pressure and salinity (H-P-X).
- **Water equation of state:** IAPS-84 (default) or IAPWS-95, as for the [Python](python/index.md) and
  [C++](cpp/index.md) interfaces.
- **Output:**
    - the phase region;
    - bulk temperature, enthalpy, density, isobaric heat capacity, viscosity, isothermal compressibility and isobaric
      expansivity;
    - per phase (liquid or single-phase fluid, vapour, halite): volume and mass fraction, salinity, density, enthalpy,
      heat capacity and viscosity.
    - In multiphase states, bulk viscosity and expansivity are shown as "not defined" or "not computed", as are heat
      capacity and compressibility where the library does not compute them (see
      [Multiphase properties](limitations.md#multiphase-properties)).
- **Errors:** a state the library does not evaluate is shown with its status (InvalidInput, OutOfRange, NotConverged)
  and the library's own message, as in [Errors and status](python/errors.md).

| Quantity | Input unit | Output unit |
|---|---|---|
| Temperature | °C | °C and K |
| Pressure | bar | bar |
| Salinity | wt % NaCl or kg/kg | wt % NaCl |
| Specific enthalpy | kJ/kg | kJ/kg |
| Density | – | kg/m³ |
| Heat capacity | – | J/(kg·K) |
| Viscosity | – | Pa·s |
| Compressibility, expansivity | – | 1/Pa, 1/K |

The evaluation runs in a background worker. A slow H-P-X state can be stopped with **Cancel**.

## Accuracy

The Explorer uses the same code as the native library, with no approximation or lookup table. Each build of the
module is compared with the native library built from the same commit, on about 14 000 states with both water
backends: phases, statuses and error messages are identical, and the properties agree within 10⁻⁹ relative (tiny
phase fractions and near-zero compressibilities within a small absolute round-off). The
[Model and limitations](limitations.md) page applies unchanged, including its notes on platform-dependent last digits
and on the narrow H-P-X regions where the library reports NotConverged.

The Explorer evaluates one state at a time. For profiles, grids, phase diagrams and lookup tables, use the
[Python interface](python/index.md) or the [command-line tool](cli.md).

## Licence and source

The WebAssembly module contains the GNU Scientific Library and is therefore distributed as a combined work under the
GNU General Public License, version 3 or later. The xThermo source code itself remains under the MIT licence. The
About section of the Explorer lists every component with its licence text and links to the corresponding source
code, including the exact xThermo revision the module was built from and the GSL source archive. See also
[About](about.md#license).
