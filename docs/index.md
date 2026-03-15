# xThermo

C++ implementation of the Driesner & Heinrich (2007) equation of state for the H₂O–NaCl system, with Python bindings via SWIG.

xThermo was created by Zhikui Guo and Lars Rüpke

<video autoplay loop muted playsinline style="width:100%;border-radius:6px;margin-bottom:1rem;">
  <source src="assets/video/PhaseChanges.mp4" type="video/mp4">
</video>

## What it does

xThermo calculates thermodynamic properties — density, enthalpy, heat capacity, viscosity, and phase state — for:

- **Pure water** — IAPWS-84 and IAPWS-95 formulations
- **Pure NaCl** — solid, liquid, and vapour phases
- **H₂O–NaCl mixtures** — full phase diagram covering liquid, vapour, halite, and three-phase (V+L+H) regions

**Valid range:** 0–1000 °C, 1–5000 bar, 0–100 wt% NaCl.

## Phase diagram

The isobaric H₂O–NaCl phase diagram at 300 bar computed with xThermo:

![H2O-NaCl phase diagram](assets/img/phase_diagram_PTX.png)

## Quick start

=== "Python"

    ```python
    from xThermal import H2ONaCl

    sw = H2ONaCl.cH2ONaCl("IAPS84")
    props = sw.UpdateState_TPX(573.15, 300e5, 0.1)  # T [K], p [Pa], X [kg/kg]
    print(props.Rho, props.H, props.phase)
    ```

=== "C++"

    ```cpp
    #include "H2ONaCl.h"
    #include <iostream>
    using namespace xThermal;

    int main() {
        H2ONaCl::cH2ONaCl eos("IAPS84");
        ThermodynamicProperties props;
        eos.UpdateState_TPX(props, 573.15, 300e5, 0.1);  // T [K], p [Pa], X [kg/kg]
        std::cout << props.Rho << std::endl;
    }
    ```

## References

Driesner, T. & Heinrich, C.A. (2007). The system H₂O–NaCl. Part I: Correlation formulae for phase relations in temperature–pressure–composition space from 0 to 1000 °C, 0 to 5000 bar, and 0 to 1 XNaCl. *Geochimica et Cosmochimica Acta*, **71**, 4880–4901.

Driesner, T. (2007). The system H₂O–NaCl. Part II: Correlations for molar volume, enthalpy, and isobaric heat capacity from 0 to 1000 °C, 1 to 5000 bar, and 0 to 1 XNaCl. *Geochimica et Cosmochimica Acta*, **71**, 4902–4919.
