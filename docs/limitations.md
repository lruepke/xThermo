# Model, verification and limitations

xThermo implements the phase relations and thermodynamic properties of the H₂O–NaCl system following the formulations of Driesner & Heinrich (2007) and Driesner (2007). The model covers temperatures from 0 to 1000 °C, pressures from 1 to 5000 bar, and compositions ranging from pure water to pure NaCl.

This page describes the scientific formulation, its numerical verification in xThermo, and relevant implementation choices and limitations.

## The model

### Phase relations

Driesner & Heinrich (2007) describe the phase relations of the H₂O–NaCl system as functions of temperature, pressure and NaCl mole fraction. Their correlations define the phase boundaries, including the critical curve, halite liquidus, vapour–liquid–halite equilibrium, and the melting, boiling and sublimation curves of NaCl. They also describe the compositions of coexisting liquid and vapour phases.

Together, these relations determine the stable phase assemblage for a given temperature, pressure and bulk composition.

The [Gallery](gallery.md#how-the-phase-diagram-changes-with-depth) illustrates the phase relations at different pressures and [in three dimensions](gallery.md#the-phase-diagram-in-three-dimensions).

### Volumetric and thermal properties

Driesner (2007) describes the molar volume, specific enthalpy and isobaric heat capacity of H₂O–NaCl fluids using a corresponding-temperature approach. In this formulation, the properties of saline fluids are related to those of pure water through separate temperature-scaling functions for volume and enthalpy. The isobaric heat capacity is obtained from the temperature derivative of the enthalpy.

The original formulation uses the IAPS-84 equation of state for pure water, which remains the reference backend in xThermo. An alternative implementation using IAPWS-95 is also available.

In regions where the scaled temperature falls outside the applicable range of the pure-water formulation, extrapolations are required. These are discussed below under implementation differences and limitations.

The properties of halite are described by separate correlations.

### Viscosity

The dynamic viscosity of H₂O–NaCl fluids is calculated following Klyukin et al. (2017), using a corresponding-temperature approach based on the viscosity of pure water.

Pure-water viscosity is evaluated using the IAPWS (2008) formulation with the IAPWS-95 backend and the NBS/NRC Steam Tables correlation implemented in PROST with the IAPS-84 backend.

### Water backends

Driesner (2007) uses the IAPS-84 formulation for the thermodynamic properties of pure water. xThermo retains IAPS-84 as its reference backend and also supports IAPWS-95 as an alternative.

| Backend | Formulation | Pure-water implementation |
|---|---|---|
| `IAPS84` | IAPS-84 (Haar et al., 1984; Kestin & Sengers, 1986) | PROST 4.1 |
| `IAPWS95` | IAPWS-95 (Wagner & Pruß, 2002) | xThermo implementation |

The choice of water backend affects the calculated thermodynamic properties of H₂O–NaCl fluids because Driesner's corresponding-temperature formulation relies on pure-water properties.

IAPS-84 provides consistency with the original formulation, while IAPWS-95 offers an alternative based on the more recent pure-water equation of state. The two backends generally give similar results, although differences can become more pronounced near phase boundaries and the critical region.

The published temperature, pressure and composition limits of the H₂O–NaCl correlations apply independently of the selected water backend.

### Model domain

The H₂O–NaCl correlations of Driesner & Heinrich (2007) and Driesner (2007) cover temperatures from **0 to 1000 °C**, pressures from **1 to 5000 bar**, and the complete composition range from pure water to pure NaCl.

While the published correlations use NaCl mole fraction, **xThermo's public API uses NaCl mass fraction (kg/kg)**. For example, seawater with 3.2 wt% NaCl is specified as `X = 0.032`.

The corresponding limits reported by the H₂O–NaCl class are `Tmin() = 273.16 K`, `Tmax() = 1273.15 K`, `pmin() = 1e5 Pa`, and `pmax() = 5e8 Pa`.

Not all inputs outside the published domain are rejected automatically. See [Input validation and error handling](#input-validation-and-error-handling) for details.

## Verification

### What was compared

The numerical implementation of xThermo has been systematically verified against the published formulations and supplementary data of Driesner & Heinrich (2007) and Driesner (2007).

The verification includes:

- **Published reference tables:** All ten supplementary tables were evaluated, comprising 29,680 values across 24 property comparisons. Of these, 13 comparisons reproduce the published values to their reported precision. Differences in the remaining comparisons have been investigated and documented, including differences arising from implementation choices, numerical behaviour near phase boundaries and inconsistencies between individual tabulated values and the published equations. A small number of marginal discrepancies remain unresolved.
- **Published equations:** Explicit analytical correlations were independently evaluated and compared with both the published tables and xThermo.
- **SoWat:** Selected results were compared with SoWat 0.2 (2009), providing an additional reference implementation of the Driesner formulation.
- **Pure-water formulations:** The IAPS-84 and IAPWS-95 implementations were checked against their respective published reference values and available independent calculations.
- **Internal consistency:** Thermodynamic derivatives and the consistency between temperature-based (T–p–X) and enthalpy-based (H–p–X) evaluations were examined.

These comparisons are maintained as automated regression tests in the development repository.

### What the agreement means

The comparisons demonstrate how closely xThermo reproduces the published Driesner formulations and reference data. They provide confidence in the numerical implementation and help identify differences associated with particular correlations, implementation choices or numerical conditions.

The agreement with reference calculations should be distinguished from the physical accuracy of the underlying thermodynamic model. Driesner & Heinrich (2007) and Driesner (2007) assessed their correlations against available experimental data and discussed the associated uncertainties and limits of applicability.

For example, Driesner (2007) reports density uncertainties of approximately 0.1–0.2% at low temperatures and pressures, and heat-capacity deviations reaching approximately 10% at elevated temperatures and salinities. These estimates characterize the published formulation rather than the numerical precision of xThermo.

### Scope of verification

The verification described above focuses primarily on the phase relations, density and enthalpy correlations of Driesner & Heinrich (2007) and Driesner (2007), together with the underlying pure-water formulations.

Some aspects of the implementation have not yet been evaluated to the same extent:

- The viscosity correlation of Klyukin et al. (2017) for saline fluids and the IAPS-84 pure-water viscosity implementation.
- The thermodynamic properties of halite and liquid NaCl.
- Independent comparisons of heat capacities and density derivatives with experimental measurements.
- Systematic verification of phase identification throughout the full temperature–pressure–composition domain, beyond the states represented in the published reference tables.
- Additional reference comparisons for IAPS-84 beyond its published check values.

The IAPWS-95 backend represents an alternative implementation of the corresponding-temperature formulation. Since Driesner's reference tables were generated using IAPS-84, they do not provide a direct reference for results obtained with IAPWS-95.

These distinctions define the present scope of the verification effort and identify opportunities for further assessment.

## Implementation choices and formulation differences

xThermo follows the formulations of Driesner & Heinrich (2007) and Driesner (2007), with a small number of documented adaptations inherited from earlier GEOMAR implementations or introduced during development.

Comparisons with the published equations, supplementary tables and other implementations, including SoWat and ProBrine (Klyukin et al., 2020), have identified differences in the treatment of selected thermodynamic properties and phase relations.

The following sections describe three cases that have been examined in detail: isobaric heat capacity, liquid enthalpy at high temperatures and low pressures, and vapour composition above the melting point of NaCl. For each case, the implementation choice and its numerical consequences are documented.

Some differences concern regions where the published formulation requires extrapolation or leaves implementation details unspecified. The current treatment is retained pending further scientific evidence.

### Heat capacity

In xThermo, the isobaric heat capacity is calculated as the temperature derivative of the specific enthalpy at constant pressure and composition, following Eq. 27 of Driesner (2007):

<p align="center"><i>c</i><sub>P</sub> = (∂<i>h</i>/∂<i>T</i>)<sub>p,x</sub></p>

SoWat 0.2 and ProBrine apply an additional factor of 1/(1 + <i>x</i>) to the published heat-capacity expression, where <i>x</i> denotes the NaCl mole fraction. This correction is not included in Eq. 27 of Driesner (2007), although SoWat identifies it as an omission in the original publication.

Applying this factor to heat capacity without a corresponding modification of enthalpy changes the relationship between the two properties. xThermo therefore retains the published equation to preserve consistency between enthalpy and its temperature derivative.

The resulting heat capacities differ from SoWat and ProBrine by approximately 1% at seawater salinity, 9% at 26 wt% NaCl and 23% at 58 wt% NaCl. The physical justification for the additional factor remains unresolved and would require further assessment against experimental data.

### Liquid enthalpy of hot brines at low pressure

At temperatures above approximately 600 °C and pressures below 390 bar, the corresponding-temperature approach of Driesner (2007) may map saline liquid states into the vapour region of the pure-water formulation, where the enthalpy correlation cannot be applied directly.

To evaluate liquid enthalpy in this region, xThermo uses a smooth extrapolation developed at GEOMAR, following the approach used for density in Driesner (2007). The extrapolation is continuous at its upper pressure boundary, and the resulting heat capacity remains consistent with the temperature derivative of enthalpy.

Comparisons with the supplementary tables of Driesner (2007) show differences of up to approximately 3%, with a median difference of 0.8%. The treatment used to generate the published tabulated values in this region is not specified in the paper.

Because these conditions require extrapolation beyond the underlying correlations, the available comparisons do not establish which treatment more accurately represents the physical properties of saline liquids.

### Vapour composition above the melting point of NaCl

Above the melting point of NaCl (800.7 °C), Eq. 15 of Driesner & Heinrich (2007) contains a parameter whose value is not explicitly defined for this temperature range.

xThermo assigns a value of unity to this parameter, consistent with the published phase-boundary conditions. SoWat uses a different, temperature-dependent treatment that is not specified in the published paper.

As a result, the calculated NaCl concentration in the vapour phase differs between the two implementations. Relative to SoWat, xThermo predicts approximately 1.3% lower concentrations at 810 °C, 6–13% lower concentrations at 850–900 °C, and differences of up to approximately 24% at 1000 °C. At pressures below 300 bar, the corresponding vapour concentrations are generally small, approximately 0.01–0.2 wt% NaCl.

Below 800.7 °C, xThermo reproduces the vapour compositions reported in the supplementary tables of Driesner & Heinrich (2007).

Available experimental data do not presently provide a basis for selecting between these treatments.

### Other implementation details

Several additional implementation choices originate from earlier GEOMAR developments or were introduced during the development of xThermo.

| Property or correlation | Implementation detail |
|---|---|
| **Low-temperature density correction** | The correction term <i>D</i>(<i>T</i>) in Eq. 13 of Driesner (2007) is omitted because it can shift the corresponding temperature below the freezing point of pure water. The resulting density differences are generally below 0.01% above approximately 50 °C and reach 0.18% at lower temperatures near halite saturation. |
| **Molar mass of water** | xThermo uses 18.015268 g/mol instead of 18.015 g/mol in Driesner (2007). The resulting differences in density and composition conversions are approximately 0.0015% or less. |
| **Low-pressure density extrapolation** | The implementation of Eq. 17 of Driesner (2007) differs slightly from the published expression, whose notation is ambiguous. Density differences reach approximately 0.16–0.46% near boiling at 1–15 bar. |
| **High-temperature density extrapolation** | The implementation of Eq. 18 of Driesner (2007) differs slightly from the published expression, with density differences of approximately 0.06% or less in the investigated high-temperature, low-pressure region. |
| **High-temperature liquid composition** | Above 800.7 °C, the liquid branch of the vapour–liquid coexistence surface uses a slightly different pressure reference, resulting in composition differences of approximately 10<sup>−5</sup> or less. |
| **Pure-water viscosity extrapolation** | With the IAPS-84 backend, the PROST water-viscosity correlation is evaluated beyond its nominal upper temperature limit of 900 °C. These results should be considered extrapolations. |

These differences are documented to distinguish implementation choices from the published formulation. Their significance depends on the thermodynamic conditions and properties of interest.

## Multiphase properties

In multiphase regions, xThermo provides the thermodynamic properties of the individual coexisting phases as well as selected properties of the bulk fluid mixture.

Phase-specific quantities, including density (`Rho_l`, `Rho_v`), enthalpy (`H_l`, `H_v`), heat capacity (`Cp_l`, `Cp_v`), viscosity (`Mu_l`, `Mu_v`) and composition (`X_l`, `X_v`), are evaluated at the equilibrium composition of each phase.

Bulk quantities require additional definitions because their values may depend on phase proportions, phase equilibrium and the treatment of phase transitions. The following table summarizes the conventions currently used in xThermo.

| Property | Single-phase state | Multiphase state |
|---|---|---|
| **Phase fractions** (`S_l`, `S_v`, `S_h`) | 1 for the phase present | Volume fractions of the coexisting phases in the H₂O–NaCl model. |
| **Bulk density** (`Rho`) | Density of the phase | Volume-weighted average of the phase densities: ρ = ∑<sub>i</sub> <i>S</i><sub>i</sub> ρ<sub>i</sub>. |
| **Bulk specific enthalpy** (`H`) | Specific enthalpy of the phase | Mass-weighted average of the phase enthalpies: <i>h</i> = (∑<sub>i</sub> <i>S</i><sub>i</sub> ρ<sub>i</sub> <i>h</i><sub>i</sub>)/ρ. |
| **Heat capacity** (`Cp`) | Isobaric heat capacity of the phase, <i>c</i><sub>P</sub> = (∂<i>h</i>/∂<i>T</i>)<sub>p,x</sub>. | In two-phase regions, the mass-weighted average of the individual phase heat capacities, evaluated at fixed phase proportions and compositions. In V+L+H, not computed (`0` is a placeholder). |
| **Isothermal compressibility** (`IsothermalCompressibility`) | Isothermal compressibility of the phase. | In two-phase H₂O–NaCl regions, calculated from the pressure dependence of bulk density at constant temperature and bulk composition, allowing the phases to re-equilibrate. In V+L+H, not computed (`0` is a placeholder). |
| **Thermal expansivity** (`IsobaricExpansivity`) | Isobaric thermal expansivity of the phase. | Not computed (`0` is a placeholder). |
| **Density derivatives** (`dRhodT`, `dRhodP`) | Temperature and pressure derivatives of phase density. | Not computed (`0` values are placeholders); use phase-specific properties where available. |
| **Viscosity** (`Mu`) | Dynamic viscosity of the phase. | No general effective multiphase viscosity model is implemented. The returned value depends on phase assemblage and evaluation method (see below). |

**Phase-fraction convention:** In the H₂O–NaCl model, `S_l`, `S_v` and `S_h` represent volume fractions. The pure-water classes use a different convention: `S_v` represents vapour mass fraction (quality). These quantities must not be interchanged without conversion.

### Heat capacity in multiphase regions

In two-phase regions, xThermo reports a phase-weighted heat capacity:

<p align="center"><i>c</i><sub>P,bulk</sub> = (∑<sub>i</sub> <i>S</i><sub>i</sub> ρ<sub>i</sub> <i>c</i><sub>P,i</sub>)/ρ</p>

This quantity represents the heat capacity at fixed phase proportions and compositions, often referred to as a *frozen* heat capacity. It does not include latent heat associated with phase transitions or changes in the equilibrium compositions of the coexisting phases.

Consequently, it generally differs from the equilibrium temperature derivative of the bulk enthalpy,

<p align="center">(∂<i>H</i>/∂<i>T</i>)<sub>p,X</sub></p>

which includes changes in phase proportions and compositions.

At fixed pressure, the temperature of a V+L+H equilibrium state is constrained by the three-phase equilibrium. xThermo does not compute a bulk heat capacity for these states; the returned value of zero is a placeholder, not a physical heat capacity.

For example, at 400 °C, 250 bar and 10 wt% NaCl, xThermo returns a phase-weighted `Cp` of approximately 9.2 kJ/(kg K), whereas the temperature derivative of bulk enthalpy is approximately 31 kJ/(kg K):

```python
import numpy as np
import xthermo
from xthermo import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")
dT = 0.01                                            # K
s = xthermo.tpx(sw, 673.15 + np.array([-dT, 0.0, dT]), 250e5, 0.10)
print(s.phase_name[1])                               # V+L
print(s.Cp[1] / 1e3)                                 # frozen: 9.2 kJ/(kg K)
print((s.H[2] - s.H[0]) / (2 * dT) / 1e3)            # equilibrium (dH/dT): 31 kJ/(kg K)
```

### Isothermal compressibility in multiphase regions

In two-phase regions, xThermo calculates the isothermal compressibility from the bulk density:

<p align="center">κ<sub>T</sub> = −(1/ρ)(∂ρ/∂<i>p</i>)<sub>T,X</sub></p>

The calculation accounts for changes in phase proportions and equilibrium compositions as pressure varies.

For three-phase vapour–liquid–halite states, this quantity is not evaluated. A returned value of zero indicates that the property has not been computed, rather than zero physical compressibility.

The pure-water classes also leave the bulk compressibility uncomputed for two-phase states obtained through enthalpy-based evaluation (H–p–X).

### Viscosity in multiphase regions

The effective viscosity of a multiphase fluid depends on its flow characteristics and cannot generally be determined from the thermodynamic properties of the individual phases alone.

xThermo does not implement a constitutive model for effective multiphase viscosity. In two-phase H₂O–NaCl regions, the bulk `Mu` field currently contains a placeholder value of `1`. In three-phase states, and in two-phase pure-water states evaluated from enthalpy-based input, it contains a volume-fraction-weighted average of liquid and vapour viscosities, excluding halite.

Neither convention should be interpreted as an effective viscosity for multiphase flow calculations. Applications should use the phase-specific viscosities (`Mu_l`, `Mu_v`) together with an appropriate multiphase-flow model.

**Absent phases:** Fields belonging to a phase that is not present retain zero in the C++ structure and classic interface; `xthermo.tpx` and `xthermo.hpx` can return NaN instead (`absent="nan"`). See [Evaluating states](python/evaluate.md#absent-phases).

## Input validation and error handling

The published model domain and the input range enforced at runtime are not identical:

| Input condition | Current behaviour |
|---|---|
| Temperature above 1000 °C | `OutOfRangeError` |
| Enthalpy above the upper model limit | `OutOfRangeError` |
| Non-finite inputs, non-positive pressure or negative salinity | `InvalidInputError` |
| Pressure below 1 bar | Not rejected automatically; results are extrapolations outside the published domain. |
| Pressure above 5000 bar | Behaviour depends on the water backend; IAPS-84 may return extrapolated values, while IAPWS-95 raises `OutOfRangeError`. |
| NaCl mass fraction above 1 | Not rejected automatically; returned results have no physical interpretation. |

**Accepting an input numerically does not imply that its result lies within the published validity range.** Check input units and model applicability before evaluation.

When evaluation fails, the library raises a typed exception rather than returning a previous or partial state. C++ exceptions derive from `xThermal::xThermalBaseError`; Python exceptions from `xthermo.EvaluationError`. These include `InvalidInputError`, `OutOfRangeError`, `NotConvergedError` and `NotAvailableError`. C++ callers can use `tryUpdateState_TPX` and `tryUpdateState_HPX` to obtain an `EvaluationStatus` without throwing.

For arrays, `xthermo.tpx` and `xthermo.hpx` raise the first error by default (`errors="raise"`). With `errors="nan"`, failed points are set to NaN and per-point outcomes are reported in a `Status` array, while valid points are preserved. See [Errors and status](python/errors.md).

When an individual property cannot be evaluated at an otherwise valid state, that property may be NaN while the remaining properties are valid. For example, the IAPS-84 pure-water viscosity has its own validity restrictions (above 5000 bar; above 3500 bar for T > 150 °C; above 3000 bar for T > 600 °C), which also affect brine viscosity derived from it.

Some phase-boundary helper functions have a different convention: outside their supported temperature ranges, they may print a warning and return zero (for example, `X_Critical`).

## Applicability and known limitations

### Limitations of the published correlations

The formulations of Driesner & Heinrich (2007) and Driesner (2007) provide a broad description of the H₂O–NaCl system. As with any empirical formulation, their applicability depends on the experimental data and assumptions underlying the correlations.

Several regions require particular consideration:

- **Extrapolated conditions:** Brine densities at low pressures near boiling, and the density and enthalpy of hot saline liquids above approximately 600 °C and below 390 bar, involve extrapolations of the published correlations. Conditions below 1 bar lie outside the published model domain.
- **Highly saline liquids at low pressure:** In a restricted region near halite saturation, at pressures of approximately 1–5 bar and temperatures of 10–100 °C, the published density correlation can yield negative isothermal compressibility. Although the associated variation in calculated density is small, negative compressibility is not thermodynamically consistent with a stable homogeneous phase. xThermo retains the published density correlation and calculates derivatives consistently with it.
- **Critical region:** Predictions near the critical point of water require particular care because of the limitations and increasing sensitivity of the underlying pure-water formulations. These effects can also influence dilute saline fluids through the corresponding-temperature approach.
- **High-temperature phase compositions and enthalpy:** The vapour-composition treatment above 800.7 °C and the extrapolation of liquid enthalpy at high temperature and low pressure are discussed under [Implementation choices and formulation differences](#implementation-choices-and-formulation-differences).

### Known numerical issues

The following numerical issues have been identified in specific regions of the model domain. They are documented here to help users assess their relevance to particular applications.

**Phase-diagram geometry**

Some functions used to construct two-dimensional phase diagrams (`Slice_constP` and `Slice_constT`) have known geometric inaccuracies:

- At pressures below approximately 221 bar, a narrow part of the vapour–liquid region may be omitted from enthalpy–salinity diagrams at very low salinity (below approximately 0.13 wt% NaCl).
- In isothermal diagrams, the liquid–halite region may extend below the three-phase pressure into the vapour–halite field. The current plotting routines generally cover this overlap with the correctly drawn vapour–halite region.
- Isothermal sections below 374 °C may contain a spurious critical point at zero salinity.

These issues concern the construction of phase-region polygons and do not affect direct thermodynamic state evaluations through `tpx` or `hpx`.

**Thermodynamic property evaluation**

- **Pure-water saturation conditions:** Very close to the saturation pressure of pure water, the IAPS-84 backend may produce inconsistent phase labels and densities or small discontinuities in calculated properties. The IAPWS-95 backend may fail to converge at the saturation boundary.
- **IAPWS-95 at low temperatures and high pressures:** For certain saline compositions, the corresponding-temperature transformation reaches conditions outside the supported range of the pure-water formulation. This can prevent evaluation at low temperatures or near the upper pressure boundary of 5000 bar. The IAPS-84 backend can evaluate these conditions.
- **Low-pressure density derivatives:** Near boiling at approximately 1–5 bar, calculated pressure derivatives of brine density can differ significantly from numerical derivatives of the density correlation. The calculated density itself is unaffected.
- **High-temperature density with IAPWS-95:** In the low-pressure extrapolation region, small numerical discontinuities may occur in density and its temperature derivative.
- **Enthalpy extrapolation below 1 bar:** A known implementation issue can affect calculated liquid enthalpies at pressures below the published model domain, or when liquid properties are requested at conditions where no stable liquid phase exists.

### Under investigation

Two related questions concerning pure NaCl near its melting point remain under investigation.

**Enthalpy-based evaluation near melting**

For pure NaCl, enthalpies between those of solid halite and liquid NaCl at the melting temperature correspond to coexistence of the two phases. In this region, enthalpy-based evaluation (`hpx`) may return a successful result whose calculated enthalpy differs from the requested value.

For example, at 250 bar, an input enthalpy of 1.031 MJ/kg produces a liquid–halite state at approximately 806.9 °C with a returned enthalpy of 0.807 MJ/kg.

Until this behaviour is resolved, enthalpy-based evaluations of pure NaCl within the melting interval should be treated with caution.

**Enthalpy of solid and liquid NaCl**

The enthalpy of solid halite is calculated using a separate correlation, while the enthalpy of liquid NaCl is obtained from the Driesner (2007) fluid formulation at pure-NaCl composition.

At 1 bar, the calculated enthalpy difference at melting is approximately 0.25 MJ/kg. The consistency of the reference states and the agreement of this difference with published enthalpy-of-fusion data have not yet been established.

Further examination is required to determine whether the observed behaviour reflects a reference-state inconsistency, a numerical issue or another aspect of the formulation.

## References

Driesner, T. & Heinrich, C.A. (2007). The system H₂O–NaCl. Part I: Correlation formulae for phase relations in
temperature–pressure–composition space from 0 to 1000 °C, 0 to 5000 bar, and 0 to 1 X<sub>NaCl</sub>. *Geochim. Cosmochim.
Acta* **71**, 4880–4901. [doi:10.1016/j.gca.2006.01.033](https://doi.org/10.1016/j.gca.2006.01.033)

Driesner, T. (2007). The system H₂O–NaCl. Part II: Correlations for molar volume, enthalpy, and isobaric heat
capacity from 0 to 1000 °C, 1 to 5000 bar, and 0 to 1 X<sub>NaCl</sub>. *Geochim. Cosmochim. Acta* **71**, 4902–4919.
[doi:10.1016/j.gca.2007.05.026](https://doi.org/10.1016/j.gca.2007.05.026)

Haar, L., Gallagher, J.S. & Kell, G.S. (1984). *NBS/NRC Steam Tables*. Hemisphere, Washington.

Kestin, J. & Sengers, J.V. (1986). New international formulations for the thermodynamic properties of light and heavy
water. *J. Phys. Chem. Ref. Data* **15**, 305. [doi:10.1063/1.555772](https://doi.org/10.1063/1.555772)

Klyukin, Yu.I., Lowell, R.P. & Bodnar, R.J. (2017). A revised empirical model to calculate the dynamic viscosity of
H₂O–NaCl fluids at elevated temperatures and pressures (≤ 1000 °C, ≤ 500 MPa, 0–100 wt% NaCl). *Fluid Phase
Equilib.* **433**, 193–205. [doi:10.1016/j.fluid.2016.11.002](https://doi.org/10.1016/j.fluid.2016.11.002)

Klyukin, Y.I., Haroldson, E.L. & Steele-MacInnis, M. (2020). A comprehensive numerical model for the thermodynamic and transport properties of H₂O–NaCl fluids. *Chemical Geology* **557**, 119840. [doi:10.1016/j.chemgeo.2020.119840](https://doi.org/10.1016/j.chemgeo.2020.119840)

IAPWS (2008). *Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary Water Substance* (R12-08).
International Association for the Properties of Water and Steam.

Wagner, W. & Pruß, A. (2002). The IAPWS formulation 1995 for the thermodynamic properties of ordinary water substance
for general and scientific use. *J. Phys. Chem. Ref. Data* **31**, 387–535.
[doi:10.1063/1.1461829](https://doi.org/10.1063/1.1461829)
