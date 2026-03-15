

# Class xThermal::H2ONaCl::cH2ONaCl



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**H2ONaCl**](namespacexThermal_1_1H2ONaCl.md) **>** [**cH2ONaCl**](classxThermal_1_1H2ONaCl_1_1cH2ONaCl.md)



_Class of_ \(H_2O-NaCl\) _EOS._[More...](#detailed-description)

* `#include <H2ONaCl.h>`



Inherits the following classes: [xThermal::cxThermal](classxThermal_1_1cxThermal.md)
























## Public Attributes inherited from xThermal::cxThermal

See [xThermal::cxThermal](classxThermal_1_1cxThermal.md)

| Type | Name |
| ---: | :--- |
|  [**bool**](classxThermal_1_1xThermalError.md) | [**m\_isShowProgressBar**](classxThermal_1_1cxThermal.md#variable-m_isshowprogressbar)  <br> |






























## Public Functions

| Type | Name |
| ---: | :--- |
|  std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; | [**D\_Tstar\_V**](#function-d_tstar_v-22) (std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; X) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**H\_phase**](#function-h_phase-14) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X, [**double**](classxThermal_1_1xThermalError.md) & H, [**PhaseType**](namespacexThermal.md#enum-phasetype) phase) <br>_Calculate specific enthalpy of liquid or vapor phase of H2O-NaCl._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**H\_phase**](#function-h_phase-24) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X, [**double**](classxThermal_1_1xThermalError.md) & H, [**double**](classxThermal_1_1xThermalError.md) & Cp, [**PhaseType**](namespacexThermal.md#enum-phasetype) phase) <br>_Calculate specific enthalpy of liquid or vapor phase of H2O-NaCl._  |
|  std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; | [**H\_phase**](#function-h_phase-34) (std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; X, [**PhaseType**](namespacexThermal.md#enum-phasetype) phase) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**H\_phase**](#function-h_phase-44) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X, [**PhaseType**](namespacexThermal.md#enum-phasetype) phase) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**HminHmax\_VLH\_triangle**](#function-hminhmax_vlh_triangle) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & H\_v0, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & H\_l0, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & H\_h0, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X\_v0, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X\_l0, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X0, [**double**](classxThermal_1_1xThermalError.md) & Hmin0, [**double**](classxThermal_1_1xThermalError.md) & Hmax0) <br>_Calculate minimum and maximum specific enthalpy of a specific VLH "triangle" for given bulk salinity X. The "VLH triangle" is boundary (edge) of a three-phase region at constant P slice. The max and min H will be calculated for the latent heat of evaporation, see also grant1979compressibility._  |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Mol2Wt**](#function-mol2wt-12) ([**double**](classxThermal_1_1xThermalError.md) X\_mol) <br>_Convert molar fraction of NaCl to mass fraction._  |
|  std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; | [**Mol2Wt**](#function-mol2wt-22) (std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; X\_mol) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**Mu\_phase**](#function-mu_phase) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X, [**double**](classxThermal_1_1xThermalError.md) & Mu, [**PhaseType**](namespacexThermal.md#enum-phasetype) phase) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**P\_Critical**](#function-p_critical-12) ([**double**](classxThermal_1_1xThermalError.md) T, [**double**](classxThermal_1_1xThermalError.md) & P\_crit) <br>_Critical pressure as a function of temperature, see Eq. 5a-c of Driesner2007Part1._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**P\_Critical**](#function-p_critical-22) (std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & P\_crit) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**P\_VLH**](#function-p_vlh-13) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Vapor-Liquid-Halite coexistence pressure as a function of temperature._  |
|  [**double**](classxThermal_1_1xThermalError.md) | [**P\_VLH**](#function-p_vlh-23) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**P\_VLH**](#function-p_vlh-33) (std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & P) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**P\_X\_Critical**](#function-p_x_critical-12) ([**double**](classxThermal_1_1xThermalError.md) T, [**double**](classxThermal_1_1xThermalError.md) & P\_crit, [**double**](classxThermal_1_1xThermalError.md) & X\_crit) <br>_Calculate critical pressure(equation 5a-c of Driesner2007Part1) and salinity(equation 7a,b of Driesner2007Part1) by given temperature._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**P\_X\_Critical**](#function-p_x_critical-22) (std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & P\_crit, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & X\_crit) <br> |
|  [**TriMesh**](structxThermal_1_1TriMesh.md) | [**PhaseBoundary\_HaliteLiquidus**](#function-phaseboundary_haliteliquidus) (std::string fmt\_out="", [**double**](classxThermal_1_1xThermalError.md) p\_max=2500E5, [**double**](classxThermal_1_1xThermalError.md) dT=10, [**double**](classxThermal_1_1xThermalError.md) dp=50E5) <br> |
|  [**DeformLinearMesh**](structxThermal_1_1DeformLinearMesh.md) | [**PhaseBoundary\_HaliteLiquidus\_DeformLinear**](#function-phaseboundary_haliteliquidus_deformlinear) ([**double**](classxThermal_1_1xThermalError.md) p\_max=2500E5, [**double**](classxThermal_1_1xThermalError.md) dT=10, [**double**](classxThermal_1_1xThermalError.md) dp=50E5) <br> |
|  [**DeformLinearMesh**](structxThermal_1_1DeformLinearMesh.md) | [**PhaseBoundary\_VH\_DeformLinear**](#function-phaseboundary_vh_deformlinear) ([**int**](classxThermal_1_1xThermalError.md) nT=100, [**int**](classxThermal_1_1xThermalError.md) np=100) <br> |
|  [**DeformLinearMesh**](structxThermal_1_1DeformLinearMesh.md) | [**PhaseBoundary\_VLH\_DeformLinear**](#function-phaseboundary_vlh_deformlinear) ([**int**](classxThermal_1_1xThermalError.md) nT=100, [**int**](classxThermal_1_1xThermalError.md) nX=60) <br> |
|  [**DeformLinearMesh**](structxThermal_1_1DeformLinearMesh.md) | [**PhaseBoundary\_VL\_DeformLinear**](#function-phaseboundary_vl_deformlinear) ([**PhaseType**](namespacexThermal.md#enum-phasetype) VaporOrLiquid, [**int**](classxThermal_1_1xThermalError.md) nT=100, [**int**](classxThermal_1_1xThermalError.md) np=200) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**Rho\_phase**](#function-rho_phase-13) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X, [**double**](classxThermal_1_1xThermalError.md) & rho, [**double**](classxThermal_1_1xThermalError.md) & dRhodP, [**double**](classxThermal_1_1xThermalError.md) & dRhodT, [**PhaseType**](namespacexThermal.md#enum-phasetype) phase) <br>_Calculate density of liquid or vapor phase._  |
|  std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; | [**Rho\_phase**](#function-rho_phase-23) (std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; X, [**PhaseType**](namespacexThermal.md#enum-phasetype) phase) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Rho\_phase**](#function-rho_phase-33) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X, [**PhaseType**](namespacexThermal.md#enum-phasetype) phase) <br> |
|  [**PhaseRegion\_Slice**](structxThermal_1_1PhaseRegion__Slice.md) | [**Slice\_constP**](#function-slice_constp) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) P0, [**size\_t**](classxThermal_1_1xThermalError.md) nPoints=500) <br> |
|  [**PhaseRegion\_Slice**](structxThermal_1_1PhaseRegion__Slice.md) | [**Slice\_constT**](#function-slice_constt) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) T0, [**size\_t**](classxThermal_1_1xThermalError.md) nPoints=500) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**T\_Critical**](#function-t_critical-12) ([**double**](classxThermal_1_1xThermalError.md) P, [**double**](classxThermal_1_1xThermalError.md) & T\_crit) <br>_Solve critical T and X using bisection method for given Temperature._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**T\_Critical**](#function-t_critical-22) (std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & T\_crit) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**T\_HPX**](#function-t_hpx) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & H, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T\_low=[**T\_MIN**](group__PhysicalConstants__H2ONaCl.md#variable-t_min), [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T\_high=[**T\_MAX**](group__PhysicalConstants__H2ONaCl.md#variable-t_max)) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**T\_VL**](#function-t_vl) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X, [**PhaseType**](namespacexThermal.md#enum-phasetype) phase) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**T\_VLH\_P0**](#function-t_vlh_p0) ([**double**](classxThermal_1_1xThermalError.md) P0, [**double**](classxThermal_1_1xThermalError.md) & T\_min, [**double**](classxThermal_1_1xThermalError.md) & T\_max) <br>_Calculate temperature boundary of VLH zone by given pressure P0._  |
|  [**double**](classxThermal_1_1xThermalError.md) | [**T\_VL\_L**](#function-t_vl_l) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T\_low, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T\_high) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**T\_VL\_V**](#function-t_vl_v) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T\_low, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T\_high) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**T\_X\_Critical**](#function-t_x_critical-12) ([**double**](classxThermal_1_1xThermalError.md) P, [**double**](classxThermal_1_1xThermalError.md) & T\_crit, [**double**](classxThermal_1_1xThermalError.md) & X\_crit) <br>_Calculate critical temperature and composition by given pressure._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**T\_X\_Critical**](#function-t_x_critical-22) (std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & T\_crit, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & X\_crit) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**T\_critical**](#function-t_critical) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Tmax**](#function-tmax) () <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Tmax\_VLH**](#function-tmax_vlh) () <br>_Calculate the maximum temperature of the V+L+H zone._  |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Tmin**](#function-tmin) () <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Tmin\_VLH**](#function-tmin_vlh) () <br> |
|  [**TriMesh**](structxThermal_1_1TriMesh.md) | [**Triangulation**](#function-triangulation) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & x\_poly, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & y\_poly, [**double**](classxThermal_1_1xThermalError.md) xIn, [**double**](classxThermal_1_1xThermalError.md) yIn, [**double**](classxThermal_1_1xThermalError.md) dx, [**double**](classxThermal_1_1xThermalError.md) dy) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Ttriple**](#function-ttriple) () <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**UpdateState\_HPX**](#function-updatestate_hpx) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & H, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X) <br> |
|  [**bool**](classxThermal_1_1xThermalError.md) | [**UpdateState\_HPX\_vl\_lowXlowP**](#function-updatestate_hpx_vl_lowxlowp) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & H, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X) <br> |
|  [**bool**](classxThermal_1_1xThermalError.md) | [**UpdateState\_HPX\_vlh**](#function-updatestate_hpx_vlh) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & H, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & T1\_T2) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**UpdateState\_TPX**](#function-updatestate_tpx) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X) <br>_Update state using T, p as independent variables._  |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Wt2Mol**](#function-wt2mol-13) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X\_wt) <br>_Convert mass fraction of NaCl to molar fraction._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**Wt2Mol**](#function-wt2mol-23) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X\_wt, [**double**](classxThermal_1_1xThermalError.md) & X\_mol) <br> |
|  std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; | [**Wt2Mol**](#function-wt2mol-33) (std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; X\_wt) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**XL\_VL**](#function-xl_vl-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Calculate liquid composition on liquid branch of V+L coexistence._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**XL\_VL**](#function-xl_vl-22) (std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**XV\_VL**](#function-xv_vl-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Calculate vapor composition on vapor branch of V+L coexistence._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**XV\_VL**](#function-xv_vl-22) (std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**X\_Critical**](#function-x_critical) ([**double**](classxThermal_1_1xThermalError.md) T, [**double**](classxThermal_1_1xThermalError.md) & X\_crit) <br>_Critical composition as a function of temperature, see Eq. 7a,b of Driesner2007Part1._  |
|  [**double**](classxThermal_1_1xThermalError.md) | [**X\_HaliteLiquidus**](#function-x_haliteliquidus-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Halite-saturated liquid composition, see equation 8 of Driesner2007Part1._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**X\_HaliteLiquidus**](#function-x_haliteliquidus-22) (std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**X\_VH**](#function-x_vh-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Halite-saturated vapor composition at V+H coexistence, see equation 9 and Fig. 8 of Driesner2007Part1._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**X\_VH**](#function-x_vh-22) (std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**X\_VL**](#function-x_vl-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**double**](classxThermal_1_1xThermalError.md) & X\_L, [**double**](classxThermal_1_1xThermalError.md) & X\_V) <br>_Calculate_ \(X_{V}, X_{L}\) _on V+L surface in one step._ |
|  [**void**](classxThermal_1_1xThermalError.md) | [**X\_VL**](#function-x_vl-22) (std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & X\_L, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & X\_V) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**X\_VLH**](#function-x_vlh-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**double**](classxThermal_1_1xThermalError.md) & X\_L, [**double**](classxThermal_1_1xThermalError.md) & X\_V) <br>_Calculate salinity of vapor and liquid at VLH coexistence._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**X\_VLH**](#function-x_vlh-22) (std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & X\_L, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & X\_V) <br> |
|   | [**cH2ONaCl**](#function-ch2onacl-12) (std::string backend\_H2O) <br> |
|   | [**cH2ONaCl**](#function-ch2onacl-22) ([**const**](classxThermal_1_1xThermalError.md) [**cH2ONaCl**](classxThermal_1_1H2ONaCl_1_1cH2ONaCl.md) & sw) <br> |
|  [**PhaseBoundaries**](structxThermal_1_1PhaseBoundaries.md) | [**calc\_PhaseBoundaries**](#function-calc_phaseboundaries) (std::string scale\_X="linear", [**double**](classxThermal_1_1xThermalError.md) ratio\_log\_to\_linear=1, [**double**](classxThermal_1_1xThermalError.md) Xcenter=0.01) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**compressibility\_LH**](#function-compressibility_lh) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T0, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p0, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X0, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & bulkRho, [**double**](classxThermal_1_1xThermalError.md) & compressibility, [**double**](classxThermal_1_1xThermalError.md) dp=1) <br>_Calculate isothermal compressibility of L+H phase._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**compressibility\_VH**](#function-compressibility_vh) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T0, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p0, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X0, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & bulkRho, [**double**](classxThermal_1_1xThermalError.md) & compressibility, [**double**](classxThermal_1_1xThermalError.md) dp=1) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**compressibility\_VL**](#function-compressibility_vl) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T0, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p0, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X0, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & bulkRho, [**double**](classxThermal_1_1xThermalError.md) & compressibility, [**double**](classxThermal_1_1xThermalError.md) dp=1) <br>_Calculate isothermal compressibility of two phase._  |
|  [**double**](classxThermal_1_1xThermalError.md) | [**dPdT\_VLH**](#function-dpdt_vlh) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T) <br> |
| virtual [**PhaseRegion**](namespacexThermal.md#enum-phaseregion) | [**findPhaseRegion\_TPX**](#function-findphaseregion_tpx) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**getPhaseRegion\_HPX**](#function-getphaseregion_hpx) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & H, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X, [**xThermal::PhaseRegion**](namespacexThermal.md#enum-phaseregion) & phase\_region, [**double**](classxThermal_1_1xThermalError.md) & T, [**double**](classxThermal_1_1xThermalError.md) S\_lvh) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**getPhaseRegion\_TPX**](#function-getphaseregion_tpx-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X, [**xThermal::PhaseRegion**](namespacexThermal.md#enum-phaseregion) & phase\_region, [**double**](classxThermal_1_1xThermalError.md) & X\_V, [**double**](classxThermal_1_1xThermalError.md) & X\_L) <br>_Calculate phase region index and composition of vapor and liquid phase if the phase exists._  |
|  std::vector&lt; [**int**](classxThermal_1_1xThermalError.md) &gt; | [**getPhaseRegion\_TPX**](#function-getphaseregion_tpx-22) (std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; X, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & X\_V, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & X\_L) <br> |
|  [**cxThermal**](classxThermal_1_1cxThermal.md) \* | [**get\_pWater**](#function-get_pwater) () const<br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**molar\_mass**](#function-molar_mass) () <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**n1n2\_Tstar\_V**](#function-n1n2_tstar_v-22) (std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; X, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & n1, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & n2) <br> |
| virtual std::string | [**name**](#function-name) () <br> |
|  std::string | [**name\_backend**](#function-name_backend) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**p\_critical**](#function-p_critical) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**pmax**](#function-pmax) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**pmin**](#function-pmin) () <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**prop\_VL**](#function-prop_vl) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X, [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & prop) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**q1q2\_Tstar\_H**](#function-q1q2_tstar_h-22) (std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; X, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & q1, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & q2) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**rhomass\_critical**](#function-rhomass_critical) () <br> |
|   | [**~cH2ONaCl**](#function-ch2onacl) () <br> |


## Public Functions inherited from xThermal::cxThermal

See [xThermal::cxThermal](classxThermal_1_1cxThermal.md)

| Type | Name |
| ---: | :--- |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_T**](classxThermal_1_1cxThermal.md#function-boiling_t-13) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_T**](classxThermal_1_1cxThermal.md#function-boiling_t-23) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**double**](classxThermal_1_1xThermalError.md) & rho\_l, [**double**](classxThermal_1_1xThermalError.md) & rho\_v) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_T**](classxThermal_1_1cxThermal.md#function-boiling_t-33) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props) <br> |
|  [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) | [**Boiling\_T\_props**](classxThermal_1_1cxThermal.md#function-boiling_t_props) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_p**](classxThermal_1_1cxThermal.md#function-boiling_p-13) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_p**](classxThermal_1_1cxThermal.md#function-boiling_p-23) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**double**](classxThermal_1_1xThermalError.md) & rho\_l, [**double**](classxThermal_1_1xThermalError.md) & rho\_v) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_p**](classxThermal_1_1cxThermal.md#function-boiling_p-33) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props) <br> |
|  [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) | [**Boiling\_p\_props**](classxThermal_1_1cxThermal.md#function-boiling_p_props) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T) <br> |
|  [**PhaseRegion**](namespacexThermal.md#enum-phaseregion) | [**Rho\_lookup**](classxThermal_1_1cxThermal.md#function-rho_lookup) ([**double**](classxThermal_1_1xThermalError.md) & Rho\_estimate, [**double**](classxThermal_1_1xThermalError.md) & Rho\_min, [**double**](classxThermal_1_1xThermalError.md) & Rho\_max, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Only lookup density for given T,P,X from loaded LUT. This information can be a good guess of Rho for IAPWS95 solution of Rho._  |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**T\_critical**](classxThermal_1_1cxThermal.md#function-t_critical) () = 0<br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Tmax**](classxThermal_1_1cxThermal.md#function-tmax) () = 0<br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Tmin**](classxThermal_1_1cxThermal.md#function-tmin) () = 0<br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**Triangulation**](classxThermal_1_1cxThermal.md#function-triangulation) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & x\_poly, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & y\_poly, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) pointInMesh, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) dxdy, [**TriMesh**](structxThermal_1_1TriMesh.md) & trimesh, [**bool**](classxThermal_1_1xThermalError.md) normalize=[**true**](classxThermal_1_1xThermalError.md)) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Ttriple**](classxThermal_1_1cxThermal.md#function-ttriple) () = 0<br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**UpdateState\_HPX**](classxThermal_1_1cxThermal.md#function-updatestate_hpx-14) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & H, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br> |
| virtual [**ThermodynamicPropertiesVector**](structxThermal_1_1ThermodynamicPropertiesVector.md) | [**UpdateState\_HPX**](classxThermal_1_1cxThermal.md#function-updatestate_hpx-24) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & H, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & p, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & X, [**bool**](classxThermal_1_1xThermalError.md) isMeshGrid=[**false**](classxThermal_1_1xThermalError.md)) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**UpdateState\_HPX**](classxThermal_1_1cxThermal.md#function-updatestate_hpx-34) ([**ThermodynamicPropertiesArray**](structxThermal_1_1ThermodynamicPropertiesArray.md) & stateArray, [**const**](classxThermal_1_1xThermalError.md) [**size\_t**](classxThermal_1_1xThermalError.md) & N, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) \* H, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) \* p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) \* X=[**NULL**](classxThermal_1_1xThermalError.md)) <br> |
|  [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) | [**UpdateState\_HPX**](classxThermal_1_1cxThermal.md#function-updatestate_hpx-44) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & H, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**UpdateState\_TPX**](classxThermal_1_1cxThermal.md#function-updatestate_tpx-14) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br>_Update state using T, p as independent variables._  |
| virtual [**ThermodynamicPropertiesVector**](structxThermal_1_1ThermodynamicPropertiesVector.md) | [**UpdateState\_TPX**](classxThermal_1_1cxThermal.md#function-updatestate_tpx-24) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & T, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & p, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & X, [**bool**](classxThermal_1_1xThermalError.md) isMeshGrid=[**false**](classxThermal_1_1xThermalError.md)) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**UpdateState\_TPX**](classxThermal_1_1cxThermal.md#function-updatestate_tpx-34) ([**ThermodynamicPropertiesArray**](structxThermal_1_1ThermodynamicPropertiesArray.md) & stateArray, [**const**](classxThermal_1_1xThermalError.md) [**size\_t**](classxThermal_1_1xThermalError.md) & N, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) \* T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) \* p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) \* X=[**NULL**](classxThermal_1_1xThermalError.md)) <br> |
|  [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) | [**UpdateState\_TPX**](classxThermal_1_1cxThermal.md#function-updatestate_tpx-44) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**createLUT\_2D**](classxThermal_1_1cxThermal.md#function-createlut_2d-12) ([**double**](classxThermal_1_1xThermalError.md) xy\_min, [**double**](classxThermal_1_1xThermalError.md) xy\_max, [**double**](classxThermal_1_1xThermalError.md) constZ, LOOKUPTABLE\_FOREST::CONST\_WHICH\_VAR const\_which\_var, LOOKUPTABLE\_FOREST::EOS\_ENERGY TorH, [**int**](classxThermal_1_1xThermalError.md) min\_level=4, [**int**](classxThermal_1_1xThermalError.md) max\_level=6, [**int**](classxThermal_1_1xThermalError.md) update\_which\_props=0) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**createLUT\_2D**](classxThermal_1_1cxThermal.md#function-createlut_2d-22) ([**double**](classxThermal_1_1xThermalError.md) xmin, [**double**](classxThermal_1_1xThermalError.md) xmax, [**double**](classxThermal_1_1xThermalError.md) ymin, [**double**](classxThermal_1_1xThermalError.md) ymax, [**double**](classxThermal_1_1xThermalError.md) constZ, LOOKUPTABLE\_FOREST::CONST\_WHICH\_VAR const\_which\_var, LOOKUPTABLE\_FOREST::EOS\_ENERGY TorH, [**int**](classxThermal_1_1xThermalError.md) min\_level=4, [**int**](classxThermal_1_1xThermalError.md) max\_level=6, [**int**](classxThermal_1_1xThermalError.md) update\_which\_props=0) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**createLUT\_3D**](classxThermal_1_1cxThermal.md#function-createlut_3d-12) ([**double**](classxThermal_1_1xThermalError.md) xyz\_min, [**double**](classxThermal_1_1xThermalError.md) xyz\_max, LOOKUPTABLE\_FOREST::EOS\_ENERGY TorH, [**int**](classxThermal_1_1xThermalError.md) min\_level=4, [**int**](classxThermal_1_1xThermalError.md) max\_level=6, [**int**](classxThermal_1_1xThermalError.md) update\_which\_props=0) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**createLUT\_3D**](classxThermal_1_1cxThermal.md#function-createlut_3d-22) ([**double**](classxThermal_1_1xThermalError.md) xmin, [**double**](classxThermal_1_1xThermalError.md) xmax, [**double**](classxThermal_1_1xThermalError.md) ymin, [**double**](classxThermal_1_1xThermalError.md) ymax, [**double**](classxThermal_1_1xThermalError.md) zmin, [**double**](classxThermal_1_1xThermalError.md) zmax, LOOKUPTABLE\_FOREST::EOS\_ENERGY TorH, [**int**](classxThermal_1_1xThermalError.md) min\_level=4, [**int**](classxThermal_1_1xThermalError.md) max\_level=6, [**int**](classxThermal_1_1xThermalError.md) update\_which\_props=0) <br> |
|   | [**cxThermal**](classxThermal_1_1cxThermal.md#function-cxthermal) () <br> |
| virtual [**PhaseRegion**](namespacexThermal.md#enum-phaseregion) | [**findPhaseRegion\_TPX**](classxThermal_1_1cxThermal.md#function-findphaseregion_tpx) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br> |
|  [**Head\_AMR\_LUT**](structxThermal_1_1Head__AMR__LUT.md) | [**getLutInfo**](classxThermal_1_1cxThermal.md#function-getlutinfo) (std::string file\_lut, [**bool**](classxThermal_1_1xThermalError.md) printSummary=[**true**](classxThermal_1_1xThermalError.md)) <br> |
| virtual [**const**](classxThermal_1_1xThermalError.md) std::map&lt; [**int**](classxThermal_1_1xThermalError.md), [**propInfo**](structxThermal_1_1propInfo.md) &gt; & | [**get\_UpdateWhichProps**](classxThermal_1_1cxThermal.md#function-get_updatewhichprops) () <br> |
|  [**int**](classxThermal_1_1xThermalError.md) | [**get\_dim\_lut**](classxThermal_1_1cxThermal.md#function-get_dim_lut) () <br> |
|  [**int**](classxThermal_1_1xThermalError.md) | [**get\_dim\_lut\_lookup**](classxThermal_1_1cxThermal.md#function-get_dim_lut_lookup) () <br> |
| virtual [**int**](classxThermal_1_1xThermalError.md) | [**get\_num\_threads**](classxThermal_1_1cxThermal.md#function-get_num_threads) () <br> |
|  [**void**](classxThermal_1_1xThermalError.md) \* | [**get\_pLUT**](classxThermal_1_1cxThermal.md#function-get_plut) () <br> |
|  [**void**](classxThermal_1_1xThermalError.md) \* | [**get\_pLUT\_lookup**](classxThermal_1_1cxThermal.md#function-get_plut_lookup) () <br> |
|  std::vector&lt; T &gt; | [**linspace**](classxThermal_1_1cxThermal.md#function-linspace-12) (T xmin, T xmax, std::size\_t n, [**bool**](classxThermal_1_1xThermalError.md) isLogScale=[**false**](classxThermal_1_1xThermalError.md)) <br>_Make a linearly spaced vector of points._  |
|  std::vector&lt; T &gt; | [**linspace**](classxThermal_1_1cxThermal.md#function-linspace-22) (T xmindxxmax) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**loadLUT**](classxThermal_1_1cxThermal.md#function-loadlut) ([**const**](classxThermal_1_1xThermalError.md) std::string & filename, [**bool**](classxThermal_1_1xThermalError.md) printStatus=[**true**](classxThermal_1_1xThermalError.md)) <br> |
|  [**LOOKUPTABLE\_FOREST::Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; 2, [**LOOKUPTABLE\_FOREST::FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md)&lt; 2 &gt; &gt; \* | [**lookup**](classxThermal_1_1cxThermal.md#function-lookup-14) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & prop, [**double**](classxThermal_1_1xThermalError.md) x, [**double**](classxThermal_1_1xThermalError.md) y) <br> |
|  [**LOOKUPTABLE\_FOREST::Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; 2, [**LOOKUPTABLE\_FOREST::FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md)&lt; 2 &gt; &gt; \* | [**lookup**](classxThermal_1_1cxThermal.md#function-lookup-24) ([**double**](classxThermal_1_1xThermalError.md) \* props, [**double**](classxThermal_1_1xThermalError.md) \* xyz\_min\_target, [**double**](classxThermal_1_1xThermalError.md) x, [**double**](classxThermal_1_1xThermalError.md) y, [**bool**](classxThermal_1_1xThermalError.md) is\_cal=[**true**](classxThermal_1_1xThermalError.md)) <br> |
|  [**LOOKUPTABLE\_FOREST::Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; 3, [**LOOKUPTABLE\_FOREST::FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md)&lt; 3 &gt; &gt; \* | [**lookup**](classxThermal_1_1cxThermal.md#function-lookup-34) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & prop, [**double**](classxThermal_1_1xThermalError.md) x, [**double**](classxThermal_1_1xThermalError.md) y, [**double**](classxThermal_1_1xThermalError.md) z) <br> |
|  [**LOOKUPTABLE\_FOREST::Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; 3, [**LOOKUPTABLE\_FOREST::FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md)&lt; 3 &gt; &gt; \* | [**lookup**](classxThermal_1_1cxThermal.md#function-lookup-44) ([**double**](classxThermal_1_1xThermalError.md) \* props, [**double**](classxThermal_1_1xThermalError.md) \* xyz\_min\_target, [**double**](classxThermal_1_1xThermalError.md) x, [**double**](classxThermal_1_1xThermalError.md) y, [**double**](classxThermal_1_1xThermalError.md) z, [**bool**](classxThermal_1_1xThermalError.md) is\_cal=[**true**](classxThermal_1_1xThermalError.md)) <br> |
|  [**LOOKUPTABLE\_FOREST::Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; 2, [**LOOKUPTABLE\_FOREST::FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md)&lt; 2 &gt; &gt; \* | [**lookup\_only**](classxThermal_1_1cxThermal.md#function-lookup_only-12) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & prop, [**double**](classxThermal_1_1xThermalError.md) x, [**double**](classxThermal_1_1xThermalError.md) y) <br> |
|  [**LOOKUPTABLE\_FOREST::Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; 3, [**LOOKUPTABLE\_FOREST::FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md)&lt; 3 &gt; &gt; \* | [**lookup\_only**](classxThermal_1_1cxThermal.md#function-lookup_only-22) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & prop, [**double**](classxThermal_1_1xThermalError.md) x, [**double**](classxThermal_1_1xThermalError.md) y, [**double**](classxThermal_1_1xThermalError.md) z) <br> |
|  T | [**max\_vector**](classxThermal_1_1cxThermal.md#function-max_vector) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; T &gt; & x) <br> |
|  T | [**mean\_vector**](classxThermal_1_1cxThermal.md#function-mean_vector) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; T &gt; & x) <br> |
|  std::vector&lt; [**size\_t**](classxThermal_1_1xThermalError.md) &gt; | [**meshgrid**](classxThermal_1_1cxThermal.md#function-meshgrid) (T xmindxxmax, T ymindyymax, std::vector&lt; T &gt; & XX, std::vector&lt; T &gt; & YY) <br> |
|  T | [**min\_vector**](classxThermal_1_1cxThermal.md#function-min_vector) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; T &gt; & x) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**molar\_mass**](classxThermal_1_1cxThermal.md#function-molar_mass) () = 0<br> |
| virtual std::string | [**name**](classxThermal_1_1cxThermal.md#function-name) () = 0<br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**normalizePhaseBoundaries**](classxThermal_1_1cxThermal.md#function-normalizephaseboundaries) ([**PhaseBoundaries**](structxThermal_1_1PhaseBoundaries.md) & phaseBoundaries) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**p\_critical**](classxThermal_1_1cxThermal.md#function-p_critical) () = 0<br> |
| virtual std::string | [**phase\_name**](classxThermal_1_1cxThermal.md#function-phase_name) ([**PhaseRegion**](namespacexThermal.md#enum-phaseregion) phase\_index) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**pmax**](classxThermal_1_1cxThermal.md#function-pmax) () = 0<br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**pmin**](classxThermal_1_1cxThermal.md#function-pmin) () = 0<br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**rhomass\_critical**](classxThermal_1_1cxThermal.md#function-rhomass_critical) () = 0<br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**rhomolar\_critical**](classxThermal_1_1cxThermal.md#function-rhomolar_critical) () <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**save\_lut\_to\_binary**](classxThermal_1_1cxThermal.md#function-save_lut_to_binary) (std::string filename) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**save\_lut\_to\_vtk**](classxThermal_1_1cxThermal.md#function-save_lut_to_vtk) (std::string filename, [**bool**](classxThermal_1_1xThermalError.md) isNormalizeXYZ=[**true**](classxThermal_1_1xThermalError.md)) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**set\_num\_threads**](classxThermal_1_1cxThermal.md#function-set_num_threads) ([**int**](classxThermal_1_1xThermalError.md) num\_threads) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**showProgressBar**](classxThermal_1_1cxThermal.md#function-showprogressbar) ([**bool**](classxThermal_1_1xThermalError.md) isShow) <br> |
|  T | [**sum\_vector**](classxThermal_1_1cxThermal.md#function-sum_vector) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; T &gt; & x) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**writeLine2VTU**](classxThermal_1_1cxThermal.md#function-writeline2vtu) (std::string vtuFile, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & X, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & Y, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & Z, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) scale\_x=1.0, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) scale\_y=1.0, [**double**](classxThermal_1_1xThermalError.md) scale\_z=1.0) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**writeMeshGrid2VTK**](classxThermal_1_1cxThermal.md#function-writemeshgrid2vtk) ([**const**](classxThermal_1_1xThermalError.md) std::string & vtkFile, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & x, [**const**](classxThermal_1_1xThermalError.md) std::string & xTitle, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & y, [**const**](classxThermal_1_1xThermalError.md) std::string & yTitle, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & z, [**const**](classxThermal_1_1xThermalError.md) std::string & zTitle, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; &gt; & props, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**propInfo**](structxThermal_1_1propInfo.md) &gt; & propsInfo, [**bool**](classxThermal_1_1xThermalError.md) isNormalize=[**false**](classxThermal_1_1xThermalError.md)) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**writePhaseBoundaries2VTU**](classxThermal_1_1cxThermal.md#function-writephaseboundaries2vtu) (std::string outputPath, [**const**](classxThermal_1_1xThermalError.md) [**PhaseBoundaries**](structxThermal_1_1PhaseBoundaries.md) & phaseBoundaries, [**double**](classxThermal_1_1xThermalError.md) scale\_T=1, [**double**](classxThermal_1_1xThermalError.md) scale\_p=1, [**double**](classxThermal_1_1xThermalError.md) scale\_X=1) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**writeTriMesh2Txt**](classxThermal_1_1cxThermal.md#function-writetrimesh2txt) ([**const**](classxThermal_1_1xThermalError.md) [**TriMesh**](structxThermal_1_1TriMesh.md) & mesh, std::string path\_out=".") <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**writeXXYYZZ2VTU**](classxThermal_1_1cxThermal.md#function-writexxyyzz2vtu) (std::string vtuFile, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; &gt; & XX, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; &gt; & YY, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; &gt; & ZZ, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) scale\_x=1.0, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) scale\_y=1.0, [**double**](classxThermal_1_1xThermalError.md) scale\_z=1.0) <br> |
| virtual  | [**~cxThermal**](classxThermal_1_1cxThermal.md#function-cxthermal) () <br> |




## Public Static Functions inherited from xThermal::cxThermal

See [xThermal::cxThermal](classxThermal_1_1cxThermal.md)

| Type | Name |
| ---: | :--- |
|  INDEX\_FLUID | [**validateFluid**](classxThermal_1_1cxThermal.md#function-validatefluid) ([**const**](classxThermal_1_1xThermalError.md) std::string & fluidName) <br> |


















































## Detailed Description


![Image](H2ONaCl/H2ONaCl_isothermal_PX.svg)
 


    
## Public Functions Documentation




### function D\_Tstar\_V [2/2]

```C++
std::vector< double > xThermal::H2ONaCl::cH2ONaCl::D_Tstar_V (
    std::vector< double > T,
    std::vector< double > P,
    std::vector< double > X
) 
```




<hr>



### function H\_phase [1/4]

_Calculate specific enthalpy of liquid or vapor phase of H2O-NaCl._ 
```C++
void xThermal::H2ONaCl::cH2ONaCl::H_phase (
    const  double & T,
    const  double & P,
    const  double & X,
    double & H,
    PhaseType phase
) 
```



For high-temperature(T&gt;800 deg.C) and low pressure (P&lt;P\_vlh), implement the similar extrapolation method with that in \like Rho\_phase .




**Warning:**

The low-T low-P(P&lt;1bar) extrapolation doesn't work for VLH surface, actually we don't need that because this part of VLH surface already below the minimum valid pressure (1bar).




**Parameters:**


* `T` [K] 
* `P` [Pa] 
* `X` [kg/kg] 
* `H` [J/kg] 
* `phase` [Liquid, Vapor] 




        

<hr>



### function H\_phase [2/4]

_Calculate specific enthalpy of liquid or vapor phase of H2O-NaCl._ 
```C++
void xThermal::H2ONaCl::cH2ONaCl::H_phase (
    const  double & T,
    const  double & P,
    const  double & X,
    double & H,
    double & Cp,
    PhaseType phase
) 
```



For high-temperature(T&gt;800 deg.C) and low pressure (P&lt;P\_vlh), implement the similar extrapolation method with that in \like Rho\_phase .



\[\begin{align}
\begin{cases}
H_{solution} (T,P,X) &= H_{H_2O}(T^*_H, P)\\
H_{extrapol}^{lowT} (T,P,X) &= o_0 + o_1 T^*_V + o_2 T^{*3}_{H},& ~~  (P<P_{H_2O}^{crit}, T<=T_{H_2O}^{crit})\\
H_{extrapol}^{highT} (T,P,X) &= o_3 + o_4 ln^{P + 1000} + o_5 P,& ~~ (P<P_{vlh}^{peak}, T>600 ^{\circ}C)
\end{cases}
\end{align}\]





**Parameters:**


* `T` 
* `P` 
* `X` 
* `H` 
* `Cp` 
* `phase` 




        

<hr>



### function H\_phase [3/4]

```C++
std::vector< double > xThermal::H2ONaCl::cH2ONaCl::H_phase (
    std::vector< double > T,
    std::vector< double > P,
    std::vector< double > X,
    PhaseType phase
) 
```




<hr>



### function H\_phase [4/4]

```C++
inline double xThermal::H2ONaCl::cH2ONaCl::H_phase (
    const  double & T,
    const  double & P,
    const  double & X,
    PhaseType phase
) 
```




<hr>



### function HminHmax\_VLH\_triangle 

_Calculate minimum and maximum specific enthalpy of a specific VLH "triangle" for given bulk salinity X. The "VLH triangle" is boundary (edge) of a three-phase region at constant P slice. The max and min H will be calculated for the latent heat of evaporation, see also grant1979compressibility._ 
```C++
void xThermal::H2ONaCl::cH2ONaCl::HminHmax_VLH_triangle (
    const  double & H_v0,
    const  double & H_l0,
    const  double & H_h0,
    const  double & X_v0,
    const  double & X_l0,
    const  double & X0,
    double & Hmin0,
    double & Hmax0
) 
```





**Parameters:**


* `H_v` 
* `H_l` 
* `H_h` 
* `X_v` 
* `X_l` 




        

<hr>



### function Mol2Wt [1/2]

_Convert molar fraction of NaCl to mass fraction._ 
```C++
double xThermal::H2ONaCl::cH2ONaCl::Mol2Wt (
    double X_mol
) 
```




\[\begin{align}
X_{wt} = \frac{X_{mol}M_{NaCl}}{X_{mol}M_{NaCl} + (1-X_{mol})M_{H_2O}}
\end{align}\]





**Parameters:**


* `X_mol` [0,1] 



**Returns:**

double [0,1] 





        

<hr>



### function Mol2Wt [2/2]

```C++
std::vector< double > xThermal::H2ONaCl::cH2ONaCl::Mol2Wt (
    std::vector< double > X_mol
) 
```




<hr>



### function Mu\_phase 

```C++
void xThermal::H2ONaCl::cH2ONaCl::Mu_phase (
    const  double & T,
    const  double & P,
    const  double & X,
    double & Mu,
    PhaseType phase
) 
```



Calculate dynamic viscosity of H2O-NaCl for given T,P,X. See klyukin2017revised equation (3,4,5,6).




**Warning:**

For low temperature case, the T\_star will be negative, e.g. T = 1 + 273.15, P = 200E5, X=0.1; see also Fig. 7a of klyukin2017revised. Use minimum T ?




**Parameters:**


* `T` [K] 
* `P` [Pa] 
* `X` [kg/kg] 
* `Mu` [Pa s] 
* `phase` [Liquid\|Vapor] 




        

<hr>



### function P\_Critical [1/2]

_Critical pressure as a function of temperature, see Eq. 5a-c of Driesner2007Part1._ 
```C++
void xThermal::H2ONaCl::cH2ONaCl::P_Critical (
    double T,
    double & P_crit
) 
```





**Note:**

Change 500 deg.C to 773.15 K to make temperature with unit [K].




**Warning:**

The critical pressure ( \(2.2054915\times 10^2\) bar) and temperature(373.976 \(^{\circ}\)C = 647.126K) of H2O used in Driesner & Heinrich(2007) formula comes from IAPS-84 EOS(see PROST code prost), which is different from that of IAPWS95\_CoolProp (m\_constants\_Water.T\_critical, H2O::P\_c). The difference between them are \(\Delta T=0.0299 ^{\circ}\)C, \(\Delta P = 9085\) Pa.




**Parameters:**


* `T` [K] 
* `P_crit` [Pa] 




        

<hr>



### function P\_Critical [2/2]

```C++
void xThermal::H2ONaCl::cH2ONaCl::P_Critical (
    std::vector< double > T,
    std::vector< double > & P_crit
) 
```




<hr>



### function P\_VLH [1/3]

_Vapor-Liquid-Halite coexistence pressure as a function of temperature._ 
```C++
void xThermal::H2ONaCl::cH2ONaCl::P_VLH (
    const  double & T,
    double & P
) 
```



See Eq. 10 of Driesner2007Part1.


![Image](H2ONaCl/H2ONaCl_P_VLH.svg)





**See also:** Fig. 9,10 of Driesner2007Part1. [**X\_HaliteLiquidus**](classxThermal_1_1H2ONaCl_1_1cH2ONaCl.md#function-x_haliteliquidus-12) and [**X\_VLH**](classxThermal_1_1H2ONaCl_1_1cH2ONaCl.md#function-x_vlh-12)


**Parameters:**


* `T` [K] 
* `P` [Pa] 




        

<hr>



### function P\_VLH [2/3]

```C++
double xThermal::H2ONaCl::cH2ONaCl::P_VLH (
    const  double & T
) 
```




<hr>



### function P\_VLH [3/3]

```C++
void xThermal::H2ONaCl::cH2ONaCl::P_VLH (
    std::vector< double > T,
    std::vector< double > & P
) 
```




<hr>



### function P\_X\_Critical [1/2]

_Calculate critical pressure(equation 5a-c of Driesner2007Part1) and salinity(equation 7a,b of Driesner2007Part1) by given temperature._ 
```C++
void xThermal::H2ONaCl::cH2ONaCl::P_X_Critical (
    double T,
    double & P_crit,
    double & X_crit
) 
```




* Benchmark test




If use the critical point of H2O in Driesner(2007a)Driesner2007Part1, this function will get completely the same result with [Electronic Annex EA-2](https://www.sciencedirect.com/science/article/pii/S0016703707002955?via%3Dihub#aep-e-component-id41) published with Driesner2007Part1.


![Image](H2ONaCl/H2ONaCl_CriticalCurve.svg)





**Parameters:**


* `T` [K] 
* `P_crit` [Pa] 
* `X_crit` [kg/kg] 




        

<hr>



### function P\_X\_Critical [2/2]

```C++
void xThermal::H2ONaCl::cH2ONaCl::P_X_Critical (
    std::vector< double > T,
    std::vector< double > & P_crit,
    std::vector< double > & X_crit
) 
```




<hr>



### function PhaseBoundary\_HaliteLiquidus 

```C++
TriMesh xThermal::H2ONaCl::cH2ONaCl::PhaseBoundary_HaliteLiquidus (
    std::string fmt_out="",
    double p_max=2500E5,
    double dT=10,
    double dp=50E5
) 
```



Calculate and return the phase boundary triangle mesh.




**Parameters:**


* `fmt_out` [txt\|vtk] 



**Returns:**







        

<hr>



### function PhaseBoundary\_HaliteLiquidus\_DeformLinear 

```C++
DeformLinearMesh xThermal::H2ONaCl::cH2ONaCl::PhaseBoundary_HaliteLiquidus_DeformLinear (
    double p_max=2500E5,
    double dT=10,
    double dp=50E5
) 
```



Create phase boundary of halite liquidus in a deformed linear "structured" mesh. 

**Parameters:**


* `fmt_out` 
* `p_max` 
* `dT` 
* `dp` 



**Returns:**







        

<hr>



### function PhaseBoundary\_VH\_DeformLinear 

```C++
DeformLinearMesh xThermal::H2ONaCl::cH2ONaCl::PhaseBoundary_VH_DeformLinear (
    int nT=100,
    int np=100
) 
```




<hr>



### function PhaseBoundary\_VLH\_DeformLinear 

```C++
DeformLinearMesh xThermal::H2ONaCl::cH2ONaCl::PhaseBoundary_VLH_DeformLinear (
    int nT=100,
    int nX=60
) 
```




<hr>



### function PhaseBoundary\_VL\_DeformLinear 

```C++
DeformLinearMesh xThermal::H2ONaCl::cH2ONaCl::PhaseBoundary_VL_DeformLinear (
    PhaseType VaporOrLiquid,
    int nT=100,
    int np=200
) 
```




<hr>



### function Rho\_phase [1/3]

_Calculate density of liquid or vapor phase._ 
```C++
void xThermal::H2ONaCl::cH2ONaCl::Rho_phase (
    const  double & T,
    const  double & P,
    const  double & X,
    double & rho,
    double & dRhodP,
    double & dRhodT,
    PhaseType phase
) 
```



The density of H2ONaCl can be calculated from water density, see Eq. 7 of Driesner2007Part2. \(V_{solution}(P,T,X) = V_{water}(T^*,P)\), so \(\frac{M_{solution}}{\rho_{solution}} = \frac{M_{water}}{\rho_{water}}\).


Therefore, \(\rho_{solution} = \frac{M_{solution}}{M_{H_2O}} \rho_{H_2O}\), where \(M_{solution} = M_{H_2O}(1-X_{mol}) + M_{NaCl}X_{mol}\)




**Note:**

Density difference calculated from IAPWS95\_CoolProp and IAPS84 will cause H2ONaCl density difference up to 32 kg/m^3, test the following points:



|T[c]   |P[bar]   |X[mole fraction]   |rho[Driesner2007Part2]   |rho[[**xThermal**](namespacexThermal.md)]    |
|-----|-----|-----|-----|-----|
|8.500000e+02   |2.000000e+03   |9.996030e-01   |1.608149e+03   |1.602652e+03    |
|8.600000e+02   |2.500000e+03   |9.937537e-01   |1.615147e+03   |1.605205e+03    |
|9.200000e+02   |5.000000e+03   |9.915095e-01   |1.651829e+03   |1.619722e+03    |
|9.243300e+02   |5.000000e+03   |1.000000e+00   |1.656120e+03   |1.623869e+03   |






For example, T=920 deg.C, P = 5000 bar. The calculated \(T^*\) for density is 1279.92 deg.C. Density of water at (T=1279.92 deg.C, P=5000 bar) is 502.235 kg/m^3 (IAPWS95\_CoolProp), 512.153 kg/m^3 (IAPS84). The coefficient \(\frac{M_{H_2O}(1-X_{mol}) + M_{NaCl}X_{mol}}{M_{H_2O}} = 3.225\) at this P-T-X condition, therefore the final density different will be \(3.225\times(512.153-502.235)=32\) (This P-T condition exceeds the valid range of IAPWS-IF97).




**Todo**

Estimate density difference range based on the density difference between IAPWS95\_CoolProp and IAPS84.






**Todo**

How to deal with properties at high-T low-P region, e.g. sw.Rho\_phase(1.700000e+02+273.15,7.914706e+00\*1E5,0,Rhol, H2ONaCl::Liquid); see mmc4 of Driesner(2007b).






**Todo**

Merge Rho\_phase of liquid and vapor together to calculate them simultaneously, because property of saturated vapor and liquid always needed simultaneous, it could save computing time in this way. 



**Bug**

T=273.16, P=20E5, X=0.1, Tstar = n1 + n2\*T\_C + D will be negative.






**Parameters:**


* `T` [K] 
* `P` [Pa] 
* `X` [kg/kg] salinity of liquid or salinity of vapor, not bulk salinity 
* `Tstar` [K] 
* `rho` [ \(kg/m^3\)] 




        

<hr>



### function Rho\_phase [2/3]

```C++
std::vector< double > xThermal::H2ONaCl::cH2ONaCl::Rho_phase (
    std::vector< double > T,
    std::vector< double > P,
    std::vector< double > X,
    PhaseType phase
) 
```




<hr>



### function Rho\_phase [3/3]

```C++
inline double xThermal::H2ONaCl::cH2ONaCl::Rho_phase (
    const  double & T,
    const  double & P,
    const  double & X,
    PhaseType phase
) 
```




<hr>



### function Slice\_constP 

```C++
PhaseRegion_Slice xThermal::H2ONaCl::cH2ONaCl::Slice_constP (
    const  double P0,
    size_t nPoints=500
) 
```



Calculate polygons of phase region at constant P. 

**Parameters:**


* `P0` [Pa] 



**Returns:**







        

<hr>



### function Slice\_constT 

```C++
PhaseRegion_Slice xThermal::H2ONaCl::cH2ONaCl::Slice_constT (
    const  double T0,
    size_t nPoints=500
) 
```




<hr>



### function T\_Critical [1/2]

_Solve critical T and X using bisection method for given Temperature._ 
```C++
void xThermal::H2ONaCl::cH2ONaCl::T_Critical (
    double P,
    double & T_crit
) 
```





**Parameters:**


* `P` [Pa] 
* `T_crit` [K] 




        

<hr>



### function T\_Critical [2/2]

```C++
void xThermal::H2ONaCl::cH2ONaCl::T_Critical (
    std::vector< double > P,
    std::vector< double > & T_crit
) 
```




<hr>



### function T\_HPX 

```C++
double xThermal::H2ONaCl::cH2ONaCl::T_HPX (
    const  double & H,
    const  double & p,
    const  double & X,
    const  double & T_low=T_MIN,
    const  double & T_high=T_MAX
) 
```




<hr>



### function T\_VL 

```C++
double xThermal::H2ONaCl::cH2ONaCl::T_VL (
    const  double & P,
    const  double & X,
    PhaseType phase
) 
```



Calculate temperature of liquid branch for given P and X. It is used for HPX calculation in the low salinity region and P&lt;Pcrit\_H2O 

**Parameters:**


* `P` [Pa] 
* `X` [kg/kg] 



**Returns:**







        

<hr>



### function T\_VLH\_P0 

_Calculate temperature boundary of VLH zone by given pressure P0._ 
```C++
void xThermal::H2ONaCl::cH2ONaCl::T_VLH_P0 (
    double P0,
    double & T_min,
    double & T_max
) 
```





**Warning:**

Only valid in pressure range [[**H2ONaCl::P\_MIN**](group__PhysicalConstants__H2ONaCl.md#variable-p_min), [**H2ONaCl::P\_Peak\_VLH**](group__PhysicalConstants__H2ONaCl.md#variable-p_peak_vlh)]




**Note:**

The pressure value should be checked in the function, for example to avoid program crash happening when P&gt;P\_Peak\_VLH, while the if statement check could decrease program performance if this function is called in large for loops. Therefore, please be careful the valid pressure range!




**Parameters:**


* `P0` [Pa] 
* `Tmin` [K] 
* `Tmax` [K] 




        

<hr>



### function T\_VL\_L 

```C++
double xThermal::H2ONaCl::cH2ONaCl::T_VL_L (
    const  double & p,
    const  double & X,
    const  double & T_low,
    const  double & T_high
) 
```



Calculate temperature on surface VL-liquid branch for given p, X. Used in HPX calculation. 

**Parameters:**


* `p` [Pa] 
* `X` [kg/kg] 



**Returns:**







        

<hr>



### function T\_VL\_V 

```C++
double xThermal::H2ONaCl::cH2ONaCl::T_VL_V (
    const  double & p,
    const  double & X,
    const  double & T_low,
    const  double & T_high
) 
```



Inverse T at VL-&gt;Vapor branch surface for give p,X 

**Parameters:**


* `p` [Pa] 
* `X` [kg/kg] 
* `Tmin` [K] 
* `Tmax` [K] 



**Returns:**







        

<hr>



### function T\_X\_Critical [1/2]

_Calculate critical temperature and composition by given pressure._ 
```C++
void xThermal::H2ONaCl::cH2ONaCl::T_X_Critical (
    double P,
    double & T_crit,
    double & X_crit
) 
```



This is the inversion of [**P\_X\_Critical**](classxThermal_1_1H2ONaCl_1_1cH2ONaCl.md#function-p_x_critical-12)




**Parameters:**


* `P` [Pa] 
* `T_crit` [K] 
* `X_crit` [kg/kg] 




        

<hr>



### function T\_X\_Critical [2/2]

```C++
void xThermal::H2ONaCl::cH2ONaCl::T_X_Critical (
    std::vector< double > P,
    std::vector< double > & T_crit,
    std::vector< double > & X_crit
) 
```




<hr>



### function T\_critical 

```C++
inline virtual double xThermal::H2ONaCl::cH2ONaCl::T_critical () 
```



Get the triple point temperature in K 


        
Implements [*xThermal::cxThermal::T\_critical*](classxThermal_1_1cxThermal.md#function-t_critical)


<hr>



### function Tmax 

```C++
inline virtual double xThermal::H2ONaCl::cH2ONaCl::Tmax () 
```



Get the minimum temperature in K. Same as H2O 


        
Implements [*xThermal::cxThermal::Tmax*](classxThermal_1_1cxThermal.md#function-tmax)


<hr>



### function Tmax\_VLH 

_Calculate the maximum temperature of the V+L+H zone._ 
```C++
double xThermal::H2ONaCl::cH2ONaCl::Tmax_VLH () 
```



It should be the root of \(P_{LVH} = P_{min}\), see also Eq. 10 of Driesner2007Part1 and [**P\_VLH**](classxThermal_1_1H2ONaCl_1_1cH2ONaCl.md#function-p_vlh-13)




**Note:**

This function is only used to calculate the constant value.




**Returns:**

double [K] 





        

<hr>



### function Tmin 

```C++
inline virtual double xThermal::H2ONaCl::cH2ONaCl::Tmin () 
```



Get the minimum temperature in K 


        
Implements [*xThermal::cxThermal::Tmin*](classxThermal_1_1cxThermal.md#function-tmin)


<hr>



### function Tmin\_VLH 

```C++
double xThermal::H2ONaCl::cH2ONaCl::Tmin_VLH () 
```




<hr>



### function Triangulation 

```C++
TriMesh xThermal::H2ONaCl::cH2ONaCl::Triangulation (
    const std::vector< double > & x_poly,
    const std::vector< double > & y_poly,
    double xIn,
    double yIn,
    double dx,
    double dy
) 
```




<hr>



### function Ttriple 

```C++
inline virtual double xThermal::H2ONaCl::cH2ONaCl::Ttriple () 
```



Get the maximum pressure in Pa 


        
Implements [*xThermal::cxThermal::Ttriple*](classxThermal_1_1cxThermal.md#function-ttriple)


<hr>



### function UpdateState\_HPX 

```C++
virtual void xThermal::H2ONaCl::cH2ONaCl::UpdateState_HPX (
    ThermodynamicProperties & props,
    const  double & H,
    const  double & p,
    const  double & X
) 
```



Implements [*xThermal::cxThermal::UpdateState\_HPX*](classxThermal_1_1cxThermal.md#function-updatestate_hpx-14)


<hr>



### function UpdateState\_HPX\_vl\_lowXlowP 

```C++
bool xThermal::H2ONaCl::cH2ONaCl::UpdateState_HPX_vl_lowXlowP (
    ThermodynamicProperties & props,
    const  double & H,
    const  double & p,
    const  double & X
) 
```



Calculate T and thermodynamic properties in V+L zone in the low salinity (&lt;0.01 kg/kg) and low pressure (&lt;P\_crit\_h2o) zone.


There are some issues in low salinity region and close VL boundary for T inversion, therefore we need to do some special process to accurately invert T. Note that below T\_crit\_h2o, the boundary of VL is an isothermal boundary, i.e., for given P, the T is not changed with X, and the T is the boiling T of H2O for given P




**Parameters:**


* `props` 
* `H` [J/kg] 
* `p` [Pa] 
* `X` [kg/kg] 



**Returns:**







        

<hr>



### function UpdateState\_HPX\_vlh 

```C++
bool xThermal::H2ONaCl::cH2ONaCl::UpdateState_HPX_vlh (
    ThermodynamicProperties & props,
    const  double & H,
    const  double & p,
    const  double & X,
    const std::vector< double > & T1_T2
) 
```



Calculate thermodynamic properties and T for given H,P,X in VLH zone. 

**Parameters:**


* `props` 
* `H` [J/kg] 
* `p` [Pa] 
* `X` [kg/kg] 
* `T1_T2` [K] Temperature of V+L+H surface for given P 



**Returns:**







        

<hr>



### function UpdateState\_TPX 

_Update state using T, p as independent variables._ 
```C++
virtual void xThermal::H2ONaCl::cH2ONaCl::UpdateState_TPX (
    ThermodynamicProperties & props,
    const  double & T,
    const  double & p,
    const  double & X
) 
```



Calculate the critical molar density in \(mol/m^3\) 




**Parameters:**


* `T` [K] 
* `p` [Pa] 




        
Implements [*xThermal::cxThermal::UpdateState\_TPX*](classxThermal_1_1cxThermal.md#function-updatestate_tpx-14)


<hr>



### function Wt2Mol [1/3]

_Convert mass fraction of NaCl to molar fraction._ 
```C++
double xThermal::H2ONaCl::cH2ONaCl::Wt2Mol (
    const  double & X_wt
) 
```




\[\begin{align}
X_{mol} = \frac{X_{wt}/M_{NaCl}}{X_{wt}/M_{NaCl} + (1-X_{wt})/M_{H_2O}}
\end{align}\]





**Parameters:**


* `X_wt` [0,1] 



**Returns:**

double [0,1] 





        

<hr>



### function Wt2Mol [2/3]

```C++
void xThermal::H2ONaCl::cH2ONaCl::Wt2Mol (
    const  double & X_wt,
    double & X_mol
) 
```




<hr>



### function Wt2Mol [3/3]

```C++
std::vector< double > xThermal::H2ONaCl::cH2ONaCl::Wt2Mol (
    std::vector< double > X_wt
) 
```




<hr>



### function XL\_VL [1/2]

_Calculate liquid composition on liquid branch of V+L coexistence._ 
```C++
double xThermal::H2ONaCl::cH2ONaCl::XL_VL (
    const  double & T,
    const  double & P
) 
```





**See also:** XL\_VL\_mol
![Image](H2ONaCl/H2ONaCl_X_VL.svg)





**See also:** Fig. 12 of Driesner2007Part1.


**Parameters:**


* `T` [K] 
* `P` [Pa] 



**Returns:**

double [kg/kg] 





        

<hr>



### function XL\_VL [2/2]

```C++
void xThermal::H2ONaCl::cH2ONaCl::XL_VL (
    std::vector< double > T,
    std::vector< double > P,
    std::vector< double > & res
) 
```




<hr>



### function XV\_VL [1/2]

_Calculate vapor composition on vapor branch of V+L coexistence._ 
```C++
double xThermal::H2ONaCl::cH2ONaCl::XV_VL (
    const  double & T,
    const  double & P
) 
```



![Image](H2ONaCl/H2ONaCl_X_VL.svg)





**See also:** Fig. 12 of Driesner2007Part1.


**Parameters:**


* `T` [K] 
* `P` [Pa] 



**Returns:**

double [kg/kg] 





        

<hr>



### function XV\_VL [2/2]

```C++
void xThermal::H2ONaCl::cH2ONaCl::XV_VL (
    std::vector< double > T,
    std::vector< double > P,
    std::vector< double > & res
) 
```




<hr>



### function X\_Critical 

_Critical composition as a function of temperature, see Eq. 7a,b of Driesner2007Part1._ 
```C++
void xThermal::H2ONaCl::cH2ONaCl::X_Critical (
    double T,
    double & X_crit
) 
```





**Parameters:**


* `T` [K] 
* `X_crit` [kg/kg] 




        

<hr>



### function X\_HaliteLiquidus [1/2]

_Halite-saturated liquid composition, see equation 8 of Driesner2007Part1._ 
```C++
double xThermal::H2ONaCl::cH2ONaCl::X_HaliteLiquidus (
    const  double & T,
    const  double & P
) 
```





**Note:**

Change coefficient e to make pressure with unit [Pa]


![Image](H2ONaCl/HaliteLiquidus3D.svg)





**See also:** Fig.7a,b of Driesner2007Part1. [**P\_VLH**](classxThermal_1_1H2ONaCl_1_1cH2ONaCl.md#function-p_vlh-13)


**Parameters:**


* `T` [K] 
* `P` [Pa] 



**Returns:**

double [kg/kg] 





        

<hr>



### function X\_HaliteLiquidus [2/2]

```C++
void xThermal::H2ONaCl::cH2ONaCl::X_HaliteLiquidus (
    std::vector< double > T,
    std::vector< double > P,
    std::vector< double > & res
) 
```




<hr>



### function X\_VH [1/2]

_Halite-saturated vapor composition at V+H coexistence, see equation 9 and Fig. 8 of Driesner2007Part1._ 
```C++
double xThermal::H2ONaCl::cH2ONaCl::X_VH (
    const  double & T,
    const  double & P
) 
```



![Image](H2ONaCl/H2ONaCl_X_HaliteSatVapor.svg)





**Warning:**

The Fig.8a of Driesner(2007) looks shift 1E-1 to the left.




**Parameters:**


* `T` [K] 
* `P` [Pa] 



**Returns:**

double [kg/kg] 





        

<hr>



### function X\_VH [2/2]

```C++
void xThermal::H2ONaCl::cH2ONaCl::X_VH (
    std::vector< double > T,
    std::vector< double > P,
    std::vector< double > & res
) 
```




<hr>



### function X\_VL [1/2]

_Calculate_ \(X_{V}, X_{L}\) _on V+L surface in one step._
```C++
void xThermal::H2ONaCl::cH2ONaCl::X_VL (
    const  double & T,
    const  double & P,
    double & X_L,
    double & X_V
) 
```



Because [**XV\_VL**](classxThermal_1_1H2ONaCl_1_1cH2ONaCl.md#function-xv_vl-12) always need to call [**XL\_VL**](classxThermal_1_1H2ONaCl_1_1cH2ONaCl.md#function-xl_vl-12) first, so this function can save computing time if both \(X_{V}\) and \(X_{L}\) are needed.


![Image](H2ONaCl/H2ONaCl_X_VL.svg)





**See also:** Fig. 12 of Driesner2007Part1.


**Parameters:**


* `T` [K] 
* `P` [Pa] 
* `X_L` [kg/kg] 
* `X_V` [kg/kg] 




        

<hr>



### function X\_VL [2/2]

```C++
void xThermal::H2ONaCl::cH2ONaCl::X_VL (
    std::vector< double > T,
    std::vector< double > P,
    std::vector< double > & X_L,
    std::vector< double > & X_V
) 
```




<hr>



### function X\_VLH [1/2]

_Calculate salinity of vapor and liquid at VLH coexistence._ 
```C++
void xThermal::H2ONaCl::cH2ONaCl::X_VLH (
    const  double & T,
    const  double & P,
    double & X_L,
    double & X_V
) 
```



Use the same formula with halite-saturated vapor composition at VH coexistence, see [**X\_VH**](classxThermal_1_1H2ONaCl_1_1cH2ONaCl.md#function-x_vh-12) and X\_HaliteSat\_mol. Where X\_HaliteSat\_mol can calculate both \(X_l\) and \(X_v\) at halite-saturated domain.


![Image](H2ONaCl/H2ONaCl_X_VLH.svg)





**See also:** Fig. 10 and 11 of Driesner2007Part1.


**Parameters:**


* `T` [K] 
* `P` [Pa] 
* `X_L` [kg/kg] 
* `X_V` [kg/kg] 




        

<hr>



### function X\_VLH [2/2]

```C++
void xThermal::H2ONaCl::cH2ONaCl::X_VLH (
    std::vector< double > T,
    std::vector< double > P,
    std::vector< double > & X_L,
    std::vector< double > & X_V
) 
```




<hr>



### function cH2ONaCl [1/2]

```C++
xThermal::H2ONaCl::cH2ONaCl::cH2ONaCl (
    std::string backend_H2O
) 
```




<hr>



### function cH2ONaCl [2/2]

```C++
xThermal::H2ONaCl::cH2ONaCl::cH2ONaCl (
    const  cH2ONaCl & sw
) 
```




<hr>



### function calc\_PhaseBoundaries 

```C++
PhaseBoundaries xThermal::H2ONaCl::cH2ONaCl::calc_PhaseBoundaries (
    std::string scale_X="linear",
    double ratio_log_to_linear=1,
    double Xcenter=0.01
) 
```



Calculate all phase boundaries in the valid p-T-X range. 

**Parameters:**


* `scale_X` ["linear" \| "log" \| "loglinear"] 
* `ratio_log_to_linear` [length of log axis]/[length of linear axis] 
* `Xcenter` [kg/kg], X&lt;Xcenter use log scale; X&gt;=Xcenter use linear scale. It is only used with scale\_X="loglinear" 



**Returns:**







        

<hr>



### function compressibility\_LH 

_Calculate isothermal compressibility of L+H phase._ 
```C++
void xThermal::H2ONaCl::cH2ONaCl::compressibility_LH (
    const  double & T0,
    const  double & p0,
    const  double & X0,
    const  double & bulkRho,
    double & compressibility,
    double dp=1
) 
```





**Parameters:**


* `T0` 
* `p0` 
* `X0` 
* `compressibility` 
* `dp` 




        

<hr>



### function compressibility\_VH 

```C++
void xThermal::H2ONaCl::cH2ONaCl::compressibility_VH (
    const  double & T0,
    const  double & p0,
    const  double & X0,
    const  double & bulkRho,
    double & compressibility,
    double dp=1
) 
```




<hr>



### function compressibility\_VL 

_Calculate isothermal compressibility of two phase._ 
```C++
void xThermal::H2ONaCl::cH2ONaCl::compressibility_VL (
    const  double & T0,
    const  double & p0,
    const  double & X0,
    const  double & bulkRho,
    double & compressibility,
    double dp=1
) 
```



According to the definition of compressibility \(\beta = - \frac{1}{V} \left( \frac{\partial V}{\partial p} \right)_T\), the compressibility of two phase can be expressed as,  \(\begin{equation}
\beta = \frac{1}{\rho} \left( \frac{\partial p}{\partial p} \right)_{T,X}
\end{equation}\) where \(\rho = S_1 \rho_1 + S_2\rho_2\) is the bulk density. The analytical expression is a bit complicated, so we can calculated it using numerical difference.  \(\begin{equation}
\beta (T,P,X) = \frac{\rho(T,P + 0.5\delta p, X) - \rho (T, P - 0.5\delta P, X)}{\delta p}
\end{equation}\)




**Parameters:**


* `T` 
* `P` 
* `X` 
* `bulkRho` 
* `compressibility` 
* `dp` 




        

<hr>



### function dPdT\_VLH 

```C++
double xThermal::H2ONaCl::cH2ONaCl::dPdT_VLH (
    const  double & T
) 
```



Calculate derivative \(\frac{dP}{dT}\) along VLH saturation curve. 

**Parameters:**


* `T` 



**Returns:**







        

<hr>



### function findPhaseRegion\_TPX 

```C++
virtual PhaseRegion xThermal::H2ONaCl::cH2ONaCl::findPhaseRegion_TPX (
    const  double & T,
    const  double & p,
    const  double & X
) 
```



Return the molar mass in kg/mol 


        
Implements [*xThermal::cxThermal::findPhaseRegion\_TPX*](classxThermal_1_1cxThermal.md#function-findphaseregion_tpx)


<hr>



### function getPhaseRegion\_HPX 

```C++
void xThermal::H2ONaCl::cH2ONaCl::getPhaseRegion_HPX (
    const  double & H,
    const  double & p,
    const  double & X,
    xThermal::PhaseRegion & phase_region,
    double & T,
    double S_lvh
) 
```



Get phase region for given H,P,X. 

**Parameters:**


* `H` 
* `p` 
* `X` 
* `phase_region` 




        

<hr>



### function getPhaseRegion\_TPX [1/2]

_Calculate phase region index and composition of vapor and liquid phase if the phase exists._ 
```C++
void xThermal::H2ONaCl::cH2ONaCl::getPhaseRegion_TPX (
    const  double & T,
    const  double & P,
    const  double & X,
    xThermal::PhaseRegion & phase_region,
    double & X_V,
    double & X_L
) 
```





**Parameters:**


* `T` [K] 
* `P` [Pa] 
* `X` Bulk salinity. [kg/kg] 
* `phase_region` Phase region index 
* `X_V` Salinity of vapor phase if it exists. [kg/kg] 
* `X_L` Salinity of liquid phase if it exists. [kg/kg]

![Image](H2ONaCl/H2ONaCl_isothermal_PhaseRegion.svg)
 


        

<hr>



### function getPhaseRegion\_TPX [2/2]

```C++
std::vector< int > xThermal::H2ONaCl::cH2ONaCl::getPhaseRegion_TPX (
    std::vector< double > T,
    std::vector< double > P,
    std::vector< double > X,
    std::vector< double > & X_V,
    std::vector< double > & X_L
) 
```




<hr>



### function get\_pWater 

```C++
inline cxThermal * xThermal::H2ONaCl::cH2ONaCl::get_pWater () const
```




<hr>



### function molar\_mass 

```C++
inline virtual double xThermal::H2ONaCl::cH2ONaCl::molar_mass () 
```



Return the critical mass density in \(kg/m^3\) 


        
Implements [*xThermal::cxThermal::molar\_mass*](classxThermal_1_1cxThermal.md#function-molar_mass)


<hr>



### function n1n2\_Tstar\_V [2/2]

```C++
void xThermal::H2ONaCl::cH2ONaCl::n1n2_Tstar_V (
    std::vector< double > P,
    std::vector< double > X,
    std::vector< double > & n1,
    std::vector< double > & n2
) 
```




<hr>



### function name 

```C++
inline virtual std::string xThermal::H2ONaCl::cH2ONaCl::name () 
```



Name of the model 


        
Implements [*xThermal::cxThermal::name*](classxThermal_1_1cxThermal.md#function-name)


<hr>



### function name\_backend 

```C++
std::string xThermal::H2ONaCl::cH2ONaCl::name_backend () 
```




<hr>



### function p\_critical 

```C++
inline virtual double xThermal::H2ONaCl::cH2ONaCl::p_critical () 
```



Return the critical temperature in K 


        
Implements [*xThermal::cxThermal::p\_critical*](classxThermal_1_1cxThermal.md#function-p_critical)


<hr>



### function pmax 

```C++
inline virtual double xThermal::H2ONaCl::cH2ONaCl::pmax () 
```



Get the minimum pressure in Pa, see Driesner2007a 


        
Implements [*xThermal::cxThermal::pmax*](classxThermal_1_1cxThermal.md#function-pmax)


<hr>



### function pmin 

```C++
inline virtual double xThermal::H2ONaCl::cH2ONaCl::pmin () 
```



Get the maximum temperature in K, see Driesner2007a 


        
Implements [*xThermal::cxThermal::pmin*](classxThermal_1_1cxThermal.md#function-pmin)


<hr>



### function prop\_VL 

```C++
void xThermal::H2ONaCl::cH2ONaCl::prop_VL (
    const  double & T,
    const  double & P,
    const  double & X,
    ThermodynamicProperties & prop
) 
```




<hr>



### function q1q2\_Tstar\_H [2/2]

```C++
void xThermal::H2ONaCl::cH2ONaCl::q1q2_Tstar_H (
    std::vector< double > P,
    std::vector< double > X,
    std::vector< double > & q1,
    std::vector< double > & q2
) 
```




<hr>



### function rhomass\_critical 

```C++
inline virtual double xThermal::H2ONaCl::cH2ONaCl::rhomass_critical () 
```



Return the critical pressure in Pa 


        
Implements [*xThermal::cxThermal::rhomass\_critical*](classxThermal_1_1cxThermal.md#function-rhomass_critical)


<hr>



### function ~cH2ONaCl 

```C++
xThermal::H2ONaCl::cH2ONaCl::~cH2ONaCl () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2ONaCl/H2ONaCl.h`

