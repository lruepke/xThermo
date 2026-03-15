

# Class xThermal::IAPWS\_IF97::cIAPWS\_IF97



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**IAPWS\_IF97**](namespacexThermal_1_1IAPWS__IF97.md) **>** [**cIAPWS\_IF97**](classxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97.md)








Inherits the following classes: [xThermal::cxThermal](classxThermal_1_1cxThermal.md)
























## Public Attributes inherited from xThermal::cxThermal

See [xThermal::cxThermal](classxThermal_1_1cxThermal.md)

| Type | Name |
| ---: | :--- |
|  [**bool**](classxThermal_1_1xThermalError.md) | [**m\_isShowProgressBar**](classxThermal_1_1cxThermal.md#variable-m_isshowprogressbar)  <br> |






























## Public Functions

| Type | Name |
| ---: | :--- |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Backward\_T\_PH\_region1**](#function-backward_t_ph_region1) ([**double**](classxThermal_1_1xThermalError.md) P, [**double**](classxThermal_1_1xThermalError.md) H) <br>_Calculate temperature in region 1 by given P and H. See equation 11 of IF97._  |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Backward\_T\_PH\_region2a**](#function-backward_t_ph_region2a) ([**double**](classxThermal_1_1xThermalError.md) P, [**double**](classxThermal_1_1xThermalError.md) H) <br>_Calculate saturated temperature by given pressure. Valid pressure rante if [H2O::m\_constants.pmin, H2O::m\_constants.p\_critical] See Equation 31 and Table 34 of IF97._  |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Backward\_T\_PH\_region2b**](#function-backward_t_ph_region2b) ([**double**](classxThermal_1_1xThermalError.md) P, [**double**](classxThermal_1_1xThermalError.md) H) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Backward\_T\_PH\_region2c**](#function-backward_t_ph_region2c) ([**double**](classxThermal_1_1xThermalError.md) P, [**double**](classxThermal_1_1xThermalError.md) H) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Backward\_T\_PH\_region3a**](#function-backward_t_ph_region3a) ([**double**](classxThermal_1_1xThermalError.md) P, [**double**](classxThermal_1_1xThermalError.md) H) <br>_The backward equation_ \(T_{3a}(p,h)\) _for subregion 3a. See equation 2 of IF97-Region3._ |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Backward\_T\_PH\_region3b**](#function-backward_t_ph_region3b) ([**double**](classxThermal_1_1xThermalError.md) P, [**double**](classxThermal_1_1xThermalError.md) H) <br>_The backward equation_ \(T_{3b}(p,h)\) _for subregion 3a. See equation 3 of IF97-Region3._ |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_p**](#function-boiling_p) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Boundary\_region23\_P2T**](#function-boundary_region23_p2t) ([**double**](classxThermal_1_1xThermalError.md) P) <br>_The boundary between regions 2 and 3 (see Fig. 1) is defined by the following simple quadratic pressure-temperature relation, the B23-equation. See equation 5 of IF97._  |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Boundary\_region2b2c\_H2P**](#function-boundary_region2b2c_h2p) ([**double**](classxThermal_1_1xThermalError.md) H) <br>_Boundary between region 2b and 2c. See also_ [_**Boundary\_region2b2c\_P2H**_](classxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97.md#function-boundary_region2b2c_p2h) _._ |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Boundary\_region2b2c\_P2H**](#function-boundary_region2b2c_p2h) ([**double**](classxThermal_1_1xThermalError.md) P) <br>_Boundary between region 2b and 2c. In order to know whether the T(p,h) equation for subregion 2b or for subregion 2c has to be used for given values of p and h, a special correlation equation for the boundary between subregions 2b and 2c (which approximates s = 5.85 kJ/kg/K_ [_**CONST\_IF97\_S\_Region2b2c**_](group__CONSTS__IF97.md#define-const_if97_s_region2b2c) _) is needed; see Fig. 2. This boundary equation, called the B2bc-equation, is a simple quadratic pressure-enthalpy relation which reads Equation 20 of IF97._ |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Boundary\_region3ab\_H2P**](#function-boundary_region3ab_h2p-12) ([**double**](classxThermal_1_1xThermalError.md) H) <br>_Boundary between region3a and 3b. See Equation 10 and Table 17 of IF97-Region3._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**Boundary\_region3ab\_H2P**](#function-boundary_region3ab_h2p-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; H, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Boundary\_region3ab\_P2H**](#function-boundary_region3ab_p2h-12) ([**double**](classxThermal_1_1xThermalError.md) P) <br>_Boundary between region3a and 3b. See Equation 1 and Table 2 of IF97-Region3._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**Boundary\_region3ab\_P2H**](#function-boundary_region3ab_p2h-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**int**](classxThermal_1_1xThermalError.md) | [**GetRegion\_PH**](#function-getregion_ph-12) ([**double**](classxThermal_1_1xThermalError.md) P, [**double**](classxThermal_1_1xThermalError.md) H) <br>_Calculate region index by given pressure and specific enthalpy._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**GetRegion\_PH**](#function-getregion_ph-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; H, std::vector&lt; [**int**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**int**](classxThermal_1_1xThermalError.md) | [**GetRegion\_PT**](#function-getregion_pt-12) ([**double**](classxThermal_1_1xThermalError.md) P, [**double**](classxThermal_1_1xThermalError.md) T) <br>_Calculate region index by given pressure and temperature._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**GetRegion\_PT**](#function-getregion_pt-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T, std::vector&lt; [**int**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**H\_sat\_P**](#function-h_sat_p-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**double**](classxThermal_1_1xThermalError.md) & H\_l, [**double**](classxThermal_1_1xThermalError.md) & H\_v) <br>_Calculate saturated liquid enthalpy and saturated vapor enthalpy for a given P._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**H\_sat\_P**](#function-h_sat_p-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & H\_l, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & H\_v) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**P\_sat\_T**](#function-p_sat_t-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T) <br>_Calculate saturated temperature by given pressure. Valid pressure rante if [H2O::T\_MIN, H2O::T\_c] See Equation 30 and Table 34 of IF97._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**P\_sat\_T**](#function-p_sat_t-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Prop\_Region1**](#function-prop_region1) ([**const**](classxThermal_1_1xThermalError.md) [**State\_Region1**](structxThermal_1_1IAPWS__IF97_1_1State__Region1.md) & state, BasicThermodynamicProperties which) <br>_Calculate thermodynamic properties by given dimensionless Gibbs free energy and its derivatives for the region 1._  |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Prop\_Region2**](#function-prop_region2) ([**State\_Region2**](structxThermal_1_1IAPWS__IF97_1_1State__Region2.md) state, BasicThermodynamicProperties which) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Prop\_Region5**](#function-prop_region5) ([**double**](classxThermal_1_1xThermalError.md) P, [**double**](classxThermal_1_1xThermalError.md) T, [**State\_Region5**](structxThermal_1_1IAPWS__IF97_1_1State__Region5.md) state, BasicThermodynamicProperties which) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**T\_PH**](#function-t_ph-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & H) <br>_Calculate temperature T by given pressure and specific enthalpy._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**T\_PH**](#function-t_ph-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; H, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**T\_critical**](#function-t_critical) () <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**T\_sat\_P**](#function-t_sat_p-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Calculate saturated temperature by given pressure. Valid pressure rante if [H2O::P\_MIN, H2O::P\_c] See Equation 31 and Table 34 of IF97._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**T\_sat\_P**](#function-t_sat_p-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Tmax**](#function-tmax) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Tmin**](#function-tmin) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Ttriple**](#function-ttriple) () <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**UpdateState\_TP**](#function-updatestate_tp) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p) <br> |
|   | [**cIAPWS\_IF97**](#function-ciapws_if97) () <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**getState\_Region1**](#function-getstate_region1) ([**double**](classxThermal_1_1xThermalError.md) P, [**double**](classxThermal_1_1xThermalError.md) T, [**State\_Region1**](structxThermal_1_1IAPWS__IF97_1_1State__Region1.md) & state) <br>_Calculate Gibbs free energy and its derivatives by given pressure and temperature._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**getState\_Region2**](#function-getstate_region2) ([**double**](classxThermal_1_1xThermalError.md) P, [**double**](classxThermal_1_1xThermalError.md) T, [**State\_Region2**](structxThermal_1_1IAPWS__IF97_1_1State__Region2.md) & state) <br>_Calculate dimensionless Gibbs free energy in the region 2. ideal-gas part_ \(\gamma^o\) _and the residual part_\(\gamma^r\) _. See equation 15-17 of IF97._ |
|  [**void**](classxThermal_1_1xThermalError.md) | [**getState\_Region5**](#function-getstate_region5) ([**double**](classxThermal_1_1xThermalError.md) P, [**double**](classxThermal_1_1xThermalError.md) T, [**State\_Region5**](structxThermal_1_1IAPWS__IF97_1_1State__Region5.md) & state) <br>_Calculate dimensionless Gibbs free energy in the region 5. ideal-gas part_ \(\gamma^o\) _and the residual part_\(\gamma^r\) _. See equation 32-34 of IF97._ |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**molar\_mass**](#function-molar_mass) () <br> |
| virtual std::string | [**name**](#function-name) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**p\_critical**](#function-p_critical) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**pmax**](#function-pmax) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**pmin**](#function-pmin) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**rhomass\_critical**](#function-rhomass_critical) () <br> |
|   | [**~cIAPWS\_IF97**](#function-ciapws_if97) () <br> |


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


















































## Public Functions Documentation




### function Backward\_T\_PH\_region1 

_Calculate temperature in region 1 by given P and H. See equation 11 of IF97._ 
```C++
double xThermal::IAPWS_IF97::cIAPWS_IF97::Backward_T_PH_region1 (
    double P,
    double H
) 
```





**Parameters:**


* `P` [Pa] 
* `H` [J/kg] 



**Returns:**

double 





        

<hr>



### function Backward\_T\_PH\_region2a 

_Calculate saturated temperature by given pressure. Valid pressure rante if [H2O::m\_constants.pmin, H2O::m\_constants.p\_critical] See Equation 31 and Table 34 of IF97._ 
```C++
double xThermal::IAPWS_IF97::cIAPWS_IF97::Backward_T_PH_region2a (
    double P,
    double H
) 
```





**Parameters:**


* `P` [Pa] 



**Returns:**

double [K]


Calculate temperature in region 2a by given P and H




**Parameters:**


* `P` [Pa] 
* `H` [J/kg] 



**Returns:**

double 





        

<hr>



### function Backward\_T\_PH\_region2b 

```C++
double xThermal::IAPWS_IF97::cIAPWS_IF97::Backward_T_PH_region2b (
    double P,
    double H
) 
```




<hr>



### function Backward\_T\_PH\_region2c 

```C++
double xThermal::IAPWS_IF97::cIAPWS_IF97::Backward_T_PH_region2c (
    double P,
    double H
) 
```




<hr>



### function Backward\_T\_PH\_region3a 

_The backward equation_ \(T_{3a}(p,h)\) _for subregion 3a. See equation 2 of IF97-Region3._
```C++
double xThermal::IAPWS_IF97::cIAPWS_IF97::Backward_T_PH_region3a (
    double P,
    double H
) 
```





**Parameters:**


* `P` 
* `H` 



**Returns:**

double 





        

<hr>



### function Backward\_T\_PH\_region3b 

_The backward equation_ \(T_{3b}(p,h)\) _for subregion 3a. See equation 3 of IF97-Region3._
```C++
double xThermal::IAPWS_IF97::cIAPWS_IF97::Backward_T_PH_region3b (
    double P,
    double H
) 
```





**Parameters:**


* `P` 
* `H` 



**Returns:**

double 





        

<hr>



### function Boiling\_p 

```C++
inline virtual double xThermal::IAPWS_IF97::cIAPWS_IF97::Boiling_p (
    const  double & T
) 
```



Calculate boiling pressure [Pa] of water for a given T [K] 

**Parameters:**


* `T` 



**Returns:**







        
Implements [*xThermal::cxThermal::Boiling\_p*](classxThermal_1_1cxThermal.md#function-boiling_p-13)


<hr>



### function Boundary\_region23\_P2T 

_The boundary between regions 2 and 3 (see Fig. 1) is defined by the following simple quadratic pressure-temperature relation, the B23-equation. See equation 5 of IF97._ 
```C++
double xThermal::IAPWS_IF97::cIAPWS_IF97::Boundary_region23_P2T (
    double P
) 
```





**Parameters:**


* `P` [Pa] 



**Returns:**

double 





        

<hr>



### function Boundary\_region2b2c\_H2P 

_Boundary between region 2b and 2c. See also_ [_**Boundary\_region2b2c\_P2H**_](classxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97.md#function-boundary_region2b2c_p2h) _._
```C++
double xThermal::IAPWS_IF97::cIAPWS_IF97::Boundary_region2b2c_H2P (
    double H
) 
```





**Parameters:**


* `H` [J/kg] 



**Returns:**

double [Pa] 





        

<hr>



### function Boundary\_region2b2c\_P2H 

_Boundary between region 2b and 2c. In order to know whether the T(p,h) equation for subregion 2b or for subregion 2c has to be used for given values of p and h, a special correlation equation for the boundary between subregions 2b and 2c (which approximates s = 5.85 kJ/kg/K_ [_**CONST\_IF97\_S\_Region2b2c**_](group__CONSTS__IF97.md#define-const_if97_s_region2b2c) _) is needed; see Fig. 2. This boundary equation, called the B2bc-equation, is a simple quadratic pressure-enthalpy relation which reads Equation 20 of IF97._
```C++
double xThermal::IAPWS_IF97::cIAPWS_IF97::Boundary_region2b2c_P2H (
    double P
) 
```





**Parameters:**


* `P` [Pa] 



**Returns:**

double [J/kg] 





        

<hr>



### function Boundary\_region3ab\_H2P [1/2]

_Boundary between region3a and 3b. See Equation 10 and Table 17 of IF97-Region3._ 
```C++
double xThermal::IAPWS_IF97::cIAPWS_IF97::Boundary_region3ab_H2P (
    double H
) 
```





**Note:**

The valid H range is \([h^{\prime}(#CONST_IF97_Tmin_Region3), h^{\prime\prime}(#CONST_IF97_Tmin_Region3)]\). See also pp.17 of IF97-Region3 . 




**Parameters:**


* `H` [J/kg] 



**Returns:**

double [Pa] 





        

<hr>



### function Boundary\_region3ab\_H2P [2/2]

```C++
void xThermal::IAPWS_IF97::cIAPWS_IF97::Boundary_region3ab_H2P (
    const std::vector< double > H,
    std::vector< double > & res
) 
```




<hr>



### function Boundary\_region3ab\_P2H [1/2]

_Boundary between region3a and 3b. See Equation 1 and Table 2 of IF97-Region3._ 
```C++
double xThermal::IAPWS_IF97::cIAPWS_IF97::Boundary_region3ab_P2H (
    double P
) 
```





**Note:**

This equation is only valid in P range of [#m\_constants.p\_critical, [**CONST\_IF97\_Pmax\_Region1**](group__CONSTS__IF97.md#define-const_if97_pmax_region1)]. 




**Parameters:**


* `P` [Pa] 



**Returns:**

double 





        

<hr>



### function Boundary\_region3ab\_P2H [2/2]

```C++
void xThermal::IAPWS_IF97::cIAPWS_IF97::Boundary_region3ab_P2H (
    const std::vector< double > P,
    std::vector< double > & res
) 
```




<hr>



### function GetRegion\_PH [1/2]

_Calculate region index by given pressure and specific enthalpy._ 
```C++
int xThermal::IAPWS_IF97::cIAPWS_IF97::GetRegion_PH (
    double P,
    double H
) 
```





**Parameters:**


* `P` [Pa] 
* `H` [J/kg] 



**Returns:**

int 





        

<hr>



### function GetRegion\_PH [2/2]

```C++
void xThermal::IAPWS_IF97::cIAPWS_IF97::GetRegion_PH (
    const std::vector< double > P,
    const std::vector< double > H,
    std::vector< int > & res
) 
```




<hr>



### function GetRegion\_PT [1/2]

_Calculate region index by given pressure and temperature._ 
```C++
int xThermal::IAPWS_IF97::cIAPWS_IF97::GetRegion_PT (
    double P,
    double T
) 
```





**Note:**

In P-T space, there is no need to introduce subregions, e.g. 2a, 2b, 2c and 3a, 3b. 




**Parameters:**


* `P` [Pa] 
* `T` [K] 



**Returns:**

int 





        

<hr>



### function GetRegion\_PT [2/2]

```C++
void xThermal::IAPWS_IF97::cIAPWS_IF97::GetRegion_PT (
    const std::vector< double > P,
    const std::vector< double > T,
    std::vector< int > & res
) 
```




<hr>



### function H\_sat\_P [1/2]

_Calculate saturated liquid enthalpy and saturated vapor enthalpy for a given P._ 
```C++
void xThermal::IAPWS_IF97::cIAPWS_IF97::H_sat_P (
    const  double & P,
    double & H_l,
    double & H_v
) 
```





**Note:**

This function is only valid in P range [#P\_MIN, [**CONST\_IF97\_Pmin\_Region3**](group__CONSTS__IF97.md#define-const_if97_pmin_region3)]. 




**Parameters:**


* `P` 
* `H_l` 
* `H_v` 




        

<hr>



### function H\_sat\_P [2/2]

```C++
void xThermal::IAPWS_IF97::cIAPWS_IF97::H_sat_P (
    const std::vector< double > P,
    std::vector< double > & H_l,
    std::vector< double > & H_v
) 
```




<hr>



### function P\_sat\_T [1/2]

_Calculate saturated temperature by given pressure. Valid pressure rante if [H2O::T\_MIN, H2O::T\_c] See Equation 30 and Table 34 of IF97._ 
```C++
double xThermal::IAPWS_IF97::cIAPWS_IF97::P_sat_T (
    const  double & T
) 
```





**Parameters:**


* `T` [K] 



**Returns:**

double [Pa] 





        

<hr>



### function P\_sat\_T [2/2]

```C++
void xThermal::IAPWS_IF97::cIAPWS_IF97::P_sat_T (
    const std::vector< double > T,
    std::vector< double > & res
) 
```




<hr>



### function Prop\_Region1 

_Calculate thermodynamic properties by given dimensionless Gibbs free energy and its derivatives for the region 1._ 
```C++
double xThermal::IAPWS_IF97::cIAPWS_IF97::Prop_Region1 (
    const  State_Region1 & state,
    BasicThermodynamicProperties which
) 
```





**Parameters:**


* `gamma` 
* `which` 



**Returns:**

double 





        

<hr>



### function Prop\_Region2 

```C++
double xThermal::IAPWS_IF97::cIAPWS_IF97::Prop_Region2 (
    State_Region2 state,
    BasicThermodynamicProperties which
) 
```




<hr>



### function Prop\_Region5 

```C++
double xThermal::IAPWS_IF97::cIAPWS_IF97::Prop_Region5 (
    double P,
    double T,
    State_Region5 state,
    BasicThermodynamicProperties which
) 
```




<hr>



### function T\_PH [1/2]

_Calculate temperature T by given pressure and specific enthalpy._ 
```C++
double xThermal::IAPWS_IF97::cIAPWS_IF97::T_PH (
    const  double & P,
    const  double & H
) 
```





**Parameters:**


* `P` 
* `H` 



**Returns:**

double 





        

<hr>



### function T\_PH [2/2]

```C++
void xThermal::IAPWS_IF97::cIAPWS_IF97::T_PH (
    const std::vector< double > P,
    const std::vector< double > H,
    std::vector< double > & res
) 
```




<hr>



### function T\_critical 

```C++
inline virtual double xThermal::IAPWS_IF97::cIAPWS_IF97::T_critical () 
```



Get the triple point temperature in K 


        
Implements [*xThermal::cxThermal::T\_critical*](classxThermal_1_1cxThermal.md#function-t_critical)


<hr>



### function T\_sat\_P [1/2]

_Calculate saturated temperature by given pressure. Valid pressure rante if [H2O::P\_MIN, H2O::P\_c] See Equation 31 and Table 34 of IF97._ 
```C++
double xThermal::IAPWS_IF97::cIAPWS_IF97::T_sat_P (
    const  double & P
) 
```





**Parameters:**


* `P` [Pa] 



**Returns:**

double [K] 





        

<hr>



### function T\_sat\_P [2/2]

```C++
void xThermal::IAPWS_IF97::cIAPWS_IF97::T_sat_P (
    const std::vector< double > P,
    std::vector< double > & res
) 
```




<hr>



### function Tmax 

```C++
inline virtual double xThermal::IAPWS_IF97::cIAPWS_IF97::Tmax () 
```



Get the minimum temperature in K 


        
Implements [*xThermal::cxThermal::Tmax*](classxThermal_1_1cxThermal.md#function-tmax)


<hr>



### function Tmin 

```C++
inline virtual double xThermal::IAPWS_IF97::cIAPWS_IF97::Tmin () 
```



Get the minimum temperature in K 


        
Implements [*xThermal::cxThermal::Tmin*](classxThermal_1_1cxThermal.md#function-tmin)


<hr>



### function Ttriple 

```C++
inline virtual double xThermal::IAPWS_IF97::cIAPWS_IF97::Ttriple () 
```



Get the maximum pressure in Pa 


        
Implements [*xThermal::cxThermal::Ttriple*](classxThermal_1_1cxThermal.md#function-ttriple)


<hr>



### function UpdateState\_TP 

```C++
inline void xThermal::IAPWS_IF97::cIAPWS_IF97::UpdateState_TP (
    ThermodynamicProperties & props,
    const  double & T,
    const  double & p
) 
```



Return the molar mass in kg/mol 


        

<hr>



### function cIAPWS\_IF97 

```C++
xThermal::IAPWS_IF97::cIAPWS_IF97::cIAPWS_IF97 () 
```




<hr>



### function getState\_Region1 

_Calculate Gibbs free energy and its derivatives by given pressure and temperature._ 
```C++
void xThermal::IAPWS_IF97::cIAPWS_IF97::getState_Region1 (
    double P,
    double T,
    State_Region1 & state
) 
```





**Parameters:**


* `P` [Pa] 
* `T` [K] 
* `state` 




        

<hr>



### function getState\_Region2 

_Calculate dimensionless Gibbs free energy in the region 2. ideal-gas part_ \(\gamma^o\) _and the residual part_\(\gamma^r\) _. See equation 15-17 of IF97._
```C++
void xThermal::IAPWS_IF97::cIAPWS_IF97::getState_Region2 (
    double P,
    double T,
    State_Region2 & state
) 
```





**Parameters:**


* `P` 
* `T` 
* `gammao` 
* `gammar` 




        

<hr>



### function getState\_Region5 

_Calculate dimensionless Gibbs free energy in the region 5. ideal-gas part_ \(\gamma^o\) _and the residual part_\(\gamma^r\) _. See equation 32-34 of IF97._
```C++
void xThermal::IAPWS_IF97::cIAPWS_IF97::getState_Region5 (
    double P,
    double T,
    State_Region5 & state
) 
```





**Parameters:**


* `P` 
* `T` 
* `state` 




        

<hr>



### function molar\_mass 

```C++
inline virtual double xThermal::IAPWS_IF97::cIAPWS_IF97::molar_mass () 
```



Return the critical mass density in \(kg/m^3\) 


        
Implements [*xThermal::cxThermal::molar\_mass*](classxThermal_1_1cxThermal.md#function-molar_mass)


<hr>



### function name 

```C++
inline virtual std::string xThermal::IAPWS_IF97::cIAPWS_IF97::name () 
```



Name of the model 


        
Implements [*xThermal::cxThermal::name*](classxThermal_1_1cxThermal.md#function-name)


<hr>



### function p\_critical 

```C++
inline virtual double xThermal::IAPWS_IF97::cIAPWS_IF97::p_critical () 
```



Return the critical temperature in K 


        
Implements [*xThermal::cxThermal::p\_critical*](classxThermal_1_1cxThermal.md#function-p_critical)


<hr>



### function pmax 

```C++
inline virtual double xThermal::IAPWS_IF97::cIAPWS_IF97::pmax () 
```



Get the minimum pressure in Pa 


        
Implements [*xThermal::cxThermal::pmax*](classxThermal_1_1cxThermal.md#function-pmax)


<hr>



### function pmin 

```C++
inline virtual double xThermal::IAPWS_IF97::cIAPWS_IF97::pmin () 
```



Get the maximum temperature in K 


        
Implements [*xThermal::cxThermal::pmin*](classxThermal_1_1cxThermal.md#function-pmin)


<hr>



### function rhomass\_critical 

```C++
inline virtual double xThermal::IAPWS_IF97::cIAPWS_IF97::rhomass_critical () 
```



Return the critical pressure in Pa 


        
Implements [*xThermal::cxThermal::rhomass\_critical*](classxThermal_1_1cxThermal.md#function-rhomass_critical)


<hr>



### function ~cIAPWS\_IF97 

```C++
xThermal::IAPWS_IF97::cIAPWS_IF97::~cIAPWS_IF97 () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2O/IAPWS-IF97.h`

