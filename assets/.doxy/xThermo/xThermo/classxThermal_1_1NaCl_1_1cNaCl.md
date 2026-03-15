

# Class xThermal::NaCl::cNaCl



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**NaCl**](namespacexThermal_1_1NaCl.md) **>** [**cNaCl**](classxThermal_1_1NaCl_1_1cNaCl.md)



_Class of NaCl EOS._ [More...](#detailed-description)

* `#include <NaCl.h>`



Inherits the following classes: [xThermal::cxThermal](classxThermal_1_1cxThermal.md)
























## Public Attributes inherited from xThermal::cxThermal

See [xThermal::cxThermal](classxThermal_1_1cxThermal.md)

| Type | Name |
| ---: | :--- |
|  [**bool**](classxThermal_1_1xThermalError.md) | [**m\_isShowProgressBar**](classxThermal_1_1cxThermal.md#variable-m_isshowprogressbar)  <br> |






























## Public Functions

| Type | Name |
| ---: | :--- |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_T**](#function-boiling_t-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Boiling temperature as a function of pressure. See also_ [_**Boiling\_p**_](classxThermal_1_1NaCl_1_1cNaCl.md#function-boiling_p-12) _._ |
|  [**void**](classxThermal_1_1xThermalError.md) | [**Boiling\_T**](#function-boiling_t-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_p**](#function-boiling_p-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T) <br>_The T-P relations of the boiling curve of NaCl. See equation 3, table 3 and figure 4 of Driesner2007Part1._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**Boiling\_p**](#function-boiling_p-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Cp\_Liquid**](#function-cp_liquid) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Specific heat capacity of molten NaCl (liquid phase), which is calculated from Eq. 27 of Driesner2007Part2 when_ \(X_{NaCl} =1\) _._ |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Cp\_Solid**](#function-cp_solid-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Specific heat capacity of halite (solid NaCl). See equation 30 of Driesner2007Part2._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**Cp\_Solid**](#function-cp_solid-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**DeltaH\_fus**](#function-deltah_fus-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & T, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**H\_Liquid**](#function-h_liquid-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Calculate specific enthalpy of liquid NaCl. See equation 21 of Driesner2007Part2._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**H\_Liquid**](#function-h_liquid-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & T, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**H\_Solid**](#function-h_solid-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Calculate specific enthalpy of halite (solid NaCl)._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**H\_Solid**](#function-h_solid-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & T, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Melting\_T**](#function-melting_t-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Melting curve. See equation 1, table 3 and figure 3 of Driesner2007Part1._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**Melting\_T**](#function-melting_t-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Melting\_T\_C**](#function-melting_t_c) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Melting temperature, the same as_ [_**Melting\_T**_](classxThermal_1_1NaCl_1_1cNaCl.md#function-melting_t-12) _except temperature unit._ |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Melting\_p**](#function-melting_p-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T) <br>_Melting pressure as a function of temperature. See also_ [_**Melting\_T**_](classxThermal_1_1NaCl_1_1cNaCl.md#function-melting_t-12) _._ |
|  [**void**](classxThermal_1_1xThermalError.md) | [**Melting\_p**](#function-melting_p-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**P\_Vapor**](#function-p_vapor-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T) <br>_Calculate vapor pressure of NaCl. See also_ Sublimation\_P _and_[_**Boiling\_p**_](classxThermal_1_1NaCl_1_1cNaCl.md#function-boiling_p-12) _._ |
|  [**void**](classxThermal_1_1xThermalError.md) | [**P\_Vapor**](#function-p_vapor-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Rho\_Liquid**](#function-rho_liquid-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Density of liquid (molten) NaCl as function of (T, P). See equation (4-6) of Driesner2007Part2._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**Rho\_Liquid**](#function-rho_liquid-22) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**double**](classxThermal_1_1xThermalError.md) & rho, [**double**](classxThermal_1_1xThermalError.md) & dRhodP, [**double**](classxThermal_1_1xThermalError.md) & dRhodT, [**double**](classxThermal_1_1xThermalError.md) & kappa, [**double**](classxThermal_1_1xThermalError.md) & beta) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Rho\_Solid**](#function-rho_solid-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Halite density as function of (T, P). See equation (1-3) of Driesner2007Part2._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**Rho\_Solid**](#function-rho_solid-22) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**double**](classxThermal_1_1xThermalError.md) & rho, [**double**](classxThermal_1_1xThermalError.md) & dRhodP, [**double**](classxThermal_1_1xThermalError.md) & dRhodT, [**double**](classxThermal_1_1xThermalError.md) & kappa, [**double**](classxThermal_1_1xThermalError.md) & beta) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Sublimation\_T**](#function-sublimation_t-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Sublimation temperature as a function of pressure. See also_ Sublimation\_P _._ |
|  [**void**](classxThermal_1_1xThermalError.md) | [**Sublimation\_T**](#function-sublimation_t-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Sublimation\_p**](#function-sublimation_p-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T) <br>_Sublimation curve. See equation 2, table 3 and figure 4 of Driesner2007Part1._  |
|  [**void**](classxThermal_1_1xThermalError.md) | [**Sublimation\_p**](#function-sublimation_p-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; T, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**T\_Vapor**](#function-t_vapor-12) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Calculate vapor temperature of NaCl. See also_ [_**Sublimation\_T**_](classxThermal_1_1NaCl_1_1cNaCl.md#function-sublimation_t-12) _and_[_**Boiling\_T**_](classxThermal_1_1NaCl_1_1cNaCl.md#function-boiling_t-12) _._ |
|  [**void**](classxThermal_1_1xThermalError.md) | [**T\_Vapor**](#function-t_vapor-22) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; P, std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & res) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**T\_critical**](#function-t_critical) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Tmax**](#function-tmax) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Tmin**](#function-tmin) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Ttriple**](#function-ttriple) () <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**UpdateState\_TPX**](#function-updatestate_tpx) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br>_Update state by given T,p as input pair._  |
|   | [**cNaCl**](#function-cnacl-12) (std::string backend\_H2O) <br> |
|   | [**cNaCl**](#function-cnacl-22) ([**const**](classxThermal_1_1xThermalError.md) [**cNaCl**](classxThermal_1_1NaCl_1_1cNaCl.md) & salt) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**molar\_mass**](#function-molar_mass) () <br> |
| virtual std::string | [**name**](#function-name) () <br> |
|  std::string | [**name\_backend**](#function-name_backend) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**p\_critical**](#function-p_critical) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**pmax**](#function-pmax) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**pmin**](#function-pmin) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**rhomass\_critical**](#function-rhomass_critical) () <br> |
|   | [**~cNaCl**](#function-cnacl) () <br> |


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


![Image](Melting_Sublimation_Boiling_NaCl.svg)
 


    
## Public Functions Documentation




### function Boiling\_T [1/2]

_Boiling temperature as a function of pressure. See also_ [_**Boiling\_p**_](classxThermal_1_1NaCl_1_1cNaCl.md#function-boiling_p-12) _._
```C++
virtual double xThermal::NaCl::cNaCl::Boiling_T (
    const  double & P
) 
```





**Parameters:**


* `P` [Pa] 



**Returns:**

double [K] 





        
Implements [*xThermal::cxThermal::Boiling\_T*](classxThermal_1_1cxThermal.md#function-boiling_t-13)


<hr>



### function Boiling\_T [2/2]

```C++
void xThermal::NaCl::cNaCl::Boiling_T (
    const std::vector< double > P,
    std::vector< double > & res
) 
```




<hr>



### function Boiling\_p [1/2]

_The T-P relations of the boiling curve of NaCl. See equation 3, table 3 and figure 4 of Driesner2007Part1._ 
```C++
virtual double xThermal::NaCl::cNaCl::Boiling_p (
    const  double & T
) 
```





**Note:**

The unit of result is Pa. Unit of [**T\_Triple**](group__CONSTANTS__NACL.md#variable-t_triple) and [**P\_Triple**](group__CONSTANTS__NACL.md#variable-p_triple) is [K] and [Pa], respectively.


![Image](NaCl/Melting_Sublimation_Boiling_NaCl.svg)





**Parameters:**


* `T` [K] 



**Returns:**

double [Pa] 





        
Implements [*xThermal::cxThermal::Boiling\_p*](classxThermal_1_1cxThermal.md#function-boiling_p-13)


<hr>



### function Boiling\_p [2/2]

```C++
void xThermal::NaCl::cNaCl::Boiling_p (
    const std::vector< double > T,
    std::vector< double > & res
) 
```




<hr>



### function Cp\_Liquid 

_Specific heat capacity of molten NaCl (liquid phase), which is calculated from Eq. 27 of Driesner2007Part2 when_ \(X_{NaCl} =1\) _._
```C++
double xThermal::NaCl::cNaCl::Cp_Liquid (
    const  double & T,
    const  double & P
) 
```





**Parameters:**


* `T` [K] 
* `P` [Pa] 



**Returns:**

double [J/kg] 





        

<hr>



### function Cp\_Solid [1/2]

_Specific heat capacity of halite (solid NaCl). See equation 30 of Driesner2007Part2._ 
```C++
double xThermal::NaCl::cNaCl::Cp_Solid (
    const  double & T,
    const  double & P
) 
```





**Parameters:**


* `T` [K] 
* `P` [Pa] 



**Returns:**

double [J/kg/K] 





        

<hr>



### function Cp\_Solid [2/2]

```C++
void xThermal::NaCl::cNaCl::Cp_Solid (
    const std::vector< double > T,
    const std::vector< double > P,
    std::vector< double > & res
) 
```




<hr>



### function DeltaH\_fus [2/2]

```C++
void xThermal::NaCl::cNaCl::DeltaH_fus (
    const std::vector< double > & T,
    const std::vector< double > & P,
    std::vector< double > & res
) 
```




<hr>



### function H\_Liquid [1/2]

_Calculate specific enthalpy of liquid NaCl. See equation 21 of Driesner2007Part2._ 
```C++
double xThermal::NaCl::cNaCl::H_Liquid (
    const  double & T,
    const  double & P
) 
```





**Parameters:**


* `T` [K] 
* `P` [Pa] 



**Returns:**

double [J/kg] 





        

<hr>



### function H\_Liquid [2/2]

```C++
void xThermal::NaCl::cNaCl::H_Liquid (
    const std::vector< double > & T,
    const std::vector< double > & P,
    std::vector< double > & res
) 
```




<hr>



### function H\_Solid [1/2]

_Calculate specific enthalpy of halite (solid NaCl)._ 
```C++
double xThermal::NaCl::cNaCl::H_Solid (
    const  double & T,
    const  double & P
) 
```




* Basic formula:  
  \[\begin{align}
H_{halite}(T, P) = \int_{(T_0,P_0)}^{(T,P)} dH = \int_{(T_0,P_0)}^{(T,P)} \left( \frac{\partial H}{\partial T}\right)_P dT + \left( \frac{\partial H}{\partial P} \right)_T dP = \int_{T_0}^{T} \left( \frac{\partial H}{\partial T}\right)_{P=P_0} dT + \int_{P_0}^{P} \left( \frac{\partial H}{\partial P} \right)_T dP
\end{align}\]

* Integration along T (See Eq. 30 of Driesner2007Part2 for \(C_p\) of halite)  \(\begin{align}
\int_{T_0}^{T} \left( \frac{\partial H}{\partial T}\right)_{P=P_0} dT &= \int_{T_0}^{T} C_p(T, P_0) dT =
\left[
\underbrace{\left(r_0 - 2r_1 T_{triple, NaCl} + 3r_2T_{triple, NaCl}^2 + \color{red}{r_{3a}P} + r_4P^2 \right)}_{R_1}\color{blue}{T}
+ \underbrace{\left(r_1 - 3r_2T_{triple, NaCl} + \frac{1}{2}\color{red}{r_{3b}P} \right)}_{R_2}\color{blue}{T^2}
+ \underbrace{\left(r_2 + \frac{1}{3}\color{red}{r_{3c}}P \right)}_{R_3}\color{blue}{T^3}
\right] \bigg|_{T_0}^{T} = H_1(T, P_0) - H_1(T_0, P_0) \\
& = \sum_{i=1}^{3}R_i T^i - R_0, ~~(R_0=\sum_{i=1}^{3}R_i T_0^i = H_1(T_0, P_0))
\end{align}\)
* Integration along P (See Eq. 29 of Driesner2007Part2 for \(\left(\frac{\partial H}{\partial P}\right)_T\) of halite, and Eq. 1-3 for density \(\rho_h\) of halite)  \(\begin{align}
\int_{P_0}^P \left( \frac{\partial H}{\partial P} \right)_T dP = \int_{P_0}^P \left[ V - {\color{orange}{T}}\left(\frac{\partial V}{\partial T} \right)_P \right] dP = \color{green}{\int_{P_0}^P V dP} -\color{orange}{T}\frac{\partial }{\partial T} \left( \color{green}{\int_{P_0}^P V dP} \right) = H_2(T, P) - H_2(T, P_0)
\end{align}\)






**Warning:**

The unit of \(\color{orange}{T}\) in front of \(\left( \frac{\partial V}{\partial T} \right)_P\) is [K]


Where  \(\begin{align}
H_2(T, P) = H_2^*(T, P) - \color{orange}{T} H_2^{**}(T, P)
\end{align}\)



\[\begin{align}
\color{green}{\int_{P_0}^P V dP} = \int_{P_0}^P \frac{1}{lP + \rho^0_h} dP = \frac{1}{l}ln\rho_h\bigg|_{P_0}^P = H_2^*(T, P) - H_2^*(T, P_0)
\end{align}\]




\[\begin{align}
\frac{\partial }{\partial T} \left( \color{green}{\int_{P_0}^P V dP} \right) = \frac{\partial H_2(T, P)}{\partial T} - \frac{\partial H_2(T, P_0)}{\partial T} = H_2^{**}(T, P) - H_2^{**}(T, P_0)
\end{align}\]




\[\begin{align}
H_2^{**}(T, P) = \frac{\partial H_2(T, P)}{\partial T} = \frac{\partial }{\partial T}\left(\frac{ln\rho_h}{l} \right) = -\frac{ln\rho_h}{l^2}\frac{\partial l}{\partial T} + \frac{1}{l\rho_h}\frac{\partial \rho_h}{\partial T} = -\frac{ln\rho_h}{l^2} \frac{l_4}{l_5}e^{T/l_5} + \frac{1}{l\rho_h} \left( l_1 + 2l_2T + P\frac{l_4}{l_5}e^{T/l_5} \right)
\end{align}\]



Therefore,  
\[\begin{align}
\boxed{ H_{halite}(T, P) = H_1(T, P) - H_1(T_0, P) + \left[ H_2(T, P) - H_2(T, P_0)\right] + C }, \\
\text{where } C=T_{halite}(T_0, P_0) \text{ is the reference value, we select a reference point and value as } \boxed{T_{halite}(T=100^{\circ}, P=100bar)=9.415867359\times10^4 ~ J/kg}
\end{align}\]





**Note:**

The unit of \(T\) in the above equations is \(^{\circ}C\). Unlike the pressure \(P\), it is not possible to simply change the coefficient to make the unit of T to [K], so just keep it as \(^{\circ} C\).


**Some original ideas as following**




**Note:**

**How to calculate halite (solid phase) enthalpy at the triple point ?**
* The enthalpy of molten NaCl (liquid phase) \(H_l\) can be calculated from Eq. 21 of Driesner2007Part2 by subsituting \(X_{NaCl} = 1\).
* The enthalpy difference \(\Delta H_{fus} = H_l - H_h\) of phase changing (solid phase to liquid phase) can be calculated from Eq. 28 of Driesner2007Part2.(1). Clapeyron slope of the NaCl melting curve is \(\frac{dP}{dT} = \frac{1}{a} = 4044325.810887325~ Pa/K\); (**notice the unit!**)(2). Triple temperature \(T_{Triple} = 1073.85~ K\);(3). Volume of melting at triple point \(\Delta V_{fus} = \frac{1}{\rho_l} - \frac{1}{\rho_h} = \frac{1}{1561.224722497394} - \frac{1}{1912.0183915214473} = 0.00011751525885297855  ~ kg/m^3 = 6.867944273144626\times 10^{-6} ~ m^3/mol\) (see [**Rho\_Liquid**](classxThermal_1_1NaCl_1_1cNaCl.md#function-rho_liquid-12) and [**Rho\_Solid**](classxThermal_1_1NaCl_1_1cNaCl.md#function-rho_solid-12), and [**M**](group__CONSTANTS__NACL.md#variable-m)).Therefore, \(\Delta H_{fus} = \frac{dP}{dT}\times T_{Triple}\times \Delta V_{fus} = 4044325.810887325 \times 1073.85 \times 6.867944273144626\times 10^{-6} = 29827.476978550334 kJ/mol\) which is constent with value \(29.8 ~ kJ/mol\) in Driesner2007Part2.
* Enthalpy of halite at triple point temperature \(H_h(T_{Triple}, P) = H_l(T_{Triple}, P) - \Delta H_{fus}(T_{Triple}, P)\).






**Todo**

Confirm that \(\Delta H_{fus} (T_{melting}, P) = H_l(T_{melting}, P) - H_s(T_{melting}, P)\) ?






**Note:**

The constant \(C(P)\) in the above equation


We easily calculate the constant C(P) once get \(H_h\), so \(C(P) = H_h(T_{Triple}, P) - H_{halite}(T_{Triple}, P)\).




**Warning:**

Just like the description in 2.2.4 of Driesner2007Part2 : "Pure NaCl vapor as a separate phase exists only at pressures substantially lower than those encountered in any geological system", no worries that the (T, P) point occurs in pure vapor region of NaCl (See figure at [**Boiling\_p(const double& T)**](classxThermal_1_1NaCl_1_1cNaCl.md#function-boiling_p-12)), i.e. NaCl always present as liquid or solid state in the geological T-P conditions.




**Parameters:**


* `T` 
* `P` 



**Returns:**

double 





        

<hr>



### function H\_Solid [2/2]

```C++
void xThermal::NaCl::cNaCl::H_Solid (
    const std::vector< double > & T,
    const std::vector< double > & P,
    std::vector< double > & res
) 
```




<hr>



### function Melting\_T [1/2]

_Melting curve. See equation 1, table 3 and figure 3 of Driesner2007Part1._ 
```C++
double xThermal::NaCl::cNaCl::Melting_T (
    const  double & P
) 
```





**Note:**

Because the unit of P\_Triple and P is [Pa], so the parameter \(a=2.4726E-7\).


![Image](NaCl/Melting_Sublimation_Boiling_NaCl.svg)





**Parameters:**


* `P` [Pa] 



**Returns:**

double [K] 





        

<hr>



### function Melting\_T [2/2]

```C++
void xThermal::NaCl::cNaCl::Melting_T (
    const std::vector< double > P,
    std::vector< double > & res
) 
```




<hr>



### function Melting\_T\_C 

_Melting temperature, the same as_ [_**Melting\_T**_](classxThermal_1_1NaCl_1_1cNaCl.md#function-melting_t-12) _except temperature unit._
```C++
double xThermal::NaCl::cNaCl::Melting_T_C (
    const  double & P
) 
```





**Parameters:**


* `P` [Pa] 



**Returns:**

double [deg.C] 





        

<hr>



### function Melting\_p [1/2]

_Melting pressure as a function of temperature. See also_ [_**Melting\_T**_](classxThermal_1_1NaCl_1_1cNaCl.md#function-melting_t-12) _._
```C++
double xThermal::NaCl::cNaCl::Melting_p (
    const  double & T
) 
```





**Parameters:**


* `T` [K] 



**Returns:**

double [Pa] 





        

<hr>



### function Melting\_p [2/2]

```C++
void xThermal::NaCl::cNaCl::Melting_p (
    const std::vector< double > T,
    std::vector< double > & res
) 
```




<hr>



### function P\_Vapor [1/2]

_Calculate vapor pressure of NaCl. See also_ Sublimation\_P _and_[_**Boiling\_p**_](classxThermal_1_1NaCl_1_1cNaCl.md#function-boiling_p-12) _._
```C++
double xThermal::NaCl::cNaCl::P_Vapor (
    const  double & T
) 
```



![Image](NaCl/Melting_Sublimation_Boiling_NaCl.svg)





**Parameters:**


* `T` [K] 



**Returns:**

double [Pa] 





        

<hr>



### function P\_Vapor [2/2]

```C++
void xThermal::NaCl::cNaCl::P_Vapor (
    const std::vector< double > T,
    std::vector< double > & res
) 
```




<hr>



### function Rho\_Liquid [1/2]

_Density of liquid (molten) NaCl as function of (T, P). See equation (4-6) of Driesner2007Part2._ 
```C++
double xThermal::NaCl::cNaCl::Rho_Liquid (
    const  double & T,
    const  double & P
) 
```





**Parameters:**


* `T` 
* `P` 



**Returns:**

double 





        

<hr>



### function Rho\_Liquid [2/2]

```C++
void xThermal::NaCl::cNaCl::Rho_Liquid (
    const  double & T,
    const  double & P,
    double & rho,
    double & dRhodP,
    double & dRhodT,
    double & kappa,
    double & beta
) 
```



Calculate melting NaCl density \(\rho\) and derivatives \(\left( \frac{\partial \rho}{\partial P} \right)_T\), \(\left( \frac{\partial \rho}{\partial T} \right)_P\)



* \(\rho = \frac{\rho^0}{1 - 0.1ln(1 + 10P \kappa)}\): Density of liquid (molten) NaCl as function of (T, P). See equation (4-6) of Driesner2007Part2.
* 
  \[v = \frac{1}{\rho} = \frac{1}{\rho^0} \left[ 1 - 0.1ln(1 + 10P \kappa) \right]\]

* 
  \[\frac{\partial v}{ \partial P} = - \frac{\kappa}{\rho^0 (1 + 10\kappa P)}\]

* 
  \[\left( \frac{\partial \rho}{\partial P} \right)_T = -\frac{1}{v^2} \frac{\partial v}{\partial P} = \frac{\kappa \rho^2}{\rho^0 (1 + 10\kappa P)}\]

* 
  \[v = \frac{1}{\rho} = \frac{1}{\rho^0} \left[ 1 - 0.1ln(1 + 10P \kappa) \right] =  \frac{1 - 0.1ln(1 + 10P \kappa)}{m_0} (m_1 + m_2T + m_3 T^2)\]

* 
  \[\left( \frac{\partial \rho}{\partial T} \right)_P = -\frac{1}{v^2} \frac{\partial v}{\partial T} =  -\rho^2 \frac{1 - 0.1ln(1 + 10P \kappa)}{m_0} (m_2 + 2m_3T)\]







**Parameters:**


* `T` [K] 
* `P` [Pa] 
* `dRhodP` [kg/m3/Pa] 
* `dRhodT` [kg/m3/K] 
* `kappa` [1/Pa] Isothermal Compressibility 
* `beta` [1/K] Isobaric Expansivity 



**Returns:**







        

<hr>



### function Rho\_Solid [1/2]

_Halite density as function of (T, P). See equation (1-3) of Driesner2007Part2._ 
```C++
double xThermal::NaCl::cNaCl::Rho_Solid (
    const  double & T,
    const  double & P
) 
```





**Parameters:**


* `T` [K] 
* `P` [Pa] 



**Returns:**

double [kg/m3] 





        

<hr>



### function Rho\_Solid [2/2]

```C++
void xThermal::NaCl::cNaCl::Rho_Solid (
    const  double & T,
    const  double & P,
    double & rho,
    double & dRhodP,
    double & dRhodT,
    double & kappa,
    double & beta
) 
```



Calculate halite density \(\rho\) and derivatives \(\left( \frac{\partial \rho}{\partial P} \right)_T\), \(\left( \frac{\partial \rho}{\partial T} \right)_P\)



* \(\rho = l_0 + l_1T + l_2 T^2 + (l_3 + l_4 e^{T/l_5}) P\): Halite density as function of (T, P). See equation (1-3) of Driesner2007Part2.
* 
  \[\left( \frac{\partial \rho}{\partial P} \right)_T = l = l_3 + l_4 e^{T/l_5}\]

* 
  \[\left( \frac{\partial \rho}{\partial T} \right)_P = l_1+ 2l_2T+ \frac{l_4}{l_5}e^{T/l_5}P\]







**Parameters:**


* `T` [K] 
* `P` [Pa] 
* `dRhodP` [kg/m3/Pa] 
* `dRhodT` [kg/m3/K] 
* `kappa` [1/Pa] Isothermal Compressibility 
* `beta` [1/K] Isobaric Expansivity 



**Returns:**







        

<hr>



### function Sublimation\_T [1/2]

_Sublimation temperature as a function of pressure. See also_ Sublimation\_P _._
```C++
double xThermal::NaCl::cNaCl::Sublimation_T (
    const  double & P
) 
```





**Parameters:**


* `P` [Pa] 



**Returns:**

double [K] 





        

<hr>



### function Sublimation\_T [2/2]

```C++
void xThermal::NaCl::cNaCl::Sublimation_T (
    const std::vector< double > P,
    std::vector< double > & res
) 
```




<hr>



### function Sublimation\_p [1/2]

_Sublimation curve. See equation 2, table 3 and figure 4 of Driesner2007Part1._ 
```C++
double xThermal::NaCl::cNaCl::Sublimation_p (
    const  double & T
) 
```



![Image](NaCl/Melting_Sublimation_Boiling_NaCl.svg)





**Parameters:**


* `T` [K] 



**Returns:**

double [Pa] 





        

<hr>



### function Sublimation\_p [2/2]

```C++
void xThermal::NaCl::cNaCl::Sublimation_p (
    const std::vector< double > T,
    std::vector< double > & res
) 
```




<hr>



### function T\_Vapor [1/2]

_Calculate vapor temperature of NaCl. See also_ [_**Sublimation\_T**_](classxThermal_1_1NaCl_1_1cNaCl.md#function-sublimation_t-12) _and_[_**Boiling\_T**_](classxThermal_1_1NaCl_1_1cNaCl.md#function-boiling_t-12) _._
```C++
double xThermal::NaCl::cNaCl::T_Vapor (
    const  double & P
) 
```





**Parameters:**


* `P` [Pa] 



**Returns:**

double [K] 





        

<hr>



### function T\_Vapor [2/2]

```C++
void xThermal::NaCl::cNaCl::T_Vapor (
    const std::vector< double > P,
    std::vector< double > & res
) 
```




<hr>



### function T\_critical 

```C++
inline virtual double xThermal::NaCl::cNaCl::T_critical () 
```



Get the triple point temperature in K 


        
Implements [*xThermal::cxThermal::T\_critical*](classxThermal_1_1cxThermal.md#function-t_critical)


<hr>



### function Tmax 

```C++
inline virtual double xThermal::NaCl::cNaCl::Tmax () 
```



Get the minimum temperature in K 


        
Implements [*xThermal::cxThermal::Tmax*](classxThermal_1_1cxThermal.md#function-tmax)


<hr>



### function Tmin 

```C++
inline virtual double xThermal::NaCl::cNaCl::Tmin () 
```



Get the minimum temperature in K 


        
Implements [*xThermal::cxThermal::Tmin*](classxThermal_1_1cxThermal.md#function-tmin)


<hr>



### function Ttriple 

```C++
inline virtual double xThermal::NaCl::cNaCl::Ttriple () 
```



Get the maximum pressure in Pa 


        
Implements [*xThermal::cxThermal::Ttriple*](classxThermal_1_1cxThermal.md#function-ttriple)


<hr>



### function UpdateState\_TPX 

_Update state by given T,p as input pair._ 
```C++
virtual void xThermal::NaCl::cNaCl::UpdateState_TPX (
    ThermodynamicProperties & props,
    const  double & T,
    const  double & p,
    const  double & X=0
) 
```



Return the molar mass in kg/mol




**Warning:**

Here we only consider liquid/solid region, because for geological conditions, the vapor region will never be reached. The minimum valid pressure [**pmin**](classxThermal_1_1NaCl_1_1cNaCl.md#function-pmin) is limited to 1E5.




**Parameters:**


* `T` [K] 
* `p` [Pa] 




        
Implements [*xThermal::cxThermal::UpdateState\_TPX*](classxThermal_1_1cxThermal.md#function-updatestate_tpx-14)


<hr>



### function cNaCl [1/2]

```C++
xThermal::NaCl::cNaCl::cNaCl (
    std::string backend_H2O
) 
```




<hr>



### function cNaCl [2/2]

```C++
xThermal::NaCl::cNaCl::cNaCl (
    const  cNaCl & salt
) 
```




<hr>



### function molar\_mass 

```C++
inline virtual double xThermal::NaCl::cNaCl::molar_mass () 
```



Return the critical mass density in \(kg/m^3\) 


        
Implements [*xThermal::cxThermal::molar\_mass*](classxThermal_1_1cxThermal.md#function-molar_mass)


<hr>



### function name 

```C++
inline virtual std::string xThermal::NaCl::cNaCl::name () 
```



Name of the model 


        
Implements [*xThermal::cxThermal::name*](classxThermal_1_1cxThermal.md#function-name)


<hr>



### function name\_backend 

```C++
std::string xThermal::NaCl::cNaCl::name_backend () 
```




<hr>



### function p\_critical 

```C++
inline virtual double xThermal::NaCl::cNaCl::p_critical () 
```



Return the critical temperature in K 


        
Implements [*xThermal::cxThermal::p\_critical*](classxThermal_1_1cxThermal.md#function-p_critical)


<hr>



### function pmax 

```C++
inline virtual double xThermal::NaCl::cNaCl::pmax () 
```



Get the minimum pressure in Pa. Same as NaCl-H2O 


        
Implements [*xThermal::cxThermal::pmax*](classxThermal_1_1cxThermal.md#function-pmax)


<hr>



### function pmin 

```C++
inline virtual double xThermal::NaCl::cNaCl::pmin () 
```



Get the maximum temperature in K.Same as NaCl-H2O 


        
Implements [*xThermal::cxThermal::pmin*](classxThermal_1_1cxThermal.md#function-pmin)


<hr>



### function rhomass\_critical 

```C++
inline virtual double xThermal::NaCl::cNaCl::rhomass_critical () 
```



Return the critical pressure in Pa 


        
Implements [*xThermal::cxThermal::rhomass\_critical*](classxThermal_1_1cxThermal.md#function-rhomass_critical)


<hr>



### function ~cNaCl 

```C++
xThermal::NaCl::cNaCl::~cNaCl () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/NaCl/NaCl.h`

