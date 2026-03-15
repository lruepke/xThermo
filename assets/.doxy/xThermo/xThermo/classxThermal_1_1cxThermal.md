

# Class xThermal::cxThermal



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**cxThermal**](classxThermal_1_1cxThermal.md)



_Top abstract class of the thermodynamic model._ 

* `#include <thermo.h>`





Inherited by the following classes: [xThermal::H2ONaCl::cH2ONaCl](classxThermal_1_1H2ONaCl_1_1cH2ONaCl.md),  [xThermal::IAPWS95::cIAPWS95](classxThermal_1_1IAPWS95_1_1cIAPWS95.md),  [xThermal::IAPWS\_IF97::cIAPWS\_IF97](classxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97.md),  [xThermal::NaCl::cNaCl](classxThermal_1_1NaCl_1_1cNaCl.md),  [xThermal::PROST::cIAPS84](classxThermal_1_1PROST_1_1cIAPS84.md)
















## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**bool**](classxThermal_1_1xThermalError.md) | [**m\_isShowProgressBar**](#variable-m_isshowprogressbar)  <br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_T**](#function-boiling_t-13) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_T**](#function-boiling_t-23) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**double**](classxThermal_1_1xThermalError.md) & rho\_l, [**double**](classxThermal_1_1xThermalError.md) & rho\_v) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_T**](#function-boiling_t-33) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props) <br> |
|  [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) | [**Boiling\_T\_props**](#function-boiling_t_props) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_p**](#function-boiling_p-13) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_p**](#function-boiling_p-23) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**double**](classxThermal_1_1xThermalError.md) & rho\_l, [**double**](classxThermal_1_1xThermalError.md) & rho\_v) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_p**](#function-boiling_p-33) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props) <br> |
|  [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) | [**Boiling\_p\_props**](#function-boiling_p_props) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T) <br> |
|  [**PhaseRegion**](namespacexThermal.md#enum-phaseregion) | [**Rho\_lookup**](#function-rho_lookup) ([**double**](classxThermal_1_1xThermalError.md) & Rho\_estimate, [**double**](classxThermal_1_1xThermalError.md) & Rho\_min, [**double**](classxThermal_1_1xThermalError.md) & Rho\_max, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P) <br>_Only lookup density for given T,P,X from loaded LUT. This information can be a good guess of Rho for IAPWS95 solution of Rho._  |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**T\_critical**](#function-t_critical) () = 0<br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Tmax**](#function-tmax) () = 0<br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Tmin**](#function-tmin) () = 0<br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**Triangulation**](#function-triangulation) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & x\_poly, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & y\_poly, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) pointInMesh, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) dxdy, [**TriMesh**](structxThermal_1_1TriMesh.md) & trimesh, [**bool**](classxThermal_1_1xThermalError.md) normalize=[**true**](classxThermal_1_1xThermalError.md)) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Ttriple**](#function-ttriple) () = 0<br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**UpdateState\_HPX**](#function-updatestate_hpx-14) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & H, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br> |
| virtual [**ThermodynamicPropertiesVector**](structxThermal_1_1ThermodynamicPropertiesVector.md) | [**UpdateState\_HPX**](#function-updatestate_hpx-24) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & H, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & p, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & X, [**bool**](classxThermal_1_1xThermalError.md) isMeshGrid=[**false**](classxThermal_1_1xThermalError.md)) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**UpdateState\_HPX**](#function-updatestate_hpx-34) ([**ThermodynamicPropertiesArray**](structxThermal_1_1ThermodynamicPropertiesArray.md) & stateArray, [**const**](classxThermal_1_1xThermalError.md) [**size\_t**](classxThermal_1_1xThermalError.md) & N, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) \* H, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) \* p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) \* X=[**NULL**](classxThermal_1_1xThermalError.md)) <br> |
|  [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) | [**UpdateState\_HPX**](#function-updatestate_hpx-44) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & H, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**UpdateState\_TPX**](#function-updatestate_tpx-14) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br>_Update state using T, p as independent variables._  |
| virtual [**ThermodynamicPropertiesVector**](structxThermal_1_1ThermodynamicPropertiesVector.md) | [**UpdateState\_TPX**](#function-updatestate_tpx-24) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & T, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & p, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & X, [**bool**](classxThermal_1_1xThermalError.md) isMeshGrid=[**false**](classxThermal_1_1xThermalError.md)) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**UpdateState\_TPX**](#function-updatestate_tpx-34) ([**ThermodynamicPropertiesArray**](structxThermal_1_1ThermodynamicPropertiesArray.md) & stateArray, [**const**](classxThermal_1_1xThermalError.md) [**size\_t**](classxThermal_1_1xThermalError.md) & N, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) \* T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) \* p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) \* X=[**NULL**](classxThermal_1_1xThermalError.md)) <br> |
|  [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) | [**UpdateState\_TPX**](#function-updatestate_tpx-44) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**createLUT\_2D**](#function-createlut_2d-12) ([**double**](classxThermal_1_1xThermalError.md) xy\_min, [**double**](classxThermal_1_1xThermalError.md) xy\_max, [**double**](classxThermal_1_1xThermalError.md) constZ, LOOKUPTABLE\_FOREST::CONST\_WHICH\_VAR const\_which\_var, LOOKUPTABLE\_FOREST::EOS\_ENERGY TorH, [**int**](classxThermal_1_1xThermalError.md) min\_level=4, [**int**](classxThermal_1_1xThermalError.md) max\_level=6, [**int**](classxThermal_1_1xThermalError.md) update\_which\_props=0) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**createLUT\_2D**](#function-createlut_2d-22) ([**double**](classxThermal_1_1xThermalError.md) xmin, [**double**](classxThermal_1_1xThermalError.md) xmax, [**double**](classxThermal_1_1xThermalError.md) ymin, [**double**](classxThermal_1_1xThermalError.md) ymax, [**double**](classxThermal_1_1xThermalError.md) constZ, LOOKUPTABLE\_FOREST::CONST\_WHICH\_VAR const\_which\_var, LOOKUPTABLE\_FOREST::EOS\_ENERGY TorH, [**int**](classxThermal_1_1xThermalError.md) min\_level=4, [**int**](classxThermal_1_1xThermalError.md) max\_level=6, [**int**](classxThermal_1_1xThermalError.md) update\_which\_props=0) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**createLUT\_3D**](#function-createlut_3d-12) ([**double**](classxThermal_1_1xThermalError.md) xyz\_min, [**double**](classxThermal_1_1xThermalError.md) xyz\_max, LOOKUPTABLE\_FOREST::EOS\_ENERGY TorH, [**int**](classxThermal_1_1xThermalError.md) min\_level=4, [**int**](classxThermal_1_1xThermalError.md) max\_level=6, [**int**](classxThermal_1_1xThermalError.md) update\_which\_props=0) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**createLUT\_3D**](#function-createlut_3d-22) ([**double**](classxThermal_1_1xThermalError.md) xmin, [**double**](classxThermal_1_1xThermalError.md) xmax, [**double**](classxThermal_1_1xThermalError.md) ymin, [**double**](classxThermal_1_1xThermalError.md) ymax, [**double**](classxThermal_1_1xThermalError.md) zmin, [**double**](classxThermal_1_1xThermalError.md) zmax, LOOKUPTABLE\_FOREST::EOS\_ENERGY TorH, [**int**](classxThermal_1_1xThermalError.md) min\_level=4, [**int**](classxThermal_1_1xThermalError.md) max\_level=6, [**int**](classxThermal_1_1xThermalError.md) update\_which\_props=0) <br> |
|   | [**cxThermal**](#function-cxthermal) () <br> |
| virtual [**PhaseRegion**](namespacexThermal.md#enum-phaseregion) | [**findPhaseRegion\_TPX**](#function-findphaseregion_tpx) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br> |
|  [**Head\_AMR\_LUT**](structxThermal_1_1Head__AMR__LUT.md) | [**getLutInfo**](#function-getlutinfo) (std::string file\_lut, [**bool**](classxThermal_1_1xThermalError.md) printSummary=[**true**](classxThermal_1_1xThermalError.md)) <br> |
| virtual [**const**](classxThermal_1_1xThermalError.md) std::map&lt; [**int**](classxThermal_1_1xThermalError.md), [**propInfo**](structxThermal_1_1propInfo.md) &gt; & | [**get\_UpdateWhichProps**](#function-get_updatewhichprops) () <br> |
|  [**int**](classxThermal_1_1xThermalError.md) | [**get\_dim\_lut**](#function-get_dim_lut) () <br> |
|  [**int**](classxThermal_1_1xThermalError.md) | [**get\_dim\_lut\_lookup**](#function-get_dim_lut_lookup) () <br> |
| virtual [**int**](classxThermal_1_1xThermalError.md) | [**get\_num\_threads**](#function-get_num_threads) () <br> |
|  [**void**](classxThermal_1_1xThermalError.md) \* | [**get\_pLUT**](#function-get_plut) () <br> |
|  [**void**](classxThermal_1_1xThermalError.md) \* | [**get\_pLUT\_lookup**](#function-get_plut_lookup) () <br> |
|  std::vector&lt; T &gt; | [**linspace**](#function-linspace-12) (T xmin, T xmax, std::size\_t n, [**bool**](classxThermal_1_1xThermalError.md) isLogScale=[**false**](classxThermal_1_1xThermalError.md)) <br>_Make a linearly spaced vector of points._  |
|  std::vector&lt; T &gt; | [**linspace**](#function-linspace-22) (T xmindxxmax) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**loadLUT**](#function-loadlut) ([**const**](classxThermal_1_1xThermalError.md) std::string & filename, [**bool**](classxThermal_1_1xThermalError.md) printStatus=[**true**](classxThermal_1_1xThermalError.md)) <br> |
|  [**LOOKUPTABLE\_FOREST::Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; 2, [**LOOKUPTABLE\_FOREST::FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md)&lt; 2 &gt; &gt; \* | [**lookup**](#function-lookup-14) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & prop, [**double**](classxThermal_1_1xThermalError.md) x, [**double**](classxThermal_1_1xThermalError.md) y) <br> |
|  [**LOOKUPTABLE\_FOREST::Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; 2, [**LOOKUPTABLE\_FOREST::FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md)&lt; 2 &gt; &gt; \* | [**lookup**](#function-lookup-24) ([**double**](classxThermal_1_1xThermalError.md) \* props, [**double**](classxThermal_1_1xThermalError.md) \* xyz\_min\_target, [**double**](classxThermal_1_1xThermalError.md) x, [**double**](classxThermal_1_1xThermalError.md) y, [**bool**](classxThermal_1_1xThermalError.md) is\_cal=[**true**](classxThermal_1_1xThermalError.md)) <br> |
|  [**LOOKUPTABLE\_FOREST::Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; 3, [**LOOKUPTABLE\_FOREST::FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md)&lt; 3 &gt; &gt; \* | [**lookup**](#function-lookup-34) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & prop, [**double**](classxThermal_1_1xThermalError.md) x, [**double**](classxThermal_1_1xThermalError.md) y, [**double**](classxThermal_1_1xThermalError.md) z) <br> |
|  [**LOOKUPTABLE\_FOREST::Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; 3, [**LOOKUPTABLE\_FOREST::FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md)&lt; 3 &gt; &gt; \* | [**lookup**](#function-lookup-44) ([**double**](classxThermal_1_1xThermalError.md) \* props, [**double**](classxThermal_1_1xThermalError.md) \* xyz\_min\_target, [**double**](classxThermal_1_1xThermalError.md) x, [**double**](classxThermal_1_1xThermalError.md) y, [**double**](classxThermal_1_1xThermalError.md) z, [**bool**](classxThermal_1_1xThermalError.md) is\_cal=[**true**](classxThermal_1_1xThermalError.md)) <br> |
|  [**LOOKUPTABLE\_FOREST::Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; 2, [**LOOKUPTABLE\_FOREST::FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md)&lt; 2 &gt; &gt; \* | [**lookup\_only**](#function-lookup_only-12) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & prop, [**double**](classxThermal_1_1xThermalError.md) x, [**double**](classxThermal_1_1xThermalError.md) y) <br> |
|  [**LOOKUPTABLE\_FOREST::Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; 3, [**LOOKUPTABLE\_FOREST::FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md)&lt; 3 &gt; &gt; \* | [**lookup\_only**](#function-lookup_only-22) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & prop, [**double**](classxThermal_1_1xThermalError.md) x, [**double**](classxThermal_1_1xThermalError.md) y, [**double**](classxThermal_1_1xThermalError.md) z) <br> |
|  T | [**max\_vector**](#function-max_vector) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; T &gt; & x) <br> |
|  T | [**mean\_vector**](#function-mean_vector) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; T &gt; & x) <br> |
|  std::vector&lt; [**size\_t**](classxThermal_1_1xThermalError.md) &gt; | [**meshgrid**](#function-meshgrid) (T xmindxxmax, T ymindyymax, std::vector&lt; T &gt; & XX, std::vector&lt; T &gt; & YY) <br> |
|  T | [**min\_vector**](#function-min_vector) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; T &gt; & x) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**molar\_mass**](#function-molar_mass) () = 0<br> |
| virtual std::string | [**name**](#function-name) () = 0<br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**normalizePhaseBoundaries**](#function-normalizephaseboundaries) ([**PhaseBoundaries**](structxThermal_1_1PhaseBoundaries.md) & phaseBoundaries) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**p\_critical**](#function-p_critical) () = 0<br> |
| virtual std::string | [**phase\_name**](#function-phase_name) ([**PhaseRegion**](namespacexThermal.md#enum-phaseregion) phase\_index) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**pmax**](#function-pmax) () = 0<br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**pmin**](#function-pmin) () = 0<br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**rhomass\_critical**](#function-rhomass_critical) () = 0<br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**rhomolar\_critical**](#function-rhomolar_critical) () <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**save\_lut\_to\_binary**](#function-save_lut_to_binary) (std::string filename) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**save\_lut\_to\_vtk**](#function-save_lut_to_vtk) (std::string filename, [**bool**](classxThermal_1_1xThermalError.md) isNormalizeXYZ=[**true**](classxThermal_1_1xThermalError.md)) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**set\_num\_threads**](#function-set_num_threads) ([**int**](classxThermal_1_1xThermalError.md) num\_threads) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**showProgressBar**](#function-showprogressbar) ([**bool**](classxThermal_1_1xThermalError.md) isShow) <br> |
|  T | [**sum\_vector**](#function-sum_vector) ([**const**](classxThermal_1_1xThermalError.md) std::vector&lt; T &gt; & x) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**writeLine2VTU**](#function-writeline2vtu) (std::string vtuFile, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & X, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & Y, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & Z, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) scale\_x=1.0, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) scale\_y=1.0, [**double**](classxThermal_1_1xThermalError.md) scale\_z=1.0) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**writeMeshGrid2VTK**](#function-writemeshgrid2vtk) ([**const**](classxThermal_1_1xThermalError.md) std::string & vtkFile, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & x, [**const**](classxThermal_1_1xThermalError.md) std::string & xTitle, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & y, [**const**](classxThermal_1_1xThermalError.md) std::string & yTitle, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; & z, [**const**](classxThermal_1_1xThermalError.md) std::string & zTitle, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; &gt; & props, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; [**propInfo**](structxThermal_1_1propInfo.md) &gt; & propsInfo, [**bool**](classxThermal_1_1xThermalError.md) isNormalize=[**false**](classxThermal_1_1xThermalError.md)) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**writePhaseBoundaries2VTU**](#function-writephaseboundaries2vtu) (std::string outputPath, [**const**](classxThermal_1_1xThermalError.md) [**PhaseBoundaries**](structxThermal_1_1PhaseBoundaries.md) & phaseBoundaries, [**double**](classxThermal_1_1xThermalError.md) scale\_T=1, [**double**](classxThermal_1_1xThermalError.md) scale\_p=1, [**double**](classxThermal_1_1xThermalError.md) scale\_X=1) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**writeTriMesh2Txt**](#function-writetrimesh2txt) ([**const**](classxThermal_1_1xThermalError.md) [**TriMesh**](structxThermal_1_1TriMesh.md) & mesh, std::string path\_out=".") <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**writeXXYYZZ2VTU**](#function-writexxyyzz2vtu) (std::string vtuFile, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; &gt; & XX, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; &gt; & YY, [**const**](classxThermal_1_1xThermalError.md) std::vector&lt; std::vector&lt; [**double**](classxThermal_1_1xThermalError.md) &gt; &gt; & ZZ, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) scale\_x=1.0, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) scale\_y=1.0, [**double**](classxThermal_1_1xThermalError.md) scale\_z=1.0) <br> |
| virtual  | [**~cxThermal**](#function-cxthermal) () <br> |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  INDEX\_FLUID | [**validateFluid**](#function-validatefluid) ([**const**](classxThermal_1_1xThermalError.md) std::string & fluidName) <br> |


























## Public Attributes Documentation




### variable m\_isShowProgressBar 

```C++
bool xThermal::cxThermal::m_isShowProgressBar;
```




<hr>
## Public Functions Documentation




### function Boiling\_T [1/3]

```C++
inline virtual double xThermal::cxThermal::Boiling_T (
    const  double & p
) 
```



Calculate boiling temperature [K] of water for a given p [Pa] 

**Parameters:**


* `p` 



**Returns:**







        

<hr>



### function Boiling\_T [2/3]

```C++
inline virtual double xThermal::cxThermal::Boiling_T (
    const  double & p,
    double & rho_l,
    double & rho_v
) 
```



Calculate both boiling temperature and density of liquid and vapor phase. 

**Parameters:**


* `p` [Pa] 
* `rho_l` [kg/m3] 
* `rho_v` [kg/m3] 



**Returns:**







        

<hr>



### function Boiling\_T [3/3]

```C++
inline virtual double xThermal::cxThermal::Boiling_T (
    const  double & p,
    ThermodynamicProperties & props
) 
```




<hr>



### function Boiling\_T\_props 

```C++
ThermodynamicProperties xThermal::cxThermal::Boiling_T_props (
    const  double & p
) 
```




<hr>



### function Boiling\_p [1/3]

```C++
inline virtual double xThermal::cxThermal::Boiling_p (
    const  double & T
) 
```



Calculate boiling pressure [Pa] of water for a given T [K] 

**Parameters:**


* `T` 



**Returns:**







        

<hr>



### function Boiling\_p [2/3]

```C++
inline virtual double xThermal::cxThermal::Boiling_p (
    const  double & T,
    double & rho_l,
    double & rho_v
) 
```



Calculate both boiling temperature and density of liquid and vapor phase. 

**Parameters:**


* `T` [K] 
* `rho_l` [kg/m3] 
* `rho_v` [kg/m3] 



**Returns:**







        

<hr>



### function Boiling\_p [3/3]

```C++
inline virtual double xThermal::cxThermal::Boiling_p (
    const  double & T,
    ThermodynamicProperties & props
) 
```




<hr>



### function Boiling\_p\_props 

```C++
ThermodynamicProperties xThermal::cxThermal::Boiling_p_props (
    const  double & T
) 
```




<hr>



### function Rho\_lookup 

_Only lookup density for given T,P,X from loaded LUT. This information can be a good guess of Rho for IAPWS95 solution of Rho._ 
```C++
PhaseRegion xThermal::cxThermal::Rho_lookup (
    double & Rho_estimate,
    double & Rho_min,
    double & Rho_max,
    const  double & T,
    const  double & P
) 
```





**Parameters:**


* `Rho_estimate` Interpolated density in the quad. 
* `Rho_min` Minimum rho in the target (searched) quad 
* `Rho_max` Maximum rho in the target (searched) quad 
* `T` [K] 
* `P` [Pa] 




        

<hr>



### function T\_critical 

```C++
virtual double xThermal::cxThermal::T_critical () = 0
```



Return the critical temperature in K 


        

<hr>



### function Tmax 

```C++
virtual double xThermal::cxThermal::Tmax () = 0
```



Get the maximum temperature in K 


        

<hr>



### function Tmin 

```C++
virtual double xThermal::cxThermal::Tmin () = 0
```



Get the minimum temperature in K 


        

<hr>



### function Triangulation 

```C++
void xThermal::cxThermal::Triangulation (
    const std::vector< double > & x_poly,
    const std::vector< double > & y_poly,
    const  double pointInMesh,
    const  double dxdy,
    TriMesh & trimesh,
    bool normalize=true
) 
```



Triangulate a polygon using the Triangle code, see [https://www.cs.cmu.edu/~quake/triangle.html](https://www.cs.cmu.edu/~quake/triangle.html) 

**Parameters:**


* `x_poly` 
* `y_poly` 
* `trimesh` 




        

<hr>



### function Ttriple 

```C++
virtual double xThermal::cxThermal::Ttriple () = 0
```



Get the triple point temperature in K 


        

<hr>



### function UpdateState\_HPX [1/4]

```C++
inline virtual void xThermal::cxThermal::UpdateState_HPX (
    ThermodynamicProperties & props,
    const  double & H,
    const  double & p,
    const  double & X=0
) 
```




<hr>



### function UpdateState\_HPX [2/4]

```C++
virtual ThermodynamicPropertiesVector xThermal::cxThermal::UpdateState_HPX (
    const std::vector< double > & H,
    const std::vector< double > & p,
    const std::vector< double > & X,
    bool isMeshGrid=false
) 
```




<hr>



### function UpdateState\_HPX [3/4]

```C++
virtual void xThermal::cxThermal::UpdateState_HPX (
    ThermodynamicPropertiesArray & stateArray,
    const  size_t & N,
    const  double * H,
    const  double * p,
    const  double * X=NULL
) 
```




<hr>



### function UpdateState\_HPX [4/4]

```C++
ThermodynamicProperties xThermal::cxThermal::UpdateState_HPX (
    const  double & H,
    const  double & p,
    const  double & X=0
) 
```




<hr>



### function UpdateState\_TPX [1/4]

_Update state using T, p as independent variables._ 
```C++
inline virtual void xThermal::cxThermal::UpdateState_TPX (
    ThermodynamicProperties & props,
    const  double & T,
    const  double & p,
    const  double & X=0
) 
```



Calculate the critical molar density in \(mol/m^3\) 




**Parameters:**


* `T` [K] 
* `p` [Pa] 




        

<hr>



### function UpdateState\_TPX [2/4]

```C++
virtual ThermodynamicPropertiesVector xThermal::cxThermal::UpdateState_TPX (
    const std::vector< double > & T,
    const std::vector< double > & p,
    const std::vector< double > & X,
    bool isMeshGrid=false
) 
```




<hr>



### function UpdateState\_TPX [3/4]

```C++
virtual void xThermal::cxThermal::UpdateState_TPX (
    ThermodynamicPropertiesArray & stateArray,
    const  size_t & N,
    const  double * T,
    const  double * p,
    const  double * X=NULL
) 
```




<hr>



### function UpdateState\_TPX [4/4]

```C++
ThermodynamicProperties xThermal::cxThermal::UpdateState_TPX (
    const  double & T,
    const  double & p,
    const  double & X=0
) 
```




<hr>



### function createLUT\_2D [1/2]

```C++
virtual void xThermal::cxThermal::createLUT_2D (
    double xy_min,
    double xy_max,
    double constZ,
    LOOKUPTABLE_FOREST::CONST_WHICH_VAR const_which_var,
    LOOKUPTABLE_FOREST::EOS_ENERGY TorH,
    int min_level=4,
    int max_level=6,
    int update_which_props=0
) 
```




<hr>



### function createLUT\_2D [2/2]

```C++
virtual void xThermal::cxThermal::createLUT_2D (
    double xmin,
    double xmax,
    double ymin,
    double ymax,
    double constZ,
    LOOKUPTABLE_FOREST::CONST_WHICH_VAR const_which_var,
    LOOKUPTABLE_FOREST::EOS_ENERGY TorH,
    int min_level=4,
    int max_level=6,
    int update_which_props=0
) 
```




<hr>



### function createLUT\_3D [1/2]

```C++
virtual void xThermal::cxThermal::createLUT_3D (
    double xyz_min,
    double xyz_max,
    LOOKUPTABLE_FOREST::EOS_ENERGY TorH,
    int min_level=4,
    int max_level=6,
    int update_which_props=0
) 
```




<hr>



### function createLUT\_3D [2/2]

```C++
virtual void xThermal::cxThermal::createLUT_3D (
    double xmin,
    double xmax,
    double ymin,
    double ymax,
    double zmin,
    double zmax,
    LOOKUPTABLE_FOREST::EOS_ENERGY TorH,
    int min_level=4,
    int max_level=6,
    int update_which_props=0
) 
```




<hr>



### function cxThermal 

```C++
xThermal::cxThermal::cxThermal () 
```




<hr>



### function findPhaseRegion\_TPX 

```C++
inline virtual PhaseRegion xThermal::cxThermal::findPhaseRegion_TPX (
    const  double & T,
    const  double & p,
    const  double & X=0
) 
```




<hr>



### function getLutInfo 

```C++
Head_AMR_LUT xThermal::cxThermal::getLutInfo (
    std::string file_lut,
    bool printSummary=true
) 
```





**Parameters:**


* `file_lut` The .bin file. 




        

<hr>



### function get\_UpdateWhichProps 

```C++
virtual const std::map< int , propInfo > & xThermal::cxThermal::get_UpdateWhichProps () 
```




<hr>



### function get\_dim\_lut 

```C++
inline int xThermal::cxThermal::get_dim_lut () 
```




<hr>



### function get\_dim\_lut\_lookup 

```C++
inline int xThermal::cxThermal::get_dim_lut_lookup () 
```




<hr>



### function get\_num\_threads 

```C++
virtual int xThermal::cxThermal::get_num_threads () 
```




<hr>



### function get\_pLUT 

```C++
void * xThermal::cxThermal::get_pLUT () 
```




<hr>



### function get\_pLUT\_lookup 

```C++
void * xThermal::cxThermal::get_pLUT_lookup () 
```




<hr>



### function linspace [1/2]

_Make a linearly spaced vector of points._ 
```C++
template<typename T>
inline std::vector< T > xThermal::cxThermal::linspace (
    T xmin,
    T xmax,
    std::size_t n,
    bool isLogScale=false
) 
```




<hr>



### function linspace [2/2]

```C++
template<typename T>
inline std::vector< T > xThermal::cxThermal::linspace (
    T xmindxxmax
) 
```




<hr>



### function loadLUT 

```C++
void xThermal::cxThermal::loadLUT (
    const std::string & filename,
    bool printStatus=true
) 
```




<hr>



### function lookup [1/4]

```C++
LOOKUPTABLE_FOREST::Quadrant < 2, LOOKUPTABLE_FOREST::FIELD_DATA < 2 > > * xThermal::cxThermal::lookup (
    ThermodynamicProperties & prop,
    double x,
    double y
) 
```




<hr>



### function lookup [2/4]

```C++
LOOKUPTABLE_FOREST::Quadrant < 2, LOOKUPTABLE_FOREST::FIELD_DATA < 2 > > * xThermal::cxThermal::lookup (
    double * props,
    double * xyz_min_target,
    double x,
    double y,
    bool is_cal=true
) 
```




<hr>



### function lookup [3/4]

```C++
LOOKUPTABLE_FOREST::Quadrant < 3, LOOKUPTABLE_FOREST::FIELD_DATA < 3 > > * xThermal::cxThermal::lookup (
    ThermodynamicProperties & prop,
    double x,
    double y,
    double z
) 
```




<hr>



### function lookup [4/4]

```C++
LOOKUPTABLE_FOREST::Quadrant < 3, LOOKUPTABLE_FOREST::FIELD_DATA < 3 > > * xThermal::cxThermal::lookup (
    double * props,
    double * xyz_min_target,
    double x,
    double y,
    double z,
    bool is_cal=true
) 
```




<hr>



### function lookup\_only [1/2]

```C++
LOOKUPTABLE_FOREST::Quadrant < 2, LOOKUPTABLE_FOREST::FIELD_DATA < 2 > > * xThermal::cxThermal::lookup_only (
    ThermodynamicProperties & prop,
    double x,
    double y
) 
```




<hr>



### function lookup\_only [2/2]

```C++
LOOKUPTABLE_FOREST::Quadrant < 3, LOOKUPTABLE_FOREST::FIELD_DATA < 3 > > * xThermal::cxThermal::lookup_only (
    ThermodynamicProperties & prop,
    double x,
    double y,
    double z
) 
```




<hr>



### function max\_vector 

```C++
template<class T>
inline T xThermal::cxThermal::max_vector (
    const std::vector< T > & x
) 
```




<hr>



### function mean\_vector 

```C++
template<class T>
inline T xThermal::cxThermal::mean_vector (
    const std::vector< T > & x
) 
```




<hr>



### function meshgrid 

```C++
template<typename T>
inline std::vector< size_t > xThermal::cxThermal::meshgrid (
    T xmindxxmax,
    T ymindyymax,
    std::vector< T > & XX,
    std::vector< T > & YY
) 
```



Create a regular mesh grid. 

**Template parameters:**


* `T` Data type 



**Parameters:**


* `xmindxxmax` [xmin, dx, xmax] 
* `ymindyymax` [ymin, dy, ymax] 
* `XX` 1D vector of X in order of [x0, x1, ..., xn; x0, x1, ..., xn; ... ; x0, x1, ..., xn] 
* `YY` 1D vector of Y in order of [y0, y0, ..., y0; y1, y1, ..., y1; ...; yn, yn, ..., yn] 



**Returns:**

[nx, ny] 





        

<hr>



### function min\_vector 

```C++
template<class T>
inline T xThermal::cxThermal::min_vector (
    const std::vector< T > & x
) 
```




<hr>



### function molar\_mass 

```C++
virtual double xThermal::cxThermal::molar_mass () = 0
```



Return the molar mass in kg/mol 


        

<hr>



### function name 

```C++
virtual std::string xThermal::cxThermal::name () = 0
```



Name of the model 


        

<hr>



### function normalizePhaseBoundaries 

```C++
virtual void xThermal::cxThermal::normalizePhaseBoundaries (
    PhaseBoundaries & phaseBoundaries
) 
```



Normalize phase boundaries p,T,X to [0,1] range. This function is designed for VTU format output, because the data visualization in ParaView is default supported x,y,z in 1:1:1 scale. 

**Parameters:**


* `phaseBoundaries` 




        

<hr>



### function p\_critical 

```C++
virtual double xThermal::cxThermal::p_critical () = 0
```



Return the critical pressure in Pa 


        

<hr>



### function phase\_name 

```C++
virtual std::string xThermal::cxThermal::phase_name (
    PhaseRegion phase_index
) 
```




<hr>



### function pmax 

```C++
virtual double xThermal::cxThermal::pmax () = 0
```



Get the maximum pressure in Pa 


        

<hr>



### function pmin 

```C++
virtual double xThermal::cxThermal::pmin () = 0
```



Get the minimum pressure in Pa 


        

<hr>



### function rhomass\_critical 

```C++
virtual double xThermal::cxThermal::rhomass_critical () = 0
```



Return the critical mass density in \(kg/m^3\) 


        

<hr>



### function rhomolar\_critical 

```C++
inline virtual double xThermal::cxThermal::rhomolar_critical () 
```




<hr>



### function save\_lut\_to\_binary 

```C++
virtual void xThermal::cxThermal::save_lut_to_binary (
    std::string filename
) 
```




<hr>



### function save\_lut\_to\_vtk 

```C++
virtual void xThermal::cxThermal::save_lut_to_vtk (
    std::string filename,
    bool isNormalizeXYZ=true
) 
```




<hr>



### function set\_num\_threads 

```C++
virtual void xThermal::cxThermal::set_num_threads (
    int num_threads
) 
```




<hr>



### function showProgressBar 

```C++
inline void xThermal::cxThermal::showProgressBar (
    bool isShow
) 
```




<hr>



### function sum\_vector 

```C++
template<class T>
inline T xThermal::cxThermal::sum_vector (
    const std::vector< T > & x
) 
```




<hr>



### function writeLine2VTU 

```C++
virtual void xThermal::cxThermal::writeLine2VTU (
    std::string vtuFile,
    const std::vector< double > & X,
    const std::vector< double > & Y,
    const std::vector< double > & Z,
    const  double scale_x=1.0,
    const  double scale_y=1.0,
    double scale_z=1.0
) 
```




<hr>



### function writeMeshGrid2VTK 

```C++
virtual void xThermal::cxThermal::writeMeshGrid2VTK (
    const std::string & vtkFile,
    const std::vector< double > & x,
    const std::string & xTitle,
    const std::vector< double > & y,
    const std::string & yTitle,
    const std::vector< double > & z,
    const std::string & zTitle,
    const std::vector< std::vector< double > > & props,
    const std::vector< propInfo > & propsInfo,
    bool isNormalize=false
) 
```




<hr>



### function writePhaseBoundaries2VTU 

```C++
virtual void xThermal::cxThermal::writePhaseBoundaries2VTU (
    std::string outputPath,
    const  PhaseBoundaries & phaseBoundaries,
    double scale_T=1,
    double scale_p=1,
    double scale_X=1
) 
```




<hr>



### function writeTriMesh2Txt 

```C++
virtual void xThermal::cxThermal::writeTriMesh2Txt (
    const  TriMesh & mesh,
    std::string path_out="."
) 
```




<hr>



### function writeXXYYZZ2VTU 

```C++
virtual void xThermal::cxThermal::writeXXYYZZ2VTU (
    std::string vtuFile,
    const std::vector< std::vector< double > > & XX,
    const std::vector< std::vector< double > > & YY,
    const std::vector< std::vector< double > > & ZZ,
    const  double scale_x=1.0,
    const  double scale_y=1.0,
    double scale_z=1.0
) 
```



Write a grid mesh to vtu file. The grid is described by three/to 2-D array. 

**Parameters:**


* `vtuFile` 
* `XX` 
* `YY` 
* `ZZ` 
* `scale_x` 
* `scale_y` 
* `scale_z` 




        

<hr>



### function ~cxThermal 

```C++
virtual xThermal::cxThermal::~cxThermal () 
```




<hr>
## Public Static Functions Documentation




### function validateFluid 

```C++
static INDEX_FLUID xThermal::cxThermal::validateFluid (
    const std::string & fluidName
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/thermo/thermo.h`

