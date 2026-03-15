

# Class xThermal::PROST::cIAPS84



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**PROST**](namespacexThermal_1_1PROST.md) **>** [**cIAPS84**](classxThermal_1_1PROST_1_1cIAPS84.md)








Inherits the following classes: [xThermal::cxThermal](classxThermal_1_1cxThermal.md)
























## Public Attributes inherited from xThermal::cxThermal

See [xThermal::cxThermal](classxThermal_1_1cxThermal.md)

| Type | Name |
| ---: | :--- |
|  [**bool**](classxThermal_1_1xThermalError.md) | [**m\_isShowProgressBar**](classxThermal_1_1cxThermal.md#variable-m_isshowprogressbar)  <br> |






























## Public Functions

| Type | Name |
| ---: | :--- |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_T**](#function-boiling_t-13) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_T**](#function-boiling_t-23) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**double**](classxThermal_1_1xThermalError.md) & rho\_l, [**double**](classxThermal_1_1xThermalError.md) & rho\_v) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_T**](#function-boiling_t-33) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_p**](#function-boiling_p-13) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_p**](#function-boiling_p-23) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**double**](classxThermal_1_1xThermalError.md) & rho\_l, [**double**](classxThermal_1_1xThermalError.md) & rho\_v) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Boiling\_p**](#function-boiling_p-33) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**T\_critical**](#function-t_critical) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Tmax**](#function-tmax) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Tmin**](#function-tmin) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**Ttriple**](#function-ttriple) () <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**UpdateState\_HPX**](#function-updatestate_hpx) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & H, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br> |
| virtual [**void**](classxThermal_1_1xThermalError.md) | [**UpdateState\_TPX**](#function-updatestate_tpx-12) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br> |
|  [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) | [**UpdateState\_TPX**](#function-updatestate_tpx-22) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br>_Calculate water properties for given T, p._  |
|   | [**cIAPS84**](#function-ciaps84-12) () <br> |
|   | [**cIAPS84**](#function-ciaps84-22) ([**const**](classxThermal_1_1xThermalError.md) [**cIAPS84**](classxThermal_1_1PROST_1_1cIAPS84.md) & water) <br> |
| virtual [**PhaseRegion**](namespacexThermal.md#enum-phaseregion) | [**findPhaseRegion\_TPX**](#function-findphaseregion_tpx) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X=0) <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**molar\_mass**](#function-molar_mass) () <br> |
| virtual std::string | [**name**](#function-name) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**p\_critical**](#function-p_critical) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**pmax**](#function-pmax) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**pmin**](#function-pmin) () <br> |
| virtual [**double**](classxThermal_1_1xThermalError.md) | [**rhomass\_critical**](#function-rhomass_critical) () <br> |
|   | [**~cIAPS84**](#function-ciaps84) () <br> |


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




### function Boiling\_T [1/3]

```C++
virtual double xThermal::PROST::cIAPS84::Boiling_T (
    const  double & p
) 
```



Calculate boiling temperature [K] of water for a given p [Pa] 

**Parameters:**


* `p` 



**Returns:**







        
Implements [*xThermal::cxThermal::Boiling\_T*](classxThermal_1_1cxThermal.md#function-boiling_t-13)


<hr>



### function Boiling\_T [2/3]

```C++
virtual double xThermal::PROST::cIAPS84::Boiling_T (
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







        
Implements [*xThermal::cxThermal::Boiling\_T*](classxThermal_1_1cxThermal.md#function-boiling_t-23)


<hr>



### function Boiling\_T [3/3]

```C++
virtual double xThermal::PROST::cIAPS84::Boiling_T (
    const  double & p,
    ThermodynamicProperties & props
) 
```



Implements [*xThermal::cxThermal::Boiling\_T*](classxThermal_1_1cxThermal.md#function-boiling_t-33)


<hr>



### function Boiling\_p [1/3]

```C++
virtual double xThermal::PROST::cIAPS84::Boiling_p (
    const  double & T
) 
```



Calculate boiling pressure [Pa] of water for a given T [K] 

**Parameters:**


* `T` 



**Returns:**







        
Implements [*xThermal::cxThermal::Boiling\_p*](classxThermal_1_1cxThermal.md#function-boiling_p-13)


<hr>



### function Boiling\_p [2/3]

```C++
virtual double xThermal::PROST::cIAPS84::Boiling_p (
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







        
Implements [*xThermal::cxThermal::Boiling\_p*](classxThermal_1_1cxThermal.md#function-boiling_p-23)


<hr>



### function Boiling\_p [3/3]

```C++
virtual double xThermal::PROST::cIAPS84::Boiling_p (
    const  double & T,
    ThermodynamicProperties & props
) 
```



Implements [*xThermal::cxThermal::Boiling\_p*](classxThermal_1_1cxThermal.md#function-boiling_p-33)


<hr>



### function T\_critical 

```C++
inline virtual double xThermal::PROST::cIAPS84::T_critical () 
```



Get the triple point temperature in K 


        
Implements [*xThermal::cxThermal::T\_critical*](classxThermal_1_1cxThermal.md#function-t_critical)


<hr>



### function Tmax 

```C++
inline virtual double xThermal::PROST::cIAPS84::Tmax () 
```



Get the minimum temperature in K 


        
Implements [*xThermal::cxThermal::Tmax*](classxThermal_1_1cxThermal.md#function-tmax)


<hr>



### function Tmin 

```C++
inline virtual double xThermal::PROST::cIAPS84::Tmin () 
```



Get the minimum temperature in K 


        
Implements [*xThermal::cxThermal::Tmin*](classxThermal_1_1cxThermal.md#function-tmin)


<hr>



### function Ttriple 

```C++
inline virtual double xThermal::PROST::cIAPS84::Ttriple () 
```



Get the maximum pressure in Pa 


        
Implements [*xThermal::cxThermal::Ttriple*](classxThermal_1_1cxThermal.md#function-ttriple)


<hr>



### function UpdateState\_HPX 

```C++
virtual void xThermal::PROST::cIAPS84::UpdateState_HPX (
    ThermodynamicProperties & props,
    const  double & H,
    const  double & p,
    const  double & X=0
) 
```



Implements [*xThermal::cxThermal::UpdateState\_HPX*](classxThermal_1_1cxThermal.md#function-updatestate_hpx-14)


<hr>



### function UpdateState\_TPX [1/2]

```C++
virtual void xThermal::PROST::cIAPS84::UpdateState_TPX (
    ThermodynamicProperties & props,
    const  double & T,
    const  double & p,
    const  double & X=0
) 
```



Calculate water properties based on IAPS84 and implementation of PROST.




**Warning:**

The viscosity will be zero when T&gt;900 deg.C, see line 332 of iaps.c in PROST.




**Warning:**

The T,p validation in the matlab code (see [**Line**](structxThermal_1_1Line.md) 39-43 of water\_tp\_IAPS84.m) is not complete, see also PROST code ([**Line**](structxThermal_1_1Line.md) 1235 of iaps.c), so for e.g. T = -11 deg.C(T\_star of viscosity correction for T=1 deg.C and X=0.7 kg/kg) p=100E5 Pa, the matlab code still give output, but it is in an invalid region.




**Parameters:**


* `props` 
* `T` [K] 
* `P` [Pa] 
* `X` [default 0 for H2O] 




        
Implements [*xThermal::cxThermal::UpdateState\_TPX*](classxThermal_1_1cxThermal.md#function-updatestate_tpx-14)


<hr>



### function UpdateState\_TPX [2/2]

_Calculate water properties for given T, p._ 
```C++
ThermodynamicProperties xThermal::PROST::cIAPS84::UpdateState_TPX (
    const  double & T,
    const  double & p,
    const  double & X=0
) 
```





**Parameters:**


* `stateArray` 
* `T` 
* `p` 
* `N` 




        

<hr>



### function cIAPS84 [1/2]

```C++
xThermal::PROST::cIAPS84::cIAPS84 () 
```




<hr>



### function cIAPS84 [2/2]

```C++
xThermal::PROST::cIAPS84::cIAPS84 (
    const  cIAPS84 & water
) 
```




<hr>



### function findPhaseRegion\_TPX 

```C++
virtual PhaseRegion xThermal::PROST::cIAPS84::findPhaseRegion_TPX (
    const  double & T,
    const  double & p,
    const  double & X=0
) 
```



Return the molar mass in kg/mol 


        
Implements [*xThermal::cxThermal::findPhaseRegion\_TPX*](classxThermal_1_1cxThermal.md#function-findphaseregion_tpx)


<hr>



### function molar\_mass 

```C++
inline virtual double xThermal::PROST::cIAPS84::molar_mass () 
```



Return the critical mass density in \(kg/m^3\) 


        
Implements [*xThermal::cxThermal::molar\_mass*](classxThermal_1_1cxThermal.md#function-molar_mass)


<hr>



### function name 

```C++
inline virtual std::string xThermal::PROST::cIAPS84::name () 
```



Name of the model 


        
Implements [*xThermal::cxThermal::name*](classxThermal_1_1cxThermal.md#function-name)


<hr>



### function p\_critical 

```C++
inline virtual double xThermal::PROST::cIAPS84::p_critical () 
```



Return the critical temperature in K 


        
Implements [*xThermal::cxThermal::p\_critical*](classxThermal_1_1cxThermal.md#function-p_critical)


<hr>



### function pmax 

```C++
inline virtual double xThermal::PROST::cIAPS84::pmax () 
```



Get the minimum pressure in Pa 


        
Implements [*xThermal::cxThermal::pmax*](classxThermal_1_1cxThermal.md#function-pmax)


<hr>



### function pmin 

```C++
inline virtual double xThermal::PROST::cIAPS84::pmin () 
```



Get the maximum temperature in K 


        
Implements [*xThermal::cxThermal::pmin*](classxThermal_1_1cxThermal.md#function-pmin)


<hr>



### function rhomass\_critical 

```C++
inline virtual double xThermal::PROST::cIAPS84::rhomass_critical () 
```



Return the critical pressure in Pa 


        
Implements [*xThermal::cxThermal::rhomass\_critical*](classxThermal_1_1cxThermal.md#function-rhomass_critical)


<hr>



### function ~cIAPS84 

```C++
xThermal::PROST::cIAPS84::~cIAPS84 () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2O/IAPS84.h`

