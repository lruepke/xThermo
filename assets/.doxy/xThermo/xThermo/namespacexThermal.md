

# Namespace xThermal



[**Namespace List**](namespaces.md) **>** [**xThermal**](namespacexThermal.md)



_Namespace of_ [_**xThermal**_](namespacexThermal.md) _library._[More...](#detailed-description)














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**H2ONaCl**](namespacexThermal_1_1H2ONaCl.md) <br> |
| namespace | [**IAPWS95**](namespacexThermal_1_1IAPWS95.md) <br> |
| namespace | [**IAPWS\_IF97**](namespacexThermal_1_1IAPWS__IF97.md) <br> |
| namespace | [**NaCl**](namespacexThermal_1_1NaCl.md) <br> |
| namespace | [**PROST**](namespacexThermal_1_1PROST.md) <br> |
| namespace | [**gsl\_guard**](namespacexThermal_1_1gsl__guard.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| struct | [**COLOR**](structxThermal_1_1COLOR.md) <br> |
| struct | [**CONSTENTS\_Thermo**](structxThermal_1_1CONSTENTS__Thermo.md) <br>_为了实现多个H2O EOS的backend，必须使用虚函数调用相应的参数，比如T\_critical()，但是在子类中频繁调用函数的性能肯定很低，所以将所有热力学常数打包为一个struct类型，作为子类的成员变量，然后在构造函数中进行初始化，后面需要常数的地方直接访问成员变量即可，可提高性能。_  |
| struct | [**DeformLinearMesh**](structxThermal_1_1DeformLinearMesh.md) <br> |
| struct | [**Head\_AMR\_LUT**](structxThermal_1_1Head__AMR__LUT.md) <br> |
| struct | [**Line**](structxThermal_1_1Line.md) <br> |
| struct | [**Line\_slice**](structxThermal_1_1Line__slice.md) <br> |
| class | [**OuterEvaluation**](classxThermal_1_1OuterEvaluation.md) <br>_Reset of a caller's result object at the outermost public evaluation entry (C7)._  |
| struct | [**PhaseBoundaries**](structxThermal_1_1PhaseBoundaries.md) <br> |
| struct | [**PhaseRegion\_Slice**](structxThermal_1_1PhaseRegion__Slice.md) <br> |
| struct | [**PhysicalDimension**](structxThermal_1_1PhysicalDimension.md) <br>_Physical dimension of a quantity can be expressed in terms of 7 basic SI unit, e.g. dimension of density is_ \(kg/m^3\) _, can be expressed in vector of [1, -3, 0, 0, 0, 0, 0]. This idea comes from OpenFOAM._ |
| struct | [**Point**](structxThermal_1_1Point.md) <br> |
| struct | [**Point\_slice**](structxThermal_1_1Point__slice.md) <br> |
| struct | [**Polygon\_slice**](structxThermal_1_1Polygon__slice.md) <br> |
| struct | [**Surface**](structxThermal_1_1Surface.md) <br> |
| struct | [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) <br>_Result of a state evaluation (UpdateState\_TPX / UpdateState\_HPX)._  |
| struct | [**ThermodynamicPropertiesArray**](structxThermal_1_1ThermodynamicPropertiesArray.md) <br> |
| struct | [**ThermodynamicPropertiesVector**](structxThermal_1_1ThermodynamicPropertiesVector.md) <br> |
| struct | [**ThermodynamicProperty**](structxThermal_1_1ThermodynamicProperty.md) <br>_Data struct of a thermodynamic property._  |
| struct | [**TriMesh**](structxThermal_1_1TriMesh.md) <br> |
| class | [**ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md) &lt;errcode&gt;<br> |
| class | [**cxThermal**](classxThermal_1_1cxThermal.md) <br>_Top abstract class of the thermodynamic model._  |
| struct | [**propInfo**](structxThermal_1_1propInfo.md) <br>_Information of a thermodynamic property._  |
| class | [**xThermalBaseError**](classxThermal_1_1xThermalBaseError.md) <br> |
| class | [**xThermalError**](classxThermal_1_1xThermalError.md) &lt;errcode&gt;<br> |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**xThermalError**](classxThermal_1_1xThermalError.md)&lt; xThermalBaseError::eAttribute &gt; | [**AttributeError**](#typedef-attributeerror)  <br> |
| typedef [**ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md)&lt; xThermalBaseError::eComposition &gt; | [**CompositionError**](#typedef-compositionerror)  <br> |
| typedef [**xThermalError**](classxThermal_1_1xThermalError.md)&lt; xThermalBaseError::eDirectorySize &gt; | [**DirectorySizeError**](#typedef-directorysizeerror)  <br> |
| enum  | [**EvaluationStatus**](#enum-evaluationstatus)  <br> |
| typedef [**xThermalError**](classxThermal_1_1xThermalError.md)&lt; xThermalBaseError::eHandle &gt; | [**HandleError**](#typedef-handleerror)  <br> |
| enum  | [**INDEX\_FLUID**](#enum-index_fluid)  <br> |
| enum  | [**INPUT\_PAIR**](#enum-input_pair)  <br> |
| typedef [**ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md)&lt; xThermalBaseError::eInput &gt; | [**InputError**](#typedef-inputerror)  <br> |
| typedef [**xThermalError**](classxThermal_1_1xThermalError.md)&lt; xThermalBaseError::eKey &gt; | [**KeyError**](#typedef-keyerror)  <br> |
| typedef std::map&lt; [**PhaseRegion**](namespacexThermal.md#enum-phaseregion), std::string &gt; | [**MAP\_PhaseRegion2Name**](#typedef-map_phaseregion2name)  <br>_Map of phase region index (enum) to phase region name (string)._  |
| typedef [**ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md)&lt; xThermalBaseError::eNotAvailable &gt; | [**NotAvailableError**](#typedef-notavailableerror)  <br> |
| typedef [**xThermalError**](classxThermal_1_1xThermalError.md)&lt; xThermalBaseError::eNotImplemented &gt; | [**NotImplementedError**](#typedef-notimplementederror)  <br> |
| typedef [**xThermalError**](classxThermal_1_1xThermalError.md)&lt; xThermalBaseError::eOutOfRange &gt; | [**OutOfRangeError**](#typedef-outofrangeerror)  <br> |
| enum  | [**PhaseRegion**](#enum-phaseregion)  <br>_Phase region index of H2O, NaCl, H2ONaCl system._  |
| enum  | [**PhaseType**](#enum-phasetype)  <br> |
| enum  | [**SOLVE\_SATURATED\_PorT**](#enum-solve_saturated_port)  <br>_Define flag of solving saturated pressure or saturated temperature._  |
| typedef [**xThermalError**](classxThermal_1_1xThermalError.md)&lt; xThermalBaseError::eSolution &gt; | [**SolutionError**](#typedef-solutionerror)  <br> |
| enum  | [**Space\_EOS**](#enum-space_eos)  <br>_Which space of the EOS is calculated or defined. For H2O, the available space are Space\_TP and Space\_HP. For H2ONaCl, the available space are Space\_TPX and Space\_HPX._  |
| typedef [**xThermalError**](classxThermal_1_1xThermalError.md)&lt; xThermalBaseError::eUnableToLoad &gt; | [**UnableToLoadError**](#typedef-unabletoloaderror)  <br> |
| typedef [**xThermalError**](classxThermal_1_1xThermalError.md)&lt; xThermalBaseError::eValue &gt; | [**ValueError**](#typedef-valueerror)  <br> |
| typedef [**ValueErrorSpec**](classxThermal_1_1ValueErrorSpec.md)&lt; xThermalBaseError::eWrongFluid &gt; | [**WrongFluidError**](#typedef-wrongfluiderror)  <br> |
| enum  | [**dimensionType**](#enum-dimensiontype)  <br>_index of dimension vector._  |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) | [**R**](#variable-r)   = `461.51805`<br> |
|  [**const**](classxThermal_1_1xThermalError.md) [**int**](classxThermal_1_1xThermalError.md) | [**nDimensions**](#variable-ndimensions)   = `7`<br> |


## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  [**MAP\_PhaseRegion2Name**](namespacexThermal.md#typedef-map_phaseregion2name) | [**map\_phase2name**](#variable-map_phase2name)   = `/* multi line expression */`<br> |














## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**\_\_getattr\_\_**](#function-__getattr__) (name name) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**check\_state\_input**](#function-check_state_input) ([**const**](classxThermal_1_1xThermalError.md) [**char**](classxThermal_1_1xThermalError.md) \* where, [**const**](classxThermal_1_1xThermalError.md) [**char**](classxThermal_1_1xThermalError.md) \* name1, [**double**](classxThermal_1_1xThermalError.md) value1, [**double**](classxThermal_1_1xThermalError.md) p, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) \* X=[**nullptr**](classxThermal_1_1xThermalError.md)) <br>_Input check at the public state-evaluation entry points (C2)._  |
|  [**EvaluationStatus**](namespacexThermal.md#enum-evaluationstatus) | [**evaluation\_status**](#function-evaluation_status) (xThermalBaseError::ErrCode code) <br> |
|  [**const**](classxThermal_1_1xThermalError.md) [**char**](classxThermal_1_1xThermalError.md) \* | [**evaluation\_status\_name**](#function-evaluation_status_name) ([**EvaluationStatus**](namespacexThermal.md#enum-evaluationstatus) status) <br> |
|  std::string | [**extname\_file**](#function-extname_file) ([**const**](classxThermal_1_1xThermalError.md) std::string & filepath) <br> |
|  std::string | [**filename\_without\_ext**](#function-filename_without_ext) ([**const**](classxThermal_1_1xThermalError.md) std::string & filepath) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**fill\_prop2data**](#function-fill_prop2data) ([**cxThermal**](classxThermal_1_1cxThermal.md) \* pEOS, [**const**](classxThermal_1_1xThermalError.md) [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) \* prop, [**const**](classxThermal_1_1xThermalError.md) std::map&lt; [**int**](classxThermal_1_1xThermalError.md), [**propInfo**](structxThermal_1_1propInfo.md) &gt; & update\_which\_props, [**double**](classxThermal_1_1xThermalError.md) \* data) <br> |
|  std::vector&lt; std::string &gt; | [**string\_split**](#function-string_split) ([**const**](classxThermal_1_1xThermalError.md) std::string & s, [**const**](classxThermal_1_1xThermalError.md) std::string & delimiter) <br> |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  [**int**](classxThermal_1_1xThermalError.md) & | [**evaluation\_depth**](#function-evaluation_depth) () <br> |
|  std::string | [**plain\_message**](#function-plain_message) ([**const**](classxThermal_1_1xThermalError.md) [**char**](classxThermal_1_1xThermalError.md) \* what) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**set\_failed\_state**](#function-set_failed_state) ([**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) & props, [**const**](classxThermal_1_1xThermalError.md) std::string & fluidName, [**bool**](classxThermal_1_1xThermalError.md) hpx, [**double**](classxThermal_1_1xThermalError.md) v1, [**double**](classxThermal_1_1xThermalError.md) p, [**double**](classxThermal_1_1xThermalError.md) X) <br> |


























## Detailed Description


 


    
## Public Types Documentation




### typedef AttributeError 

```Python
typedef xThermalError< xThermalBaseError::eAttribute > xThermal::AttributeError;
```




<hr>



### typedef CompositionError 

```Python
typedef ValueErrorSpec< xThermalBaseError::eComposition > xThermal::CompositionError;
```




<hr>



### typedef DirectorySizeError 

```Python
typedef xThermalError< xThermalBaseError::eDirectorySize > xThermal::DirectorySizeError;
```




<hr>



### enum EvaluationStatus 

```Python
enum xThermal::EvaluationStatus {
    Ok = 0,
    InvalidInput,
    OutOfRange,
    NotConverged,
    NotAvailable,
    Failed
};
```



Outcome of a non-throwing state evaluation ([**cxThermal::tryUpdateState\_TPX**](classxThermal_1_1cxThermal.md#function-tryupdatestate_tpx), [**cxThermal::tryUpdateState\_HPX**](classxThermal_1_1cxThermal.md#function-tryupdatestate_hpx)). Each failure corresponds to the [**xThermal**](namespacexThermal.md) exception that the throwing call (UpdateState\_TPX / UpdateState\_HPX) raises for the same input. 


        

<hr>



### typedef HandleError 

```Python
typedef xThermalError< xThermalBaseError::eHandle > xThermal::HandleError;
```




<hr>



### enum INDEX\_FLUID 

```Python
enum xThermal::INDEX_FLUID {
    Fluid_Unknown = -1,
    FLUID_Water,
    FLUID_H2O_NaCl
};
```




<hr>



### enum INPUT\_PAIR 

```Python
enum xThermal::INPUT_PAIR {
    TPX,
    HPX
};
```




<hr>



### typedef InputError 

```Python
typedef ValueErrorSpec< xThermalBaseError::eInput > xThermal::InputError;
```




<hr>



### typedef KeyError 

```Python
typedef xThermalError< xThermalBaseError::eKey > xThermal::KeyError;
```




<hr>



### typedef MAP\_PhaseRegion2Name 

_Map of phase region index (enum) to phase region name (string)._ 
```Python
typedef std::map<PhaseRegion,std::string> xThermal.MAP_PhaseRegion2Name;
```




<hr>



### typedef NotAvailableError 

```Python
typedef ValueErrorSpec< xThermalBaseError::eNotAvailable > xThermal::NotAvailableError;
```




<hr>



### typedef NotImplementedError 

```Python
typedef xThermalError< xThermalBaseError::eNotImplemented > xThermal::NotImplementedError;
```




<hr>



### typedef OutOfRangeError 

```Python
typedef xThermalError< xThermalBaseError::eOutOfRange > xThermal::OutOfRangeError;
```




<hr>



### enum PhaseRegion 

_Phase region index of H2O, NaCl, H2ONaCl system._ 
```Python
enum xThermal::PhaseRegion {
    MixPhaseRegion =-1,
    SinglePhase_L,
    SinglePhase_V,
    SinglePhase_S,
    Supercritical,
    Supercritical_vapor,
    Supercritical_liquid,
    Critical,
    TwoPhase_VL_Water,
    TwoPhase_LH,
    TwoPhase_VH,
    TwoPhase_VL,
    ThreePhase_VLH,
    Unknown,
    NotImposed
};
```




<hr>



### enum PhaseType 

```Python
enum xThermal::PhaseType {
    Vapor,
    Liquid
};
```



only used as input argument for Rho\_phase; 


        

<hr>



### enum SOLVE\_SATURATED\_PorT 

_Define flag of solving saturated pressure or saturated temperature._ 
```Python
enum xThermal::SOLVE_SATURATED_PorT {
    SOLVE_SATURATED_P,
    SOLVE_SATURATED_T
};
```




<hr>



### typedef SolutionError 

```Python
typedef xThermalError< xThermalBaseError::eSolution > xThermal::SolutionError;
```




<hr>



### enum Space\_EOS 

_Which space of the EOS is calculated or defined. For H2O, the available space are Space\_TP and Space\_HP. For H2ONaCl, the available space are Space\_TPX and Space\_HPX._ 
```Python
enum xThermal::Space_EOS {
    Space_TP,
    Space_HP,
    Space_TPX,
    Space_HPX
};
```




<hr>



### typedef UnableToLoadError 

```Python
typedef xThermalError< xThermalBaseError::eUnableToLoad > xThermal::UnableToLoadError;
```




<hr>



### typedef ValueError 

```Python
typedef xThermalError< xThermalBaseError::eValue > xThermal::ValueError;
```




<hr>



### typedef WrongFluidError 

```Python
typedef ValueErrorSpec< xThermalBaseError::eWrongFluid > xThermal::WrongFluidError;
```




<hr>



### enum dimensionType 

_index of dimension vector._ 
```Python
enum xThermal::dimensionType {
    MASS,
    LENGTH,
    TIME,
    TEMPERATURE,
    MOLES,
    CURRENT,
    LUMINOUS_INTENSITY
};
```




<hr>
## Public Attributes Documentation




### variable R 

```Python
const double xThermal.R;
```



Specific gass constant: \(J/kg/K\) 


        

<hr>



### variable nDimensions 

```Python
const int xThermal.nDimensions;
```




<hr>
## Public Static Attributes Documentation




### variable map\_phase2name 

```Python
MAP_PhaseRegion2Name xThermal.map_phase2name;
```




<hr>
## Public Functions Documentation




### function \_\_getattr\_\_ 

```Python
xThermal::__getattr__ (
    name name
) 
```




<hr>



### function check\_state\_input 

_Input check at the public state-evaluation entry points (C2)._ 
```Python
void xThermal::check_state_input (
    const  char * where,
    const  char * name1,
    double value1,
    double p,
    const  double * X=nullptr
) 
```



Throws InputError if an input is not finite, if p &lt;= 0, or (when X is given) if X &lt; 0. These inputs cannot be evaluated by any model; no model-specific range is checked here (each model keeps its own range handling). Called only for the outermost evaluation ([**OuterEvaluation::is\_outer**](classxThermal_1_1OuterEvaluation.md#function-is_outer)), so internal calls are unaffected.




**Parameters:**


* `where` Function name for the message, e.g. "cH2ONaCl::UpdateState\_TPX" 
* `name1` Name of the first input ("T" or "H"), value1 its value 
* `p` Pressure [Pa] 
* `X` Bulk salinity [kg/kg], or nullptr for models without composition (water, NaCl) 




        

<hr>



### function evaluation\_status 

```Python
inline EvaluationStatus xThermal::evaluation_status (
    xThermalBaseError::ErrCode code
) 
```



Status corresponding to an [**xThermal**](namespacexThermal.md) error code 


        

<hr>



### function evaluation\_status\_name 

```Python
inline const  char * xThermal::evaluation_status_name (
    EvaluationStatus status
) 
```



Name of a status, e.g. "OutOfRange" (for messages and logs) 


        

<hr>



### function extname\_file 

```Python
std::string xThermal::extname_file (
    const std::string & filepath
) 
```




<hr>



### function filename\_without\_ext 

```Python
std::string xThermal::filename_without_ext (
    const std::string & filepath
) 
```




<hr>



### function fill\_prop2data 

```Python
void xThermal::fill_prop2data (
    cxThermal * pEOS,
    const  ThermodynamicProperties * prop,
    const std::map< int , propInfo > & update_which_props,
    double * data
) 
```




<hr>



### function string\_split 

```Python
std::vector< std::string > xThermal::string_split (
    const std::string & s,
    const std::string & delimiter
) 
```




<hr>
## Public Static Functions Documentation




### function evaluation\_depth 

```Python
static int & xThermal::evaluation_depth () 
```




<hr>



### function plain\_message 

```Python
static std::string xThermal::plain_message (
    const  char * what
) 
```




<hr>



### function set\_failed\_state 

```Python
static void xThermal::set_failed_state (
    ThermodynamicProperties & props,
    const std::string & fluidName,
    bool hpx,
    double v1,
    double p,
    double X
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/API/python/xThermal/__init__.py`

