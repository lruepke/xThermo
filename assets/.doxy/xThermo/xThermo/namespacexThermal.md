

# Namespace xThermal



[**Namespace List**](namespaces.md) **>** [**xThermal**](namespacexThermal.md)



_Namespace of_ [_**xThermal**_](namespacexThermal.md) _library._














## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**H2ONaCl**](namespacexThermal_1_1H2ONaCl.md) <br> |
| namespace | [**IAPWS95**](namespacexThermal_1_1IAPWS95.md) <br> |
| namespace | [**IAPWS\_IF97**](namespacexThermal_1_1IAPWS__IF97.md) <br> |
| namespace | [**NaCl**](namespacexThermal_1_1NaCl.md) <br> |
| namespace | [**PROST**](namespacexThermal_1_1PROST.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| struct | [**COLOR**](structxThermal_1_1COLOR.md) <br> |
| struct | [**CONSTENTS\_Thermo**](structxThermal_1_1CONSTENTS__Thermo.md) <br>_为了实现多个H2O EOS的backend，必须使用虚函数调用相应的参数，比如T\_critical()，但是在子类中频繁调用函数的性能肯定很低，所以将所有热力学常数打包为一个struct类型，作为子类的成员变量，然后在构造函数中进行初始化，后面需要常数的地方直接访问成员变量即可，可提高性能。_  |
| struct | [**DeformLinearMesh**](structxThermal_1_1DeformLinearMesh.md) <br> |
| struct | [**Head\_AMR\_LUT**](structxThermal_1_1Head__AMR__LUT.md) <br> |
| struct | [**Line**](structxThermal_1_1Line.md) <br> |
| struct | [**Line\_slice**](structxThermal_1_1Line__slice.md) <br> |
| struct | [**PhaseBoundaries**](structxThermal_1_1PhaseBoundaries.md) <br> |
| struct | [**PhaseRegion\_Slice**](structxThermal_1_1PhaseRegion__Slice.md) <br> |
| struct | [**PhysicalDimension**](structxThermal_1_1PhysicalDimension.md) <br>_Physical dimension of a quantity can be expressed in terms of 7 basic SI unit, e.g. dimension of density is_ \(kg/m^3\) _, can be expressed in vector of [1, -3, 0, 0, 0, 0, 0]. This idea comes from OpenFOAM._ |
| struct | [**Point**](structxThermal_1_1Point.md) <br> |
| struct | [**Point\_slice**](structxThermal_1_1Point__slice.md) <br> |
| struct | [**Polygon\_slice**](structxThermal_1_1Polygon__slice.md) <br> |
| struct | [**Surface**](structxThermal_1_1Surface.md) <br> |
| struct | [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) <br> |
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
|  std::string | [**extname\_file**](#function-extname_file) ([**const**](classxThermal_1_1xThermalError.md) std::string & filepath) <br> |
|  std::string | [**filename\_without\_ext**](#function-filename_without_ext) ([**const**](classxThermal_1_1xThermalError.md) std::string & filepath) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**fill\_prop2data**](#function-fill_prop2data) ([**cxThermal**](classxThermal_1_1cxThermal.md) \* pEOS, [**const**](classxThermal_1_1xThermalError.md) [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) \* prop, [**const**](classxThermal_1_1xThermalError.md) std::map&lt; [**int**](classxThermal_1_1xThermalError.md), [**propInfo**](structxThermal_1_1propInfo.md) &gt; & update\_which\_props, [**double**](classxThermal_1_1xThermalError.md) \* data) <br> |
|  std::vector&lt; std::string &gt; | [**string\_split**](#function-string_split) ([**const**](classxThermal_1_1xThermalError.md) std::string & s, [**const**](classxThermal_1_1xThermalError.md) std::string & delimiter) <br> |




























## Public Types Documentation




### typedef AttributeError 

```C++
typedef xThermalError< xThermalBaseError::eAttribute > xThermal::AttributeError;
```




<hr>



### typedef CompositionError 

```C++
typedef ValueErrorSpec< xThermalBaseError::eComposition > xThermal::CompositionError;
```




<hr>



### typedef DirectorySizeError 

```C++
typedef xThermalError< xThermalBaseError::eDirectorySize > xThermal::DirectorySizeError;
```




<hr>



### typedef HandleError 

```C++
typedef xThermalError< xThermalBaseError::eHandle > xThermal::HandleError;
```




<hr>



### enum INDEX\_FLUID 

```C++
enum xThermal::INDEX_FLUID {
    Fluid_Unknown = -1,
    FLUID_Water,
    FLUID_H2O_NaCl
};
```




<hr>



### enum INPUT\_PAIR 

```C++
enum xThermal::INPUT_PAIR {
    TPX,
    HPX
};
```




<hr>



### typedef InputError 

```C++
typedef ValueErrorSpec< xThermalBaseError::eInput > xThermal::InputError;
```




<hr>



### typedef KeyError 

```C++
typedef xThermalError< xThermalBaseError::eKey > xThermal::KeyError;
```




<hr>



### typedef MAP\_PhaseRegion2Name 

_Map of phase region index (enum) to phase region name (string)._ 
```C++
typedef std::map<PhaseRegion,std::string> xThermal::MAP_PhaseRegion2Name;
```




<hr>



### typedef NotAvailableError 

```C++
typedef ValueErrorSpec< xThermalBaseError::eNotAvailable > xThermal::NotAvailableError;
```




<hr>



### typedef NotImplementedError 

```C++
typedef xThermalError< xThermalBaseError::eNotImplemented > xThermal::NotImplementedError;
```




<hr>



### typedef OutOfRangeError 

```C++
typedef xThermalError< xThermalBaseError::eOutOfRange > xThermal::OutOfRangeError;
```




<hr>



### enum PhaseRegion 

_Phase region index of H2O, NaCl, H2ONaCl system._ 
```C++
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

```C++
enum xThermal::PhaseType {
    Vapor,
    Liquid
};
```



only used as input argument for Rho\_phase; 


        

<hr>



### enum SOLVE\_SATURATED\_PorT 

_Define flag of solving saturated pressure or saturated temperature._ 
```C++
enum xThermal::SOLVE_SATURATED_PorT {
    SOLVE_SATURATED_P,
    SOLVE_SATURATED_T
};
```




<hr>



### typedef SolutionError 

```C++
typedef xThermalError< xThermalBaseError::eSolution > xThermal::SolutionError;
```




<hr>



### enum Space\_EOS 

_Which space of the EOS is calculated or defined. For H2O, the available space are Space\_TP and Space\_HP. For H2ONaCl, the available space are Space\_TPX and Space\_HPX._ 
```C++
enum xThermal::Space_EOS {
    Space_TP,
    Space_HP,
    Space_TPX,
    Space_HPX
};
```




<hr>



### typedef UnableToLoadError 

```C++
typedef xThermalError< xThermalBaseError::eUnableToLoad > xThermal::UnableToLoadError;
```




<hr>



### typedef ValueError 

```C++
typedef xThermalError< xThermalBaseError::eValue > xThermal::ValueError;
```




<hr>



### typedef WrongFluidError 

```C++
typedef ValueErrorSpec< xThermalBaseError::eWrongFluid > xThermal::WrongFluidError;
```




<hr>



### enum dimensionType 

_index of dimension vector._ 
```C++
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

```C++
const double xThermal::R;
```



Specific gass constant: \(J/kg/K\) 


        

<hr>



### variable nDimensions 

```C++
const int xThermal::nDimensions;
```




<hr>
## Public Static Attributes Documentation




### variable map\_phase2name 

```C++
MAP_PhaseRegion2Name xThermal::map_phase2name;
```




<hr>
## Public Functions Documentation




### function extname\_file 

```C++
std::string xThermal::extname_file (
    const std::string & filepath
) 
```




<hr>



### function filename\_without\_ext 

```C++
std::string xThermal::filename_without_ext (
    const std::string & filepath
) 
```




<hr>



### function fill\_prop2data 

```C++
void xThermal::fill_prop2data (
    cxThermal * pEOS,
    const  ThermodynamicProperties * prop,
    const std::map< int , propInfo > & update_which_props,
    double * data
) 
```




<hr>



### function string\_split 

```C++
std::vector< std::string > xThermal::string_split (
    const std::string & s,
    const std::string & delimiter
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2O/IAPS84.cpp`

