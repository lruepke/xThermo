

# File test.cxx



[**FileList**](files.md) **>** [**H2ONaCl**](dir_0a7d68cdfa0215953309325142e4c674.md) **>** [**test.cxx**](H2ONaCl_2test_8cxx.md)

[Go to the source code of this file](H2ONaCl_2test_8cxx_source.md)

_Test H2ONaCl model._ [More...](#detailed-description)

* `#include "H2ONaCl.h"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  [**void**](classxThermal_1_1xThermalError.md) | [**debug\_HaliteLiquidus\_RhoL**](#function-debug_haliteliquidus_rhol) (H2ONaCl::cH2ONaCl \* sw) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**debug\_VL\_props**](#function-debug_vl_props) (H2ONaCl::cH2ONaCl \* sw) <br> |
|  [**int**](classxThermal_1_1xThermalError.md) | [**lutGen\_2D**](#function-lutgen_2d) ([**int**](classxThermal_1_1xThermalError.md) argc, [**char**](classxThermal_1_1xThermalError.md) \*\* argv) <br> |
|  [**int**](classxThermal_1_1xThermalError.md) | [**lutInfo**](#function-lutinfo) ([**int**](classxThermal_1_1xThermalError.md) argc, [**char**](classxThermal_1_1xThermalError.md) \*\* argv) <br> |
|  [**int**](classxThermal_1_1xThermalError.md) | [**main**](#function-main) ([**int**](classxThermal_1_1xThermalError.md) argc, [**char**](classxThermal_1_1xThermalError.md) \*\* argv) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**test\_HaliteLiquidus**](#function-test_haliteliquidus) (H2ONaCl::cH2ONaCl \* sw) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**test\_LUT\_lookup**](#function-test_lut_lookup) ([**int**](classxThermal_1_1xThermalError.md) argc, [**char**](classxThermal_1_1xThermalError.md) \*\* argv, [**xThermal::cxThermal**](classxThermal_1_1cxThermal.md) & thermo) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**test\_Props\_TPX**](#function-test_props_tpx) (H2ONaCl::cH2ONaCl \* sw) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**test\_T\_HPX**](#function-test_t_hpx) (H2ONaCl::cH2ONaCl \* sw) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**test\_T\_VL**](#function-test_t_vl) (H2ONaCl::cH2ONaCl \* sw) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**test\_compressibility\_multiphase**](#function-test_compressibility_multiphase) (H2ONaCl::cH2ONaCl \* sw, [**double**](classxThermal_1_1xThermalError.md) S=0.5, [**double**](classxThermal_1_1xThermalError.md) porosity=0.1, [**double**](classxThermal_1_1xThermalError.md) rho\_r=2000, [**double**](classxThermal_1_1xThermalError.md) cp\_r=1000) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**test\_constants**](#function-test_constants) (H2ONaCl::cH2ONaCl \* sw) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**test\_criticalCurve**](#function-test_criticalcurve) (H2ONaCl::cH2ONaCl \* sw) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**test\_phaseBoundaries2VTU**](#function-test_phaseboundaries2vtu) (H2ONaCl::cH2ONaCl \* sw) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**test\_phaseRegion\_HPX**](#function-test_phaseregion_hpx) ([**xThermal::cxThermal**](classxThermal_1_1cxThermal.md) \* sw) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**test\_sliceP**](#function-test_slicep) (H2ONaCl::cH2ONaCl \* sw) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**test\_sliceT**](#function-test_slicet) (H2ONaCl::cH2ONaCl \* sw) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**test\_viscosity**](#function-test_viscosity) ([**xThermal::cxThermal**](classxThermal_1_1cxThermal.md) \* sw) <br> |




























## Detailed Description




**Author:**

Zhikui Guo ([zguo@geomar.de](mailto:zguo@geomar.de)) 




**Version:**

0.1 




**Date:**

2022-03-27




**Copyright:**

Copyright (c) 2022 





    
## Public Functions Documentation




### function debug\_HaliteLiquidus\_RhoL 

```C++
void debug_HaliteLiquidus_RhoL (
    H2ONaCl::cH2ONaCl * sw
) 
```




<hr>



### function debug\_VL\_props 

```C++
void debug_VL_props (
    H2ONaCl::cH2ONaCl * sw
) 
```




<hr>



### function lutGen\_2D 

```C++
int lutGen_2D (
    int argc,
    char ** argv
) 
```




<hr>



### function lutInfo 

```C++
int lutInfo (
    int argc,
    char ** argv
) 
```




<hr>



### function main 

```C++
int main (
    int argc,
    char ** argv
) 
```




<hr>



### function test\_HaliteLiquidus 

```C++
void test_HaliteLiquidus (
    H2ONaCl::cH2ONaCl * sw
) 
```




<hr>



### function test\_LUT\_lookup 

```C++
void test_LUT_lookup (
    int argc,
    char ** argv,
    xThermal::cxThermal & thermo
) 
```




<hr>



### function test\_Props\_TPX 

```C++
void test_Props_TPX (
    H2ONaCl::cH2ONaCl * sw
) 
```




<hr>



### function test\_T\_HPX 

```C++
void test_T_HPX (
    H2ONaCl::cH2ONaCl * sw
) 
```




<hr>



### function test\_T\_VL 

```C++
void test_T_VL (
    H2ONaCl::cH2ONaCl * sw
) 
```




<hr>



### function test\_compressibility\_multiphase 

```C++
void test_compressibility_multiphase (
    H2ONaCl::cH2ONaCl * sw,
    double S=0.5,
    double porosity=0.1,
    double rho_r=2000,
    double cp_r=1000
) 
```




<hr>



### function test\_constants 

```C++
void test_constants (
    H2ONaCl::cH2ONaCl * sw
) 
```




<hr>



### function test\_criticalCurve 

```C++
void test_criticalCurve (
    H2ONaCl::cH2ONaCl * sw
) 
```




<hr>



### function test\_phaseBoundaries2VTU 

```C++
void test_phaseBoundaries2VTU (
    H2ONaCl::cH2ONaCl * sw
) 
```




<hr>



### function test\_phaseRegion\_HPX 

```C++
void test_phaseRegion_HPX (
    xThermal::cxThermal * sw
) 
```




<hr>



### function test\_sliceP 

```C++
void test_sliceP (
    H2ONaCl::cH2ONaCl * sw
) 
```




<hr>



### function test\_sliceT 

```C++
void test_sliceT (
    H2ONaCl::cH2ONaCl * sw
) 
```




<hr>



### function test\_viscosity 

```C++
void test_viscosity (
    xThermal::cxThermal * sw
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2ONaCl/test.cxx`

