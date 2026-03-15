

# File test.cxx



[**FileList**](files.md) **>** [**H2O**](dir_22b94ec4f3699d3ab17d67318090f0eb.md) **>** [**test.cxx**](H2O_2test_8cxx.md)

[Go to the source code of this file](H2O_2test_8cxx_source.md)



* `#include "IAPS84.h"`
* `#include "IAPWS95_CoolProp.h"`
* `#include "IAPWS95.h"`
* `#include "IAPWS-IF97.h"`
* `#include "LookUpTableForestI.H"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  [**int**](classxThermal_1_1xThermalError.md) | [**lutinfo**](#function-lutinfo) ([**int**](classxThermal_1_1xThermalError.md) argc, [**char**](classxThermal_1_1xThermalError.md) \*\* argv) <br> |
|  [**int**](classxThermal_1_1xThermalError.md) | [**main**](#function-main) ([**int**](classxThermal_1_1xThermalError.md) argc, [**char**](classxThermal_1_1xThermalError.md) \*\* argv) <br> |
|  [**int**](classxThermal_1_1xThermalError.md) | [**test\_LUT\_lookup**](#function-test_lut_lookup) ([**int**](classxThermal_1_1xThermalError.md) argc, [**char**](classxThermal_1_1xThermalError.md) \*\* argv, [**xThermal::cxThermal**](classxThermal_1_1cxThermal.md) & thermo) <br> |
|  [**int**](classxThermal_1_1xThermalError.md) | [**test\_LUTgen**](#function-test_lutgen) ([**int**](classxThermal_1_1xThermalError.md) argc, [**char**](classxThermal_1_1xThermalError.md) \*\* argv, [**xThermal::cxThermal**](classxThermal_1_1cxThermal.md) & thermo) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**test\_PH**](#function-test_ph) ([**xThermal::cxThermal**](classxThermal_1_1cxThermal.md) \* thermo) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**test\_Performance**](#function-test_performance) () <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**test\_PhaseDiagram**](#function-test_phasediagram) ([**xThermal::cxThermal**](classxThermal_1_1cxThermal.md) \* thermo, [**double**](classxThermal_1_1xThermalError.md) T, [**double**](classxThermal_1_1xThermalError.md) p) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**test\_Rho**](#function-test_rho) ([**xThermal::cxThermal**](classxThermal_1_1cxThermal.md) & thermo, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**double**](classxThermal_1_1xThermalError.md) & Rho) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**test\_compressibility\_twophase\_in\_porousmedium**](#function-test_compressibility_twophase_in_porousmedium) ([**xThermal::cxThermal**](classxThermal_1_1cxThermal.md) \* thermo, [**double**](classxThermal_1_1xThermalError.md) S=0.5, [**double**](classxThermal_1_1xThermalError.md) porosity=0.1, [**double**](classxThermal_1_1xThermalError.md) rho\_r=2000, [**double**](classxThermal_1_1xThermalError.md) cp\_r=1000) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**test\_constants**](#function-test_constants) ([**xThermal::cxThermal**](classxThermal_1_1cxThermal.md) \* thermo) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**test\_viscosity**](#function-test_viscosity) ([**xThermal::cxThermal**](classxThermal_1_1cxThermal.md) \* thermo) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**test\_water\_profile**](#function-test_water_profile) ([**xThermal::cxThermal**](classxThermal_1_1cxThermal.md) \* thermo, [**int**](classxThermal_1_1xThermalError.md) nT=100, [**int**](classxThermal_1_1xThermalError.md) nP=100) <br> |




























## Public Functions Documentation




### function lutinfo 

```C++
int lutinfo (
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



### function test\_LUT\_lookup 

```C++
int test_LUT_lookup (
    int argc,
    char ** argv,
    xThermal::cxThermal & thermo
) 
```




<hr>



### function test\_LUTgen 

```C++
int test_LUTgen (
    int argc,
    char ** argv,
    xThermal::cxThermal & thermo
) 
```




<hr>



### function test\_PH 

```C++
void test_PH (
    xThermal::cxThermal * thermo
) 
```




<hr>



### function test\_Performance 

```C++
void test_Performance () 
```




<hr>



### function test\_PhaseDiagram 

```C++
void test_PhaseDiagram (
    xThermal::cxThermal * thermo,
    double T,
    double p
) 
```




<hr>



### function test\_Rho 

```C++
void test_Rho (
    xThermal::cxThermal & thermo,
    const  double & T,
    const  double & P,
    double & Rho
) 
```




<hr>



### function test\_compressibility\_twophase\_in\_porousmedium 

```C++
void test_compressibility_twophase_in_porousmedium (
    xThermal::cxThermal * thermo,
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
    xThermal::cxThermal * thermo
) 
```




<hr>



### function test\_viscosity 

```C++
void test_viscosity (
    xThermal::cxThermal * thermo
) 
```




<hr>



### function test\_water\_profile 

```C++
void test_water_profile (
    xThermal::cxThermal * thermo,
    int nT=100,
    int nP=100
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2O/test.cxx`

