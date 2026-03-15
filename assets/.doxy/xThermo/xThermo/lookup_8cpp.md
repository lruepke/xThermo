

# File lookup.cpp



[**FileList**](files.md) **>** [**H2ONaCl**](dir_0a7d68cdfa0215953309325142e4c674.md) **>** [**Matlab**](dir_a83d64805f95e67d3b37610e6eb407bf.md) **>** [**lookup.cpp**](lookup_8cpp.md)

[Go to the source code of this file](lookup_8cpp_source.md)



* `#include "mex.h"`
* `#include <matrix.h>`
* `#include <math.h>`
* `#include <string.h>`
* `#include <iostream>`
* `#include <vector>`
* `#include "H2ONaCl.h"`
* `#include "IAPS84.h"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  [**xThermal::PhaseRegion**](namespacexThermal.md#enum-phaseregion) | [**lookup\_TorHPX**](#function-lookup_torhpx) ([**xThermal::cxThermal**](classxThermal_1_1cxThermal.md) \* m\_pEOS, [**string**](classxThermal_1_1xThermalError.md) m\_fileName\_LUT, [**const**](classxThermal_1_1xThermalError.md) [**int**](classxThermal_1_1xThermalError.md) & dim, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & TorH, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X, [**double**](classxThermal_1_1xThermalError.md) \* props, [**bool**](classxThermal_1_1xThermalError.md) isCal, [**bool**](classxThermal_1_1xThermalError.md) isCout=[**false**](classxThermal_1_1xThermalError.md)) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**mGetMatrix**](#function-mgetmatrix) ([**const**](classxThermal_1_1xThermalError.md) [**mxArray**](classxThermal_1_1xThermalError.md) \* prhs, [**double**](classxThermal_1_1xThermalError.md) \*\* out, [**const**](classxThermal_1_1xThermalError.md) [**char**](classxThermal_1_1xThermalError.md) \* varname, [**mwSize**](classxThermal_1_1xThermalError.md) \* out\_m, [**mwSize**](classxThermal_1_1xThermalError.md) \* out\_n) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**mexFunction**](#function-mexfunction) ([**int**](classxThermal_1_1xThermalError.md) nlhs, [**mxArray**](classxThermal_1_1xThermalError.md) \* plhs, [**int**](classxThermal_1_1xThermalError.md) nrhs, [**const**](classxThermal_1_1xThermalError.md) [**mxArray**](classxThermal_1_1xThermalError.md) \* prhs) <br> |



























## Macros

| Type | Name |
| ---: | :--- |
| define  | [**Kelvin**](lookup_8cpp.md#define-kelvin)  `273.15`<br> |

## Public Functions Documentation




### function lookup\_TorHPX 

```C++
xThermal::PhaseRegion lookup_TorHPX (
    xThermal::cxThermal * m_pEOS,
    string m_fileName_LUT,
    const  int & dim,
    const  double & TorH,
    const  double & P,
    const  double & X,
    double * props,
    bool isCal,
    bool isCout=false
) 
```




<hr>



### function mGetMatrix 

```C++
void mGetMatrix (
    const  mxArray * prhs,
    double ** out,
    const  char * varname,
    mwSize * out_m,
    mwSize * out_n
) 
```




<hr>



### function mexFunction 

```C++
void mexFunction (
    int nlhs,
    mxArray * plhs,
    int nrhs,
    const  mxArray * prhs
) 
```




<hr>
## Macro Definition Documentation





### define Kelvin 

```C++
#define Kelvin `273.15`
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2ONaCl/Matlab/lookup.cpp`

