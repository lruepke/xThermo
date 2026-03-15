

# File H2ONaCl.h



[**FileList**](files.md) **>** [**H2ONaCl**](dir_0a7d68cdfa0215953309325142e4c674.md) **>** [**H2ONaCl.h**](H2ONaCl_8h.md)

[Go to the source code of this file](H2ONaCl_8h_source.md)

_Head file of H2ONaCl._ [More...](#detailed-description)

* `#include "IAPS84.h"`
* `#include "IAPWS95_CoolProp.h"`
* `#include "IAPWS95.h"`
* `#include "NaCl.h"`
* `#include <vector>`
* `#include <cstdlib>`
* `#include <cstdio>`
* `#include <gsl/gsl_vector.h>`
* `#include <gsl/gsl_multiroots.h>`
* `#include <gsl/gsl_errno.h>`
* `#include <gsl/gsl_math.h>`
* `#include <gsl/gsl_roots.h>`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**xThermal**](namespacexThermal.md) <br>_Namespace of_ [_**xThermal**_](namespacexThermal.md) _library._ |
| namespace | [**H2ONaCl**](namespacexThermal_1_1H2ONaCl.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| struct | [**Coeffs\_Estimate\_CriticalT**](structxThermal_1_1H2ONaCl_1_1Coeffs__Estimate__CriticalT.md) <br>_Poly fit coefficient for the estimate critical temperature._  |
| struct | [**Coeffs\_Viscosity**](structxThermal_1_1H2ONaCl_1_1Coeffs__Viscosity.md) <br> |
| struct | [**Params\_Inversion\_PTX**](structxThermal_1_1H2ONaCl_1_1Params__Inversion__PTX.md) <br>_Parameters for pressure and temperature calculation in VLH zone._  |
| struct | [**Params\_P2CriticalT**](structxThermal_1_1H2ONaCl_1_1Params__P2CriticalT.md) <br>_Parameters for inversion of critical temperature from pressure._  |
| struct | [**Params\_PX2T\_VL**](structxThermal_1_1H2ONaCl_1_1Params__PX2T__VL.md) <br> |
| struct | [**Table4\_Driesner2007a\_CriticalCurve**](structxThermal_1_1H2ONaCl_1_1Table4__Driesner2007a__CriticalCurve.md) <br> |
| struct | [**Table6\_Pressure\_VLH**](structxThermal_1_1H2ONaCl_1_1Table6__Pressure__VLH.md) <br> |
| struct | [**Table7\_XL\_VL**](structxThermal_1_1H2ONaCl_1_1Table7__XL__VL.md) <br>_Parameters for liquid composition,_ \(X_{NaCl}^{VL,liq}\) _, on V+L coexistence surface. V+L liquid branch. See Table 7 of Driesner2007Part1._ |
| struct | [**Table8\_VaporComposition**](structxThermal_1_1H2ONaCl_1_1Table8__VaporComposition.md) <br> |
| class | [**cH2ONaCl**](classxThermal_1_1H2ONaCl_1_1cH2ONaCl.md) <br>_Class of_ \(H_2O-NaCl\) _EOS._ |


















































## Detailed Description




**Author:**

Zhikui Guo ([zguo@geomar.de](mailto:zguo@geomar.de)) 




**Version:**

0.1 




**Date:**

2022-03-27




**Copyright:**

Copyright (c) 2022 





    

------------------------------
The documentation for this class was generated from the following file `Library/H2ONaCl/H2ONaCl.h`

