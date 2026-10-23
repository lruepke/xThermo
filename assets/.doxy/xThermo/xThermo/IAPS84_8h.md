

# File IAPS84.h



[**FileList**](files.md) **>** [**H2O**](dir_22b94ec4f3699d3ab17d67318090f0eb.md) **>** [**IAPS84.h**](IAPS84_8h.md)

[Go to the source code of this file](IAPS84_8h_source.md)

_xThermo backend "IAPS84": a C++ wrapper around PROST 4.1, which implements the IAPS-84 equation of state for water (Haar, Gallagher & Kell 1984; prost)._ [More...](#detailed-description)

* `#include "thermo.h"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**xThermal**](namespacexThermal.md) <br>_Namespace of_ [_**xThermal**_](namespacexThermal.md) _library._ |
| namespace | [**PROST**](namespacexThermal_1_1PROST.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| class | [**cIAPS84**](classxThermal_1_1PROST_1_1cIAPS84.md) <br> |

















































## Macros

| Type | Name |
| ---: | :--- |
| define  | [**Name\_Backend\_IAPS84**](IAPS84_8h.md#define-name_backend_iaps84)  `"IAPS84"`<br> |

## Detailed Description




**Author:**

Zhikui Guo ([zguo@geomar.de](mailto:zguo@geomar.de)) 


The equation of state, its derivatives, the saturation routines and the viscosity correlation are computed entirely by PROST 4.1 (O. Bauer and O. Engel, TU Hamburg-Harburg, 1998; ThirdParties/PROST), which is distributed under the GNU Library General Public License v2.0 (ThirdParties/PROST/COPYING.LIB-2.0). This wrapper contains only the xThermo interface: state set-up, unit and derivative conversion, phase classification and result mapping. The changes xThermo made to PROST itself are listed in ThirdParties/PROST/XTHERMO\_MODIFICATIONS.txt. 

**Version:**

0.1 




**Date:**

2022-03-27




**Copyright:**

Copyright (c) 2022. This wrapper is part of xThermo (MIT licence, see LICENSE). PROST is not covered by xThermo's licence (see THIRD\_PARTY\_NOTICES.md). 





    
## Macro Definition Documentation





### define Name\_Backend\_IAPS84 

```C++
#define Name_Backend_IAPS84 `"IAPS84"`
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2O/IAPS84.h`

