

# File IAPS84.cpp



[**FileList**](files.md) **>** [**H2O**](dir_22b94ec4f3699d3ab17d67318090f0eb.md) **>** [**IAPS84.cpp**](IAPS84_8cpp.md)

[Go to the source code of this file](IAPS84_8cpp_source.md)

_Implementation of the xThermo backend "IAPS84", a C++ wrapper around PROST 4.1 (see_ [_**IAPS84.h**_](IAPS84_8h.md) _)._[More...](#detailed-description)

* `#include "IAPS84.h"`
* `#include "steam4.h"`
* `#include "iaps.h"`
* `#include "LookUpTableForest.h"`
* `#include "LookUpTableForestI.H"`
* `#include "interpolationI.H"`
* `#include "AMR_LUT_RefineFuncI.H"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**xThermal**](namespacexThermal.md) <br>_Namespace of_ [_**xThermal**_](namespacexThermal.md) _library._ |
| namespace | [**PROST**](namespacexThermal_1_1PROST.md) <br> |




















































## Detailed Description




**Author:**

Zhikui Guo ([zguo@geomar.de](mailto:zguo@geomar.de)) 


All IAPS-84 computations (equation of state, derivatives, saturation, viscosity) are done by PROST 4.1 (O. Bauer and O. Engel, TU Hamburg-Harburg, 1998; ThirdParties/PROST; GNU LGPL v2.0). This file contains only the xThermo interface around it. 

**Version:**

0.1 




**Date:**

2022-03-27




**Copyright:**

Copyright (c) 2022. This wrapper is part of xThermo (MIT licence, see LICENSE). PROST is not covered by xThermo's licence (see THIRD\_PARTY\_NOTICES.md). 





    

------------------------------
The documentation for this class was generated from the following file `Library/H2O/IAPS84.cpp`

