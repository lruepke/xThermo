

# File IAPWS95.h



[**FileList**](files.md) **>** [**H2O**](dir_22b94ec4f3699d3ab17d67318090f0eb.md) **>** [**IAPWS95.h**](IAPWS95_8h.md)

[Go to the source code of this file](IAPWS95_8h_source.md)

_Implementation of IAPWS-95 EOS._ [More...](#detailed-description)

* `#include "thermo.h"`
* `#include "IAPWS-IF97.h"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**xThermal**](namespacexThermal.md) <br>_Namespace of_ [_**xThermal**_](namespacexThermal.md) _library._ |
| namespace | [**IAPWS95**](namespacexThermal_1_1IAPWS95.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| struct | [**Coeff\_phi\_o**](structxThermal_1_1IAPWS95_1_1Coeff__phi__o.md) <br> |
| struct | [**Coeff\_phi\_r**](structxThermal_1_1IAPWS95_1_1Coeff__phi__r.md) <br>_Numerical values of the coefficients and parameters of the residual part of the dimensionless Helmholtz free energy. Eq.6 in IAPWS-95_ [_**cIAPWS95::phi\_r**_](classxThermal_1_1IAPWS95_1_1cIAPWS95.md#function-phi_r-13) _._ |
| struct | [**Constants\_Viscosity2008\_Water**](structxThermal_1_1IAPWS95_1_1Constants__Viscosity2008__Water.md) <br>_Constants for water viscosity calculation, see mu2008._  |
| struct | [**HelmholtzEnergy\_dimensionless**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless.md) <br>_Dimensionless Helmholtz free energy_ \(\phi = f/(RT)\) _and its partial derivatives._ |
| struct | [**HelmholtzEnergy\_dimensionless\_SinglePhase**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless__SinglePhase.md) <br> |
| struct | [**HelmholtzEnergy\_dimensionless\_TwoPhase**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless__TwoPhase.md) <br> |
| struct | [**Params\_HP2RhoT**](structxThermal_1_1IAPWS95_1_1Params__HP2RhoT.md) <br> |
| struct | [**Params\_SolvePhaseEquilibrium**](structxThermal_1_1IAPWS95_1_1Params__SolvePhaseEquilibrium.md) <br>_Data struct for phase equilibrium solving. It is used in_ [_**func\_PhaseEquilibrium**_](namespacexThermal_1_1IAPWS95.md#function-func_phaseequilibrium) _,_cIAPWS95::Boiling\_P _and_cIAPWS95::Boiling\_T _._ |
| struct | [**Params\_TP2Rho**](structxThermal_1_1IAPWS95_1_1Params__TP2Rho.md) <br> |
| struct | [**Params\_T\_Sat\_estimate**](structxThermal_1_1IAPWS95_1_1Params__T__Sat__estimate.md) <br> |
| struct | [**STRUCT\_delta\_TwoPhase**](structxThermal_1_1IAPWS95_1_1STRUCT__delta__TwoPhase.md) <br> |
| struct | [**State**](structxThermal_1_1IAPWS95_1_1State.md) <br> |
| class | [**cIAPWS95**](classxThermal_1_1IAPWS95_1_1cIAPWS95.md) <br>_Class of IAPWS-95 formula of H2O, which inherits from_ [_**IAPWS\_IF97::cIAPWS\_IF97**_](classxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97.md) _._ |

















































## Macros

| Type | Name |
| ---: | :--- |
| define  | [**Name\_Backend\_IAPWS95**](IAPWS95_8h.md#define-name_backend_iapws95)  `"IAPWS95"`<br> |

## Detailed Description




**Author:**

Zhikui Guo ([zguo@geomar.de](mailto:zguo@geomar.de)) 




**Version:**

0.1 




**Date:**

2022-04-11




**Copyright:**

Copyright (c) 2022 





    
## Macro Definition Documentation





### define Name\_Backend\_IAPWS95 

```C++
#define Name_Backend_IAPWS95 `"IAPWS95"`
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2O/IAPWS95.h`

