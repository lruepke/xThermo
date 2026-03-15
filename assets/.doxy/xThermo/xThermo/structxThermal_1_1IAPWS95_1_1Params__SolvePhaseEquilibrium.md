

# Struct xThermal::IAPWS95::Params\_SolvePhaseEquilibrium



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**IAPWS95**](namespacexThermal_1_1IAPWS95.md) **>** [**Params\_SolvePhaseEquilibrium**](structxThermal_1_1IAPWS95_1_1Params__SolvePhaseEquilibrium.md)



_Data struct for phase equilibrium solving. It is used in_ [_**func\_PhaseEquilibrium**_](namespacexThermal_1_1IAPWS95.md#function-func_phaseequilibrium) _,_cIAPWS95::Boiling\_P _and_cIAPWS95::Boiling\_T _._

* `#include <IAPWS95.h>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**double**](classxThermal_1_1xThermalError.md) | [**P**](#variable-p)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**RT**](#variable-rt)  <br> |
|  [**SOLVE\_SATURATED\_PorT**](namespacexThermal.md#enum-solve_saturated_port) | [**Solve\_PorT**](#variable-solve_port)  <br> |
|  [**union**](classxThermal_1_1xThermalError.md) [**xThermal::IAPWS95::Params\_SolvePhaseEquilibrium**](structxThermal_1_1IAPWS95_1_1Params__SolvePhaseEquilibrium.md) | [**TorP**](#variable-torp)  <br>_One of_ \(\tau\) _and_\(p\) _, which is used for solving saturated pressure and saturated temperature, respectively._ |
|  [**cIAPWS95**](classxThermal_1_1IAPWS95_1_1cIAPWS95.md) \* | [**iapws**](#variable-iapws)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**tau**](#variable-tau)  <br> |












































## Public Attributes Documentation




### variable P 

```C++
double xThermal::IAPWS95::Params_SolvePhaseEquilibrium::P;
```



Pressure [Pa] 


        

<hr>



### variable RT 

```C++
double xThermal::IAPWS95::Params_SolvePhaseEquilibrium::RT;
```



[**xThermal::R**](namespacexThermal.md#variable-r) times T: \(RT\) [J/kg] 


        

<hr>



### variable Solve\_PorT 

```C++
SOLVE_SATURATED_PorT xThermal::IAPWS95::Params_SolvePhaseEquilibrium::Solve_PorT;
```



Flag of solving \(P_{sat}\) or \(T_{sat}\), the value will be one of SOLVE\_SATURATED\_P, SOLVE\_SATURATED\_T 


        

<hr>



### variable TorP 

_One of_ \(\tau\) _and_\(p\) _, which is used for solving saturated pressure and saturated temperature, respectively._
```C++
union xThermal::IAPWS95::Params_SolvePhaseEquilibrium xThermal::IAPWS95::Params_SolvePhaseEquilibrium::TorP;
```





**Note:**

Because for solving phase equilibrium, only P or T is given, so here I use an union to store it. 





        

<hr>



### variable iapws 

```C++
cIAPWS95* xThermal::IAPWS95::Params_SolvePhaseEquilibrium::iapws;
```




<hr>



### variable tau 

```C++
double xThermal::IAPWS95::Params_SolvePhaseEquilibrium::tau;
```



\(T_c/T\) 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2O/IAPWS95.h`

