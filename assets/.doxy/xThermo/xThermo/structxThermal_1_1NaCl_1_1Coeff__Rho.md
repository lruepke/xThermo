

# Struct xThermal::NaCl::Coeff\_Rho



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**NaCl**](namespacexThermal_1_1NaCl.md) **>** [**Coeff\_Rho**](structxThermal_1_1NaCl_1_1Coeff__Rho.md)



_Parameters for halite and liquid NaCl volumetric properties. See table 3 of Driesner2007Part2._ [More...](#detailed-description)

* `#include <NaCl.h>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) | [**l**](#variable-l)   = `{2170.4, -0.24599, -9.5797[**E**](classxThermal_1_1xThermalError.md)-05, 0.005727[**E**](classxThermal_1_1xThermalError.md)-5, 0.002715[**E**](classxThermal_1_1xThermalError.md)-5, 733.4}`<br> |
|  [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) | [**m**](#variable-m)   = `{58443, 23.772, 0.018639, -1.9687[**E**](classxThermal_1_1xThermalError.md)-06, -1.5259[**E**](classxThermal_1_1xThermalError.md)-10, 5.5058[**E**](classxThermal_1_1xThermalError.md)-13}`<br> |












































## Detailed Description




**Note:**

In order to make the input pressure with unit [Pa], change \(l_3, l_4\) to be 0.005727E-5 and 0.002715E-5, respectively. Because the pressure unit in equation (1-3) is [bar]. In addition, change \(m_4, m_5\) to be -1.5259E-10, 5.5058E-13 for the same reason. (See equation 5,6 of Driesner2007Part2) 





    
## Public Attributes Documentation




### variable l 

```C++
const double xThermal::NaCl::Coeff_Rho::l[6];
```




<hr>



### variable m 

```C++
const double xThermal::NaCl::Coeff_Rho::m[6];
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/NaCl/NaCl.h`

