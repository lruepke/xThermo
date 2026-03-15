

# Struct xThermal::H2ONaCl::Coeffs\_Estimate\_CriticalT



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**H2ONaCl**](namespacexThermal_1_1H2ONaCl.md) **>** [**Coeffs\_Estimate\_CriticalT**](structxThermal_1_1H2ONaCl_1_1Coeffs__Estimate__CriticalT.md)



_Poly fit coefficient for the estimate critical temperature._ [More...](#detailed-description)

* `#include <H2ONaCl.h>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) | [**coeffs**](#variable-coeffs)   = `/* multi line expression */`<br> |
|  [**const**](classxThermal_1_1xThermalError.md) [**int**](classxThermal_1_1xThermalError.md) | [**num**](#variable-num)   = `9`<br> |












































## Detailed Description


\(T_c = \sum_{i=0}^{9} a_i y^i\), where \(y=log_{10}P[Pa]\) 


    
## Public Attributes Documentation




### variable coeffs 

```C++
const double xThermal::H2ONaCl::Coeffs_Estimate_CriticalT::coeffs[9];
```




<hr>



### variable num 

```C++
const int xThermal::H2ONaCl::Coeffs_Estimate_CriticalT::num;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2ONaCl/H2ONaCl.h`

