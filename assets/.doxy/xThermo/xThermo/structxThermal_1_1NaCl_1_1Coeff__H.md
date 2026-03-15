

# Struct xThermal::NaCl::Coeff\_H



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**NaCl**](namespacexThermal_1_1NaCl.md) **>** [**Coeff\_H**](structxThermal_1_1NaCl_1_1Coeff__H.md)



_Parameters for enthalpies. See table 5 of Driesner2007Part2._ [More...](#detailed-description)

* `#include <NaCl.h>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**double**](classxThermal_1_1xThermalError.md) | [**r**](#variable-r)   = `/* multi line expression */`<br> |
|  [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) | [**r3**](#variable-r3)   = `{-0.0017099[**E**](classxThermal_1_1xThermalError.md)-5, - 3.82734[**E**](classxThermal_1_1xThermalError.md)-11, - 8.65455[**E**](classxThermal_1_1xThermalError.md)-14}`<br> |












































## Detailed Description




**Note:**

r[3] is a function of T, need to updated inside function Cp.




**Note:**

Again, in order to make input P with unit [Pa], multiply 1E-5 and 1E-10 to \(r_3, r_4\), respectively, see equation 30 of Driesner2007Part2.




**Warning:**

Why the r4 is written as summation of some tiny values in the Tab. 5 ?




**Todo**

Need to confirm r4 in Tab.5 is correct




    
## Public Attributes Documentation




### variable r 

```C++
double xThermal::NaCl::Coeff_H::r[5];
```




<hr>



### variable r3 

```C++
const double xThermal::NaCl::Coeff_H::r3[3];
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/NaCl/NaCl.h`

