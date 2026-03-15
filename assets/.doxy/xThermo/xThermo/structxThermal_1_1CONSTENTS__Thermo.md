

# Struct xThermal::CONSTENTS\_Thermo



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**CONSTENTS\_Thermo**](structxThermal_1_1CONSTENTS__Thermo.md)



_为了实现多个H2O EOS的backend，必须使用虚函数调用相应的参数，比如T\_critical()，但是在子类中频繁调用函数的性能肯定很低，所以将所有热力学常数打包为一个struct类型，作为子类的成员变量，然后在构造函数中进行初始化，后面需要常数的地方直接访问成员变量即可，可提高性能。_ 

* `#include <DataStructures.h>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**double**](classxThermal_1_1xThermalError.md) | [**R**](#variable-r)   = `461.51805`<br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**T\_critical**](#variable-t_critical)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Tmax**](#variable-tmax)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Tmin**](#variable-tmin)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**Ttriple**](#variable-ttriple)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**molar\_mass**](#variable-molar_mass)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**p\_critical**](#variable-p_critical)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**pmax**](#variable-pmax)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**pmin**](#variable-pmin)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**rhomass\_critical**](#variable-rhomass_critical)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**rhomolar\_critical**](#variable-rhomolar_critical)  <br> |












































## Public Attributes Documentation




### variable R 

```C++
double xThermal::CONSTENTS_Thermo::R;
```




<hr>



### variable T\_critical 

```C++
double xThermal::CONSTENTS_Thermo::T_critical;
```




<hr>



### variable Tmax 

```C++
double xThermal::CONSTENTS_Thermo::Tmax;
```




<hr>



### variable Tmin 

```C++
double xThermal::CONSTENTS_Thermo::Tmin;
```




<hr>



### variable Ttriple 

```C++
double xThermal::CONSTENTS_Thermo::Ttriple;
```




<hr>



### variable molar\_mass 

```C++
double xThermal::CONSTENTS_Thermo::molar_mass;
```




<hr>



### variable p\_critical 

```C++
double xThermal::CONSTENTS_Thermo::p_critical;
```




<hr>



### variable pmax 

```C++
double xThermal::CONSTENTS_Thermo::pmax;
```




<hr>



### variable pmin 

```C++
double xThermal::CONSTENTS_Thermo::pmin;
```




<hr>



### variable rhomass\_critical 

```C++
double xThermal::CONSTENTS_Thermo::rhomass_critical;
```




<hr>



### variable rhomolar\_critical 

```C++
double xThermal::CONSTENTS_Thermo::rhomolar_critical;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/thermo/DataStructures.h`

