

# Struct xThermal::NaCl::H\_halite\_ref



[**ClassList**](annotated.md) **>** [**NaCl**](namespacexThermal_1_1NaCl.md) **>** [**H\_halite\_ref**](structxThermal_1_1NaCl_1_1H__halite__ref.md)



[More...](#detailed-description)

* `#include <NaCl.h>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) | [**H**](#variable-h)   = `9.415867359[**e4**](classxThermal_1_1xThermalError.md)`<br> |
|  [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) | [**P**](#variable-p)   = `100E5`<br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**R0**](#variable-r0)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**R1**](#variable-r1)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**R2**](#variable-r2)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**R3**](#variable-r3)  <br> |
|  [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) | [**T\_C**](#variable-t_c)   = `100`<br> |












































## Detailed Description


Reference point and enthalpy of halite. The JANAF table can be found at [https://janaf.nist.gov](https://janaf.nist.gov). While this reference point data comes from Falko Vehling(personal communication: he said this table comes from an attached table to a paper) 


    
## Public Attributes Documentation




### variable H 

```C++
const double xThermal::NaCl::H_halite_ref::H;
```



Enthalpy of halite at reference point (T0, P0). 

**Todo**

need to confirm the reference point of halite enthalpy. 




        

<hr>



### variable P 

```C++
const double xThermal::NaCl::H_halite_ref::P;
```



Pressure [Pa] at reference point 


        

<hr>



### variable R0 

```C++
double xThermal::NaCl::H_halite_ref::R0;
```



H1 = R0 + R1\*T + R2\*T^2 + R3\*T^3 


        

<hr>



### variable R1 

```C++
double xThermal::NaCl::H_halite_ref::R1;
```



Coefficient R1 as function of P in equation of H1: R1 = R1[0] + R1[1]\*P0 + R1[2]\*P0^2 


        

<hr>



### variable R2 

```C++
double xThermal::NaCl::H_halite_ref::R2;
```



Coefficient R2 as function of P in equation of H1: R2 = R2[0] + R2[1]\*P0 


        

<hr>



### variable R3 

```C++
double xThermal::NaCl::H_halite_ref::R3;
```



Coefficient R2 as function of P in equation of H1: R2 = R1[0] + R2[1]\*P0 


        

<hr>



### variable T\_C 

```C++
const double xThermal::NaCl::H_halite_ref::T_C;
```



Temperature [deg.C] at reference point 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `Library/NaCl/NaCl.h`

