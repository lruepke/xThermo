

# Struct xThermal::IAPWS\_IF97::GibbsEnergy\_dimensionless



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**IAPWS\_IF97**](namespacexThermal_1_1IAPWS__IF97.md) **>** [**GibbsEnergy\_dimensionless**](structxThermal_1_1IAPWS__IF97_1_1GibbsEnergy__dimensionless.md)



_Dimensionless Gibbs free energy_ \(\gamma = g/(RT)\) _and its partial derivatives._

* `#include <IAPWS-IF97.h>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**double**](classxThermal_1_1xThermalError.md) | [**p**](#variable-p)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**pp**](#variable-pp)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**pt**](#variable-pt)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**t**](#variable-t)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**tt**](#variable-tt)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**value**](#variable-value)  <br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  [**GibbsEnergy\_dimensionless**](structxThermal_1_1IAPWS__IF97_1_1GibbsEnergy__dimensionless.md) | [**operator+**](#function-operator) ([**const**](classxThermal_1_1xThermalError.md) [**GibbsEnergy\_dimensionless**](structxThermal_1_1IAPWS__IF97_1_1GibbsEnergy__dimensionless.md) & gamma) const<br> |
|  [**GibbsEnergy\_dimensionless**](structxThermal_1_1IAPWS__IF97_1_1GibbsEnergy__dimensionless.md) | [**operator-**](#function-operator-) ([**const**](classxThermal_1_1xThermalError.md) [**GibbsEnergy\_dimensionless**](structxThermal_1_1IAPWS__IF97_1_1GibbsEnergy__dimensionless.md) & gamma) const<br> |
|  [**bool**](classxThermal_1_1xThermalError.md) | [**operator==**](#function-operator_1) ([**const**](classxThermal_1_1xThermalError.md) [**GibbsEnergy\_dimensionless**](structxThermal_1_1IAPWS__IF97_1_1GibbsEnergy__dimensionless.md) & gamma) const<br> |




























## Public Attributes Documentation




### variable p 

```C++
double xThermal::IAPWS_IF97::GibbsEnergy_dimensionless::p;
```



\(\left[ \frac{\partial \gamma}{\partial \pi} \right]_{\tau}\) 


        

<hr>



### variable pp 

```C++
double xThermal::IAPWS_IF97::GibbsEnergy_dimensionless::pp;
```



\(\left[ \frac{\partial^2 \gamma}{\partial \pi^2} \right]_{\tau}\) 


        

<hr>



### variable pt 

```C++
double xThermal::IAPWS_IF97::GibbsEnergy_dimensionless::pt;
```



\(\frac{\partial^2 \gamma}{\partial \pi \partial \tau}\) 


        

<hr>



### variable t 

```C++
double xThermal::IAPWS_IF97::GibbsEnergy_dimensionless::t;
```



\(\left[ \frac{\partial \gamma}{\partial \tau} \right]_{\pi}\) 


        

<hr>



### variable tt 

```C++
double xThermal::IAPWS_IF97::GibbsEnergy_dimensionless::tt;
```



\(\left[ \frac{\partial^2 \gamma}{\partial \tau^2} \right]_{\pi}\) 


        

<hr>



### variable value 

```C++
double xThermal::IAPWS_IF97::GibbsEnergy_dimensionless::value;
```



\(\gamma\) 


        

<hr>
## Public Functions Documentation




### function operator+ 

```C++
inline GibbsEnergy_dimensionless xThermal::IAPWS_IF97::GibbsEnergy_dimensionless::operator+ (
    const  GibbsEnergy_dimensionless & gamma
) const
```




<hr>



### function operator- 

```C++
inline GibbsEnergy_dimensionless xThermal::IAPWS_IF97::GibbsEnergy_dimensionless::operator- (
    const  GibbsEnergy_dimensionless & gamma
) const
```




<hr>



### function operator== 

```C++
inline bool xThermal::IAPWS_IF97::GibbsEnergy_dimensionless::operator== (
    const  GibbsEnergy_dimensionless & gamma
) const
```




<hr>## Friends Documentation





### friend operator&lt;&lt; 

```C++
inline std::ostream & xThermal::IAPWS_IF97::GibbsEnergy_dimensionless::operator<< (
    std::ostream & os,
    const  GibbsEnergy_dimensionless & gamma
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2O/IAPWS-IF97.h`

