

# Struct xThermal::IAPWS95::HelmholtzEnergy\_dimensionless



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**IAPWS95**](namespacexThermal_1_1IAPWS95.md) **>** [**HelmholtzEnergy\_dimensionless**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless.md)



_Dimensionless Helmholtz free energy_ \(\phi = f/(RT)\) _and its partial derivatives._

* `#include <IAPWS95.h>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**double**](classxThermal_1_1xThermalError.md) | [**d**](#variable-d)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**dd**](#variable-dd)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**dt**](#variable-dt)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**t**](#variable-t)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**tt**](#variable-tt)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**value**](#variable-value)  <br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  [**HelmholtzEnergy\_dimensionless**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless.md) | [**operator+**](#function-operator) ([**const**](classxThermal_1_1xThermalError.md) [**HelmholtzEnergy\_dimensionless**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless.md) & phi) const<br> |
|  [**HelmholtzEnergy\_dimensionless**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless.md) | [**operator-**](#function-operator-) ([**const**](classxThermal_1_1xThermalError.md) [**HelmholtzEnergy\_dimensionless**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless.md) & phi) const<br> |
|  [**bool**](classxThermal_1_1xThermalError.md) | [**operator==**](#function-operator_1) ([**const**](classxThermal_1_1xThermalError.md) [**HelmholtzEnergy\_dimensionless**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless.md) & phi) const<br> |




























## Public Attributes Documentation




### variable d 

```C++
double xThermal::IAPWS95::HelmholtzEnergy_dimensionless::d;
```



\(\left[ \frac{\partial \phi}{\partial \delta} \right]_{\tau}\) 


        

<hr>



### variable dd 

```C++
double xThermal::IAPWS95::HelmholtzEnergy_dimensionless::dd;
```



\(\left[ \frac{\partial^2 \phi}{\partial \delta^2} \right]_{\tau}\) 


        

<hr>



### variable dt 

```C++
double xThermal::IAPWS95::HelmholtzEnergy_dimensionless::dt;
```



\(\frac{\partial^2 \phi}{\partial \delta \partial \tau}\) 


        

<hr>



### variable t 

```C++
double xThermal::IAPWS95::HelmholtzEnergy_dimensionless::t;
```



\(\left[ \frac{\partial \phi}{\partial \tau} \right]_{\delta}\) 


        

<hr>



### variable tt 

```C++
double xThermal::IAPWS95::HelmholtzEnergy_dimensionless::tt;
```



\(\left[ \frac{\partial^2 \phi}{\partial \tau^2} \right]_{\delta}\) 


        

<hr>



### variable value 

```C++
double xThermal::IAPWS95::HelmholtzEnergy_dimensionless::value;
```



\(\phi\) 


        

<hr>
## Public Functions Documentation




### function operator+ 

```C++
inline HelmholtzEnergy_dimensionless xThermal::IAPWS95::HelmholtzEnergy_dimensionless::operator+ (
    const  HelmholtzEnergy_dimensionless & phi
) const
```




<hr>



### function operator- 

```C++
inline HelmholtzEnergy_dimensionless xThermal::IAPWS95::HelmholtzEnergy_dimensionless::operator- (
    const  HelmholtzEnergy_dimensionless & phi
) const
```




<hr>



### function operator== 

```C++
inline bool xThermal::IAPWS95::HelmholtzEnergy_dimensionless::operator== (
    const  HelmholtzEnergy_dimensionless & phi
) const
```




<hr>## Friends Documentation





### friend operator&lt;&lt; 

```C++
inline std::ostream & xThermal::IAPWS95::HelmholtzEnergy_dimensionless::operator<< (
    std::ostream & os,
    const  HelmholtzEnergy_dimensionless & phi
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2O/IAPWS95.h`

