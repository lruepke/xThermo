

# Group CONSTANTS\_NACL



[**Modules**](modules.md) **>** [**CONSTANTS\_NACL**](group__CONSTANTS__NACL.md)



[More...](#detailed-description)
















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**xThermal::NaCl::H\_halite\_ref**](structxThermal_1_1NaCl_1_1H__halite__ref.md) <br> |






## Public Attributes

| Type | Name |
| ---: | :--- |
|  double const | [**M**](#variable-m)   = `0.058443`<br> |
|  double const | [**P\_Triple**](#variable-p_triple)   = `50`<br> |
|  double const | [**T\_Triple**](#variable-t_triple)   = `T\_Triple\_C + 273.15`<br> |
|  double const | [**T\_Triple\_C**](#variable-t_triple_c)   = `800.7`<br> |
|  double const | [**T\_Triple\_sqr**](#variable-t_triple_sqr)   = `T\_Triple\*T\_Triple`<br> |
|  double const | [**T\_Triple\_sqr\_C**](#variable-t_triple_sqr_c)   = `T\_Triple\_C\*T\_Triple\_C`<br> |
|  double const | [**inv\_T\_Triple**](#variable-inv_t_triple)   = `1.0/T\_Triple`<br> |
|  double const | [**log10\_P\_Triple**](#variable-log10_p_triple)   = `log10(P\_Triple)`<br> |
|  double const | [**slope\_Clapeyron**](#variable-slope_clapeyron)   = `1.0/2.4726e-7`<br> |












































## Detailed Description


see last pargraph of section 3.2 in Driesner2007Part1. 


    
## Public Attributes Documentation




### variable M 

```
double const xThermal::NaCl::M;
```



molar mass [kg/mol] 


        

<hr>



### variable P\_Triple 

```
double const xThermal::NaCl::P_Triple;
```



Pressure at triple point of NaCl, [Pa] 


        

<hr>



### variable T\_Triple 

```
double const xThermal::NaCl::T_Triple;
```



Temperature at triple point of NaCl, [K] 


        

<hr>



### variable T\_Triple\_C 

```
double const xThermal::NaCl::T_Triple_C;
```



Temperature at triple point of NaCl [ \(^{\circ} C\)] 


        

<hr>



### variable T\_Triple\_sqr 

```
double const xThermal::NaCl::T_Triple_sqr;
```




<hr>



### variable T\_Triple\_sqr\_C 

```
double const xThermal::NaCl::T_Triple_sqr_C;
```




<hr>



### variable inv\_T\_Triple 

```
double const xThermal::NaCl::inv_T_Triple;
```



\(\frac{1}{T_{Triple}}\) 


        

<hr>



### variable log10\_P\_Triple 

```
double const xThermal::NaCl::log10_P_Triple;
```



log\_{10}^{P\_{Triple}} 


        

<hr>



### variable slope\_Clapeyron 

```
double const xThermal::NaCl::slope_Clapeyron;
```



Clapeyron slope \(\frac{dP}{dT} = \frac{1}{a} ~ Pa/K\) of the NaCl melting curve (see Tab. 3 of Driesner2007Part1), which is used in \(\Delta H_{fus}\) calculation (Eq. 28 of Driesner2007Part2) 

**Note:**

Notice the unit! 





        

<hr>

------------------------------


