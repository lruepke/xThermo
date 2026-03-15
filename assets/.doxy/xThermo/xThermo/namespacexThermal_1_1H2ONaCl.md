

# Namespace xThermal::H2ONaCl



[**Namespace List**](namespaces.md) **>** [**xThermal**](namespacexThermal.md) **>** [**H2ONaCl**](namespacexThermal_1_1H2ONaCl.md)




















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**Coeffs\_Estimate\_CriticalT**](structxThermal_1_1H2ONaCl_1_1Coeffs__Estimate__CriticalT.md) <br>_Poly fit coefficient for the estimate critical temperature._  |
| struct | [**Coeffs\_Viscosity**](structxThermal_1_1H2ONaCl_1_1Coeffs__Viscosity.md) <br> |
| struct | [**Params\_Inversion\_PTX**](structxThermal_1_1H2ONaCl_1_1Params__Inversion__PTX.md) <br>_Parameters for pressure and temperature calculation in VLH zone._  |
| struct | [**Params\_P2CriticalT**](structxThermal_1_1H2ONaCl_1_1Params__P2CriticalT.md) <br>_Parameters for inversion of critical temperature from pressure._  |
| struct | [**Params\_PX2T\_VL**](structxThermal_1_1H2ONaCl_1_1Params__PX2T__VL.md) <br> |
| struct | [**Table4\_Driesner2007a\_CriticalCurve**](structxThermal_1_1H2ONaCl_1_1Table4__Driesner2007a__CriticalCurve.md) <br> |
| struct | [**Table6\_Pressure\_VLH**](structxThermal_1_1H2ONaCl_1_1Table6__Pressure__VLH.md) <br> |
| struct | [**Table7\_XL\_VL**](structxThermal_1_1H2ONaCl_1_1Table7__XL__VL.md) <br>_Parameters for liquid composition,_ \(X_{NaCl}^{VL,liq}\) _, on V+L coexistence surface. V+L liquid branch. See Table 7 of Driesner2007Part1._ |
| struct | [**Table8\_VaporComposition**](structxThermal_1_1H2ONaCl_1_1Table8__VaporComposition.md) <br> |
| class | [**cH2ONaCl**](classxThermal_1_1H2ONaCl_1_1cH2ONaCl.md) <br>_Class of_ \(H_2O-NaCl\) _EOS._ |






















## Public Functions

| Type | Name |
| ---: | :--- |
|  [**double**](classxThermal_1_1xThermalError.md) | [**func\_P2CriticalT**](#function-func_p2criticalt) ([**double**](classxThermal_1_1xThermalError.md) T, [**void**](classxThermal_1_1xThermalError.md) \* params) <br>_Construct function of_ \(P_{critical} = f(T)\) _, see Eq. 5a-c of Driesner2007Part1._ |
|  [**double**](classxThermal_1_1xThermalError.md) | [**func\_TL\_VL**](#function-func_tl_vl) ([**double**](classxThermal_1_1xThermalError.md) T, [**void**](classxThermal_1_1xThermalError.md) \* params) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**func\_T\_HPX**](#function-func_t_hpx) ([**double**](classxThermal_1_1xThermalError.md) T, [**void**](classxThermal_1_1xThermalError.md) \* params) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**func\_T\_VLH**](#function-func_t_vlh) ([**double**](classxThermal_1_1xThermalError.md) T, [**void**](classxThermal_1_1xThermalError.md) \* params) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**func\_T\_VL\_L**](#function-func_t_vl_l) ([**double**](classxThermal_1_1xThermalError.md) T, [**void**](classxThermal_1_1xThermalError.md) \* params) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**func\_T\_VL\_V**](#function-func_t_vl_v) ([**double**](classxThermal_1_1xThermalError.md) T, [**void**](classxThermal_1_1xThermalError.md) \* params) <br> |




























## Public Functions Documentation




### function func\_P2CriticalT 

_Construct function of_ \(P_{critical} = f(T)\) _, see Eq. 5a-c of Driesner2007Part1._
```C++
double xThermal::H2ONaCl::func_P2CriticalT (
    double T,
    void * params
) 
```





**Parameters:**


* `x` Temperature [K] 
* `params` 
* `f` Pressure [Pa] 



**Returns:**

int 





        

<hr>



### function func\_TL\_VL 

```C++
double xThermal::H2ONaCl::func_TL_VL (
    double T,
    void * params
) 
```




<hr>



### function func\_T\_HPX 

```C++
double xThermal::H2ONaCl::func_T_HPX (
    double T,
    void * params
) 
```




<hr>



### function func\_T\_VLH 

```C++
double xThermal::H2ONaCl::func_T_VLH (
    double T,
    void * params
) 
```




<hr>



### function func\_T\_VL\_L 

```C++
double xThermal::H2ONaCl::func_T_VL_L (
    double T,
    void * params
) 
```




<hr>



### function func\_T\_VL\_V 

```C++
double xThermal::H2ONaCl::func_T_VL_V (
    double T,
    void * params
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2ONaCl/H2ONaCl.cpp`

