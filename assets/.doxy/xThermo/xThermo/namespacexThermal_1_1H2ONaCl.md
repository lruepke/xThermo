

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


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  [**void**](classxThermal_1_1xThermalError.md) | [**dn1n2dP\_Tstar\_V**](#function-dn1n2dp_tstar_v) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & X, [**double**](classxThermal_1_1xThermalError.md) & dn1dP, [**double**](classxThermal_1_1xThermalError.md) & dn2dP) <br>_Calculate density of liquid or vapor phase._  |


























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
## Public Static Functions Documentation




### function dn1n2dP\_Tstar\_V 

_Calculate density of liquid or vapor phase._ 
```C++
static void xThermal::H2ONaCl::dn1n2dP_Tstar_V (
    const  double & P,
    const  double & X,
    double & dn1dP,
    double & dn2dP
) 
```



The density of H2ONaCl can be calculated from water density, see Eq. 7 of Driesner2007Part2. \(V_{solution}(P,T,X) = V_{water}(T^*,P)\), so \(\frac{M_{solution}}{\rho_{solution}} = \frac{M_{water}}{\rho_{water}}\).


Therefore, \(\rho_{solution} = \frac{M_{solution}}{M_{H_2O}} \rho_{H_2O}\), where \(M_{solution} = M_{H_2O}(1-X_{mol}) + M_{NaCl}X_{mol}\)




**Note:**

Density difference calculated from IAPWS95\_CoolProp and IAPS84 will cause H2ONaCl density difference up to 32 kg/m^3, test the following points:



|T[c]   |P[bar]   |X[mole fraction]   |rho[Driesner2007Part2]   |rho[[**xThermal**](namespacexThermal.md)]    |
|-----|-----|-----|-----|-----|
|8.500000e+02   |2.000000e+03   |9.996030e-01   |1.608149e+03   |1.602652e+03    |
|8.600000e+02   |2.500000e+03   |9.937537e-01   |1.615147e+03   |1.605205e+03    |
|9.200000e+02   |5.000000e+03   |9.915095e-01   |1.651829e+03   |1.619722e+03    |
|9.243300e+02   |5.000000e+03   |1.000000e+00   |1.656120e+03   |1.623869e+03   |






For example, T=920 deg.C, P = 5000 bar. The calculated \(T^*\) for density is 1279.92 deg.C. Density of water at (T=1279.92 deg.C, P=5000 bar) is 502.235 kg/m^3 (IAPWS95\_CoolProp), 512.153 kg/m^3 (IAPS84). The coefficient \(\frac{M_{H_2O}(1-X_{mol}) + M_{NaCl}X_{mol}}{M_{H_2O}} = 3.225\) at this P-T-X condition, therefore the final density different will be \(3.225\times(512.153-502.235)=32\) (This P-T condition exceeds the valid range of IAPWS-IF97).




**Todo**

Estimate density difference range based on the density difference between IAPWS95\_CoolProp and IAPS84.






**Todo**

How to deal with properties at high-T low-P region, e.g. sw.Rho\_phase(1.700000e+02+273.15,7.914706e+00\*1E5,0,Rhol, H2ONaCl::Liquid); see mmc4 of Driesner(2007b).






**Todo**

Merge Rho\_phase of liquid and vapor together to calculate them simultaneously, because property of saturated vapor and liquid always needed simultaneous, it could save computing time in this way. 



**Bug**

T=273.16, P=20E5, X=0.1, Tstar = n1 + n2\*T\_C + D will be negative.






**Parameters:**


* `T` [K] 
* `P` [Pa] 
* `X` [kg/kg] salinity of liquid or salinity of vapor, not bulk salinity 
* `Tstar` [K] 
* `rho` [ \(kg/m^3\)]

Pressure derivatives of \(n_1, n_2\): the expressions of n1n2\_Tstar\_V, differentiated analytically.


Needed for the chain rule of the density derivatives through \(T^*_V = n_1 + n_2 T\) (C12, LIB-25).




**Parameters:**


* `P` [bar] 
* `X` [mol/mol] 
* `dn1dP` [deg.C/bar] 
* `dn2dP` [1/bar] 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2ONaCl/H2ONaCl.cpp`

