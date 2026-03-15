

# Namespace xThermal::IAPWS95



[**Namespace List**](namespaces.md) **>** [**xThermal**](namespacexThermal.md) **>** [**IAPWS95**](namespacexThermal_1_1IAPWS95.md)




















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**Coeff\_phi\_o**](structxThermal_1_1IAPWS95_1_1Coeff__phi__o.md) <br> |
| struct | [**Coeff\_phi\_r**](structxThermal_1_1IAPWS95_1_1Coeff__phi__r.md) <br>_Numerical values of the coefficients and parameters of the residual part of the dimensionless Helmholtz free energy. Eq.6 in IAPWS-95_ [_**cIAPWS95::phi\_r**_](classxThermal_1_1IAPWS95_1_1cIAPWS95.md#function-phi_r-13) _._ |
| struct | [**Constants\_Viscosity2008\_Water**](structxThermal_1_1IAPWS95_1_1Constants__Viscosity2008__Water.md) <br>_Constants for water viscosity calculation, see mu2008._  |
| struct | [**HelmholtzEnergy\_dimensionless**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless.md) <br>_Dimensionless Helmholtz free energy_ \(\phi = f/(RT)\) _and its partial derivatives._ |
| struct | [**HelmholtzEnergy\_dimensionless\_SinglePhase**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless__SinglePhase.md) <br> |
| struct | [**HelmholtzEnergy\_dimensionless\_TwoPhase**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless__TwoPhase.md) <br> |
| struct | [**Params\_HP2RhoT**](structxThermal_1_1IAPWS95_1_1Params__HP2RhoT.md) <br> |
| struct | [**Params\_SolvePhaseEquilibrium**](structxThermal_1_1IAPWS95_1_1Params__SolvePhaseEquilibrium.md) <br>_Data struct for phase equilibrium solving. It is used in_ [_**func\_PhaseEquilibrium**_](namespacexThermal_1_1IAPWS95.md#function-func_phaseequilibrium) _,_cIAPWS95::Boiling\_P _and_cIAPWS95::Boiling\_T _._ |
| struct | [**Params\_TP2Rho**](structxThermal_1_1IAPWS95_1_1Params__TP2Rho.md) <br> |
| struct | [**Params\_T\_Sat\_estimate**](structxThermal_1_1IAPWS95_1_1Params__T__Sat__estimate.md) <br> |
| struct | [**STRUCT\_delta\_TwoPhase**](structxThermal_1_1IAPWS95_1_1STRUCT__delta__TwoPhase.md) <br> |
| struct | [**State**](structxThermal_1_1IAPWS95_1_1State.md) <br> |
| class | [**cIAPWS95**](classxThermal_1_1IAPWS95_1_1cIAPWS95.md) <br>_Class of IAPWS-95 formula of H2O, which inherits from_ [_**IAPWS\_IF97::cIAPWS\_IF97**_](classxThermal_1_1IAPWS__IF97_1_1cIAPWS__IF97.md) _._ |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef std::map&lt; [**unsigned**](classxThermal_1_1xThermalError.md) [**long**](classxThermal_1_1xThermalError.md) [**int**](classxThermal_1_1xThermalError.md), [**ThermodynamicProperty**](structxThermal_1_1ThermodynamicProperty.md) &gt; | [**MAP\_INDEX2PROP**](#typedef-map_index2prop)  <br> |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  [**int**](classxThermal_1_1xThermalError.md) | [**func\_HP2RhoT**](#function-func_hp2rhot) ([**const**](classxThermal_1_1xThermalError.md) [**gsl\_vector**](classxThermal_1_1xThermalError.md) \* x, [**void**](classxThermal_1_1xThermalError.md) \* params, [**gsl\_vector**](classxThermal_1_1xThermalError.md) \* f) <br>_Implementation of_ \(p(\delta, \tau), h(\delta, \tau)\) _used in SinglePhase\_HP._ |
|  [**int**](classxThermal_1_1xThermalError.md) | [**func\_PhaseEquilibrium**](#function-func_phaseequilibrium) ([**const**](classxThermal_1_1xThermalError.md) [**gsl\_vector**](classxThermal_1_1xThermalError.md) \* x, [**void**](classxThermal_1_1xThermalError.md) \* params, [**gsl\_vector**](classxThermal_1_1xThermalError.md) \* f) <br>_Construct nonlinear functions according to phase-equilibrium condition. See Table 3 (pp.12) in IAPWS-95._  |
|  [**int**](classxThermal_1_1xThermalError.md) | [**func\_TP2Rho**](#function-func_tp2rho) ([**const**](classxThermal_1_1xThermalError.md) [**gsl\_vector**](classxThermal_1_1xThermalError.md) \* x, [**void**](classxThermal_1_1xThermalError.md) \* params, [**gsl\_vector**](classxThermal_1_1xThermalError.md) \* f) <br>_Construct nonlinear function of_ \(f(\rho, p, T)\) _, which is used to solve density by given Temperautre and pressure._ |
|  [**double**](classxThermal_1_1xThermalError.md) | [**func\_TP2Rho**](#function-func_tp2rho) ([**double**](classxThermal_1_1xThermalError.md) rho, [**void**](classxThermal_1_1xThermalError.md) \* params) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**func\_T\_Sat\_estimate**](#function-func_t_sat_estimate) ([**double**](classxThermal_1_1xThermalError.md) T\_K, [**void**](classxThermal_1_1xThermalError.md) \* params) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**print\_state\_HP2RhoT**](#function-print_state_hp2rhot) ([**size\_t**](classxThermal_1_1xThermalError.md) iter, [**gsl\_multiroot\_fsolver**](classxThermal_1_1xThermalError.md) \* s) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**print\_state\_PhaseEquilibrium**](#function-print_state_phaseequilibrium) ([**size\_t**](classxThermal_1_1xThermalError.md) iter, [**gsl\_multiroot\_fsolver**](classxThermal_1_1xThermalError.md) \* s) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**print\_state\_TP2Rho**](#function-print_state_tp2rho) ([**size\_t**](classxThermal_1_1xThermalError.md) iter, [**gsl\_multiroot\_fsolver**](classxThermal_1_1xThermalError.md) \* s) <br> |




























## Public Types Documentation




### typedef MAP\_INDEX2PROP 

```C++
typedef std::map<unsigned long int, ThermodynamicProperty> xThermal::IAPWS95::MAP_INDEX2PROP;
```




<hr>
## Public Functions Documentation




### function func\_HP2RhoT 

_Implementation of_ \(p(\delta, \tau), h(\delta, \tau)\) _used in SinglePhase\_HP._
```C++
int xThermal::IAPWS95::func_HP2RhoT (
    const  gsl_vector * x,
    void * params,
    gsl_vector * f
) 
```





**Parameters:**


* `x` 
* `params` 
* `f` 



**Returns:**

int 





        

<hr>



### function func\_PhaseEquilibrium 

_Construct nonlinear functions according to phase-equilibrium condition. See Table 3 (pp.12) in IAPWS-95._ 
```C++
int xThermal::IAPWS95::func_PhaseEquilibrium (
    const  gsl_vector * x,
    void * params,
    gsl_vector * f
) 
```





**Parameters:**


* `x` 
* `params` 
* `f` 



**Returns:**

int 





        

<hr>



### function func\_TP2Rho 

_Construct nonlinear function of_ \(f(\rho, p, T)\) _, which is used to solve density by given Temperautre and pressure._
```C++
int xThermal::IAPWS95::func_TP2Rho (
    const  gsl_vector * x,
    void * params,
    gsl_vector * f
) 
```





**Parameters:**


* `x` 
* `params` 
* `f` 



**Returns:**

int 





        

<hr>



### function func\_TP2Rho 

```C++
double xThermal::IAPWS95::func_TP2Rho (
    double rho,
    void * params
) 
```




<hr>



### function func\_T\_Sat\_estimate 

```C++
double xThermal::IAPWS95::func_T_Sat_estimate (
    double T_K,
    void * params
) 
```




<hr>



### function print\_state\_HP2RhoT 

```C++
void xThermal::IAPWS95::print_state_HP2RhoT (
    size_t iter,
    gsl_multiroot_fsolver * s
) 
```




<hr>



### function print\_state\_PhaseEquilibrium 

```C++
void xThermal::IAPWS95::print_state_PhaseEquilibrium (
    size_t iter,
    gsl_multiroot_fsolver * s
) 
```




<hr>



### function print\_state\_TP2Rho 

```C++
void xThermal::IAPWS95::print_state_TP2Rho (
    size_t iter,
    gsl_multiroot_fsolver * s
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2O/IAPWS95.cpp`

