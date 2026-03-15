

# Struct xThermal::IAPWS95::State



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**IAPWS95**](namespacexThermal_1_1IAPWS95.md) **>** [**State**](structxThermal_1_1IAPWS95_1_1State.md)


























## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**union**](classxThermal_1_1xThermalError.md) [**xThermal::IAPWS95::State**](structxThermal_1_1IAPWS95_1_1State.md) | [**delta**](#variable-delta)  <br> |
|  [**PhaseRegion**](namespacexThermal.md#enum-phaseregion) | [**phase**](#variable-phase)  <br> |
|  [**union**](classxThermal_1_1xThermalError.md) [**xThermal::IAPWS95::State**](structxThermal_1_1IAPWS95_1_1State.md) | [**phi**](#variable-phi)  <br> |
|  [**HelmholtzEnergy\_dimensionless\_SinglePhase**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless__SinglePhase.md) | [**singlePhase**](#variable-singlephase-12)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**singlePhase**](#variable-singlephase-22)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**tau**](#variable-tau)  <br> |
|  [**HelmholtzEnergy\_dimensionless\_TwoPhase**](structxThermal_1_1IAPWS95_1_1HelmholtzEnergy__dimensionless__TwoPhase.md) | [**twoPhase**](#variable-twophase-12)  <br> |
|  [**STRUCT\_delta\_TwoPhase**](structxThermal_1_1IAPWS95_1_1STRUCT__delta__TwoPhase.md) | [**twoPhase**](#variable-twophase-22)  <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**x**](#variable-x)   = `0`<br> |












































## Public Attributes Documentation




### variable delta 

```C++
union xThermal::IAPWS95::State xThermal::IAPWS95::State::delta;
```




<hr>



### variable phase 

```C++
PhaseRegion xThermal::IAPWS95::State::phase;
```




<hr>



### variable phi 

```C++
union xThermal::IAPWS95::State xThermal::IAPWS95::State::phi;
```




<hr>



### variable singlePhase [1/2]

```C++
HelmholtzEnergy_dimensionless_SinglePhase xThermal::IAPWS95::State::singlePhase;
```




<hr>



### variable singlePhase [2/2]

```C++
double xThermal::IAPWS95::State::singlePhase;
```




<hr>



### variable tau 

```C++
double xThermal::IAPWS95::State::tau;
```




<hr>



### variable twoPhase [1/2]

```C++
HelmholtzEnergy_dimensionless_TwoPhase xThermal::IAPWS95::State::twoPhase;
```




<hr>



### variable twoPhase [2/2]

```C++
STRUCT_delta_TwoPhase xThermal::IAPWS95::State::twoPhase;
```




<hr>



### variable x 

```C++
double xThermal::IAPWS95::State::x;
```



vapor mass fraction 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2O/IAPWS95.h`

