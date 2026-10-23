

# Namespace xThermal::PROST



[**Namespace List**](namespaces.md) **>** [**xThermal**](namespacexThermal.md) **>** [**PROST**](namespacexThermal_1_1PROST.md)




















## Classes

| Type | Name |
| ---: | :--- |
| class | [**cIAPS84**](classxThermal_1_1PROST_1_1cIAPS84.md) <br> |
























## Public Static Functions

| Type | Name |
| ---: | :--- |
|  [**double**](classxThermal_1_1xThermalError.md) | [**bisect\_density**](#function-bisect_density) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**double**](classxThermal_1_1xThermalError.md) lo, [**double**](classxThermal_1_1xThermalError.md) hi) <br> |
|  [**bool**](classxThermal_1_1xThermalError.md) | [**density\_residual\_ok**](#function-density_residual_ok) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**Prop**](classxThermal_1_1xThermalError.md) \* prop) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**eos\_at**](#function-eos_at) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & d, [**double**](classxThermal_1_1xThermalError.md) & p, [**double**](classxThermal_1_1xThermalError.md) & dpdrho, [**double**](classxThermal_1_1xThermalError.md) & g) <br> |
|  [**Prop**](classxThermal_1_1xThermalError.md) \* | [**near\_critical\_band\_state**](#function-near_critical_band_state) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**bool**](classxThermal_1_1xThermalError.md) & liquid) <br> |
|  [**bool**](classxThermal_1_1xThermalError.md) | [**stable\_on\_phase\_side**](#function-stable_on_phase_side) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**Prop**](classxThermal_1_1xThermalError.md) \* prop) <br> |
|  [**double**](classxThermal_1_1xThermalError.md) | [**viscosity\_or\_nan**](#function-viscosity_or_nan) ([**Prop**](classxThermal_1_1xThermalError.md) \* prop) <br> |
|  [**void**](classxThermal_1_1xThermalError.md) | [**water\_tp\_verified**](#function-water_tp_verified) ([**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & T, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & P, [**const**](classxThermal_1_1xThermalError.md) [**double**](classxThermal_1_1xThermalError.md) & dp, [**Prop**](classxThermal_1_1xThermalError.md) \* prop) <br> |


























## Public Static Functions Documentation




### function bisect\_density 

```C++
static double xThermal::PROST::bisect_density (
    const  double & T,
    const  double & P,
    double lo,
    double hi
) 
```




<hr>



### function density\_residual\_ok 

```C++
static bool xThermal::PROST::density_residual_ok (
    const  double & T,
    const  double & P,
    Prop * prop
) 
```




<hr>



### function eos\_at 

```C++
static void xThermal::PROST::eos_at (
    const  double & T,
    const  double & d,
    double & p,
    double & dpdrho,
    double & g
) 
```




<hr>



### function near\_critical\_band\_state 

```C++
static Prop * xThermal::PROST::near_critical_band_state (
    const  double & T,
    const  double & P,
    bool & liquid
) 
```




<hr>



### function stable\_on\_phase\_side 

```C++
static bool xThermal::PROST::stable_on_phase_side (
    const  double & T,
    const  double & P,
    Prop * prop
) 
```




<hr>



### function viscosity\_or\_nan 

```C++
static double xThermal::PROST::viscosity_or_nan (
    Prop * prop
) 
```




<hr>



### function water\_tp\_verified 

```C++
static void xThermal::PROST::water_tp_verified (
    const  double & T,
    const  double & P,
    const  double & dp,
    Prop * prop
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/H2O/IAPS84.cpp`

