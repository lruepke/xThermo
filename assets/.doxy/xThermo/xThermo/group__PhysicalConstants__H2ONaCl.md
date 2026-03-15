

# Group PhysicalConstants\_H2ONaCl



[**Modules**](modules.md) **>** [**PhysicalConstants\_H2ONaCl**](group__PhysicalConstants__H2ONaCl.md)



[More...](#detailed-description)






















## Public Attributes

| Type | Name |
| ---: | :--- |
|  const double | [**P\_MAX**](#variable-p_max)   = `5000E5`<br> |
|  const double | [**P\_MIN**](#variable-p_min)   = `1E5`<br> |
|  const double | [**P\_Peak\_VLH**](#variable-p_peak_vlh)   = `39014744.433797`<br> |
|  const double | [**P\_Peak\_VLH\_bar**](#variable-p_peak_vlh_bar)   = `390.14744433797`<br> |
|  const double | [**T\_MAX**](#variable-t_max)   = `1273.15`<br> |
|  const double | [**T\_MAX\_VLH**](#variable-t_max_vlh)   = `1073.5662157838103`<br> |
|  const double | [**T\_MIN**](#variable-t_min)   = `273.16`<br> |
|  const double | [**T\_MIN\_VLH**](#variable-t_min_vlh)   = `380.912102015643`<br> |
|  const double | [**T\_Peak\_VLH**](#variable-t_peak_vlh)   = `867.782443`<br> |
|  const double | [**T\_Peak\_VLH\_C**](#variable-t_peak_vlh_c)   = `867.782443 - 273.15`<br> |
|  const double | [**X\_MAX**](#variable-x_max)   = `1.0`<br> |
|  const double | [**X\_MIN**](#variable-x_min)   = `0.0`<br> |












































## Detailed Description




**Todo**

Need to discuss how to deal with the discrepancy between the critical point of H2O in the original paper (Driesner2007Part1) and that in IAPWS95\_CoolProp release.




    
## Public Attributes Documentation




### variable P\_MAX 

```
const double xThermal::H2ONaCl::P_MAX;
```



Maximum valid pressure [Pa], same as that of H2O 


        

<hr>



### variable P\_MIN 

```
const double xThermal::H2ONaCl::P_MIN;
```



Minimum valid pressure 1bar = 1E5 [Pa]. See Driesner2007Part1 


        

<hr>



### variable P\_Peak\_VLH 

```
const double xThermal::H2ONaCl::P_Peak_VLH;
```



Pressure [Pa] at peak of VLH zone. Calculate by P\_VLH function and [**T\_Peak\_VLH**](group__PhysicalConstants__H2ONaCl.md#variable-t_peak_vlh) constant 


        

<hr>



### variable P\_Peak\_VLH\_bar 

```
const double xThermal::H2ONaCl::P_Peak_VLH_bar;
```



Pressure [bar] at peak of VLH zone. Calculate by P\_VLH function and [**T\_Peak\_VLH**](group__PhysicalConstants__H2ONaCl.md#variable-t_peak_vlh) constant 


        

<hr>



### variable T\_MAX 

```
const double xThermal::H2ONaCl::T_MAX;
```



Maximum valid temperature [K] = 1000 \(^{\circ}\)C 


        

<hr>



### variable T\_MAX\_VLH 

```
const double xThermal::H2ONaCl::T_MAX_VLH;
```



Maximum temperature [K] of the VLH zone. Calculated by function [**cH2ONaCl::Tmax\_VLH**](classxThermal_1_1H2ONaCl_1_1cH2ONaCl.md#function-tmax_vlh) 


        

<hr>



### variable T\_MIN 

```
const double xThermal::H2ONaCl::T_MIN;
```



Minimum valid temperature 2 [deg.C]. See Driesner2007Part1 


        

<hr>



### variable T\_MIN\_VLH 

```
const double xThermal::H2ONaCl::T_MIN_VLH;
```



Minimum temperature [K] of the VLH zone. Calculated by function cH2ONaCl::Tmin\_VLH 


        

<hr>



### variable T\_Peak\_VLH 

```
const double xThermal::H2ONaCl::T_Peak_VLH;
```



Temperature [K] at peak of VLH zone. Calculated by function cH2ONaCl::T\_Pmax\_VLH 
 


        

<hr>



### variable T\_Peak\_VLH\_C 

```
const double xThermal::H2ONaCl::T_Peak_VLH_C;
```



Temperature [deg.C] at peak of VLH zone. Calculated by function cH2ONaCl::T\_Pmax\_VLH 
 


        

<hr>



### variable X\_MAX 

```
const double xThermal::H2ONaCl::X_MAX;
```



Maximum valid salinity [kg/kg] 


        

<hr>



### variable X\_MIN 

```
const double xThermal::H2ONaCl::X_MIN;
```



Minimum valid salinity [kg/kg] 


        

<hr>

------------------------------


