

# Group Initialization\_CommonlyUsedProps



[**Modules**](modules.md) **>** [**Initialization\_CommonlyUsedProps**](group__Initialization__CommonlyUsedProps.md)





































































## Macros

| Type | Name |
| ---: | :--- |
| define  | [**thermodynamicProperty\_H**](group__Initialization__CommonlyUsedProps.md#define-thermodynamicproperty_h)  `{0, "Specific enthalpy",    "$H$",      "J/kg",  [**dimension\_SI\_H**](group__DIMENSION__COMMONLYUSED.md#define-dimension_si_h)}`<br> |
| define  | [**thermodynamicProperty\_P**](group__Initialization__CommonlyUsedProps.md#define-thermodynamicproperty_p)  `{0, "Pressure",             "$p$",      "Pa",    [**dimension\_SI\_P**](group__DIMENSION__COMMONLYUSED.md#define-dimension_si_p)}`<br> |
| define  | [**thermodynamicProperty\_Rho**](group__Initialization__CommonlyUsedProps.md#define-thermodynamicproperty_rho)  `{0, "Density",              "$\rho$",   "kg/m3", [**dimension\_SI\_Rho**](group__DIMENSION__COMMONLYUSED.md#define-dimension_si_rho)}`<br> |
| define  | [**thermodynamicProperty\_T**](group__Initialization__CommonlyUsedProps.md#define-thermodynamicproperty_t)  `{0, "Temperature",          "$T$",      "K",     [**dimension\_SI\_T**](group__DIMENSION__COMMONLYUSED.md#define-dimension_si_t)}`<br> |
| define  | [**thermodynamicProperty\_Undefined**](group__Initialization__CommonlyUsedProps.md#define-thermodynamicproperty_undefined)  `{0, "Undefined",  "Undefined",      "Undefined",     [**dimension\_SI\_None**](group__DIMENSION__COMMONLYUSED.md#define-dimension_si_none)}`<br> |
| define  | [**thermodynamicProperty\_X**](group__Initialization__CommonlyUsedProps.md#define-thermodynamicproperty_x)  `{0, "Vapor mass fraction",  "$x$",      "1",     [**dimension\_SI\_None**](group__DIMENSION__COMMONLYUSED.md#define-dimension_si_none)}`<br> |

## Macro Definition Documentation





### define thermodynamicProperty\_H 

```
#define thermodynamicProperty_H `{0, "Specific enthalpy",    "$H$",      "J/kg", dimension_SI_H }`
```



initialization of property specific enthalpy 


        

<hr>



### define thermodynamicProperty\_P 

```
#define thermodynamicProperty_P `{0, "Pressure",             "$p$",      "Pa", dimension_SI_P }`
```



initialization of property pressure 


        

<hr>



### define thermodynamicProperty\_Rho 

```
#define thermodynamicProperty_Rho `{0, "Density",              "$\rho$",   "kg/m3", dimension_SI_Rho }`
```



initialization of property density 


        

<hr>



### define thermodynamicProperty\_T 

```
#define thermodynamicProperty_T `{0, "Temperature",          "$T$",      "K", dimension_SI_T }`
```



initialization of property temperature 


        

<hr>



### define thermodynamicProperty\_Undefined 

```
#define thermodynamicProperty_Undefined `{0, "Undefined",  "Undefined",      "Undefined", dimension_SI_None }`
```



Undefined property 


        

<hr>



### define thermodynamicProperty\_X 

```
#define thermodynamicProperty_X `{0, "Vapor mass fraction",  "$x$",      "1", dimension_SI_None }`
```



initialization of property vapor mass fraction: \(x=\frac{m^{\prime\prime}}{m^{\prime} + m^{\prime}} = \frac{h - h^{\prime}}{h^{\prime\prime} - h^{\prime}} =\frac{v - v^{\prime}}{v^{\prime\prime} - v^{\prime}} = \frac{1/\rho - 1/\rho^{\prime}}{1/\rho^{\prime\prime} - 1/\rho^{\prime}}\), where the superscript \(\prime, \prime\prime\) denote liquid and vapor, respectively. 


        

<hr>

------------------------------


