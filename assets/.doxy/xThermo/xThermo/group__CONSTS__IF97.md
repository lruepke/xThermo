

# Group CONSTS\_IF97



[**Modules**](modules.md) **>** [**CONSTS\_IF97**](group__CONSTS__IF97.md)





































































## Macros

| Type | Name |
| ---: | :--- |
| define  | [**CONST\_IF97\_H\_c**](group__CONSTS__IF97.md#define-const_if97_h_c)  `2.08754684511650027707E+06`<br> |
| define  | [**CONST\_IF97\_Hl\_Tmin\_Region3**](group__CONSTS__IF97.md#define-const_if97_hl_tmin_region3)  `1.670858218E6`<br> |
| define  | [**CONST\_IF97\_Hv\_Tmin\_Region3**](group__CONSTS__IF97.md#define-const_if97_hv_tmin_region3)  `2.563592004E6`<br> |
| define  | [**CONST\_IF97\_P\_Region2a2b**](group__CONSTS__IF97.md#define-const_if97_p_region2a2b)  `4E6`<br> |
| define  | [**CONST\_IF97\_Pmax\_Region1**](group__CONSTS__IF97.md#define-const_if97_pmax_region1)  `100E6`<br> |
| define  | [**CONST\_IF97\_Pmax\_Region5**](group__CONSTS__IF97.md#define-const_if97_pmax_region5)  `50E6`<br> |
| define  | [**CONST\_IF97\_Pmin\_Region2c**](group__CONSTS__IF97.md#define-const_if97_pmin_region2c)  `6.54670E6`<br> |
| define  | [**CONST\_IF97\_Pmin\_Region3**](group__CONSTS__IF97.md#define-const_if97_pmin_region3)  `16.5292E6`<br> |
| define  | [**CONST\_IF97\_S\_Region2b2c**](group__CONSTS__IF97.md#define-const_if97_s_region2b2c)  `5.85E3`<br> |
| define  | [**CONST\_IF97\_Tmax\_Region2**](group__CONSTS__IF97.md#define-const_if97_tmax_region2)  `1073.15`<br> |
| define  | [**CONST\_IF97\_Tmax\_Region5**](group__CONSTS__IF97.md#define-const_if97_tmax_region5)  `2273.15`<br> |
| define  | [**CONST\_IF97\_Tmin\_Region3**](group__CONSTS__IF97.md#define-const_if97_tmin_region3)  `623.15`<br> |

## Macro Definition Documentation





### define CONST\_IF97\_H\_c 

```
#define CONST_IF97_H_c `2.08754684511650027707E+06`
```



Specific enthalpy of critical point, calculate from IAPWS\_IF97::cIAPWS\_IF97::Boundary\_region3ab\_P2H (H2O::P\_c) 
 


        

<hr>



### define CONST\_IF97\_Hl\_Tmin\_Region3 

```
#define CONST_IF97_Hl_Tmin_Region3 `1.670858218E6`
```



Saturated liquid enthalpy at [**CONST\_IF97\_Tmin\_Region3**](group__CONSTS__IF97.md#define-const_if97_tmin_region3) , see section 4.3, pp.17 of IF97-Region3. 


        

<hr>



### define CONST\_IF97\_Hv\_Tmin\_Region3 

```
#define CONST_IF97_Hv_Tmin_Region3 `2.563592004E6`
```



Saturated vapor enthalpy at [**CONST\_IF97\_Tmin\_Region3**](group__CONSTS__IF97.md#define-const_if97_tmin_region3) , see section 4.3, pp.17 of IF97-Region3. 


        

<hr>



### define CONST\_IF97\_P\_Region2a2b 

```
#define CONST_IF97_P_Region2a2b `4E6`
```



The boundary between the subregions 2a and 2b is the isobar \(p = 4 MPa\), see IF97 


        

<hr>



### define CONST\_IF97\_Pmax\_Region1 

```
#define CONST_IF97_Pmax_Region1 `100E6`
```



Maximum pressure [MPa] of the region 1, 2, 3, 4 


        

<hr>



### define CONST\_IF97\_Pmax\_Region5 

```
#define CONST_IF97_Pmax_Region5 `50E6`
```



Maximum pressure [MPa] of the region 5 


        

<hr>



### define CONST\_IF97\_Pmin\_Region2c 

```
#define CONST_IF97_Pmin_Region2c `6.54670E6`
```



Equations (20) and (21) give the boundary line between subregions 2b and 2c from the saturation state at T=554.485 K and \(p_s\) =6.54670 MPa to T=1019.32 K and p=100 MPa. See IF97 


        

<hr>



### define CONST\_IF97\_Pmin\_Region3 

```
#define CONST_IF97_Pmin_Region3 `16.5292E6`
```



Minimum pressure of the region 3. See pp.6 of IF97 


        

<hr>



### define CONST\_IF97\_S\_Region2b2c 

```
#define CONST_IF97_S_Region2b2c `5.85E3`
```



the boundary between the subregions 2b and 2c corresponds to the entropy line s = 5.85 kJ/kg/K. See section 6.3 of IF97 


        

<hr>



### define CONST\_IF97\_Tmax\_Region2 

```
#define CONST_IF97_Tmax_Region2 `1073.15`
```



Maximum temperature [K] of the retion 2 


        

<hr>



### define CONST\_IF97\_Tmax\_Region5 

```
#define CONST_IF97_Tmax_Region5 `2273.15`
```



Maximum temperature [K] of the retion 5 


        

<hr>



### define CONST\_IF97\_Tmin\_Region3 

```
#define CONST_IF97_Tmin_Region3 `623.15`
```



Maximum temperature of the region 3. See pp.6 of IF97. 


        

<hr>

------------------------------


