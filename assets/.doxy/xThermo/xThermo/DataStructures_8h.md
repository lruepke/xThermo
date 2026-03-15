

# File DataStructures.h



[**FileList**](files.md) **>** [**Library**](dir_5ad7f572bbca03234e8e621e192fc099.md) **>** [**thermo**](dir_d760ccf1b5c74bc66b0c51c2e0ac61aa.md) **>** [**DataStructures.h**](DataStructures_8h.md)

[Go to the source code of this file](DataStructures_8h_source.md)

[More...](#detailed-description)

* `#include "stdfunc.h"`
* `#include <iostream>`
* `#include <fstream>`
* `#include <sstream>`
* `#include <utility>`
* `#include <vector>`
* `#include <algorithm>`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**xThermal**](namespacexThermal.md) <br>_Namespace of_ [_**xThermal**_](namespacexThermal.md) _library._ |


## Classes

| Type | Name |
| ---: | :--- |
| struct | [**COLOR**](structxThermal_1_1COLOR.md) <br> |
| struct | [**CONSTENTS\_Thermo**](structxThermal_1_1CONSTENTS__Thermo.md) <br>_为了实现多个H2O EOS的backend，必须使用虚函数调用相应的参数，比如T\_critical()，但是在子类中频繁调用函数的性能肯定很低，所以将所有热力学常数打包为一个struct类型，作为子类的成员变量，然后在构造函数中进行初始化，后面需要常数的地方直接访问成员变量即可，可提高性能。_  |
| struct | [**DeformLinearMesh**](structxThermal_1_1DeformLinearMesh.md) <br> |
| struct | [**Head\_AMR\_LUT**](structxThermal_1_1Head__AMR__LUT.md) <br> |
| struct | [**Line**](structxThermal_1_1Line.md) <br> |
| struct | [**Line\_slice**](structxThermal_1_1Line__slice.md) <br> |
| struct | [**PhaseBoundaries**](structxThermal_1_1PhaseBoundaries.md) <br> |
| struct | [**PhaseRegion\_Slice**](structxThermal_1_1PhaseRegion__Slice.md) <br> |
| struct | [**PhysicalDimension**](structxThermal_1_1PhysicalDimension.md) <br>_Physical dimension of a quantity can be expressed in terms of 7 basic SI unit, e.g. dimension of density is_ \(kg/m^3\) _, can be expressed in vector of [1, -3, 0, 0, 0, 0, 0]. This idea comes from OpenFOAM._ |
| struct | [**Point**](structxThermal_1_1Point.md) <br> |
| struct | [**Point\_slice**](structxThermal_1_1Point__slice.md) <br> |
| struct | [**Polygon\_slice**](structxThermal_1_1Polygon__slice.md) <br> |
| struct | [**Surface**](structxThermal_1_1Surface.md) <br> |
| struct | [**ThermodynamicProperties**](structxThermal_1_1ThermodynamicProperties.md) <br> |
| struct | [**ThermodynamicPropertiesArray**](structxThermal_1_1ThermodynamicPropertiesArray.md) <br> |
| struct | [**ThermodynamicPropertiesVector**](structxThermal_1_1ThermodynamicPropertiesVector.md) <br> |
| struct | [**ThermodynamicProperty**](structxThermal_1_1ThermodynamicProperty.md) <br>_Data struct of a thermodynamic property._  |
| struct | [**TriMesh**](structxThermal_1_1TriMesh.md) <br> |
| struct | [**propInfo**](structxThermal_1_1propInfo.md) <br>_Information of a thermodynamic property._  |

















































## Macros

| Type | Name |
| ---: | :--- |
| define  | [**COLOR\_LH**](DataStructures_8h.md#define-color_lh)  `COLOR\_lightblue`<br> |
| define  | [**COLOR\_VH**](DataStructures_8h.md#define-color_vh)  `COLOR\_darkgray`<br> |
| define  | [**COLOR\_VL**](DataStructures_8h.md#define-color_vl)  `COLOR\_plum`<br> |
| define  | [**COLOR\_VLH**](DataStructures_8h.md#define-color_vlh)  `{1.0, 0.7137254901960784, 0.37254901960784315}`<br> |
| define  | [**COLOR\_VLH\_highT**](DataStructures_8h.md#define-color_vlh_hight)  `COLOR\_yellow`<br> |
| define  | [**COLOR\_VLH\_lowT**](DataStructures_8h.md#define-color_vlh_lowt)  `COLOR\_navajowhite`<br> |
| define  | [**COLOR\_blue**](DataStructures_8h.md#define-color_blue)  `{0.000000, 0.000000, 1.000000}`<br> |
| define  | [**COLOR\_crimson**](DataStructures_8h.md#define-color_crimson)  `{0.862745, 0.078431, 0.235294}`<br> |
| define  | [**COLOR\_darkgray**](DataStructures_8h.md#define-color_darkgray)  `{0.662745, 0.662745, 0.662745}`<br> |
| define  | [**COLOR\_darkgray**](DataStructures_8h.md#define-color_darkgray)  `{0.662745, 0.662745, 0.662745}`<br> |
| define  | [**COLOR\_deeppink**](DataStructures_8h.md#define-color_deeppink)  `{1.000000, 0.078431, 0.576471}`<br> |
| define  | [**COLOR\_gainsboro**](DataStructures_8h.md#define-color_gainsboro)  `{0.862745, 0.862745, 0.862745}`<br> |
| define  | [**COLOR\_gold**](DataStructures_8h.md#define-color_gold)  `{1.000000, 0.843137, 0.000000}`<br> |
| define  | [**COLOR\_k**](DataStructures_8h.md#define-color_k)  `{0.000000, 0.000000, 0.000000}`<br> |
| define  | [**COLOR\_lc\_LH**](DataStructures_8h.md#define-color_lc_lh)  `COLOR\_darkgray`<br> |
| define  | [**COLOR\_lc\_VH**](DataStructures_8h.md#define-color_lc_vh)  `{1.0, 0.5137254901960784, 1.0}`<br> |
| define  | [**COLOR\_lc\_VLH**](DataStructures_8h.md#define-color_lc_vlh)  `{1.0, 0.7137254901960784, 0.37254901960784315}`<br> |
| define  | [**COLOR\_lc\_VL\_L**](DataStructures_8h.md#define-color_lc_vl_l)  `{0.0, 0.6745098039215687, 1.0}`<br> |
| define  | [**COLOR\_lc\_VL\_V**](DataStructures_8h.md#define-color_lc_vl_v)  `{0.0, 0.6784313725490196, 0.6745098039215687}`<br> |
| define  | [**COLOR\_lc\_twoPhaseWater**](DataStructures_8h.md#define-color_lc_twophasewater)  `COLOR\_yellow`<br> |
| define  | [**COLOR\_lightblue**](DataStructures_8h.md#define-color_lightblue)  `{0.678431, 0.847059, 0.901961}`<br> |
| define  | [**COLOR\_lightgray**](DataStructures_8h.md#define-color_lightgray)  `{0.827451, 0.827451, 0.827451}`<br> |
| define  | [**COLOR\_lime**](DataStructures_8h.md#define-color_lime)  `{0.000000, 1.000000, 0.000000}`<br> |
| define  | [**COLOR\_mec\_CriticalPoint**](DataStructures_8h.md#define-color_mec_criticalpoint)  `COLOR\_k`<br> |
| define  | [**COLOR\_mec\_twoPhaseWater**](DataStructures_8h.md#define-color_mec_twophasewater)  `COLOR\_crimson`<br> |
| define  | [**COLOR\_mfc\_CriticalPoint**](DataStructures_8h.md#define-color_mfc_criticalpoint)  `COLOR\_deeppink`<br> |
| define  | [**COLOR\_mfc\_twoPhaseWater**](DataStructures_8h.md#define-color_mfc_twophasewater)  `COLOR\_gold`<br> |
| define  | [**COLOR\_navajowhite**](DataStructures_8h.md#define-color_navajowhite)  `{1.0, 0.8705882352941177, 0.6784313725490196}`<br> |
| define  | [**COLOR\_orange**](DataStructures_8h.md#define-color_orange)  `{1.0, 0.6470588235294118, 0.0}`<br> |
| define  | [**COLOR\_plum**](DataStructures_8h.md#define-color_plum)  `{0.866667, 0.627451, 0.866667}`<br> |
| define  | [**COLOR\_red**](DataStructures_8h.md#define-color_red)  `{1.000000, 0.000000, 0.000000}`<br> |
| define  | [**COLOR\_yellow**](DataStructures_8h.md#define-color_yellow)  `{1.000000, 1.000000, 0.000000}`<br> |
| define  | [**SCALE\_X\_linear**](DataStructures_8h.md#define-scale_x_linear)  `1`<br> |
| define  | [**SCALE\_X\_log**](DataStructures_8h.md#define-scale_x_log)  `2`<br> |
| define  | [**SCALE\_X\_loglinear**](DataStructures_8h.md#define-scale_x_loglinear)  `3`<br> |
| define  | [**STR\_LENGTH\_PROPINFO**](DataStructures_8h.md#define-str_length_propinfo)  `30`<br> |
| define  | [**Update\_prop\_CommonlyUsedHPX**](DataStructures_8h.md#define-update_prop_commonlyusedhpx)  `(Update\_prop\_H \| Update\_prop\_T \| Update\_prop\_Rho \| Update\_prop\_Cp \| Update\_prop\_Mu \| Update\_prop\_IsothermalCompressibility \|Update\_prop\_IsobaricExpansivity )`<br> |
| define  | [**Update\_prop\_CommonlyUsedTPX**](DataStructures_8h.md#define-update_prop_commonlyusedtpx)  `(Update\_prop\_Rho \| Update\_prop\_Cp \| Update\_prop\_Mu \| Update\_prop\_IsothermalCompressibility \|Update\_prop\_IsobaricExpansivity )`<br> |
| define  | [**Update\_prop\_Cp**](DataStructures_8h.md#define-update_prop_cp)  `32`<br> |
| define  | [**Update\_prop\_H**](DataStructures_8h.md#define-update_prop_h)  `4`<br> |
| define  | [**Update\_prop\_IsobaricExpansivity**](DataStructures_8h.md#define-update_prop_isobaricexpansivity)  `256`<br> |
| define  | [**Update\_prop\_IsothermalCompressibility**](DataStructures_8h.md#define-update_prop_isothermalcompressibility)  `128`<br> |
| define  | [**Update\_prop\_Mu**](DataStructures_8h.md#define-update_prop_mu)  `64`<br> |
| define  | [**Update\_prop\_Rho**](DataStructures_8h.md#define-update_prop_rho)  `2`<br> |
| define  | [**Update\_prop\_T**](DataStructures_8h.md#define-update_prop_t)  `16`<br> |

## Detailed Description




**Author:**

Zhikui Guo ([zguo@geomar.de](mailto:zguo@geomar.de)) 




**Version:**

0.1 




**Date:**

2022-03-27




**Copyright:**

Copyright (c) 2022 





    
## Macro Definition Documentation





### define COLOR\_LH 

```C++
#define COLOR_LH `COLOR_lightblue`
```




<hr>



### define COLOR\_VH 

```C++
#define COLOR_VH `COLOR_darkgray`
```




<hr>



### define COLOR\_VL 

```C++
#define COLOR_VL `COLOR_plum`
```




<hr>



### define COLOR\_VLH 

```C++
#define COLOR_VLH `{1.0, 0.7137254901960784, 0.37254901960784315}`
```




<hr>



### define COLOR\_VLH\_highT 

```C++
#define COLOR_VLH_highT `COLOR_yellow`
```




<hr>



### define COLOR\_VLH\_lowT 

```C++
#define COLOR_VLH_lowT `COLOR_navajowhite`
```




<hr>



### define COLOR\_blue 

```C++
#define COLOR_blue `{0.000000, 0.000000, 1.000000}`
```




<hr>



### define COLOR\_crimson 

```C++
#define COLOR_crimson `{0.862745, 0.078431, 0.235294}`
```




<hr>



### define COLOR\_darkgray 

```C++
#define COLOR_darkgray `{0.662745, 0.662745, 0.662745}`
```




<hr>



### define COLOR\_darkgray 

```C++
#define COLOR_darkgray `{0.662745, 0.662745, 0.662745}`
```




<hr>



### define COLOR\_deeppink 

```C++
#define COLOR_deeppink `{1.000000, 0.078431, 0.576471}`
```




<hr>



### define COLOR\_gainsboro 

```C++
#define COLOR_gainsboro `{0.862745, 0.862745, 0.862745}`
```




<hr>



### define COLOR\_gold 

```C++
#define COLOR_gold `{1.000000, 0.843137, 0.000000}`
```




<hr>



### define COLOR\_k 

```C++
#define COLOR_k `{0.000000, 0.000000, 0.000000}`
```




<hr>



### define COLOR\_lc\_LH 

```C++
#define COLOR_lc_LH `COLOR_darkgray`
```




<hr>



### define COLOR\_lc\_VH 

```C++
#define COLOR_lc_VH `{1.0, 0.5137254901960784, 1.0}`
```




<hr>



### define COLOR\_lc\_VLH 

```C++
#define COLOR_lc_VLH `{1.0, 0.7137254901960784, 0.37254901960784315}`
```




<hr>



### define COLOR\_lc\_VL\_L 

```C++
#define COLOR_lc_VL_L `{0.0, 0.6745098039215687, 1.0}`
```




<hr>



### define COLOR\_lc\_VL\_V 

```C++
#define COLOR_lc_VL_V `{0.0, 0.6784313725490196, 0.6745098039215687}`
```




<hr>



### define COLOR\_lc\_twoPhaseWater 

```C++
#define COLOR_lc_twoPhaseWater `COLOR_yellow`
```




<hr>



### define COLOR\_lightblue 

```C++
#define COLOR_lightblue `{0.678431, 0.847059, 0.901961}`
```




<hr>



### define COLOR\_lightgray 

```C++
#define COLOR_lightgray `{0.827451, 0.827451, 0.827451}`
```




<hr>



### define COLOR\_lime 

```C++
#define COLOR_lime `{0.000000, 1.000000, 0.000000}`
```




<hr>



### define COLOR\_mec\_CriticalPoint 

```C++
#define COLOR_mec_CriticalPoint `COLOR_k`
```




<hr>



### define COLOR\_mec\_twoPhaseWater 

```C++
#define COLOR_mec_twoPhaseWater `COLOR_crimson`
```




<hr>



### define COLOR\_mfc\_CriticalPoint 

```C++
#define COLOR_mfc_CriticalPoint `COLOR_deeppink`
```




<hr>



### define COLOR\_mfc\_twoPhaseWater 

```C++
#define COLOR_mfc_twoPhaseWater `COLOR_gold`
```




<hr>



### define COLOR\_navajowhite 

```C++
#define COLOR_navajowhite `{1.0, 0.8705882352941177, 0.6784313725490196}`
```




<hr>



### define COLOR\_orange 

```C++
#define COLOR_orange `{1.0, 0.6470588235294118, 0.0}`
```




<hr>



### define COLOR\_plum 

```C++
#define COLOR_plum `{0.866667, 0.627451, 0.866667}`
```




<hr>



### define COLOR\_red 

```C++
#define COLOR_red `{1.000000, 0.000000, 0.000000}`
```




<hr>



### define COLOR\_yellow 

```C++
#define COLOR_yellow `{1.000000, 1.000000, 0.000000}`
```




<hr>



### define SCALE\_X\_linear 

```C++
#define SCALE_X_linear `1`
```




<hr>



### define SCALE\_X\_log 

```C++
#define SCALE_X_log `2`
```




<hr>



### define SCALE\_X\_loglinear 

```C++
#define SCALE_X_loglinear `3`
```




<hr>



### define STR\_LENGTH\_PROPINFO 

```C++
#define STR_LENGTH_PROPINFO `30`
```




<hr>



### define Update\_prop\_CommonlyUsedHPX 

```C++
#define Update_prop_CommonlyUsedHPX `(Update_prop_H | Update_prop_T | Update_prop_Rho | Update_prop_Cp | Update_prop_Mu | Update_prop_IsothermalCompressibility |Update_prop_IsobaricExpansivity )`
```




<hr>



### define Update\_prop\_CommonlyUsedTPX 

```C++
#define Update_prop_CommonlyUsedTPX `(Update_prop_Rho | Update_prop_Cp | Update_prop_Mu | Update_prop_IsothermalCompressibility |Update_prop_IsobaricExpansivity )`
```




<hr>



### define Update\_prop\_Cp 

```C++
#define Update_prop_Cp `32`
```




<hr>



### define Update\_prop\_H 

```C++
#define Update_prop_H `4`
```




<hr>



### define Update\_prop\_IsobaricExpansivity 

```C++
#define Update_prop_IsobaricExpansivity `256`
```




<hr>



### define Update\_prop\_IsothermalCompressibility 

```C++
#define Update_prop_IsothermalCompressibility `128`
```




<hr>



### define Update\_prop\_Mu 

```C++
#define Update_prop_Mu `64`
```




<hr>



### define Update\_prop\_Rho 

```C++
#define Update_prop_Rho `2`
```




<hr>



### define Update\_prop\_T 

```C++
#define Update_prop_T `16`
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/thermo/DataStructures.h`

