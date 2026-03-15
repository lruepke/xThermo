

# File LookUpTableForest.h



[**FileList**](files.md) **>** [**AMR\_LUT**](dir_eb711624b693511de215f59cb195e6f2.md) **>** [**LookUpTableForest.h**](LookUpTableForest_8h.md)

[Go to the source code of this file](LookUpTableForest_8h_source.md)



* `#include <vector>`
* `#include <map>`
* `#include <iostream>`
* `#include <fstream>`
* `#include "stdfunc.h"`
* `#include <cmath>`
* `#include "thermo.h"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**LOOKUPTABLE\_FOREST**](namespaceLOOKUPTABLE__FOREST.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| struct | [**FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md) &lt;dim&gt;<br> |
| struct | [**LeafQuad**](structLOOKUPTABLE__FOREST_1_1LeafQuad.md) &lt;dim, typename USER\_DATA&gt;<br> |
| class | [**LookUpTableForest**](classLOOKUPTABLE__FOREST_1_1LookUpTableForest.md) &lt;dim, typename USER\_DATA&gt;<br>_Pass dimension and data type to the class._  |
| struct | [**NonLeafQuad**](structLOOKUPTABLE__FOREST_1_1NonLeafQuad.md) &lt;dim, typename USER\_DATA&gt;<br> |
| struct | [**PropsData**](structLOOKUPTABLE__FOREST_1_1PropsData.md) <br> |
| struct | [**Quad\_index**](structLOOKUPTABLE__FOREST_1_1Quad__index.md) <br> |
| struct | [**Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md) &lt;dim, typename USER\_DATA&gt;<br> |
| struct | [**RMSD\_RefineCriterion**](structLOOKUPTABLE__FOREST_1_1RMSD__RefineCriterion.md) <br>_Property refinement criterion, minimum RMSD of a quadran, if the RMSD of a property in a quadran grater than this criterion, it will be refined._  |


## Public Types

| Type | Name |
| ---: | :--- |
| union  | [**Quad\_data**](#union-quad_data)  &lt;dim, typename USER\_DATA&gt;<br> |















































## Macros

| Type | Name |
| ---: | :--- |
| define  | [**CHECK\_REFINE\_PROP\_RMSD**](LookUpTableForest_8h.md#define-check_refine_prop_rmsd) (PROP) `/* multi line expression */`<br>_Check the criterion of a property and determine if need to refine. Should use relative error criterion, rather absolute error otherwise it is not fair for vapor region._  |
| define  | [**ExtensionName\_PointIndexFile**](LookUpTableForest_8h.md#define-extensionname_pointindexfile)  `"pi"`<br> |
| define  | [**MAX\_FOREST\_LEVEL**](LookUpTableForest_8h.md#define-max_forest_level)  `29`<br> |

## Public Types Documentation




### union Quad\_data 

```C++

```




<hr>
## Macro Definition Documentation





### define CHECK\_REFINE\_PROP\_RMSD 

_Check the criterion of a property and determine if need to refine. Should use relative error criterion, rather absolute error otherwise it is not fair for vapor region._ 
```C++
#define CHECK_REFINE_PROP_RMSD (
    PROP
) `/* multi line expression */`
```




<hr>



### define ExtensionName\_PointIndexFile 

```C++
#define ExtensionName_PointIndexFile `"pi"`
```




<hr>



### define MAX\_FOREST\_LEVEL 

```C++
#define MAX_FOREST_LEVEL `29`
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/AMR_LUT/LookUpTableForest.h`

