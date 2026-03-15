

# Namespace LOOKUPTABLE\_FOREST



[**Namespace List**](namespaces.md) **>** [**LOOKUPTABLE\_FOREST**](namespaceLOOKUPTABLE__FOREST.md)




















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
| enum  | [**CONST\_WHICH\_VAR**](#enum-const_which_var)  <br>_For 2D case, define which variable is constant and the variable order of xy._  |
| enum  | [**EOS\_ENERGY**](#enum-eos_energy)  <br>_Use which variable to express energy._  |
| union  | [**Quad\_data**](#union-quad_data)  &lt;dim, typename USER\_DATA&gt;<br> |
| typedef [**LookUpTableForest**](classLOOKUPTABLE__FOREST_1_1LookUpTableForest.md)&lt; 2, [**FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md)&lt; 2 &gt; &gt; | [**LookUpTableForest\_2D**](#typedef-lookuptableforest_2d)  <br> |
| typedef [**LookUpTableForest**](classLOOKUPTABLE__FOREST_1_1LookUpTableForest.md)&lt; 3, [**FIELD\_DATA**](structLOOKUPTABLE__FOREST_1_1FIELD__DATA.md)&lt; 3 &gt; &gt; | [**LookUpTableForest\_3D**](#typedef-lookuptableforest_3d)  <br> |
| enum  | [**NeedRefine**](#enum-needrefine)  <br> |
| typedef unsigned int | [**int\_pointIndex**](#typedef-int_pointindex)  <br> |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  int | [**get\_dim\_from\_binary**](#function-get_dim_from_binary) (std::string filename) <br> |




























## Public Types Documentation




### enum CONST\_WHICH\_VAR 

_For 2D case, define which variable is constant and the variable order of xy._ 
```C++
enum LOOKUPTABLE_FOREST::CONST_WHICH_VAR {
    CONST_NO_VAR_TorHPX,
    CONST_TorH_VAR_XP,
    CONST_P_VAR_XTorH,
    CONST_X_VAR_TorHP
};
```




<hr>



### enum EOS\_ENERGY 

_Use which variable to express energy._ 
```C++
enum LOOKUPTABLE_FOREST::EOS_ENERGY {
    EOS_ENERGY_T,
    EOS_ENERGY_H
};
```




<hr>



### union Quad\_data 

```C++

```




<hr>



### typedef LookUpTableForest\_2D 

```C++
typedef LookUpTableForest<2, FIELD_DATA<2> > LOOKUPTABLE_FOREST::LookUpTableForest_2D;
```




<hr>



### typedef LookUpTableForest\_3D 

```C++
typedef LookUpTableForest<3, FIELD_DATA<3> > LOOKUPTABLE_FOREST::LookUpTableForest_3D;
```




<hr>



### enum NeedRefine 

```C++
enum LOOKUPTABLE_FOREST::NeedRefine {
    NeedRefine_NoNeed,
    NeedRefine_PhaseBoundary,
    NeedRefine_Rho,
    NeedRefine_H,
    NeedRefine_Mu
};
```




<hr>



### typedef int\_pointIndex 

```C++
typedef unsigned int LOOKUPTABLE_FOREST::int_pointIndex;
```




<hr>
## Public Functions Documentation




### function get\_dim\_from\_binary 

```C++
inline int LOOKUPTABLE_FOREST::get_dim_from_binary (
    std::string filename
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/AMR_LUT/LookUpTableForest.h`

