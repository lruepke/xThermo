

# Class LOOKUPTABLE\_FOREST::LookUpTableForest

**template &lt;int dim, typename USER\_DATA&gt;**



[**ClassList**](annotated.md) **>** [**LOOKUPTABLE\_FOREST**](namespaceLOOKUPTABLE__FOREST.md) **>** [**LookUpTableForest**](classLOOKUPTABLE__FOREST_1_1LookUpTableForest.md)



_Pass dimension and data type to the class._ [More...](#detailed-description)

* `#include <LookUpTableForest.h>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**RMSD\_RefineCriterion**](structLOOKUPTABLE__FOREST_1_1RMSD__RefineCriterion.md) | [**m\_RMSD\_RefineCriterion**](#variable-m_rmsd_refinecriterion)  <br> |
|  EOS\_ENERGY | [**m\_TorH**](#variable-m_torh)  <br> |
|  double | [**m\_constZ**](#variable-m_constz)  <br> |
|  CONST\_WHICH\_VAR | [**m\_const\_which\_var**](#variable-m_const_which_var)  <br> |
|  void \* | [**m\_eosPointer**](#variable-m_eospointer)  <br> |
|  std::map&lt; int, int &gt; | [**m\_map\_prop2index**](#variable-m_map_prop2index)  <br> |
|  std::map&lt; int, [**xThermal::propInfo**](structxThermal_1_1propInfo.md) &gt; | [**m\_map\_props**](#variable-m_map_props)  <br> |
|  int | [**m\_max\_level**](#variable-m_max_level)  <br> |
|  int | [**m\_min\_level**](#variable-m_min_level)  <br> |
|  int | [**m\_num\_children**](#variable-m_num_children)  <br> |
|  int | [**m\_num\_node\_per\_quad**](#variable-m_num_node_per_quad)  <br> |
|  [**PropsData**](structLOOKUPTABLE__FOREST_1_1PropsData.md) | [**m\_props\_unique\_points\_leaves**](#variable-m_props_unique_points_leaves)  <br> |
|  double | [**m\_xyz\_max**](#variable-m_xyz_max)  <br> |
|  double | [**m\_xyz\_min**](#variable-m_xyz_min)  <br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**LookUpTableForest**](#function-lookuptableforest-13) (double xyz\_min, double xyz\_max, EOS\_ENERGY TorH, int max\_level, std::map&lt; int, [**xThermal::propInfo**](structxThermal_1_1propInfo.md) &gt; name\_props, void \* eosPointer=NULL) <br>_Construct a new Look Up Table Forest object. This is always used to create a 3D table xyz would be corresponding to TPX or PHX. Note that the unit of T is K, unit of P is Pa, unit of X is wt% NaCl (e.g., seawater is 0.032), unit of H is J/kg. The same as H2ONaCl::cH2ONaCl::prop\_pTX and The same as H2ONaCl::cH2ONaCl::prop\_pHX. For 3D case, the order of the variable MUST BE [T/H, p, X], T or H is specify by argument eos\_space._  |
|   | [**LookUpTableForest**](#function-lookuptableforest-23) (double xy\_min, double xy\_max, double constZ, CONST\_WHICH\_VAR const\_which\_var, EOS\_ENERGY TorH, int max\_level, std::map&lt; int, [**xThermal::propInfo**](structxThermal_1_1propInfo.md) &gt; name\_props, void \* eosPointer=NULL) <br>_Construct a new Look Up Table Forest object. This is always used to create a 2D table xyz would be corresponding to TPX or PHX. Note that the unit of T is K, unit of P is Pa, unit of X is wt% NaCl (e.g., seawater is 0.032), unit of H is J/kg. The same as H2ONaCl::cH2ONaCl::prop\_pTX and The same as H2ONaCl::cH2ONaCl::prop\_pHX. For 2D case, (1) if constant variable is X, xy order MUST BE T/H, P; (2) if constant variable is T/H, xy order MUST BE X, P; (3) if constant variable is P, xy order MUST BE X, T/H._  |
|   | [**LookUpTableForest**](#function-lookuptableforest-33) (std::string filename\_forest, void \* pointer=NULL, bool printStatus=true) <br> |
|  void | [**assemble\_data**](#function-assemble_data) (void(\*)([**LookUpTableForest**](classLOOKUPTABLE__FOREST_1_1LookUpTableForest.md)&lt; dim, USER\_DATA &gt; \*forest, std::map&lt; [**Quad\_index**](structLOOKUPTABLE__FOREST_1_1Quad__index.md), double \* &gt; &map\_ijk2data) cal\_prop) <br> |
|  std::string | [**byte2string**](#function-byte2string) (double bytes) <br> |
|  void | [**construct\_props\_leaves**](#function-construct_props_leaves) (void(\*)([**LookUpTableForest**](classLOOKUPTABLE__FOREST_1_1LookUpTableForest.md)&lt; dim, USER\_DATA &gt; \*forest, std::map&lt; [**Quad\_index**](structLOOKUPTABLE__FOREST_1_1Quad__index.md), unsigned int &gt; &map\_ijk2data, double \*\*data) cal\_prop) <br> |
|  void | [**destory**](#function-destory) () <br> |
|  void | [**get\_ijk\_nodes\_quadrant**](#function-get_ijk_nodes_quadrant) ([**Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; dim, USER\_DATA &gt; \* quad, const [**Quad\_index**](structLOOKUPTABLE__FOREST_1_1Quad__index.md) \* ijk\_quad, int num\_nodes\_per\_quad, [**Quad\_index**](structLOOKUPTABLE__FOREST_1_1Quad__index.md) \* ijk) <br> |
|  int | [**get\_num\_leaves**](#function-get_num_leaves) () <br> |
|  int | [**get\_num\_need\_refine**](#function-get_num_need_refine) () <br> |
|  int | [**get\_num\_quads**](#function-get_num_quads) () <br> |
|  void | [**get\_quadrant\_physical\_length**](#function-get_quadrant_physical_length) (int level, double physical\_length) <br> |
|  [**Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; dim, USER\_DATA &gt; \* | [**get\_root**](#function-get_root) () <br> |
|  void | [**ijk2xyz**](#function-ijk2xyz) (const [**Quad\_index**](structLOOKUPTABLE__FOREST_1_1Quad__index.md) \* ijk, double & x, double & y, double & z) <br> |
|  void | [**print\_summary**](#function-print_summary) () <br> |
|  void | [**refine**](#function-refine-22) (bool(\*)([**LookUpTableForest**](classLOOKUPTABLE__FOREST_1_1LookUpTableForest.md)&lt; dim, USER\_DATA &gt; \*forest, [**Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; dim, USER\_DATA &gt; \*quad, double xmin\_quad, double ymin\_quad, double zmin\_quad, int max\_level) is\_refine) <br> |
|  void | [**searchQuadrant**](#function-searchquadrant-22) ([**Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; dim, USER\_DATA &gt; \*& targetLeaf, double \* xyz\_min\_target, double x, double y, double z) <br> |
|  void | [**set\_min\_level**](#function-set_min_level) (int min\_level) <br> |
|  void | [**union\_ijk2xyz**](#function-union_ijk2xyz) ([**Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; dim, USER\_DATA &gt; \* quad, [**Quad\_index**](structLOOKUPTABLE__FOREST_1_1Quad__index.md) & ijk\_backup) <br> |
|  void | [**write\_point\_index**](#function-write_point_index-12) (std::string filename\_forest) <br> |
|  void | [**write\_point\_index**](#function-write_point_index-22) (FILE \* fpout\_point\_index, [**Quadrant**](structLOOKUPTABLE__FOREST_1_1Quadrant.md)&lt; dim, USER\_DATA &gt; \* quad) <br> |
|  void | [**write\_to\_binary**](#function-write_to_binary) (std::string filename, bool is\_write\_data=true) <br> |
|  void | [**write\_to\_vtk**](#function-write_to_vtk) (std::string filename, bool write\_data=true, bool isNormalizeXYZ=true) <br> |
|   | [**~LookUpTableForest**](#function-lookuptableforest) () <br> |




























## Detailed Description




**Template parameters:**


* `dim` 
* `USER_DATA` 




    
## Public Attributes Documentation




### variable m\_RMSD\_RefineCriterion 

```C++
RMSD_RefineCriterion LOOKUPTABLE_FOREST::LookUpTableForest< dim, USER_DATA >::m_RMSD_RefineCriterion;
```




<hr>



### variable m\_TorH 

```C++
EOS_ENERGY LOOKUPTABLE_FOREST::LookUpTableForest< dim, USER_DATA >::m_TorH;
```




<hr>



### variable m\_constZ 

```C++
double LOOKUPTABLE_FOREST::LookUpTableForest< dim, USER_DATA >::m_constZ;
```




<hr>



### variable m\_const\_which\_var 

```C++
CONST_WHICH_VAR LOOKUPTABLE_FOREST::LookUpTableForest< dim, USER_DATA >::m_const_which_var;
```




<hr>



### variable m\_eosPointer 

```C++
void* LOOKUPTABLE_FOREST::LookUpTableForest< dim, USER_DATA >::m_eosPointer;
```




<hr>



### variable m\_map\_prop2index 

```C++
std::map<int, int> LOOKUPTABLE_FOREST::LookUpTableForest< dim, USER_DATA >::m_map_prop2index;
```




<hr>



### variable m\_map\_props 

```C++
std::map<int, xThermal::propInfo> LOOKUPTABLE_FOREST::LookUpTableForest< dim, USER_DATA >::m_map_props;
```




<hr>



### variable m\_max\_level 

```C++
int LOOKUPTABLE_FOREST::LookUpTableForest< dim, USER_DATA >::m_max_level;
```




<hr>



### variable m\_min\_level 

```C++
int LOOKUPTABLE_FOREST::LookUpTableForest< dim, USER_DATA >::m_min_level;
```




<hr>



### variable m\_num\_children 

```C++
int LOOKUPTABLE_FOREST::LookUpTableForest< dim, USER_DATA >::m_num_children;
```




<hr>



### variable m\_num\_node\_per\_quad 

```C++
int LOOKUPTABLE_FOREST::LookUpTableForest< dim, USER_DATA >::m_num_node_per_quad;
```




<hr>



### variable m\_props\_unique\_points\_leaves 

```C++
PropsData LOOKUPTABLE_FOREST::LookUpTableForest< dim, USER_DATA >::m_props_unique_points_leaves;
```




<hr>



### variable m\_xyz\_max 

```C++
double LOOKUPTABLE_FOREST::LookUpTableForest< dim, USER_DATA >::m_xyz_max[dim];
```




<hr>



### variable m\_xyz\_min 

```C++
double LOOKUPTABLE_FOREST::LookUpTableForest< dim, USER_DATA >::m_xyz_min[dim];
```




<hr>
## Public Functions Documentation




### function LookUpTableForest [1/3]

_Construct a new Look Up Table Forest object. This is always used to create a 3D table xyz would be corresponding to TPX or PHX. Note that the unit of T is K, unit of P is Pa, unit of X is wt% NaCl (e.g., seawater is 0.032), unit of H is J/kg. The same as H2ONaCl::cH2ONaCl::prop\_pTX and The same as H2ONaCl::cH2ONaCl::prop\_pHX. For 3D case, the order of the variable MUST BE [T/H, p, X], T or H is specify by argument eos\_space._ 
```C++
LOOKUPTABLE_FOREST::LookUpTableForest::LookUpTableForest (
    double xyz_min,
    double xyz_max,
    EOS_ENERGY TorH,
    int max_level,
    std::map< int, xThermal::propInfo > name_props,
    void * eosPointer=NULL
) 
```





**Parameters:**


* `xyz_min` 
* `xyz_max` 
* `max_level` 
* `eosPointer` 




        

<hr>



### function LookUpTableForest [2/3]

_Construct a new Look Up Table Forest object. This is always used to create a 2D table xyz would be corresponding to TPX or PHX. Note that the unit of T is K, unit of P is Pa, unit of X is wt% NaCl (e.g., seawater is 0.032), unit of H is J/kg. The same as H2ONaCl::cH2ONaCl::prop\_pTX and The same as H2ONaCl::cH2ONaCl::prop\_pHX. For 2D case, (1) if constant variable is X, xy order MUST BE T/H, P; (2) if constant variable is T/H, xy order MUST BE X, P; (3) if constant variable is P, xy order MUST BE X, T/H._ 
```C++
LOOKUPTABLE_FOREST::LookUpTableForest::LookUpTableForest (
    double xy_min,
    double xy_max,
    double constZ,
    CONST_WHICH_VAR const_which_var,
    EOS_ENERGY TorH,
    int max_level,
    std::map< int, xThermal::propInfo > name_props,
    void * eosPointer=NULL
) 
```





**Parameters:**


* `xy_min` 
* `xy_max` 
* `constZ` 
* `const_which_var` 
* `max_level` 
* `eosPointer` 




        

<hr>



### function LookUpTableForest [3/3]

```C++
LOOKUPTABLE_FOREST::LookUpTableForest::LookUpTableForest (
    std::string filename_forest,
    void * pointer=NULL,
    bool printStatus=true
) 
```




<hr>



### function assemble\_data 

```C++
void LOOKUPTABLE_FOREST::LookUpTableForest::assemble_data (
    void(*)( LookUpTableForest < dim, USER_DATA > *forest, std::map< Quad_index , double * > &map_ijk2data) cal_prop
) 
```




<hr>



### function byte2string 

```C++
std::string LOOKUPTABLE_FOREST::LookUpTableForest::byte2string (
    double bytes
) 
```




<hr>



### function construct\_props\_leaves 

```C++
void LOOKUPTABLE_FOREST::LookUpTableForest::construct_props_leaves (
    void(*)( LookUpTableForest < dim, USER_DATA > *forest, std::map< Quad_index , unsigned int > &map_ijk2data, double **data) cal_prop
) 
```




<hr>



### function destory 

```C++
void LOOKUPTABLE_FOREST::LookUpTableForest::destory () 
```




<hr>



### function get\_ijk\_nodes\_quadrant 

```C++
void LOOKUPTABLE_FOREST::LookUpTableForest::get_ijk_nodes_quadrant (
    Quadrant < dim, USER_DATA > * quad,
    const Quad_index * ijk_quad,
    int num_nodes_per_quad,
    Quad_index * ijk
) 
```




<hr>



### function get\_num\_leaves 

```C++
inline int LOOKUPTABLE_FOREST::LookUpTableForest::get_num_leaves () 
```




<hr>



### function get\_num\_need\_refine 

```C++
inline int LOOKUPTABLE_FOREST::LookUpTableForest::get_num_need_refine () 
```




<hr>



### function get\_num\_quads 

```C++
inline int LOOKUPTABLE_FOREST::LookUpTableForest::get_num_quads () 
```




<hr>



### function get\_quadrant\_physical\_length 

```C++
void LOOKUPTABLE_FOREST::LookUpTableForest::get_quadrant_physical_length (
    int level,
    double physical_length
) 
```




<hr>



### function get\_root 

```C++
inline Quadrant < dim, USER_DATA > * LOOKUPTABLE_FOREST::LookUpTableForest::get_root () 
```




<hr>



### function ijk2xyz 

```C++
void LOOKUPTABLE_FOREST::LookUpTableForest::ijk2xyz (
    const Quad_index * ijk,
    double & x,
    double & y,
    double & z
) 
```




<hr>



### function print\_summary 

```C++
void LOOKUPTABLE_FOREST::LookUpTableForest::print_summary () 
```




<hr>



### function refine [2/2]

```C++
void LOOKUPTABLE_FOREST::LookUpTableForest::refine (
    bool(*)( LookUpTableForest < dim, USER_DATA > *forest, Quadrant < dim, USER_DATA > *quad, double xmin_quad, double ymin_quad, double zmin_quad, int max_level) is_refine
) 
```




<hr>



### function searchQuadrant [2/2]

```C++
void LOOKUPTABLE_FOREST::LookUpTableForest::searchQuadrant (
    Quadrant < dim, USER_DATA > *& targetLeaf,
    double * xyz_min_target,
    double x,
    double y,
    double z
) 
```




<hr>



### function set\_min\_level 

```C++
inline void LOOKUPTABLE_FOREST::LookUpTableForest::set_min_level (
    int min_level
) 
```




<hr>



### function union\_ijk2xyz 

```C++
void LOOKUPTABLE_FOREST::LookUpTableForest::union_ijk2xyz (
    Quadrant < dim, USER_DATA > * quad,
    Quad_index & ijk_backup
) 
```




<hr>



### function write\_point\_index [1/2]

```C++
void LOOKUPTABLE_FOREST::LookUpTableForest::write_point_index (
    std::string filename_forest
) 
```




<hr>



### function write\_point\_index [2/2]

```C++
void LOOKUPTABLE_FOREST::LookUpTableForest::write_point_index (
    FILE * fpout_point_index,
    Quadrant < dim, USER_DATA > * quad
) 
```




<hr>



### function write\_to\_binary 

```C++
void LOOKUPTABLE_FOREST::LookUpTableForest::write_to_binary (
    std::string filename,
    bool is_write_data=true
) 
```




<hr>



### function write\_to\_vtk 

```C++
void LOOKUPTABLE_FOREST::LookUpTableForest::write_to_vtk (
    std::string filename,
    bool write_data=true,
    bool isNormalizeXYZ=true
) 
```




<hr>



### function ~LookUpTableForest 

```C++
LOOKUPTABLE_FOREST::LookUpTableForest::~LookUpTableForest () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/AMR_LUT/LookUpTableForest.h`

