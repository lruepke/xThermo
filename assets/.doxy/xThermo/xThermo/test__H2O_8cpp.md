

# File test\_H2O.cpp



[**FileList**](files.md) **>** [**benchmark\_test**](dir_c3298c2f310225057389213c5b89b78a.md) **>** [**test\_H2O.cpp**](test__H2O_8cpp.md)

[Go to the source code of this file](test__H2O_8cpp_source.md)



* `#include "H2O.h"`
* `#include "IAPWS-95.h"`
* `#include "IAPWS-IF97.h"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**std**](namespacestd.md) <br> |
























## Public Functions

| Type | Name |
| ---: | :--- |
|  int | [**main**](#function-main) (int argc, char \*\* argv) <br> |
|  void | [**test\_Boiling\_P**](#function-test_boiling_p) (double T\_K) <br> |
|  void | [**test\_Boiling\_T**](#function-test_boiling_t) (double P) <br> |
|  void | [**test\_H**](#function-test_h) (double T\_K, double P) <br> |
|  double | [**test\_TP2Rho**](#function-test_tp2rho) (double T\_K, double P\_Pa) <br> |
|  void | [**test\_state\_HP**](#function-test_state_hp) (double H, double P) <br> |
|  void | [**test\_state\_HP\_2D**](#function-test_state_hp_2d) () <br> |
|  void | [**validate\_SinglePhaseRegion\_Table7**](#function-validate_singlephaseregion_table7) () <br>_Compare and validate result given in Table 7 of IAPWS-95._  |
|  void | [**validate\_TwoPhaseRegion\_Table8**](#function-validate_twophaseregion_table8) () <br>_Compare and validate result given in Table 8 of IAPWS-95._  |
|  void | [**validate\_phio\_phir\_Table6**](#function-validate_phio_phir_table6) () <br>_Validate ideal-gas part_ \(\phi^o\) _and residual part_\(\phi^r\) _of the dimensionless Helmholtz energy and its derivatives for T=500 K and_\(\rho = 838.025 kg/m^3\) _. See Table 6 of IAPWS-95._ |




























## Public Functions Documentation




### function main 

```C++
int main (
    int argc,
    char ** argv
) 
```




<hr>



### function test\_Boiling\_P 

```C++
void test_Boiling_P (
    double T_K
) 
```




<hr>



### function test\_Boiling\_T 

```C++
void test_Boiling_T (
    double P
) 
```




<hr>



### function test\_H 

```C++
void test_H (
    double T_K,
    double P
) 
```




<hr>



### function test\_TP2Rho 

```C++
double test_TP2Rho (
    double T_K,
    double P_Pa
) 
```




<hr>



### function test\_state\_HP 

```C++
void test_state_HP (
    double H,
    double P
) 
```




<hr>



### function test\_state\_HP\_2D 

```C++
void test_state_HP_2D () 
```




<hr>



### function validate\_SinglePhaseRegion\_Table7 

_Compare and validate result given in Table 7 of IAPWS-95._ 
```C++
void validate_SinglePhaseRegion_Table7 () 
```




<hr>



### function validate\_TwoPhaseRegion\_Table8 

_Compare and validate result given in Table 8 of IAPWS-95._ 
```C++
void validate_TwoPhaseRegion_Table8 () 
```




<hr>



### function validate\_phio\_phir\_Table6 

_Validate ideal-gas part_ \(\phi^o\) _and residual part_\(\phi^r\) _of the dimensionless Helmholtz energy and its derivatives for T=500 K and_\(\rho = 838.025 kg/m^3\) _. See Table 6 of IAPWS-95._
```C++
void validate_phio_phir_Table6 () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/benchmark_test/test_H2O.cpp`

