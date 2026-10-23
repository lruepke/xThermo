

# Namespace xThermal::gsl\_guard



[**Namespace List**](namespaces.md) **>** [**xThermal**](namespacexThermal.md) **>** [**gsl\_guard**](namespacexThermal_1_1gsl__guard.md)




















## Classes

| Type | Name |
| ---: | :--- |
| class | [**GuardedFunction**](classxThermal_1_1gsl__guard_1_1GuardedFunction.md) <br> |
| class | [**GuardedMultirootFunction**](classxThermal_1_1gsl__guard_1_1GuardedMultirootFunction.md) <br> |
| struct | [**MultirootSolverFree**](structxThermal_1_1gsl__guard_1_1MultirootSolverFree.md) <br> |
| struct | [**RootSolverFree**](structxThermal_1_1gsl__guard_1_1RootSolverFree.md) <br> |
| struct | [**VectorFree**](structxThermal_1_1gsl__guard_1_1VectorFree.md) <br> |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef std::unique\_ptr&lt; [**gsl\_multiroot\_fsolver**](classxThermal_1_1xThermalError.md), [**MultirootSolverFree**](structxThermal_1_1gsl__guard_1_1MultirootSolverFree.md) &gt; | [**MultirootSolver**](#typedef-multirootsolver)  <br> |
| typedef std::unique\_ptr&lt; [**gsl\_root\_fsolver**](classxThermal_1_1xThermalError.md), [**RootSolverFree**](structxThermal_1_1gsl__guard_1_1RootSolverFree.md) &gt; | [**RootSolver**](#typedef-rootsolver)  <br> |
| typedef std::unique\_ptr&lt; [**gsl\_vector**](classxThermal_1_1xThermalError.md), [**VectorFree**](structxThermal_1_1gsl__guard_1_1VectorFree.md) &gt; | [**Vector**](#typedef-vector)  <br> |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  std::string | [**num**](#function-num) ([**double**](classxThermal_1_1xThermalError.md) x) <br> |




























## Public Types Documentation




### typedef MultirootSolver 

```C++
typedef std::unique_ptr<gsl_multiroot_fsolver, MultirootSolverFree> xThermal::gsl_guard::MultirootSolver;
```




<hr>



### typedef RootSolver 

```C++
typedef std::unique_ptr<gsl_root_fsolver, RootSolverFree> xThermal::gsl_guard::RootSolver;
```




<hr>



### typedef Vector 

```C++
typedef std::unique_ptr<gsl_vector, VectorFree> xThermal::gsl_guard::Vector;
```




<hr>
## Public Functions Documentation




### function num 

```C++
inline std::string xThermal::gsl_guard::num (
    double x
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/thermo/gsl_guard.hpp`

