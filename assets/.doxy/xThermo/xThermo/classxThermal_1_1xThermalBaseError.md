

# Class xThermal::xThermalBaseError



[**ClassList**](annotated.md) **>** [**xThermal**](namespacexThermal.md) **>** [**xThermalBaseError**](classxThermal_1_1xThermalBaseError.md)








Inherits the following classes: std::exception,  std::exception


Inherited by the following classes: [xThermal::xThermalError](classxThermal_1_1xThermalError.md),  [xThermal::xThermalError](classxThermal_1_1xThermalError.md)












## Public Types

| Type | Name |
| ---: | :--- |
| enum  | [**ErrCode**](#enum-errcode-12)  <br> |
| enum  | [**ErrCode**](#enum-errcode-12)  <br> |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  ErrCode | [**code**](#function-code-12) () <br> |
|  ErrCode | [**code**](#function-code-12) () <br> |
| virtual [**const**](classxThermal_1_1xThermalError.md) [**char**](classxThermal_1_1xThermalError.md) \* | [**what**](#function-what-12) () const<br> |
| virtual [**const**](classxThermal_1_1xThermalError.md) [**char**](classxThermal_1_1xThermalError.md) \* | [**what**](#function-what-12) () const<br> |
|   | [**xThermalBaseError**](#function-xthermalbaseerror-12) ([**const**](classxThermal_1_1xThermalError.md) std::string & err, ErrCode code) <br> |
|   | [**xThermalBaseError**](#function-xthermalbaseerror-12) ([**const**](classxThermal_1_1xThermalError.md) std::string & err, ErrCode code) <br> |
|   | [**~xThermalBaseError**](#function-xthermalbaseerror-12) () <br> |
|   | [**~xThermalBaseError**](#function-xthermalbaseerror-12) () <br> |




























## Public Types Documentation




### enum ErrCode [1/2]

```C++
enum xThermal::xThermalBaseError::ErrCode {
    eNotImplemented,
    eSolution,
    eAttribute,
    eOutOfRange,
    eValue,
    eWrongFluid,
    eComposition,
    eInput,
    eNotAvailable,
    eHandle,
    eKey,
    eUnableToLoad,
    eDirectorySize,
    eNotImplemented,
    eSolution,
    eAttribute,
    eOutOfRange,
    eValue,
    eWrongFluid,
    eComposition,
    eInput,
    eNotAvailable,
    eHandle,
    eKey,
    eUnableToLoad,
    eDirectorySize
};
```




<hr>



### enum ErrCode [1/2]

```C++
enum xThermal::xThermalBaseError::ErrCode {
    eNotImplemented,
    eSolution,
    eAttribute,
    eOutOfRange,
    eValue,
    eWrongFluid,
    eComposition,
    eInput,
    eNotAvailable,
    eHandle,
    eKey,
    eUnableToLoad,
    eDirectorySize,
    eNotImplemented,
    eSolution,
    eAttribute,
    eOutOfRange,
    eValue,
    eWrongFluid,
    eComposition,
    eInput,
    eNotAvailable,
    eHandle,
    eKey,
    eUnableToLoad,
    eDirectorySize
};
```




<hr>
## Public Functions Documentation




### function code [1/2]

```C++
inline ErrCode xThermal::xThermalBaseError::code () 
```




<hr>



### function code [1/2]

```C++
inline ErrCode xThermal::xThermalBaseError::code () 
```




<hr>



### function what [1/2]

```C++
inline virtual const  char * xThermal::xThermalBaseError::what () const
```




<hr>



### function what [1/2]

```C++
inline virtual const  char * xThermal::xThermalBaseError::what () const
```




<hr>



### function xThermalBaseError [1/2]

```C++
inline xThermal::xThermalBaseError::xThermalBaseError (
    const std::string & err,
    ErrCode code
) 
```




<hr>



### function xThermalBaseError [1/2]

```C++
inline xThermal::xThermalBaseError::xThermalBaseError (
    const std::string & err,
    ErrCode code
) 
```




<hr>



### function ~xThermalBaseError [1/2]

```C++
inline xThermal::xThermalBaseError::~xThermalBaseError () 
```




<hr>



### function ~xThermalBaseError [1/2]

```C++
inline xThermal::xThermalBaseError::~xThermalBaseError () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `Library/test/test_templateClass/Exception.h`

