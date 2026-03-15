

# File Exception.h

[**File List**](files.md) **>** [**Library**](dir_5ad7f572bbca03234e8e621e192fc099.md) **>** [**test**](dir_fd7e605cf972a87e804527375e37bc17.md) **>** [**test\_templateClass**](dir_954facde36cda98a0027ddce7db85b3d.md) **>** [**Exception.h**](test_2test__templateClass_2Exception_8h.md)

[Go to the documentation of this file](test_2test__templateClass_2Exception_8h.md)


```C++


#ifndef CPEXCEPTIONS_H
#define CPEXCEPTIONS_H

#include <exception>
#include <iostream>
#include "stdfunc.h"
namespace xThermal
{

class xThermalBaseError: public std::exception
{
public:
    enum ErrCode { eNotImplemented, eSolution, eAttribute, eOutOfRange, eValue, eWrongFluid, eComposition, eInput, eNotAvailable, eHandle, eKey, eUnableToLoad,eDirectorySize};
    xThermalBaseError(const std::string &err, ErrCode code) throw() : m_code(code), m_err(COLOR_RED + err + COLOR_DEFAULT) {}
    ~xThermalBaseError() throw() {};
    virtual const char* what() const throw() { return m_err.c_str(); }
    ErrCode code() { return m_code; }
private:
    ErrCode m_code;
    std::string m_err;
};

template <xThermalBaseError::ErrCode errcode>
class xThermalError : public xThermalBaseError {
public:
    xThermalError(const std::string &err = "", ErrCode ecode = errcode) throw() : xThermalBaseError(err, ecode) {}
};

typedef xThermalError<xThermalBaseError::eNotImplemented> NotImplementedError;
typedef xThermalError<xThermalBaseError::eSolution> SolutionError;
typedef xThermalError<xThermalBaseError::eAttribute> AttributeError;
typedef xThermalError<xThermalBaseError::eOutOfRange> OutOfRangeError;
typedef xThermalError<xThermalBaseError::eValue> ValueError;
typedef xThermalError<xThermalBaseError::eKey> KeyError;
typedef xThermalError<xThermalBaseError::eHandle> HandleError;
typedef xThermalError<xThermalBaseError::eUnableToLoad> UnableToLoadError;
typedef xThermalError<xThermalBaseError::eDirectorySize> DirectorySizeError;

// ValueError specializations
template <xThermalBaseError::ErrCode errcode>
class ValueErrorSpec : public ValueError {
public:
    ValueErrorSpec(const std::string &err = "", ErrCode ecode = errcode) throw() : ValueError(err, ecode) {}
};

typedef ValueErrorSpec<xThermalBaseError::eWrongFluid> WrongFluidError;
typedef ValueErrorSpec<xThermalBaseError::eComposition> CompositionError;
typedef ValueErrorSpec<xThermalBaseError::eInput> InputError;
typedef ValueErrorSpec<xThermalBaseError::eNotAvailable> NotAvailableError;

}; /* namespace xThermal */
#endif
```


