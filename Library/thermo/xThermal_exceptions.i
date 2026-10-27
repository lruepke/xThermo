// Translate C++ exceptions thrown by the library into target-language exceptions (C3, PY-02).
// Without this, any C++ exception escaping a wrapped call terminates the host interpreter.
// Include it in a module's .i file after the %{ ... %} header block and before the wrapped declarations.
//
// Python mapping (xThermalBaseError::code()):
//   eOutOfRange, eValue, eWrongFluid, eComposition, eInput, eNotAvailable -> ValueError
//   eNotImplemented -> NotImplementedError, eAttribute -> AttributeError, eKey -> KeyError
//   eUnableToLoad -> OSError, eSolution and all others -> RuntimeError
//   std::out_of_range -> IndexError, std::invalid_argument -> ValueError, std::bad_alloc -> MemoryError,
//   other std::exception and unknown exceptions -> RuntimeError
// ANSI colour codes that the library adds to exception texts are removed from the message.

%include "exception.i"

%{
#include <new>
#include <stdexcept>
#include <string>
#include "Exception.h"

static std::string xThermal_strip_ansi(const char* s)
{
    std::string out;
    for (const char* c = s; *c; ++c) {
        if (*c == '\033' && *(c + 1) == '[') {
            c += 2;
            while (*c && !((*c >= '@' && *c <= '~'))) ++c;
            if (!*c) break;
            continue;
        }
        out += *c;
    }
    return out;
}
%}

#ifdef SWIGPYTHON
%{
static PyObject* xThermal_python_exception(xThermal::xThermalBaseError::ErrCode code)
{
    switch (code) {
        case xThermal::xThermalBaseError::eOutOfRange:
        case xThermal::xThermalBaseError::eValue:
        case xThermal::xThermalBaseError::eWrongFluid:
        case xThermal::xThermalBaseError::eComposition:
        case xThermal::xThermalBaseError::eInput:
        case xThermal::xThermalBaseError::eNotAvailable: return PyExc_ValueError;
        case xThermal::xThermalBaseError::eNotImplemented: return PyExc_NotImplementedError;
        case xThermal::xThermalBaseError::eAttribute: return PyExc_AttributeError;
        case xThermal::xThermalBaseError::eKey: return PyExc_KeyError;
        case xThermal::xThermalBaseError::eUnableToLoad: return PyExc_OSError;
        default: return PyExc_RuntimeError;
    }
}
%}
%exception {
    try {
        $action
    } catch (xThermal::xThermalBaseError& e) {
        PyErr_SetString(xThermal_python_exception(e.code()), xThermal_strip_ansi(e.what()).c_str());
        SWIG_fail;
    } catch (std::out_of_range& e) {
        PyErr_SetString(PyExc_IndexError, e.what());
        SWIG_fail;
    } catch (std::invalid_argument& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        SWIG_fail;
    } catch (std::bad_alloc& e) {
        PyErr_SetString(PyExc_MemoryError, e.what());
        SWIG_fail;
    } catch (std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, xThermal_strip_ansi(e.what()).c_str());
        SWIG_fail;
    } catch (...) {
        PyErr_SetString(PyExc_RuntimeError, "unknown C++ exception in xThermal");
        SWIG_fail;
    }
}
#else
%{
static int xThermal_swig_error_code(xThermal::xThermalBaseError::ErrCode code)
{
    switch (code) {
        case xThermal::xThermalBaseError::eOutOfRange:
        case xThermal::xThermalBaseError::eValue:
        case xThermal::xThermalBaseError::eWrongFluid:
        case xThermal::xThermalBaseError::eComposition:
        case xThermal::xThermalBaseError::eInput:
        case xThermal::xThermalBaseError::eNotAvailable: return SWIG_ValueError;
        case xThermal::xThermalBaseError::eAttribute: return SWIG_AttributeError;
        case xThermal::xThermalBaseError::eUnableToLoad: return SWIG_IOError;
        default: return SWIG_RuntimeError;
    }
}
%}
%exception {
    try {
        $action
    } catch (xThermal::xThermalBaseError& e) {
        SWIG_exception(xThermal_swig_error_code(e.code()), xThermal_strip_ansi(e.what()).c_str());
    } catch (std::out_of_range& e) {
        SWIG_exception(SWIG_IndexError, e.what());
    } catch (std::exception& e) {
        SWIG_exception(SWIG_RuntimeError, xThermal_strip_ansi(e.what()).c_str());
    } catch (...) {
        SWIG_exception(SWIG_UnknownError, "unknown C++ exception in xThermal");
    }
}
#endif

// C++-only API (stage 4): the non-throwing evaluation reports through a status and an output string; the Python
// batch error policy will be designed separately. The input check is internal.
%ignore xThermal::cxThermal::tryUpdateState_TPX;
%ignore xThermal::cxThermal::tryUpdateState_HPX;
%ignore xThermal::check_state_input;
