

#ifndef CPEXCEPTIONS_xThermal_H
#define CPEXCEPTIONS_xThermal_H

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
    ErrCode code() const { return m_code; }
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

/**
 * Outcome of a non-throwing state evaluation (cxThermal::tryUpdateState_TPX, cxThermal::tryUpdateState_HPX).
 * Each failure corresponds to the xThermal exception that the throwing call (UpdateState_TPX / UpdateState_HPX)
 * raises for the same input.
 */
enum class EvaluationStatus
{
    Ok = 0,        /**< the state was evaluated */
    InvalidInput,  /**< non-finite T, H, p or X, X < 0, or p <= 0 (InputError) */
    OutOfRange,    /**< outside the range in which the model evaluates (OutOfRangeError) */
    NotConverged,  /**< a numerical solution failed (SolutionError) */
    NotAvailable,  /**< not implemented or not available for this model or state (NotImplementedError, NotAvailableError) */
    Failed         /**< any other xThermal error, e.g. ValueError */
};

/** Status corresponding to an xThermal error code */
inline EvaluationStatus evaluation_status(xThermalBaseError::ErrCode code)
{
    switch (code)
    {
        case xThermalBaseError::eInput: return EvaluationStatus::InvalidInput;
        case xThermalBaseError::eOutOfRange: return EvaluationStatus::OutOfRange;
        case xThermalBaseError::eSolution: return EvaluationStatus::NotConverged;
        case xThermalBaseError::eNotImplemented:
        case xThermalBaseError::eNotAvailable: return EvaluationStatus::NotAvailable;
        default: return EvaluationStatus::Failed;
    }
}

/** Name of a status, e.g. "OutOfRange" (for messages and logs) */
inline const char* evaluation_status_name(EvaluationStatus status)
{
    switch (status)
    {
        case EvaluationStatus::Ok: return "Ok";
        case EvaluationStatus::InvalidInput: return "InvalidInput";
        case EvaluationStatus::OutOfRange: return "OutOfRange";
        case EvaluationStatus::NotConverged: return "NotConverged";
        case EvaluationStatus::NotAvailable: return "NotAvailable";
        default: return "Failed";
    }
}

}; /* namespace xThermal */
#endif
