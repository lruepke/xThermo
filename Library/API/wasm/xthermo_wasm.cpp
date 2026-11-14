/**
 * @file xthermo_wasm.cpp
 * @brief WebAssembly binding of the H2O-NaCl model (xThermo Explorer, W2a; review/W1_WASM_FEASIBILITY.md).
 *
 * A flat C ABI over the public C++ API, nothing else: cH2ONaCl with the IAPS84 or IAPWS95 water backend and the batch
 * evaluators cxThermal::evaluate_TPX / evaluate_HPX, which are also behind the Python functions xthermo.tpx and
 * xthermo.hpx. A point that cannot be evaluated has status != 0 (EvaluationStatus), phase Unknown, NaN fields (the
 * inputs kept) and the library's message; no other point is affected. No numerical code, approximation or fallback
 * lives here. Use from JavaScript through xthermo.js.
 *
 * Results are written field-major: out[k * n + i] is field k (xt_field_name(k), SI units) of point i.
 */
#include <cstddef>
#include <exception>
#include <string>
#include <vector>
#include <emscripten/emscripten.h>
#include <emscripten/version.h>
#include <gsl/gsl_version.h>
#include "H2ONaCl.h"

#ifndef XTHERMO_WASM_VERSION
#define XTHERMO_WASM_VERSION "unknown"
#endif
#ifndef XTHERMO_WASM_COMMIT
#define XTHERMO_WASM_COMMIT "unknown"
#endif
#define XT_STR2(x) #x
#define XT_STR(x) XT_STR2(x)

using namespace xThermal;

namespace
{
    const int ABI_VERSION = 1;
    // name and member of every exported field, in output order
    struct Field { const char* name; std::vector<double> ThermodynamicPropertiesVector::* member; };
    const Field FIELDS[] = {
        {"T", &ThermodynamicPropertiesVector::T}, {"p", &ThermodynamicPropertiesVector::p},
        {"X", &ThermodynamicPropertiesVector::X}, {"H", &ThermodynamicPropertiesVector::H},
        {"S_l", &ThermodynamicPropertiesVector::S_l}, {"S_v", &ThermodynamicPropertiesVector::S_v},
        {"S_h", &ThermodynamicPropertiesVector::S_h},
        {"X_l", &ThermodynamicPropertiesVector::X_l}, {"X_v", &ThermodynamicPropertiesVector::X_v},
        {"Rho", &ThermodynamicPropertiesVector::Rho}, {"Rho_l", &ThermodynamicPropertiesVector::Rho_l},
        {"Rho_v", &ThermodynamicPropertiesVector::Rho_v}, {"Rho_h", &ThermodynamicPropertiesVector::Rho_h},
        {"H_l", &ThermodynamicPropertiesVector::H_l}, {"H_v", &ThermodynamicPropertiesVector::H_v},
        {"H_h", &ThermodynamicPropertiesVector::H_h},
        {"Cp", &ThermodynamicPropertiesVector::Cp}, {"Cp_l", &ThermodynamicPropertiesVector::Cp_l},
        {"Cp_v", &ThermodynamicPropertiesVector::Cp_v}, {"Cp_h", &ThermodynamicPropertiesVector::Cp_h},
        {"Mu", &ThermodynamicPropertiesVector::Mu}, {"Mu_l", &ThermodynamicPropertiesVector::Mu_l},
        {"Mu_v", &ThermodynamicPropertiesVector::Mu_v},
        {"IsothermalCompressibility", &ThermodynamicPropertiesVector::IsothermalCompressibility},
        {"IsobaricExpansivity", &ThermodynamicPropertiesVector::IsobaricExpansivity},
    };
    const int NFIELDS = sizeof(FIELDS) / sizeof(FIELDS[0]);
    const int NSTATUS = (int)EvaluationStatus::Failed + 1;

    std::vector<std::string> g_messages;   // per point of the last xt_eval ("" for evaluated points)
    std::string g_error;                    // last error that is not a per-point evaluation failure
    std::string g_text;                     // storage for returned strings
}

extern "C" {

/** Version of the layout below (fields, codes, signatures); xthermo.js refuses a different one. */
EMSCRIPTEN_KEEPALIVE int xt_abi_version() { return ABI_VERSION; }

/** JSON: xThermo version and commit, GSL and Emscripten versions, ABI version. */
EMSCRIPTEN_KEEPALIVE const char* xt_build_info()
{
    g_text = std::string("{\"xthermo\":\"") + XTHERMO_WASM_VERSION + "\",\"commit\":\"" + XTHERMO_WASM_COMMIT +
             "\",\"gsl\":\"" + GSL_VERSION + "\",\"emscripten\":\"" XT_STR(__EMSCRIPTEN_major__) "." XT_STR(__EMSCRIPTEN_minor__)
             "." XT_STR(__EMSCRIPTEN_tiny__) "\",\"abi\":" + std::to_string(ABI_VERSION) + "}";
    return g_text.c_str();
}

/** New H2O-NaCl model. backend 0: IAPS84, 1: IAPWS95. NULL on failure (xt_last_error). */
EMSCRIPTEN_KEEPALIVE void* xt_new(int backend)
{
    g_error.clear();
    if (backend != 0 && backend != 1) { g_error = "backend must be 0 (IAPS84) or 1 (IAPWS95)"; return nullptr; }
    try {
        return new H2ONaCl::cH2ONaCl(backend == 0 ? "IAPS84" : "IAPWS95");
    } catch (const std::exception& e) {
        g_error = e.what();
    } catch (...) {
        g_error = "unknown error";
    }
    return nullptr;
}

EMSCRIPTEN_KEEPALIVE void xt_free(void* model) { delete static_cast<H2ONaCl::cH2ONaCl*>(model); }

EMSCRIPTEN_KEEPALIVE const char* xt_last_error() { return g_error.c_str(); }

EMSCRIPTEN_KEEPALIVE int xt_nfields() { return NFIELDS; }

EMSCRIPTEN_KEEPALIVE const char* xt_field_name(int k) { return k >= 0 && k < NFIELDS ? FIELDS[k].name : ""; }

/**
 * Evaluate n points. mode 0: T-P-X (a = T [K]); mode 1: H-P-X (a = H [J/kg]); p [Pa]; X [kg/kg].
 * out: n * xt_nfields() doubles (field-major); status, phase: n ints each.
 * Returns the number of failed points, or -1 if the call itself failed (xt_last_error; outputs unchanged).
 */
EMSCRIPTEN_KEEPALIVE int xt_eval(void* model, int mode, int n, const double* a, const double* p, const double* X,
                                 double* out, int* status, int* phase)
{
    g_error.clear();
    g_messages.clear();
    if (!model || (mode != 0 && mode != 1) || n < 0) { g_error = "invalid arguments"; return -1; }
    try {
        const std::vector<double> va(a, a + n), vp(p, p + n), vX(X, X + n);
        auto* sw = static_cast<H2ONaCl::cH2ONaCl*>(model);
        const ThermodynamicPropertiesVector r = mode == 0 ? sw->evaluate_TPX(va, vp, vX) : sw->evaluate_HPX(va, vp, vX);
        g_messages.assign(n, std::string());
        int nfail = 0;
        for (int i = 0; i < n; ++i)
        {
            status[i] = r.status[i];
            phase[i] = r.phase[i];
            if (r.status[i] != 0) g_messages[i] = r.failed_message[nfail++];
        }
        for (int k = 0; k < NFIELDS; ++k)
        {
            const std::vector<double>& f = r.*(FIELDS[k].member);
            for (int i = 0; i < n; ++i) out[(size_t)k * n + i] = f[i];
        }
        return nfail;
    } catch (const std::exception& e) {
        g_error = e.what();
    } catch (...) {
        g_error = "unknown error";
    }
    g_messages.clear();
    return -1;
}

/** Library message of point i of the last xt_eval ("" if evaluated). */
EMSCRIPTEN_KEEPALIVE const char* xt_message(int i)
{
    return i >= 0 && i < (int)g_messages.size() ? g_messages[i].c_str() : "";
}

/** Name of a PhaseRegion code (the library's own map, as model.phase_name in Python). */
EMSCRIPTEN_KEEPALIVE const char* xt_phase_name(int code)
{
    if (code < MixPhaseRegion || code > NotImposed) return "";
    static H2ONaCl::cH2ONaCl* names = nullptr;   // phase_name does not depend on the model; one instance for the map
    try {
        if (!names) names = new H2ONaCl::cH2ONaCl("IAPS84");
        g_text = names->phase_name(static_cast<PhaseRegion>(code));
    } catch (...) {
        g_text.clear();
    }
    return g_text.c_str();
}

/** Name of an EvaluationStatus code (as xthermo.Status in Python). */
EMSCRIPTEN_KEEPALIVE const char* xt_status_name(int status)
{
    return status >= 0 && status < NSTATUS ? evaluation_status_name(static_cast<EvaluationStatus>(status)) : "";
}

}   // extern "C"
