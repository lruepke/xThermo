/**
 * @file gsl_guard.hpp
 * @brief Internal (not installed): RAII ownership and exception-safe callbacks for the GSL solvers (C1, LIB-02/18, LIB-27).
 *
 * GSL's process-global error handler is never read or written. Instead, the conditions under which GSL would call it
 * (and, with its default handler, abort the host process) are prevented here:
 *  - Brent bracket: lower > upper (or NaN), non-finite endpoint values and endpoints that do not straddle zero are
 *    detected before gsl_root_fsolver_set, with exactly GSL's criteria, and reported as SolutionError;
 *  - a non-finite function value, or a C++ exception, inside a callback is recorded; GSL receives a finite placeholder
 *    and the caller checks the record after every set/iterate and throws before it reads any solver result.
 * No C++ exception crosses a GSL frame. All state belongs to one solve: no static or global mutable state.
 *
 * Valid results are unchanged: the bracket check evaluates f(lower), f(upper) in GSL's order, and GSL's own evaluation
 * of the same endpoints in gsl_root_fsolver_set is served from that result, so the sequence of function evaluations
 * is identical to the unguarded code. Hybrid (multiroot) solvers keep their non-finite and non-convergence behaviour;
 * only callback exceptions are contained there.
 * Not preventable without a handler: GSL's allocation failures (GSL_ENOMEM) in *_alloc.
 */
#ifndef XTHERMAL_GSL_GUARD_HPP
#define XTHERMAL_GSL_GUARD_HPP

#include <cmath>
#include <cstdio>
#include <exception>
#include <memory>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>
#include "Exception.h"

namespace xThermal
{
    namespace gsl_guard
    {
        struct RootSolverFree { void operator()(gsl_root_fsolver* s) const { gsl_root_fsolver_free(s); } };
        struct MultirootSolverFree { void operator()(gsl_multiroot_fsolver* s) const { gsl_multiroot_fsolver_free(s); } };
        struct VectorFree { void operator()(gsl_vector* v) const { gsl_vector_free(v); } };
        typedef std::unique_ptr<gsl_root_fsolver, RootSolverFree> RootSolver;
        typedef std::unique_ptr<gsl_multiroot_fsolver, MultirootSolverFree> MultirootSolver;
        typedef std::unique_ptr<gsl_vector, VectorFree> Vector;

        inline std::string num(double x)
        {
            char buf[32];
            std::snprintf(buf, sizeof(buf), "%.17g", x);
            return buf;
        }

        /** One-dimensional function for a Brent solve. Must outlive every use of the solver it is set on. */
        class GuardedFunction
        {
        public:
            GuardedFunction(const gsl_function& f, const char* where) : m_f(f), m_where(where)
            {
                m_gsl.function = &eval;
                m_gsl.params = this;
            }
            /** Check the bracket as gsl_root_fsolver_set / brent_init would, then set the solver. */
            void set(gsl_root_fsolver* s, double x_lower, double x_upper)
            {
                if (!(x_lower <= x_upper))
                    throw SolutionError(std::string("Fatal error in ") + m_where + ": invalid bracket [" + num(x_lower) + ", " + num(x_upper) + "]");
                // evaluated here, in C++ frames: an exception from the function propagates normally
                const double f_lower = GSL_FN_EVAL(&m_f, x_lower);
                const double f_upper = GSL_FN_EVAL(&m_f, x_upper);
                if (!std::isfinite(f_lower) || !std::isfinite(f_upper))
                    throw SolutionError(std::string("Fatal error in ") + m_where + ": function value is not finite at the bracket ["
                                        + num(x_lower) + ", " + num(x_upper) + "]: f = " + num(f_lower) + ", " + num(f_upper));
                if ((f_lower < 0.0 && f_upper < 0.0) || (f_lower > 0.0 && f_upper > 0.0))
                    throw SolutionError(std::string("Fatal error in ") + m_where + ": bracket [" + num(x_lower) + ", " + num(x_upper)
                                        + "] does not contain a root: f = " + num(f_lower) + ", " + num(f_upper));
                m_cache_x[0] = x_lower; m_cache_f[0] = f_lower;
                m_cache_x[1] = x_upper; m_cache_f[1] = f_upper;
                m_cache_n = 2; m_cache_next = 0;
                const int status = gsl_root_fsolver_set(s, &m_gsl, x_lower, x_upper);
                m_cache_n = 0;
                check(status, "gsl_root_fsolver_set");
            }
            /** One Brent iteration; throws if the function failed during it. Returns GSL's status (GSL_SUCCESS). */
            int iterate(gsl_root_fsolver* s)
            {
                const int status = gsl_root_fsolver_iterate(s);
                check(status, "gsl_root_fsolver_iterate");
                return status;
            }

        private:
            static double eval(double x, void* self)
            {
                GuardedFunction* g = static_cast<GuardedFunction*>(self);
                if (g->m_cache_next < g->m_cache_n && x == g->m_cache_x[g->m_cache_next])
                    return g->m_cache_f[g->m_cache_next++];
                g->m_cache_n = 0;
                if (g->m_exception || g->m_nonfinite) return 0.0;   // already failed: do not evaluate again
                try
                {
                    const double y = GSL_FN_EVAL(&g->m_f, x);
                    if (std::isfinite(y)) return y;
                    g->m_nonfinite = true;
                    g->m_x_failed = x;
                    g->m_f_failed = y;
                }
                catch (...)
                {
                    g->m_exception = std::current_exception();
                }
                return 0.0;   // finite placeholder; check() throws before any solver result is used
            }
            void check(int status, const char* call) const
            {
                if (m_exception) std::rethrow_exception(m_exception);
                if (m_nonfinite)
                    throw SolutionError(std::string("Fatal error in ") + m_where + ": function value is not finite: f("
                                        + num(m_x_failed) + ") = " + num(m_f_failed));
                if (status != GSL_SUCCESS)
                    throw SolutionError(std::string("Fatal error in ") + m_where + ": " + call + ": " + gsl_strerror(status));
            }
            gsl_function m_f, m_gsl;
            const char* m_where;
            std::exception_ptr m_exception;
            bool m_nonfinite = false;
            double m_x_failed = 0, m_f_failed = 0;
            double m_cache_x[2] = {0, 0}, m_cache_f[2] = {0, 0};
            int m_cache_n = 0, m_cache_next = 0;
        };

        /** Multidimensional function for a hybrid solve: contains callback exceptions; values pass through unchanged. */
        class GuardedMultirootFunction
        {
        public:
            explicit GuardedMultirootFunction(const gsl_multiroot_function& f) : m_f(f)
            {
                m_gsl.f = &eval;
                m_gsl.n = f.n;
                m_gsl.params = this;
            }
            /** Returns GSL's status unchanged (the caller's convergence handling is not changed). */
            int set(gsl_multiroot_fsolver* s, const gsl_vector* x)
            {
                const int status = gsl_multiroot_fsolver_set(s, &m_gsl, x);
                check();
                return status;
            }
            int iterate(gsl_multiroot_fsolver* s)
            {
                const int status = gsl_multiroot_fsolver_iterate(s);
                check();
                return status;
            }

        private:
            static int eval(const gsl_vector* x, void* self, gsl_vector* f)
            {
                GuardedMultirootFunction* g = static_cast<GuardedMultirootFunction*>(self);
                if (!g->m_exception)
                {
                    try
                    {
                        return g->m_f.f(x, g->m_f.params, f);
                    }
                    catch (...)
                    {
                        g->m_exception = std::current_exception();
                    }
                }
                gsl_vector_set_zero(f);   // finite placeholder; GSL returns the status without calling its error handler
                return GSL_EBADFUNC;
            }
            void check() const
            {
                if (m_exception) std::rethrow_exception(m_exception);
            }
            gsl_multiroot_function m_f, m_gsl;
            std::exception_ptr m_exception;
        };
    }
}
#endif
