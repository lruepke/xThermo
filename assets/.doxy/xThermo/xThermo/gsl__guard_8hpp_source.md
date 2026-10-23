

# File gsl\_guard.hpp

[**File List**](files.md) **>** [**Library**](dir_5ad7f572bbca03234e8e621e192fc099.md) **>** [**thermo**](dir_d760ccf1b5c74bc66b0c51c2e0ac61aa.md) **>** [**gsl\_guard.hpp**](gsl__guard_8hpp.md)

[Go to the documentation of this file](gsl__guard_8hpp.md)


```C++

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

        class GuardedFunction
        {
        public:
            GuardedFunction(const gsl_function& f, const char* where) : m_f(f), m_where(where)
            {
                m_gsl.function = &eval;
                m_gsl.params = this;
            }
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

        class GuardedMultirootFunction
        {
        public:
            explicit GuardedMultirootFunction(const gsl_multiroot_function& f) : m_f(f)
            {
                m_gsl.f = &eval;
                m_gsl.n = f.n;
                m_gsl.params = this;
            }
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
```


