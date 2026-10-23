

# File gsl\_guard.hpp



[**FileList**](files.md) **>** [**Library**](dir_5ad7f572bbca03234e8e621e192fc099.md) **>** [**thermo**](dir_d760ccf1b5c74bc66b0c51c2e0ac61aa.md) **>** [**gsl\_guard.hpp**](gsl__guard_8hpp.md)

[Go to the source code of this file](gsl__guard_8hpp_source.md)

_Internal (not installed): RAII ownership and exception-safe callbacks for the GSL solvers (C1, LIB-02/18, LIB-27)._ [More...](#detailed-description)

* `#include <cmath>`
* `#include <cstdio>`
* `#include <exception>`
* `#include <memory>`
* `#include <string>`
* `#include <gsl/gsl_errno.h>`
* `#include <gsl/gsl_multiroots.h>`
* `#include <gsl/gsl_roots.h>`
* `#include <gsl/gsl_vector.h>`
* `#include "Exception.h"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**xThermal**](namespacexThermal.md) <br>_Namespace of_ [_**xThermal**_](namespacexThermal.md) _library._ |
| namespace | [**gsl\_guard**](namespacexThermal_1_1gsl__guard.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| class | [**GuardedFunction**](classxThermal_1_1gsl__guard_1_1GuardedFunction.md) <br> |
| class | [**GuardedMultirootFunction**](classxThermal_1_1gsl__guard_1_1GuardedMultirootFunction.md) <br> |
| struct | [**MultirootSolverFree**](structxThermal_1_1gsl__guard_1_1MultirootSolverFree.md) <br> |
| struct | [**RootSolverFree**](structxThermal_1_1gsl__guard_1_1RootSolverFree.md) <br> |
| struct | [**VectorFree**](structxThermal_1_1gsl__guard_1_1VectorFree.md) <br> |


















































## Detailed Description


GSL's process-global error handler is never read or written. Instead, the conditions under which GSL would call it (and, with its default handler, abort the host process) are prevented here:
* Brent bracket: lower &gt; upper (or NaN), non-finite endpoint values and endpoints that do not straddle zero are detected before gsl\_root\_fsolver\_set, with exactly GSL's criteria, and reported as SolutionError;
* a non-finite function value, or a C++ exception, inside a callback is recorded; GSL receives a finite placeholder and the caller checks the record after every set/iterate and throws before it reads any solver result. No C++ exception crosses a GSL frame. All state belongs to one solve: no static or global mutable state.




Valid results are unchanged: the bracket check evaluates f(lower), f(upper) in GSL's order, and GSL's own evaluation of the same endpoints in gsl\_root\_fsolver\_set is served from that result, so the sequence of function evaluations is identical to the unguarded code. Hybrid (multiroot) solvers keep their non-finite and non-convergence behaviour; only callback exceptions are contained there. Not preventable without a handler: GSL's allocation failures (GSL\_ENOMEM) in \*\_alloc. 


    

------------------------------
The documentation for this class was generated from the following file `Library/thermo/gsl_guard.hpp`

