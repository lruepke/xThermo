# Errors and threads

## Exceptions

When a state cannot be evaluated, `UpdateState_TPX` and `UpdateState_HPX` throw an exception derived from
`xThermal::xThermalBaseError` (itself a `std::exception`; `what()` gives the message, `code()` the error code). The
library does not return a previous state, a partial result or default values instead.

| Exception | Cause | `EvaluationStatus` |
|---|---|---|
| `InputError` (a `ValueError`) | non-finite T, H, p or X; p ≤ 0; X < 0 | `InvalidInput` |
| `OutOfRangeError` | outside the range in which the model evaluates, e.g. T above 1273.15 K for H₂O–NaCl | `OutOfRange` |
| `SolutionError` | a numerical solution failed, e.g. IAPWS-95 water exactly at its saturation pressure | `NotConverged` |
| `NotImplementedError`, `NotAvailableError` | not provided for this model or state, e.g. H–p–X for pure NaCl | `NotAvailable` |
| other `xThermalBaseError` (e.g. `ValueError`) | any other error reported by the library | `Failed` |

```cpp
#include "H2ONaCl.h"
#include <iostream>

using namespace xThermal;

int main()
{
    H2ONaCl::cH2ONaCl brine("IAPS84");
    ThermodynamicProperties props;
    try {
        brine.UpdateState_TPX(props, 1100 + 273.15, 300e5, 0.1);         // above the model's 1000 °C
    } catch (const OutOfRangeError& e) {
        std::cout << "out of range\n";
    } catch (const xThermalBaseError& e) {                               // any other failure
        std::cout << evaluation_status_name(evaluation_status(e.code())) << ": " << e.what() << "\n";
    }
}
```

## Status codes

`tryUpdateState_*` and the batch functions `evaluate_*` report the same failures as an `EvaluationStatus`
([Evaluating states](tpx.md#without-exceptions)). `evaluation_status(e.code())` converts an exception's code, and
`evaluation_status_name(status)` gives its name (`"OutOfRange"` for `EvaluationStatus::OutOfRange`). In the batch
results, `status[i]` holds the value as an integer:

| Value | `EvaluationStatus` | Python (`xthermo.Status`) |
|---|---|---|
| 0 | `Ok` | `Ok` |
| 1 | `InvalidInput` | `InvalidInput` |
| 2 | `OutOfRange` | `OutOfRange` |
| 3 | `NotConverged` | `NotConverged` |
| 4 | `NotAvailable` | `NotAvailable` |
| 5 | `Failed` | `Failed` |

Errors that are not evaluation failures (for example `std::bad_alloc`) are not converted into a status; they
propagate from every function.

## Input checks

Every evaluation rejects non-finite inputs, p ≤ 0 and X < 0 (`InputError`), and H₂O–NaCl temperatures above
1273.15 K (`OutOfRangeError`). Other inputs outside the published model domain (0–1000 °C, 1–5000 bar, X 0–1) are
**not** all rejected: pressures below 1 bar, and X above 1, are evaluated as extrapolations. Check inputs against the
domain where it matters; [Domain, input checks and errors](../limitations.md#domain-input-checks-and-errors) has the
full list.

A single property that cannot be evaluated at an otherwise valid state is NaN, and the other properties remain valid
([Result fields](tpx.md#result-fields)).

## Threads

**A model object is not thread-safe.** Code that evaluates states concurrently needs **one model object per thread**.
Independent objects share no state and give the same results as a serial evaluation. Construct a new object per
thread rather than copying one:

```cpp
#include "H2ONaCl.h"
#include <iostream>
#include <thread>
#include <vector>

using namespace xThermal;

int main()
{
    std::vector<double> T = {300 + 273.15, 350 + 273.15, 400 + 273.15, 450 + 273.15};   // K
    std::vector<double> rho(T.size());
    std::vector<std::thread> threads;
    for (size_t i = 0; i < T.size(); i++)
        threads.emplace_back([i, &T, &rho]() {
            H2ONaCl::cH2ONaCl brine("IAPS84");        // one model per thread
            rho[i] = brine.UpdateState_TPX(T[i], 300e5, 0.1).Rho;
        });
    for (auto& t : threads) t.join();
    for (double r : rho) std::cout << r << " kg/m3\n";
}
```

MPI ranks have their own objects anyway. In OpenMP code of your own, create the model inside the parallel region.

**Internal OpenMP** (optional). A library built with `-DXTHERMO_OPENMP=ON` can evaluate the vector functions
(`evaluate_*` and the vector `UpdateState_*`) in parallel: `model.set_num_threads(n)` sets the number of threads for
that object. It creates one model per thread itself and gives the same results as a serial evaluation. The default
is 1 thread, and the build option is off by default. Without it, `set_num_threads` has no effect and
`get_num_threads()` returns 1. In a program that is already parallel (a simulator with MPI or its own threads), keep
the internal threading off.
