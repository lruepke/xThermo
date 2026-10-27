# Errors and status

A point cannot always be evaluated: the input can be invalid (a NaN, a negative pressure), outside the range of the
model (above 1000 °C for H₂O–NaCl), or a numerical solution can fail. `tpx` and `hpx` handle this with the `errors`
argument.

## errors="raise" (the default)

`errors="raise"` is the default. If any point cannot be evaluated, the call raises an exception for the first such
point (in C order). Its message names the point's index and inputs, the library's message and how many further points
failed:

```python
import numpy as np
import xthermo
from xthermo import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")
T = np.array([300.0, 5000.0, np.nan])                 # K: valid, above the model range, invalid
try:
    xthermo.tpx(sw, T, 250e5, 0.032)
except xthermo.EvaluationError as e:
    print(type(e).__name__)                           # OutOfRangeError
    # (1,) Status.OutOfRange {'T': 5000.0, 'p': 25000000.0, 'X': 0.032}
    print(e.index, e.status, e.inputs)
    print(e)                                          # ends in "[1 further failure]"
```

The exception classes are subclasses of `xthermo.EvaluationError` and of the matching built-in exception, so existing
`except ValueError` code keeps working:

| Exception | Also a | Status | Cause |
|---|---|---|---|
| `xthermo.InvalidInputError` | `ValueError` | `InvalidInput` | non-finite T, H or p; p ≤ 0; X < 0 or non-finite X (H₂O–NaCl) |
| `xthermo.OutOfRangeError` | `ValueError` | `OutOfRange` | outside the range in which the model evaluates, e.g. T > 1273.15 K for H₂O–NaCl |
| `xthermo.NotConvergedError` | `RuntimeError` | `NotConverged` | a numerical solution failed, e.g. IAPWS-95 water exactly at its saturation pressure |
| `xthermo.NotAvailableError` | `ValueError` | `NotAvailable` | not implemented for this model or state, e.g. `hpx` with the pure NaCl model |
| `xthermo.EvaluationError` | `Exception` | `Failed` | any other error reported by the library |

Each exception has the attributes `index` (a tuple; `()` for scalar input), `status`, `inputs` (a dict) and
`message` (the library's message).

## errors="nan"

With `errors="nan"` the call never raises for points that cannot be evaluated. Those points are NaN in every number
field (their inputs are kept), their phase is `Unknown`, and their status says why:

```python
import numpy as np
import xthermo
from xthermo import H2ONaCl

sw = H2ONaCl.cH2ONaCl("IAPS84")
T = np.array([300.0, 5000.0, np.nan])
s = xthermo.tpx(sw, T, 250e5, 0.032, errors="nan")

print(s.ok)                                           # [ True False False]
print([xthermo.Status(int(c)).name for c in s.status])   # ['Ok', 'OutOfRange', 'InvalidInput']
print(s.Rho)                                          # [value, nan, nan]
for index, status, message in s.failures:
    print(index, status.name, message[:60])
```

This is the policy for grids and maps: evaluate everything, then use `s.ok` to select or mask. The plotting helpers
draw points that were not evaluated in grey.

Errors that are not about a point, such as running out of memory or a programming error, are raised in both modes.

## Status values

`xthermo.Status` is an integer enumeration; `s.status` holds these integers (a single `Status` for scalar input).

| Status | Value | Meaning |
|---|---|---|
| `Ok` | 0 | evaluated |
| `InvalidInput` | 1 | invalid input (non-finite, p ≤ 0, X < 0) |
| `OutOfRange` | 2 | outside the range in which the model evaluates |
| `NotConverged` | 3 | a numerical solution failed |
| `NotAvailable` | 4 | not implemented or not available for this model or state |
| `Failed` | 5 | any other error of the library |

## Failed points and absent phases

A point that was not evaluated and a phase that is not present at an evaluated point are different things:

| Question | Test |
|---|---|
| Was the point evaluated? | `s.ok` |
| Is there vapour at an evaluated point? | `s.has("vapour")` (False at failed points) |
| Evaluated, but without vapour? | `s.ok & ~s.has("vapour")` |

By default the fields of an absent phase are 0, as in the library; with `absent="nan"` they are NaN (see
[Evaluating states](evaluate.md#absent-phases)). Failed points are NaN in both cases.

## What is checked, and what is not

The library rejects non-finite inputs, non-positive pressure and negative salinity, and H₂O–NaCl temperatures above
1273.15 K. It does **not** reject every input outside the published model domain (0–1000 °C, 1–5000 bar, X = 0–1):
pressures below 1 bar are evaluated as extrapolations, and a salinity above 1, which is not physical, is evaluated
too. Check your inputs against the domain; see
[Model and limitations](../limitations.md#domain-input-checks-and-errors).
