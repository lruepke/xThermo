# Command-line tool

`xThermal` evaluates the models from the shell: single states, profiles, grids and files of points. It is installed
with the library in `<prefix>/bin/`. Every example on this page runs in the documentation's test suite exactly as
printed.

```text
xThermal <program> [options]
xThermal <program> -h          help for one program
xThermal -h | -v               general help, version
```

| Program | Purpose |
|---|---|
| `thermo` | evaluate states from T–p–X or H–p–X input |
| `lutgen`, `lutinfo`, `lookup`, `lut2vtu` | build, inspect, query and convert adaptive lookup tables ([below](#lookup-tables)) |

## A first example

Seawater-like brine (10 wt% NaCl) at 300 °C and 300 bar:

```bash
xThermal thermo -D0 -VPTX -T300 -P300E5 -X0.1
```

```text
H2O-NaCl
T=300 deg.C = 573.15 K, P=300 bar, X=10 wt% NaCl, H=1.21693e+06 J/kg
  phase: Liquid
  Rho=849.034 [kg/m^3]
  H=1.21693e+06 [J/kg]
  Cp=4407.72 [J/kg/K]
  Mu=0.000114076 [Pa S]
  dRhodT=-1.3541 [kg/m3/K]
  dRhodP=1.06661e-06 [kg/m3/Pa]
  Isothermal compressibility=1.25627e-09 [1/Pa]
  Isobaric expansivity=0.00159487 [1/K]
```

The values are those of the Python and C++ interfaces (`xthermo.tpx(sw, 573.15, 300e5, 0.1)`).

## Units

!!! warning "Temperature: °C on the command line, K in files"
    `-T` and the temperature ranges of `-R` are in **degrees Celsius**. Input files (`-G`) and every output file
    (CSV and VTK) use **kelvin**.

| Quantity | Options | Unit |
|---|---|---|
| temperature | `-T`, `-R` | °C |
| temperature | `-G` files; CSV and VTK output | K |
| pressure | `-P`, `-R`, files, output | Pa (300 bar = `300E5`) |
| salinity | `-X`, `-R`, files, output | kg NaCl / kg fluid, 0–1 (seawater: 0.032) |
| specific enthalpy | `-H`, `-R`, files, output | J/kg |

## Options of `thermo`

A value follows its option letter directly or after a space: `-P300E5` and `-P 300E5` are the same.

| Option | Meaning |
|---|---|
| `-D` *n* | **required**: 0 = single points or a file of points, 1 = profile, 2 = grid, 3 = 3-D grid |
| `-V` *vars* | **required**: the varying variables, in order (see the table below) |
| `-T`, `-P`, `-X`, `-H` | fixed values of the variables not in `-V` (`-T` and `-H` exclude each other) |
| `-R` *min/step/max…* | ranges of the variables in `-V`, in the same order (`-D1` to `-D3`) |
| `-G` *file* | file of points (`-D0`) |
| `-O` *file* | output file ([Output](#output)) |
| `-F` *fluid* | `H2O-NaCl` (default), `H2O` or `NaCl` |
| `-B` *model* | water model: `IAPS84` (default) or `IAPWS95` |
| `-t` *n* | threads; effective only in a library built with `-DXTHERMO_OPENMP=ON`, otherwise ignored |
| `-p` | show a progress bar (for interactive use on a terminal) |
| `-n` | accepted for compatibility; has no effect |
| `-h`, `-v` | help, version |

| `-D` | `-V` | Needs |
|---|---|---|
| 0 | `PTX` or `PHX` (any order) | `-T`/`-H`, `-P`, `-X`, or a file `-G` |
| 1 | `T`, `P`, `X` or `H` | the two other variables fixed, `-R` with 3 values, `-O` |
| 2 | `PT`, `PX`, `TX`, `PH` or `HX` (either order) | the third variable fixed, `-R` with 6 values |
| 3 | `PTX` or `PHX` (any order) | `-R` with 9 values; H₂O–NaCl only |

**Ranges.** `-Rmin/step/max` gives evenly spaced values from min to max, both included: (max − min)/step + 1
points, rounded down, with the spacing adjusted to fit. min must be smaller than max, and the step positive and not
larger than max − min. `-R1/100/1000` gives 10 values, 111 K apart.

With `-F H2O` the salinity is 0 and with `-F NaCl` it is 1, so `-X` is not needed.

## Examples

**Single point from enthalpy**:

```bash
xThermal thermo -D0 -VPHX -H1500E3 -P316E5 -X0.032
```

**A file of points**. Three numbers per line in the order of `-V`, temperature in kelvin; blank lines are skipped:

```text title="points.txt"
300E5 573.15 0.1
250E5 673.15 0.032
```

```bash
xThermal thermo -D0 -VPTX -Gpoints.txt -Opoints_result.csv
```

**Profiles** at fixed pressure and salinity, in temperature and in enthalpy:

```bash
xThermal thermo -D1 -VT -P250E5 -X0.032 -R1/1/1000 -OT_profile.csv
xThermal thermo -D1 -VH -P250E5 -X0.032 -R100E3/10E3/3000E3 -OH_profile.csv
```

and in pressure and salinity:

```bash
xThermal thermo -D1 -VP -T350 -X0.032 -R5E5/5E5/1000E5 -OP_profile.csv
xThermal thermo -D1 -VX -T350 -P250E5 -R0/0.01/1 -OX_profile.csv
```

**Grids**: temperature–salinity at 250 bar, and pressure–temperature for seawater:

```bash
xThermal thermo -D2 -VTX -R1/5/1000/0/0.01/1 -P250E5 -OTX_250bar.vtk
xThermal thermo -D2 -VPT -R1E5/10E5/1000E5/1/5/800 -X0.032 -OPT_seawater.vtk
```

**A 3-D grid** in pressure, temperature and salinity:

```bash
xThermal thermo -D3 -VPTX -R100E5/100E5/1000E5/100/50/1000/0/0.05/1 -OPTX.vtk
```

**Pure water** with IAPWS-95, and pure NaCl:

```bash
xThermal thermo -D0 -VPTX -T380 -P250E5 -FH2O -BIAPWS95
xThermal thermo -D0 -VPTX -T900 -P300E5 -FNaCl
```

## Output

| Mode | Output |
|---|---|
| `-D0` with `-T/-H`, `-P`, `-X` | the property listing on the screen (stdout) |
| `-D0` with `-G`, `-D1` | a CSV file: `-O` name, or `out_<V>.csv` for `-G` without `-O`; `-D1` requires `-O` |
| `-D2`, `-D3` | a legacy VTK file and a ParaView script (`<name>.vtk.py`): `-O` name, or `<V>.vtk` without `-O` |

**CSV files** have one row per point. The first columns are `T[K], P[Pa], X[kg/kg]`, then the phase compositions,
the phase region (`Phase[index]`, `Phase[name]`), and the bulk and phase densities, enthalpies, heat capacities,
viscosities and saturations. For H₂O–NaCl the saturations `Sl[m3/m3], Sv, Sh` are volume fractions; for pure water,
`Sl[kg/kg], Sv` are mass fractions. The bulk heat capacity and viscosity of multiphase states follow the definitions
on [Model and limitations](limitations.md#multiphase-properties) (the viscosity column of a two-phase point holds a
placeholder 1).

A CSV name that does not end in `.csv` gets `.csv` appended (`result.txt` gives `result.txt.csv`). A grid is always
written as `.vtk`: `-Ogrid.csv` or `-Ogrid` gives `grid.vtk`, with a notice.

## Exit status and errors

| Status | Meaning |
|---|---|
| 0 | success (also `-h` and `-v`) |
| 1 | the evaluation failed: an input outside the model's range, a state the model cannot evaluate, an unreadable or malformed input file, an unwritable output file |
| 2 | the command line cannot be interpreted: an unknown program, option, fluid or water model; a missing or malformed value; a missing required option; an invalid `-D`, `-V` or `-R`; both `-T` and `-H` |

Errors and warnings go to stderr, results to stdout. A failing point in a file, profile or grid stops the run
(status 1, no output file).

```bash
# above the model's 1000 °C: exit status 1
xThermal thermo -D0 -VPTX -T1100 -P300E5 -X0.1 || echo "exit status $?"
# -X missing: exit status 2
xThermal thermo -D0 -VPTX -T300 -P300E5 || echo "exit status $?"
```

The command-line tool checks T, p and X against the model's range before evaluating. Inputs outside the published
model domain are rejected there, including pressures below 1 bar, which the library itself would evaluate (see
[Domain, input checks and errors](limitations.md#domain-input-checks-and-errors)).

## Lookup tables

`lutgen` builds an adaptive lookup table of the properties over a 2-D or 3-D region, `lutinfo` describes it, `lookup`
evaluates points from it, and `lut2vtu` converts it for ParaView:

```bash
xThermal lutgen -D2 -VXT -R0/0.05/1/100/10/500 -P250E5
xThermal lutinfo lut_constP_XT_6.bin
xThermal lookup lut_constP_XT_6.bin -D0 -VTPX -T300 -P250E5 -X0.1
```

`lutgen` writes its files into the current directory, with names from the table type and refinement level
(`lut_constP_XT_<level>.bin`, one file per property, and a `.vtu` file); it does not use `-O`. `-L` sets the
minimum refinement level (default 4). Interpolated values from a table are not covered by the verification of the
model.

## Known limitations

- `-t` has an effect only in a library built with OpenMP; `-n` has no effect.
- Error messages for a failed evaluation give the library's explanation but not the input point.
