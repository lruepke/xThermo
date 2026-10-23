# Command-line tool

After building and installing xThermo, the `xThermal` binary provides command-line access to the full EOS.

## Usage

```bash
xThermal [program] [arguments]
xThermal [program] -h          # help for a specific sub-program
```

Available sub-programs:

| Program | Description |
|---|---|
| `thermo` | Calculate thermodynamic properties from T-P-X or H-P-X input |
| `lutgen` | Generate an AMR lookup-table (binary `.bin`) |
| `lutinfo` | Print information about a lookup-table file |
| `lut2vtu` | Convert a lookup-table to VTU format for visualisation |
| `lookup` | Query an existing lookup-table at given T/H-P-X points |

!!! warning "No spaces between flag and value"
    The argument parser requires the value to follow the flag letter **immediately**, with no space:
    `-P300E5` ✓ — `-P 300E5` ✗

## Units

| Quantity | Flag | Unit |
|---|---|---|
| Temperature | `-T`, `-R` | °C |
| Pressure | `-P`, `-R` | Pa |
| Salinity | `-X`, `-R` | kg/kg (mass fraction 0–1; seawater = 0.032) |
| Enthalpy | `-H`, `-R` | J/kg |

## Arguments (`thermo`)

### Required

| Flag | Description |
|---|---|
| `-D[0\|1\|2\|3]` | Calculation **D**imension: 0 = scatter, 1 = 1-D, 2 = 2-D, 3 = 3-D |
| `-V…` | Independent **V**ariables — must match `-D` (see table below) |
| `-Pval` / `-Tval` / `-Xval` / `-Hval` | Fixed value for variables **not** in `-V` |
| `-Rmin/step/max[/…]` | Range for each variable in `-V` order (1-D/2-D/3-D only) |

**Valid `-V` options by dimension:**

| `-D` | Valid `-V` |
|---|---|
| `0` | `PTX`, `PHX` |
| `1` | `T`, `P`, `X`, `H` |
| `2` | `PT`, `PX`, `TX`, `PH`, `HX` |
| `3` | `PTX`, `PHX` |

### Optional

| Flag | Description |
|---|---|
| `-Gfile` | Input text file for scatter calculation (columns match `-V` order) |
| `-Ofile` | Output file — format from extension: `.csv`, `.txt`, `.vtk` |
| `-tN` | Number of OpenMP threads |
| `-n` | Normalise VTK output (2-D/3-D) |
| `-BIAPS84` | Water backend: `IAPS84` (default) or `IAPWS95` |
| `-FH2O-NaCl` | Fluid: `H2O-NaCl` (default), `H2O`, `NaCl` |

## Verification

The CLI and Python API give identical results. For T = 300 °C, P = 300 bar (300×10⁵ Pa), X = 0.1 kg/kg:

```bash
xThermal thermo -D0 -VPTX -P300E5 -T300 -X0.1
```

```
H2O-NaCl
T=300 deg.C = 573.15 K, P=300 bar, X=10 wt% NaCl, H=1.21693e+06 J/kg
  phase: Liquid
  Rho=849.034 [kg/m^3]
  H=1.21693e+06 [J/kg]
  Cp=4407.72 [J/kg/K]
  Mu=0.000114076 [Pa S]
```

Equivalent Python:

```python
from xThermal import H2ONaCl
sw = H2ONaCl.cH2ONaCl("IAPS84")
props = sw.UpdateState_TPX(300 + 273.15, 300e5, 0.1)
print(props.Rho)   # 849.034 kg/m³
```

## Examples

### Single point — T-P-X

```bash
xThermal thermo -D0 -VPTX -P300E5 -T300 -X0.1
```

### Single point — H-P-X (enthalpy flash)

```bash
xThermal thermo -D0 -VPHX -H1500E3 -P316E5 -X0.032
```

### Multi-point from file

```bash
# File columns must match -V order (e.g. P  T  X for -VPTX)
xThermal thermo -D0 -VPTX -G../test/PTX.txt -OPTX_result.csv
xThermal thermo -D0 -VPHX -G../test/PHX.txt -OPHX_result.csv
```

### 1-D sweep

```bash
# Temperature sweep at fixed P and X
xThermal thermo -D1 -VT -P399E5 -X0.032 -R0/1/1000 -OT_1D.csv

# Enthalpy sweep at fixed P and X
xThermal thermo -D1 -VH -P399E5 -X0.032 -R43E3/1E3/3000E3 -OH_1D.csv

# Pressure sweep at fixed T and X
xThermal thermo -D1 -VP -T100 -X0.032 -R5E5/1E5/1000E5 -OP_1D.csv

# Salinity sweep at fixed T and P
xThermal thermo -D1 -VX -T100 -P399E5 -R0/0.001/1 -OX_1D.csv
```

### 2-D grid

```bash
# P-T grid at fixed X, output as VTK
xThermal thermo -D2 -VPT -R1E5/1E5/1000E5/0/1/500 -X0.032 -OPT_2D.vtk

# P-X grid at fixed T
xThermal thermo -D2 -VPX -R100E5/1E5/800E5/0/0.01/1 -T100 -OPX_2D.vtk

# T-X grid at fixed P
xThermal thermo -D2 -VTX -R0/1/800/0/0.01/1 -P100E5 -OTX_2D.vtk
```

### 3-D grid

```bash
# Full P-T-X volume, 8 threads
xThermal thermo -D3 -VPTX -R1E5/10E5/500E5/0/10/600/0/0.01/1 -t8
```
