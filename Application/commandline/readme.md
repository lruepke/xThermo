# Command line tool: xThermal

The `xThermal` binary exposes the EOS as a command-line tool.
General usage:

```bash
xThermal [program] [arguments]
xThermal [program] -h        # help for a specific program
```

Available programs: `thermo`, `lutgen`, `lutinfo`, `lut2vtu`, `lookup`.

> **Important:** values must be written directly after the flag letter with **no space**:
> `-P316E5` ✓   `-P 316E5` ✗

## Units

| Quantity | Flag | Unit |
|---|---|---|
| Temperature | `-T`, `-R` | °C |
| Pressure | `-P`, `-R` | Pa |
| Salinity | `-X`, `-R` | kg/kg (mass fraction 0–1; seawater = 0.032) |
| Enthalpy | `-H`, `-R` | J/kg |

## Arguments (`thermo` program)

### Required

* `-D[0|1|2|3]` — calculation **D**imension: 0 = scatter points, 1 = 1-D, 2 = 2-D, 3 = 3-D
* `-V[PTX|PHX|T|P|X|H|PT|PX|TX|PH|HX]` — select independent **V**ariables (must match `-D`):
    1. `-D0`: `PTX` or `PHX`
    2. `-D1`: `T`, `P`, `X`, or `H`
    3. `-D2`: `PT`, `PX`, `TX`, `PH`, or `HX`
    4. `-D3`: `PTX` or `PHX`
* `-Rmin/step/max[/min/step/max...]` — range for each variable in `-V` order (1-D/2-D/3-D only)
* `-P` / `-T` / `-X` / `-H` — fixed value for variables **not** in `-V`

### Optional

* `-Gfile` — input text file (columns match `-V` order) for multi-point scatter calculation
* `-Ofile` — output file (`.csv`, `.txt`, `.vtk`)
* `-tN` — number of OpenMP threads
* `-n` — normalise output in VTK files (2-D/3-D)
* `-Bbackend` — water backend: `IAPS84` (default) or `IAPWS95`
* `-Ffluid` — fluid: `H2O-NaCl` (default), `H2O`, `NaCl`

## Examples

### Single point (T-P-X)
```bash
xThermal thermo -D0 -VPTX -P316E5 -T100 -X0.032
```

### Single point (H-P-X)
```bash
xThermal thermo -D0 -VPHX -H1500E3 -P316E5 -X0.032
```

### Multi-point from file
```bash
xThermal thermo -D0 -VPTX -G../test/PTX.txt -OPTX_0D.csv
xThermal thermo -D0 -VPHX -G../test/PHX.txt -OPHX_0D.csv
```

### 1-D: sweep temperature at fixed P, X
```bash
xThermal thermo -D1 -VT -X0.032 -P399E5 -R0/1/1000 -OT_1D.csv
```

### 1-D: sweep enthalpy at fixed P, X
```bash
xThermal thermo -D1 -VH -X0.032 -P399E5 -R43E3/1E3/3000E3 -OH_1D.csv
```

### 1-D: sweep pressure at fixed T, X
```bash
xThermal thermo -D1 -VP -X0.032 -T100 -R5E5/1E5/1000E5 -OP_1D.csv
```

### 1-D: sweep salinity at fixed T, P
```bash
xThermal thermo -D1 -VX -T100 -P399E5 -R0/0.001/1 -OX_1D.csv
```

### 2-D: P-T grid
```bash
xThermal thermo -D2 -VPT -R1E5/1E5/1000E5/0/1/500 -X0.032 -OPT_2D.vtk
```

### 2-D: P-X grid
```bash
xThermal thermo -D2 -VPX -R100E5/1E5/800E5/0/0.01/1 -T100 -OPX_2D.vtk
```

### 2-D: T-X grid
```bash
xThermal thermo -D2 -VTX -R0/1/800/0/0.01/1 -P100E5 -OTX_2D.vtk
```

### 3-D: P-T-X grid (parallel, 8 threads)
```bash
xThermal thermo -D3 -VPTX -R1E5/10E5/500E5/0/10/600/0/0.01/1 -t8
```
