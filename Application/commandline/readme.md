# Command-line tool `xThermal`

`xThermal` evaluates the xThermo models from the shell (single states, profiles, grids, files of points, lookup
tables). It is built and installed with the library (`<prefix>/bin/xThermal`).

The documentation, with units, options, input-file and output formats, exit status and runnable examples, is at

**https://lruepke.github.io/xThermo/cli/**

and in this repository in `docs/cli.md`. Quick start:

    xThermal thermo -D0 -VPTX -T300 -P300E5 -X0.1     # 300 °C, 300 bar, 10 wt% NaCl
    xThermal thermo -h                                # help
