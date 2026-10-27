# About

## Development

xThermo is developed and maintained by the [Seafloor Modeling group](https://www.geomar.de/en/research/fb4/fb4-muhs/research-topics/modelings) at the [GEOMAR Helmholtz Centre for Ocean Research Kiel](https://www.geomar.de).

**Zhikui Guo** developed the original C++ library, including the water and NaCl backends and Python bindings, building on earlier MATLAB implementations of the Driesner formulation developed at GEOMAR.

**Lars Rüpke** leads the project and maintains its scientific development.

xThermo supports research on hydrothermal systems and is used, among other applications, in OpenFOAM-based hydrothermal flow simulations with HydrothermalFoam.

## How to cite

xThermo implements the H₂O–NaCl formulations developed by Thomas Driesner and Christoph Heinrich. When publishing results obtained with xThermo, **please cite both foundational papers**:

Driesner, T. & Heinrich, C.A. (2007). The system H₂O–NaCl. Part I: Correlation formulae for phase relations in
temperature–pressure–composition space from 0 to 1000 °C, 0 to 5000 bar, and 0 to 1 X<sub>NaCl</sub>.
*Geochimica et Cosmochimica Acta*, **71**, 4880–4901.
[doi:10.1016/j.gca.2006.01.033](https://doi.org/10.1016/j.gca.2006.01.033)

Driesner, T. (2007). The system H₂O–NaCl. Part II: Correlations for molar volume, enthalpy, and isobaric heat
capacity from 0 to 1000 °C, 1 to 5000 bar, and 0 to 1 X<sub>NaCl</sub>. *Geochimica et Cosmochimica Acta*, **71**,
4902–4919. [doi:10.1016/j.gca.2007.05.026](https://doi.org/10.1016/j.gca.2007.05.026)

Please also cite the pure-water formulation used in your calculations:

- IAPS-84: Haar, L., Gallagher, J.S. & Kell, G.S. (1984). *NBS/NRC Steam Tables*. Hemisphere, Washington; and
  Kestin, J. & Sengers, J.V. (1986). *J. Phys. Chem. Ref. Data* **15**, 305.
  [doi:10.1063/1.555772](https://doi.org/10.1063/1.555772)
- IAPWS-95: Wagner, W. & Pruß, A. (2002). *J. Phys. Chem. Ref. Data* **31**, 387–535.
  [doi:10.1063/1.1461829](https://doi.org/10.1063/1.1461829)

For calculations involving brine viscosity, please cite Klyukin, Yu.I., Lowell, R.P. & Bodnar, R.J. (2017). *Fluid Phase Equilib.* **433**,
193–205. [doi:10.1016/j.fluid.2016.11.002](https://doi.org/10.1016/j.fluid.2016.11.002). All references are listed
on [Model and limitations](limitations.md#references).

To cite the software itself, use *xThermo, GEOMAR Helmholtz Centre for Ocean Research Kiel,
https://github.com/lruepke/xThermo*, and specify the version used (`xthermo.__version__`). A software DOI is planned.

## Acknowledgements

We gratefully acknowledge **Falko Vehling, Jörg Hasenclever and Karthik Iyer** for their contributions to earlier MATLAB implementations of the Driesner H₂O–NaCl formulation at GEOMAR, which provided a foundation for the subsequent development of xThermo.

xThermo also builds on the following third-party software:

- **PROST 4.1** (Olaf Bauer and Ole Engel, TU Hamburg-Harburg) computes the IAPS-84 water properties.
- **The GNU Scientific Library** (GSL) provides the one- and multidimensional root finders.
- **Triangle** (Jonathan R. Shewchuk) provides 2-D triangulation (`cxThermal::Triangulation`).
- **SWIG** generates the Python bindings.

## License

**xThermo's own source code is released under the MIT licence** (below).

**Builds of xThermo contain third-party components under their own licences:**

| Component | Licence | In a build |
|---|---|---|
| PROST 4.1, the IAPS-84 water backend | GNU LGPL v2.0 | compiled in |
| GSL | GNU GPL v3 or later | linked statically |
| Triangle 1.6 (J. R. Shewchuk) | free for **non-commercial** use and redistribution only; commercial distribution needs the author's permission | compiled in |

A compiled xThermo library or Python package is therefore **not covered by the MIT licence alone**. Its
redistribution is subject to these licences, in particular the GPL and Triangle's terms. If you redistribute
binaries, check these terms first. The details, including optional components that are not built by default, are in
[`THIRD_PARTY_NOTICES.md`](https://github.com/lruepke/xThermo/blob/main/THIRD_PARTY_NOTICES.md).

```text
MIT License

Copyright (c) 2024 Zhikui Guo, Lars Rüpke
GEOMAR Helmholtz Centre for Ocean Research Kiel

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
