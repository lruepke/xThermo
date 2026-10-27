# Third-party components

xThermo's own source code is licensed under the MIT licence (`LICENSE`). The repository also contains, downloads or
links third-party components that keep their own licences. **The MIT licence does not apply to them**, and binaries
built from xThermo contain some of them (see "Binary distributions" below).

| Component | Authors | Licence (as stated by the component) | Location | In an xThermo build |
|---|---|---|---|---|
| **PROST 4.1** (PROperties of water and STeam; IAPS-84) | Olaf Bauer (1998); Ole Engel (1998); TU Hamburg-Harburg, Technische Thermodynamik | GNU Library General Public License v2.0 (`ThirdParties/PROST/COPYING.LIB-2.0`; stated in every source file) | `ThirdParties/PROST/` (bundled source) | **compiled into libxThermal and the Python module**; provides the `IAPS84` water backend. Modified for xThermo: see `ThirdParties/PROST/XTHERMO_MODIFICATIONS.txt`. |
| **Triangle 1.6** (two-dimensional mesh generator and Delaunay triangulator) | Jonathan Richard Shewchuk (1993–2005) | Triangle's own licence, in the header of `Library/thermo/triangle.cpp`: free redistribution if the copyright notices are kept and **no compensation is received**; modified versions must stay under the author's copyright, be freely available and give clear notice of the modifications; **distribution as part of a commercial system only by direct arrangement with the author**; no warranty. | `Library/thermo/triangle.cpp`, `triangle.h` (bundled; file renamed from `triangle.c`) | **compiled into libxThermal** (`thermo.cpp`: triangulation of phase-boundary surfaces). The copy contains configuration definitions in the source (`TRILIBRARY`, `REAL`, `ANSI_DECLARATORS`, `VOID`) that the original distribution leaves to the build; there is no modification notice yet. |
| **GNU Scientific Library (GSL) 2.7** | The GSL authors (see the GSL `AUTHORS` file) | GNU General Public License v3 or later | not in the repository; downloaded and built by `ThirdParties/gsl.sh` | **statically linked into libxThermal** |
| **numpy.i** (SWIG interface for NumPy arrays) | NumPy Developers (2005–2015) | BSD 3-clause (header of the file) | `Library/H2ONaCl/numpy.i` | used to generate the Python bindings |
| freesteam 2.1 (IAPWS-IF97) | John Pye (2004–2009, per the source file headers) | GNU General Public License v2 (`ThirdParties/freesteam-2.1/LICENSE.txt`) | `ThirdParties/freesteam-2.1/` (bundled source) | **not built** (only the optional comparison configuration, which the default build does not include) |
| CoolProp | Ian H. Bell and other CoolProp developers | MIT (`ThirdParties/coolprop/LICENSE`) | git submodule `ThirdParties/coolprop` | **not built** by default (`USE_COOLPROP` OFF) |

## Binary distributions (libraries, Python modules, installers, packages)

Not legal advice; the points below record what the licence texts require, for the distribution-policy decision.

- **GSL (GPL v3 or later) is linked statically.** A distributed binary that contains GSL is a combined work that may
  only be distributed under the terms of the GPL v3: with the corresponding source code (or a written offer) and the
  licence text. xThermo's MIT-licensed source is compatible with this, but such binaries cannot be offered under the
  MIT licence alone.
- **PROST (LGPL v2.0)** inside a binary requires a notice that PROST is used and is covered by the LGPL, a copy of the
  licence, and the possibility for the recipient to relink with a modified PROST (satisfied by providing the complete
  source).
- **Triangle's licence** permits free, non-commercial redistribution only, and requires an arrangement with the author
  for distribution as part of a commercial system. It is not an open-source licence and adds restrictions that the
  GPL v3 does not allow for a combined work; a binary that contains both GSL and Triangle therefore cannot satisfy both
  licences without such an arrangement or a change of components.
- **numpy.i (BSD)** requires its copyright notice in the documentation of binary distributions.

The distribution policy (for example GPL-licensed binaries, dynamic linking, replacing a component, or an arrangement
with an author) is a separate decision; nothing in the build has been changed. Details:
`review/IAPS84_ATTRIBUTION_LICENSING.md`.
