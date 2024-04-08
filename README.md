# ImpactX

[![CI Status](https://github.com/ECP-WarpX/impactx/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/ECP-WarpX/impactx/actions/workflows/ubuntu.yml)
[![Documentation Status](https://readthedocs.org/projects/impactx/badge/?version=latest)](https://impactx.readthedocs.io)
[![License ImpactX](https://img.shields.io/badge/license-BSD--3--Clause--LBNL-blue.svg)](https://spdx.org/licenses/BSD-3-Clause-LBNL.html)
[![Supported Platforms](https://img.shields.io/badge/platforms-linux%20|%20osx%20|%20win-blue)](https://impactx.readthedocs.io/en/latest/install/users.html)  
[![DOI (source)](https://img.shields.io/badge/DOI%20(source)-10.5281/zenodo.6954922-blue.svg)](https://doi.org/10.5281/zenodo.6954922)
[![DOI (paper)](https://img.shields.io/badge/DOI%20(paper)-10.18429%2FJACoW--NAPAC2022--TUYE2-blue.svg)](https://doi.org/10.18429/JACoW-NAPAC2022-TUYE2)  
[![Language: C++17](https://img.shields.io/badge/language-C%2B%2B17-orange.svg)](https://isocpp.org/)
[![Language: Python](https://img.shields.io/badge/language-Python-orange.svg)](https://python.org/)

ImpactX: an s-based beam dynamics code including space charge effects.
This is the next generation of the [IMPACT-Z](https://github.com/impact-lbl/IMPACT-Z) code.

## Documentation

In order to learn how to install and run the code, please see the online documentation:
https://impactx.readthedocs.io

* ImpactX Doxygen: https://impactx.readthedocs.io/en/latest/_static/doxyhtml
* AMReX Doxygen: https://amrex-codes.github.io/amrex/doxygen
* WarpX Doxygen: https://warpx.readthedocs.io/en/latest/_static/doxyhtml

## Contributing

[![AMReX](https://img.shields.io/static/v1?label="runs%20on"&message="AMReX"&color="blueviolet")](https://amrex-codes.github.io/)

Our workflow is described in [CONTRIBUTING.rst](CONTRIBUTING.rst).

## Developer Environment

Please see our [developer installation section](https://impactx.readthedocs.io/en/latest/install/dependencies.html#install-dependencies) of the documentation for an easy install of our software dependencies.

## Get the Source Code

Before you start, you will need a copy of the ImpactX source code:

```bash
git clone git@github.com:ECP-WarpX/impactx.git
cd impactx
```

## Compile

```bash
# find dependencies & configure
cmake -S . -B build

# compile
cmake --build build -j 4
```

That's all!
ImpactX binaries are now in `build/bin/`.
Most people execute these binaries directly or copy them out.

You can inspect and modify build options after running `cmake -S . -B` build with either

```bash
ccmake build
```

or by adding arguments with `-D<OPTION>=<VALUE>` to the first CMake call, e.g.:

```bash
cmake -S . -B build -DImpactX_COMPUTE=CUDA -DImpactX_MPI=OFF
```

### Python Compile

```bash
# find dependencies & configure
cmake -S . -B build -DImpactX_PYTHON=ON

# compile & install
cmake --build build -j 4 --target pip_install
```

## Run

An executable ImpactX binary with the current compile-time options encoded in its file name will be created in `build/bin/`.

Additionally, a symbolic link named `impactx` can be found in that directory, which points to the last built ImpactX executable.

The command-line syntax for this executable is:
```console
Usage: impactx <inputs-file> [some.overwritten.option=value]...

Mandatory arguments (remove the <>):
  inputs-file     the path to an input file; can be relative to the current
                  working directory or absolute.
                  Example: input_fodo.in

Optional arguments (remove the []):
  options         this can overwrite any line in an inputs-file
                  Example: quad1.ds=0.5 sbend1.rc=1.5

Examples:
  In the current working directory, there is a file "input_fodo.in" and the
  "impactx" executable.
  The line to execute would look like this:
    ./impactx input_fodo.in

  In the current working directory, there is a file "input_fodo.in" and the
  executable "impactx" is in a directory that is listed in the "PATH"
  environment variable.
  The line to execute would look like this:
    impactx input_fodo.in

  In the current working directory, there is a file "input_fodo.in" and the
  "impactx" executable. We want to voerwrite the segment length of the beamline
  element "quad1" that is already defined in it. We also want to change the
  radius of curvature of the bending magnet "sbend1" to a different value than
  in the file "input_fodo.in".
  The line to execute would look like this:
    ./impactx input_fodo.in quad1.ds=0.5 sbend1.rc=1.5
```

## Test

In order to run our tests, you need to have a few Python packages installed:
```console
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade build packaging setuptools wheel pytest
python3 -m pip install --upgrade -r tests/python/requirements.txt
```

You can run all our tests with:

```console
ctest --test-dir build --output-on-failure
```

Further options:
* help: `ctest --test-dir build --help`
* list all tests: `ctest --test-dir build -N`
* only run tests that have "FODO" in their name: `ctest --test-dir build -R FODO`

## Acknowledgements

This work was supported by the Laboratory Directed Research and Development Program of Lawrence Berkeley National Laboratory under U.S. Department of Energy Contract No. DE-AC02-05CH11231.

ImpactX is supported by the CAMPA collaboration, a project of the U.S. Department of Energy, Office of Science, Office of Advanced Scientific Computing Research and Office of High Energy Physics, Scientific Discovery through Advanced Computing (SciDAC) program.

## Copyright Notice

ImpactX Copyright (c) 2022, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).
All rights reserved.

If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Intellectual Property Office at IPO@lbl.gov.

NOTICE. This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit others to do so.

Please see the full license agreement in [LICENSE.txt](LICENSE.txt). The SPDX license identifier is `BSD-3-Clause-LBNL`.
