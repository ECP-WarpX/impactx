# ImpactX

[![CI Status](https://github.com/ECP-WarpX/impactx/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/ECP-WarpX/impactx/actions/workflows/ubuntu.yml)
[![Documentation Status](https://readthedocs.org/projects/impactx/badge/?version=latest)](https://impactx.readthedocs.io)
[![License ImpactX](https://img.shields.io/badge/license-BSD--3--Clause--LBNL-blue.svg)](https://spdx.org/licenses/BSD-3-Clause-LBNL.html)
[![Supported Platforms](https://img.shields.io/badge/platforms-linux%20|%20osx%20|%20win-blue)](https://impactx.readthedocs.io/en/latest/install/users.html)  
[![DOI (source)](https://img.shields.io/badge/DOI%20(source)-10.5281/zenodo.6954922-blue.svg)](https://doi.org/10.5281/zenodo.6954922)
[![DOI (paper)](https://img.shields.io/badge/DOI%20(paper)-10.18429%2FJACoW--NAPAC2022--TUYE2-blue.svg)](https://doi.org/10.18429/JACoW-NAPAC2022-TUYE2)  
[![Development Status](https://img.shields.io/badge/development%20status-beta-orange.svg)](https://en.wikipedia.org/wiki/Software_release_life_cycle)
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

Please prepare you local development environment as follows.
Pick *one* of the methods below:

### Perlmutter (NERSC)

```bash
ssh perlmutter-p1.nersc.gov
```

Now ``cd`` to your ImpactX source directory.

```bash
module load cmake/3.22.0
module load cray-hdf5-parallel/1.12.2.1

# CCache for faster rebuilds
export PATH=/global/common/software/spackecp/perlmutter/e4s-22.05/78535/spack/opt/spack/cray-sles15-zen3/gcc-11.2.0/ccache-4.5.1-ybl7xefvggn6hov4dsdxxnztji74tolj/bin:$PATH

# Python
module load cray-python/3.9.13.1
if [ -d "$HOME/sw/perlmutter/venvs/impactx" ]
then
  source $HOME/sw/perlmutter/venvs/impactx/bin/activate
else
  python3 -m pip install --user --upgrade pip
  python3 -m pip install --user virtualenv
  python3 -m venv $HOME/sw/perlmutter/venvs/impactx
  source $HOME/sw/perlmutter/venvs/impactx/bin/activate

  python3 -m pip install --upgrade pip
  MPICC="cc -target-accel=nvidia80 -shared" python3 -m pip install -U --no-cache-dir -v mpi4py
  python3 -m pip install --upgrade pytest
  python3 -m pip install --upgrade -r requirements.txt
  python3 -m pip install --upgrade -r examples/requirements.txt
fi

# GPU-aware MPI
export MPICH_GPU_SUPPORT_ENABLED=1

# necessary to use CUDA-Aware MPI and run a job
export CRAY_ACCEL_TARGET=nvidia80

# optimize CUDA compilation for A100
export AMREX_CUDA_ARCH=8.0

# compiler environment hints
export CC=cc
export CXX=CC
export FC=ftn
export CUDACXX=$(which nvcc)
export CUDAHOSTCXX=CC
```

```bash
# work on an interactive node
salloc -N 1 --ntasks-per-node=4 -t 1:00:00 -C gpu --gpu-bind=single:1 -c 32 -G 4 --qos=interactive -A m3906_g

# configure
cmake -S . -B build_perlmutter -DImpactX_COMPUTE=CUDA -DImpactX_PYTHON=ON

# compile
cmake --build build_perlmutter -j 64

# test
ctest --test-dir build_perlmutter -E AMReX --output-on-failure

# run
cd build_perlmutter/bin
srun ./impactx ../../examples/fodo/input_fodo.in
```

### Cori KNL (NERSC)

```bash
ssh cori.nersc.gov
```

Now ``cd`` to your ImpactX source directory.

```bash
module swap craype-haswell craype-mic-knl
module swap PrgEnv-intel PrgEnv-gnu
module load cmake/3.22.1
module load cray-hdf5-parallel/1.10.5.2
module load cray-fftw/3.3.8.10

# Python
module load cray-python/3.9.7.1
if [ -d "$HOME/sw/knl/venvs/impactx" ]
then
  source $HOME/sw/knl/venvs/impactx/bin/activate
else
  python3 -m pip install --user --upgrade pip
  python3 -m pip install --user virtualenv
  python3 -m venv $HOME/sw/knl/venvs/impactx
  source $HOME/sw/knl/venvs/impactx/bin/activate

  python3 -m pip install --upgrade pip
  MPICC="cc -shared" python3 -m pip install -U --no-cache-dir -v mpi4py
  python3 -m pip install --upgrade pytest
  python3 -m pip install --upgrade -r requirements.txt
  python3 -m pip install --upgrade -r examples/requirements.txt
fi

# tune exactly for KNL sub-architecture
export CXXFLAGS="-march=knl"
export CFLAGS="-march=knl"
```

```bash
# configure
cmake -S . -B build_knl

# compile
cmake --build build_knl -j 8

# test
srun -C knl -N 1 -t 30 -q debug ctest --test-dir build_knl --output-on-failure

# run
cd build_knl/bin
srun -C knl -N 1 -t 30 -q debug ./impactx ../../examples/fodo/input_fodo.in
```

### Homebrew (macOS)

```bash
brew update
brew install adios2      # for openPMD
brew install ccache
brew install cmake
brew install fftw
brew install git
brew install hdf5-mpi    # for openPMD
brew install libomp      # for OpenMP
brew install pkg-config  # for fftw
brew install open-mpi
```

### Apt (Debian/Ubuntu)

```bash
sudo apt update
sudo apt install build-essential ccache cmake g++ git libfftw3-mpi-dev libfftw3-dev libhdf5-openmpi-dev libopenmpi-dev pkg-config python3 python3-matplotlib python3-numpy python3-scipy
```

### Spack (Linux)

```bash
spack env create impactx-dev
spack env activate impactx-dev
spack add adios2        # for openPMD
spack add ccache
spack add cmake
spack add fftw
spack add hdf5          # for openPMD
spack add mpi
spack add pkgconfig     # for fftw
spack add python
spack add py-pip
spack add py-setuptools
spack add py-wheel

# OpenMP support on macOS
[[ $OSTYPE == 'darwin'* ]] && spack add llvm-openmp

# optional: Linux only
#spack add cuda

spack install
python3 -m pip install matplotlib numpy openpmd-api pandas pytest scipy
```

In new terminals, re-activate the environment with `spack env activate impactx-dev` again.

### Conda (Linux/macOS/Windows)

```bash
conda create -n impactx-dev -c conda-forge adios2 ccache cmake compilers git hdf5 fftw matplotlib ninja numpy openpmd-api pandas pytest scipy
conda activate impactx-dev

# compile with -DImpactX_MPI=OFF
```

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
python3 -m pip install --upgrade pip setuptools wheel pytest
python3 -m pip install --upgrade -r examples/requirements.txt
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
