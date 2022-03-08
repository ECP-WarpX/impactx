# ImpactX

[![CI Status](https://github.com/ECP-WarpX/impactx/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/ECP-WarpX/impactx/actions/workflows/ubuntu.yml)
[![Documentation Status](https://readthedocs.org/projects/impactx/badge/?version=latest)](https://impactx.readthedocs.io)
[![License ImpactX](https://img.shields.io/badge/license-BSD--3--Clause--LBNL-blue.svg)](https://spdx.org/licenses/BSD-3-Clause-LBNL.html)  
[![Supported Platforms](https://img.shields.io/badge/platforms-linux%20|%20osx%20|%20win-blue)](https://impactx.readthedocs.io/en/latest/install/users.html)
[![Development Status](https://img.shields.io/badge/development%20status-pre--alpha-orange.svg)]()
[![Language: C++17](https://img.shields.io/badge/language-C%2B%2B17-orange.svg)](https://isocpp.org/)
[![Language: Python](https://img.shields.io/badge/language-Python-orange.svg)](https://python.org/)

ImpactX: the next generation of the [IMPACT-Z](https://github.com/impact-lbl/IMPACT-Z) code

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
```bash
module load cmake/git-20210830  # 3.22-dev
module swap PrgEnv-nvidia PrgEnv-gnu
module swap gcc gcc/9.3.0
module load cuda
module load cray-hdf5-parallel/1.12.0.7

# GPU-aware MPI
export MPICH_GPU_SUPPORT_ENABLED=1

# optimize CUDA compilation for A100
export AMREX_CUDA_ARCH=8.0

# compiler environment hints
export CC=$(which gcc)
export CXX=$(which g++)
export FC=$(which gfortran)
export CUDACXX=$(which nvcc)
export CUDAHOSTCXX=$(which g++)
```

```bash
# configure
cmake -S . -B build_perlmutter -DImpactX_COMPUTE=CUDA

# compile
cmake --build build_perlmutter -j 10

# run
cd build_perlmutter/bin
srun -N 1 --ntasks-per-node=4 -t 0:10:00 -C gpu -c 32 -G 4 --qos=debug -A m3906_g ./impactx ../../examples/input_fodo.in
```

### Cori KNL (NERSC)

```bash
ssh cori.nersc.gov
```

```bash
module swap craype-haswell craype-mic-knl
module swap PrgEnv-intel PrgEnv-gnu
module load cmake/3.21.3
module load cray-hdf5-parallel/1.10.5.2
module load cray-fftw/3.3.8.4
module load cray-python/3.7.3.2
```

```bash
# configure
cmake -S . -B build_cori

# compile
cmake --build build_cori -j 8

# run
cd build_cori/bin
srun -C knl -N 1 -t 30 -q debug ./impactx ../../examples/input_fodo.in
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

# OpenMP support on macOS
[[ $OSTYPE == 'darwin'* ]] && spack add llvm-openmp

# optional:
# spack add cuda
# spack add python
# spack add py-pip
# spack add py-pandas
# spack add py-numpy
# spack add py-scipy

spack install
```

(in new terminals, re-activate the environment with `spack env activate impactx-dev` again)

### Conda (Linux/macOS/Windows)

```bash
conda create -n impactx-dev -c conda-forge adios2 ccache cmake compilers git hdf5 fftw matplotlib ninja numpy pandas scipy
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
python3 -m pip install -U pip setuptools wheel
python3 -m pip install -r requirements.txt
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

## License

Copyright (c) 2021-2022, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.

If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Intellectual Property Office at IPO@lbl.gov.

This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit others to do so.
