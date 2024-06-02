#!/bin/bash
#
# Copyright 2023 The ImpactX Community
#
# This file is part of ImpactX.
#
# Author: Axel Huebl, Luca Fedeli
# License: BSD-3-Clause-LBNL

# Exit on first error encountered #############################################
#
set -eu -o pipefail


# Check: ######################################################################
#
#   Was lumi_cpu_impactx.profile sourced and configured correctly?
if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your lumi_cpu_impactx.profile file! Please edit its line 2 to continue!"; exit 1; fi


# Remove old dependencies #####################################################
#
SRC_DIR="${HOME}/src"
SW_DIR="${HOME}/sw/lumi/cpu"
rm -rf ${SW_DIR}
mkdir -p ${SW_DIR}
mkdir -p ${SRC_DIR}

# remove common user mistakes in python, located in .local instead of a venv
python3 -m pip uninstall -qq -y pyimpactx
python3 -m pip uninstall -qq -y impactx
python3 -m pip uninstall -qqq -y mpi4py 2>/dev/null || true


# General extra dependencies ##################################################
#

# tmpfs build directory: avoids issues often seen with $HOME and is faster
build_dir=$(mktemp -d)

# c-blosc (I/O compression, for openPMD)
if [ -d ${SRC_DIR}/c-blosc ]
then
  cd ${SRC_DIR}/c-blosc
  git fetch --prune
  git checkout v1.21.1
  cd -
else
  git clone -b v1.21.1 https://github.com/Blosc/c-blosc.git ${SRC_DIR}/c-blosc
fi
rm -rf ${SRC_DIR}/c-blosc-build
cmake -S ${SRC_DIR}/c-blosc             \
      -B ${build_dir}/c-blosc-build  \
      -DBUILD_TESTS=OFF                 \
      -DBUILD_BENCHMARKS=OFF            \
      -DDEACTIVATE_AVX2=OFF             \
      -DCMAKE_INSTALL_PREFIX=${SW_DIR}/c-blosc-1.21.1
cmake --build ${build_dir}/c-blosc-build --target install --parallel 16
rm -rf ${build_dir}/c-blosc-build

# HDF5 (for openPMD)
if [ -d ${SRC_DIR}/hdf5 ]
then
  cd ${SRC_DIR}/hdf5
  git fetch --prune
  git checkout hdf5-1_14_1-2
  cd -
else
  git clone -b hdf5-1_14_1-2 https://github.com/HDFGroup/hdf5.git ${SRC_DIR}/hdf5
fi
cmake -S ${SRC_DIR}/hdf5          \
      -B ${build_dir}/hdf5-build  \
      -DBUILD_TESTING=OFF         \
      -DHDF5_ENABLE_PARALLEL=ON   \
      -DCMAKE_INSTALL_PREFIX=${SW_DIR}/hdf5-1.14.1.2
cmake --build ${build_dir}/hdf5-build --target install --parallel 10
rm -rf ${build_dir}/hdf5-build

# ADIOS2 (for openPMD)
if [ -d ${SRC_DIR}/adios2 ]
then
  cd ${SRC_DIR}/adios2
  git fetch --prune
  git checkout v2.8.3
  cd -
else
  git clone -b v2.8.3 https://github.com/ornladios/ADIOS2.git ${SRC_DIR}/adios2
fi
rm -rf ${SRC_DIR}/adios2-build
cmake -S ${SRC_DIR}/adios2             \
      -B ${build_dir}/adios2-build  \
      -DADIOS2_USE_Blosc=ON            \
      -DADIOS2_USE_Fortran=OFF         \
      -DADIOS2_USE_HDF5=OFF            \
      -DADIOS2_USE_Python=OFF          \
      -DADIOS2_USE_ZeroMQ=OFF          \
      -DCMAKE_INSTALL_PREFIX=${HOME}/sw/lumi/cpu/adios2-2.8.3
cmake --build ${build_dir}/adios2-build --target install -j 16
rm -rf ${build_dir}/adios2-build


# Python ######################################################################
#
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade virtualenv
python3 -m pip cache purge
rm -rf ${SW_DIR}/venvs/impactx-cpu-lumi
python3 -m venv ${SW_DIR}/venvs/impactx-cpu-lumi
source ${SW_DIR}/venvs/impactx-cpu-lumi/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade build
python3 -m pip install --upgrade packaging
python3 -m pip install --upgrade wheel
python3 -m pip install --upgrade setuptools
python3 -m pip install --upgrade cython
python3 -m pip install --upgrade numpy
python3 -m pip install --upgrade pandas
python3 -m pip install --upgrade scipy
MPICC="cc -shared" python3 -m pip install --upgrade mpi4py --no-cache-dir --no-build-isolation --no-binary mpi4py
python3 -m pip install --upgrade openpmd-api
python3 -m pip install --upgrade matplotlib
python3 -m pip install --upgrade yt
# install or update ImpactX dependencies
python3 -m pip install --upgrade -r ${SRC_DIR}/impactx/requirements.txt
# ML & optimization (optional)
#python3 -m pip install --upgrade torch --index-url https://download.pytorch.org/whl/cpu
#python3 -m pip install --upgrade optimas[all]
