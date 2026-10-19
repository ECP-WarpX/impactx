#!/bin/bash
#
# Copyright 2023 The ImpactX Community
#
# This file is part of ImpactX.
#
# Author: Axel Huebl
# License: BSD-3-Clause-LBNL

# Exit on first error encountered #############################################
#
set -eu -o pipefail


# Check: ######################################################################
#
#   Was perlmutter_gpu_impactx.profile sourced and configured correctly?
if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your perlmutter_gpu_impactx.profile file! Please edit its line 2 to continue!"; exit 1; fi


# Check $proj variable is correct and has a corresponding CFS directory #######
#
if [ ! -d "${CFS}/${proj}/" ]
then
    echo "WARNING: The directory ${CFS}/${proj}/ does not exist!"
    echo "Is the \$proj environment variable of value \"$proj\" correctly set? "
    echo "Please edit line 2 of your perlmutter_gpu_impactx.profile file to continue!"
    exit
fi


# Remove old dependencies #####################################################
#
SW_DIR="${CFS}/${proj}/${USER}/sw/perlmutter/gpu"
rm -rf ${SW_DIR}
mkdir -p ${SW_DIR}

# remove common user mistakes in python, located in .local instead of a venv
python3 -m pip uninstall -qq -y impactx
python3 -m pip uninstall -qqq -y mpi4py 2>/dev/null || true

# clean out caches, e.g., depending on old system modules
python3 -m pip cache purge || true
rm -rf ${HOME}/.cupy/kernel_cache ${HOME}/.nv/ComputeCache ${HOME}/.cache/numba ${HOME}/.triton


# General extra dependencies ##################################################
#

# build parallelism
PARALLEL=16

# tmpfs build directory: avoids issues often seen with $HOME and is faster
build_dir=$(mktemp -d)

# CCache (compiler cache)
curl -Lo ${build_dir}/ccache.tar.xz https://github.com/ccache/ccache/releases/download/v4.10.2/ccache-4.10.2-linux-x86_64.tar.xz
tar -xf ${build_dir}/ccache.tar.xz -C ${build_dir}
mv ${build_dir}/ccache-4.10.2-linux-x86_64 ${SW_DIR}/ccache-4.10.2

# Boost (header-only libraries, for the synmadx MAD-X lattice parser)
curl -Lo ${build_dir}/boost.tar.gz https://archives.boost.io/release/1.82.0/source/boost_1_82_0.tar.gz
tar -xzf ${build_dir}/boost.tar.gz -C ${build_dir}
rm -rf ${build_dir}/boost.tar.gz
(
  cd ${build_dir}/boost_1_82_0
  ./bootstrap.sh --with-libraries=headers --prefix=${SW_DIR}/boost-1.82.0
  ./b2 cxxflags="-std=c++20" install -j ${PARALLEL}
)
rm -rf ${build_dir}/boost_1_82_0

# c-blosc2 (I/O compression)
if [ -d $HOME/src/c-blosc2 ]
then
  cd $HOME/src/c-blosc2
  git fetch --prune --tags
  git checkout v2.23.1
  cd -
else
  git clone -b v2.23.1 https://github.com/Blosc/c-blosc2.git $HOME/src/c-blosc2
fi
cmake -S $HOME/src/c-blosc2 -B ${build_dir}/c-blosc2-pm-gpu-build -DBUILD_STATIC=OFF -DBUILD_TESTS=OFF -DBUILD_FUZZERS=OFF -DBUILD_BENCHMARKS=OFF -DBUILD_EXAMPLES=OFF -DDEACTIVATE_AVX2=OFF -DCMAKE_INSTALL_PREFIX=${SW_DIR}/c-blosc2-2.23.1
cmake --build ${build_dir}/c-blosc2-pm-gpu-build --target install --parallel ${PARALLEL}
rm -rf ${build_dir}/c-blosc2-pm-gpu-build

# ADIOS2
if [ -d $HOME/src/adios2 ]
then
  cd $HOME/src/adios2
  git fetch --prune --tags
  git checkout v2.12.1
  cd -
else
  git clone -b v2.12.1 https://github.com/ornladios/ADIOS2.git $HOME/src/adios2
fi
cmake -S $HOME/src/adios2 -B ${build_dir}/adios2-pm-gpu-build -DADIOS2_USE_Blosc2=ON -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_Python=OFF -DADIOS2_USE_ZeroMQ=OFF -DCMAKE_INSTALL_PREFIX=${SW_DIR}/adios2-2.12.1
cmake --build ${build_dir}/adios2-pm-gpu-build --target install -j ${PARALLEL}
rm -rf ${build_dir}/adios2-pm-gpu-build


# Python ######################################################################
#
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade virtualenv
rm -rf ${SW_DIR}/venvs/impactx
python3 -m venv ${SW_DIR}/venvs/impactx
source ${SW_DIR}/venvs/impactx/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install --upgrade build
python3 -m pip install --upgrade packaging
python3 -m pip install --upgrade wheel
python3 -m pip install --upgrade setuptools[core]
python3 -m pip install --upgrade cython
python3 -m pip install --upgrade numpy
python3 -m pip install --upgrade pandas
python3 -m pip install --upgrade scipy
MPICC="cc -target-accel=nvidia80 -shared" python3 -m pip install --upgrade mpi4py --no-cache-dir --no-build-isolation --no-binary mpi4py
python3 -m pip install --upgrade openpmd-api
python3 -m pip install --upgrade matplotlib
python3 -m pip install --upgrade yt
# install or update impactx dependencies
python3 -m pip install --upgrade -r $HOME/src/impactx/requirements.txt
python3 -m pip install --upgrade cupy-cuda13x  # CUDA 13 compatible wheel
python3 -m pip install --upgrade torch --index-url https://download.pytorch.org/whl/cu132  # CUDA 13.2 compatible wheel


# remove build temporary directory
rm -rf ${build_dir}
