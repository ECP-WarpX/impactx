#!/usr/bin/env bash
#
# Copyright 2021-2023 The ImpactX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

sudo apt-get -qqq update
sudo apt-get install -y \
    build-essential     \
    ca-certificates     \
    ccache              \
    clang               \
    clang-tidy          \
    cmake               \
    gnupg               \
    libfftw3-dev        \
    libfftw3-mpi-dev    \
    libhdf5-openmpi-dev \
    libomp-dev          \
    libopenmpi-dev      \
    ninja-build         \
    pkg-config          \
    python3             \
    python3-pip         \
    wget

python3 -m pip install -U pip
python3 -m pip install -U build packaging setuptools wheel
python3 -m pip install -U cmake
python3 -m pip install -U -r requirements_mpi.txt
python3 -m pip install -U -r src/python/impactx/dashboard/requirements.txt
python3 -m pip install -U -r examples/requirements.txt
python3 -m pip install -U -r tests/python/requirements.txt

# cmake-easyinstall
#
sudo curl -L -o /usr/local/bin/cmake-easyinstall https://raw.githubusercontent.com/ax3l/cmake-easyinstall/main/cmake-easyinstall
sudo chmod a+x /usr/local/bin/cmake-easyinstall
export CEI_SUDO="sudo"
export CEI_TMP="/tmp/cei"
