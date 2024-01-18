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
    cmake               \
    gnupg               \
    libhdf5-openmpi-dev \
    libopenmpi-dev      \
    ninja-build         \
    pkg-config          \
    python3             \
    python3-pip         \
    wget

python3 -m pip install -U pip
python3 -m pip install -U build packaging setuptools wheel
python3 -m pip install -U cmake pytest
python3 -m pip install -U -r requirements_mpi.txt
python3 -m pip install -U -r examples/requirements.txt

python3 -m pip install -U openPMD-validator
