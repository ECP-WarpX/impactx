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
    libhdf5-dev         \
    ninja-build         \
    pkg-config          \
    python3             \
    python3-pip         \
    wget

python3 -m pip install -U pip setuptools wheel
python3 -m pip install -U cmake pytest
python3 -m pip install -r examples/requirements.txt
