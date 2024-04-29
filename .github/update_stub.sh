#!/usr/bin/env bash
#
# Copyright 2021-2024 The ImpactX Community
#
# This script updates the .pyi stub files for documentation and interactive use.
# To run this script, the Python bindings of ImpactX need to be installed
# and importable.
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
set -eu -o pipefail

# we are in the source directory, .github/
this_dir=$(cd $(dirname $0) && pwd)

pybind11-stubgen --exit-code -o ${this_dir}/../src/python/ impactx
