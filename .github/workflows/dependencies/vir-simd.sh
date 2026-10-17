#!/usr/bin/env bash
#
# Copyright 2022-2026 The ImpactX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

# vir-simd, needed for ImpactX_SIMD=ON
#   built outside of the source tree, so it cannot end up in a commit
#   pinned to the commit of a release tag, which cannot be moved
vir_simd_sha=f936143148fbd8b9b2967238ab4aa267175f4d31  # v0.4.4
cd "$(mktemp -d)"
wget -O vir-simd.tar.gz https://github.com/mattkretz/vir-simd/archive/${vir_simd_sha}.tar.gz
tar -xf vir-simd.tar.gz
cmake -S vir-simd-${vir_simd_sha} -B vir-simd-build
sudo cmake --build vir-simd-build --target install
