#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import numpy as np
import openpmd_api as io

# initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial = series.iterations[1].particles["beam"].to_df()
is_double = initial["position_x"].dtype == np.float64
beam_final = series.iterations[last_step].particles["beam"]
final = beam_final.to_df()

# compare number of particles
num_particles = 1
assert num_particles == len(initial)
assert num_particles == len(final)

# compute differences
error_xi = initial["position_x"].abs().max()
error_pxi = initial["momentum_x"].abs().max()
error_yi = initial["position_y"].abs().max()
error_pyi = initial["momentum_y"].abs().max()
error_ti = initial["position_t"].abs().max()
error_pti = initial["momentum_t"].abs().max()

error_xf = final["position_x"].abs().max()
error_pxf = final["momentum_x"].abs().max()
error_yf = final["position_y"].abs().max()
error_pyf = final["momentum_y"].abs().max()
error_tf = final["position_t"].abs().max()
error_ptf = final["momentum_t"].abs().max()

print("Initial beam, maximum absolute difference in each coordinate:")
print("Difference x, px, y, py, t, pt:")
print(error_xi)
print(error_pxi)
print(error_yi)
print(error_pyi)
print(error_ti)
print(error_pti)

atol = 1.0e-13 if is_double else 3.0e-7
print(f"  atol={atol}")

assert np.allclose(
    [error_xi, error_pxi, error_yi, error_pyi, error_ti, error_pti],
    [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
    atol=atol,
)

print("")
print("Final beam, maximum absolute difference in each coordinate:")
print("Difference x, px, y, py, t, pt:")
print(error_xf)
print(error_pxf)
print(error_yf)
print(error_pyf)
print(error_tf)
print(error_ptf)

atol = 1.0e-13 if is_double else 3.0e-7
print(f"  atol={atol}")

assert np.allclose(
    [error_xf, error_pxf, error_yf, error_pyf, error_tf, error_ptf],
    [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
    atol=atol,
)
