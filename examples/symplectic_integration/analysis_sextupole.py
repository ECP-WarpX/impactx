#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import numpy as np
import openpmd_api as io
import pandas as pd

# initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial = series.iterations[1].particles["beam"].to_df()
is_double = initial["position_x"].dtype == np.float64
beam_final = series.iterations[last_step].particles["beam"]
final = beam_final.to_df()

# compare number of particles
num_particles = 24
# assert num_particles == len(initial)
# assert num_particles == len(final)

# load particle data
df_initial = pd.read_csv("./initial_coords.csv", sep=" ")
df_final = pd.read_csv("./final_coords.csv", sep=" ")

# compute differences
error_xi = (df_initial["x"] - initial["position_x"]).abs().max()
error_pxi = (df_initial["px"] - initial["momentum_x"]).abs().max()
error_yi = (df_initial["y"] - initial["position_y"]).abs().max()
error_pyi = (df_initial["py"] - initial["momentum_y"]).abs().max()
error_ti = (df_initial["t"] - initial["position_t"]).abs().max()
error_pti = (df_initial["pt"] - initial["momentum_t"]).abs().max()

error_xf = (df_final["x"] - final["position_x"]).abs().max()
error_pxf = (df_final["px"] - final["momentum_x"]).abs().max()
error_yf = (df_final["y"] - final["position_y"]).abs().max()
error_pyf = (df_final["py"] - final["momentum_y"]).abs().max()
error_tf = (df_final["t"] - final["position_t"]).abs().max()
error_ptf = (df_final["pt"] - final["momentum_t"]).abs().max()

xf_max = df_final["x"].abs().max()
pxf_max = df_final["px"].abs().max()
yf_max = df_final["y"].abs().max()
pyf_max = df_final["py"].abs().max()
tf_max = df_final["t"].abs().max()
ptf_max = df_final["pt"].abs().max()

print("Initial beam, maximum absolute difference in each coordinate:")
print("Difference x, px, y, py, t, pt:")
print(error_xi)
print(error_pxi)
print(error_yi)
print(error_pyi)
print(error_ti)
print(error_pti)

atol = 1.0e-13 if is_double else 1.0e-9
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
print("Final beam, maximum relative difference in transverse coordinates:")
print("Difference x, px, y, py:")
print(error_xf / xf_max)
print(error_pxf / pxf_max)
print(error_yf / yf_max)
print(error_pyf / pyf_max)

atol = 1.0e-7 if is_double else 1.0e-6
print(f"  tol={atol}")

assert np.allclose(
    [error_xf / xf_max, error_pxf / pxf_max, error_yf / yf_max, error_pyf / pyf_max],
    [
        0.0,
        0.0,
        0.0,
        0.0,
    ],
    atol=atol,
)

print("")
print("Final beam, maximum absolute difference in longitudinal coordinates:")
print("Difference t, pt:")
print(error_tf)
print(error_ptf)

atol = 1.0e-13
if not is_double:
    # FIXME: in single precision, the final t deviates by 3.5e-8 from the reference,
    # which itself is below 7.6e-12
    # https://github.com/BLAST-ImpactX/impactx/issues/1683
    atol = 1.0e-7
print(f"  atol={atol}")

assert np.allclose(
    [error_tf, error_ptf],
    [
        0.0,
        0.0,
    ],
    atol=atol,
)
