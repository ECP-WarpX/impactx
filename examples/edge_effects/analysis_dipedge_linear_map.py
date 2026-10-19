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
num_particles = 6
assert num_particles == len(initial)
assert num_particles == len(final)

# initial data
xi = np.array([1, 0, 0, 0, 0, 0])
pxi = np.array([0, 1, 0, 0, 0, 0])
yi = np.array([0, 0, 1, 0, 0, 0])
pyi = np.array([0, 0, 0, 1, 0, 0])
ti = np.array([0, 0, 0, 0, 1, 0])
pti = np.array([0, 0, 0, 0, 0, 1])

# problem parameters
rc = 10.0
g = 1.0
h = 1.0/rc
edge = np.pi/8.0
hgap = g/2.0
fint = 0.5

# derived parameters
corr = (h + h) * hgap * fint
tanedg = np.tan(edge)
secedg = 1.0 / np.cos(edge)
psip = edge - corr * secedg * (1.0 + np.sin(edge)**2)

# exact transfer matrix
r21 = h*tanedg
r43 = -h*np.tan(psip)

# final data
xf = np.array([1, 0, 0, 0, 0, 0])
pxf = np.array([r21, 1, 0, 0, 0, 0])
yf = np.array([0, 0, 1, 0, 0, 0])
pyf = np.array([0, 0, r43, 1, 0, 0])
tf = np.array([0, 0, 0, 0, 1.0, 0])
ptf = np.array([0, 0, 0, 0, 0, 1.0])

# compute differences
error_xi = np.max(np.abs(xi - initial["position_x"].to_numpy()))
error_pxi = np.max(np.abs(pxi - initial["momentum_x"].to_numpy()))
error_yi = np.max(np.abs(yi - initial["position_y"].to_numpy()))
error_pyi = np.max(np.abs(pyi - initial["momentum_y"].to_numpy()))
error_ti = np.max(np.abs(ti - initial["position_t"].to_numpy()))
error_pti = np.max(np.abs(pti - initial["momentum_t"].to_numpy()))

error_xf = np.max(np.abs(xf - final["position_x"].to_numpy()))
error_pxf = np.max(np.abs(pxf - final["momentum_x"].to_numpy()))
error_yf = np.max(np.abs(yf - final["position_y"].to_numpy()))
error_pyf = np.max(np.abs(pyf - final["momentum_y"].to_numpy()))
error_tf = np.max(np.abs(tf - final["position_t"].to_numpy()))
error_ptf = np.max(np.abs(ptf - final["momentum_t"].to_numpy()))

xf_max = np.max(np.abs(xf))
pxf_max = np.max(np.abs(pxf))
yf_max = np.max(np.abs(yf))
pyf_max = np.max(np.abs(pyf))
tf_max = np.max(np.abs(tf))
ptf_max = np.max(np.abs(ptf))

print("Initial maximum absolute difference in each column:")
print("Difference x, px, y, py, t, pt:")
print(error_xi)
print(error_pxi)
print(error_yi)
print(error_pyi)
print(error_ti)
print(error_pti)

atol = 1.0e-13
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
print("Final maximum relative difference in each column:")
print("Difference x, px, y, py, t, pt:")
print(error_xf / xf_max)
print(error_pxf / pxf_max)
print(error_yf / yf_max)
print(error_pyf / pyf_max)
print(error_tf / tf_max)
print(error_ptf / ptf_max)

atol = 1.0e-7 if is_double else 1.0e-5
print(f"  tol={atol}")

assert np.allclose(
    [
        error_xf / xf_max,
        error_pxf / pxf_max,
        error_yf / yf_max,
        error_pyf / pyf_max,
        error_tf / tf_max,
        error_ptf / ptf_max,
    ],
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
