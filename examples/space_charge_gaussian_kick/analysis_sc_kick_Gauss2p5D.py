#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import numpy as np
import openpmd_api as io
import scipy.constants as sc
from scipy.special import expi

# initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
beam_initial = series.iterations[1].particles["beam"]
initial_sort = beam_initial.to_df().set_index("id")
is_double = initial_sort["position_x"].dtype == np.float64
beam_final = series.iterations[last_step].particles["beam"]
final_sort = beam_final.to_df().set_index("id")

# Physical constants
qe = sc.elementary_charge
me = sc.electron_mass
clite = sc.speed_of_light
eps0 = sc.epsilon_0
me_eV = me * clite**2 / qe

# Basic beam parameters
bunch_charge_C = 1.0e-9
kin_energy_MeV = 100
gamma = 1.0e6 * kin_energy_MeV / me_eV + 1.0
beta_gamma_2 = gamma**2 - 1.0
beta = np.sqrt(beta_gamma_2) / gamma

# Space charge intensity parameter (units of length, meters)
Kscale = (
    qe * bunch_charge_C / (4.0 * np.pi * eps0 * me * clite**2 * beta_gamma_2 * gamma)
)

# Beam distribution parameters
sigmax = 4.0e-5
sigmay = 4.0e-5
sigmact = 1.0e-3
sigmaz = beta * sigmact

# Simulation slice length
L = 1.0

# Longitudinal kick scale parameter
gauss_long_scale = 6.0

# Initial particle data

xi = initial_sort["position_x"]
pxi = initial_sort["momentum_x"]
yi = initial_sort["position_y"]
pyi = initial_sort["momentum_y"]
ti = initial_sort["position_t"]
pti = initial_sort["momentum_t"]
zi = -beta * ti

# Predicted momentum kick, from J. Qiang, Phys. Rev. Accel. Beams 28, 114602 (2025), eqs. (31-32)

ri_2 = xi**2 + yi**2
gauss_exp = np.exp(-ri_2 / (2.0 * sigmax**2))
gauss_xy_factor = 2.0 * (1.0 - gauss_exp) / ri_2
gauss_pdf_z = np.exp(-(zi**2) / (2.0 * sigmaz**2)) / (np.sqrt(2.0 * np.pi) * sigmaz)
px_predicted = L * Kscale * gauss_pdf_z * xi * gauss_xy_factor
py_predicted = L * Kscale * gauss_pdf_z * yi * gauss_xy_factor

potential_xy_factor = expi(-ri_2 / (2.0 * sigmax**2)) + np.log(
    gauss_long_scale**2 / ri_2
)
d_gauss_pdf_z = -zi / sigmaz**2 * gauss_pdf_z

pz_predicted = -L * Kscale * potential_xy_factor * d_gauss_pdf_z
pt_predicted = -beta * pz_predicted

# Maximum momentum kick
px_max = px_predicted.abs().max()
py_max = py_predicted.abs().max()
pt_max = pt_predicted.abs().max()

# Final particle data

xf = final_sort["position_x"]
pxf = final_sort["momentum_x"]
yf = final_sort["position_y"]
pyf = final_sort["momentum_y"]
tf = final_sort["position_t"]
ptf = final_sort["momentum_t"]

# Difference between value and prediction

dpx_max = (pxf - px_predicted).abs().max()
dpy_max = (pyf - py_predicted).abs().max()
dpt_max = (ptf - pt_predicted).abs().max()

dpx_rms = np.sqrt(np.mean(np.square(pxf - px_predicted)))
dpy_rms = np.sqrt(np.mean(np.square(pyf - py_predicted)))
dpt_rms = np.sqrt(np.mean(np.square(ptf - pt_predicted)))

print()
print("Difference between predicted and computed final momentum, absolute rms:")
print("dpx_rms", dpx_rms)
print("dpy_rms", dpy_rms)
print("dpt_rms", dpt_rms)

print()
print("Difference between predicted and computed final momentum, absolute max:")
print("dpx_max", dpx_max)
print("dpy_max", dpy_max)
print("dpt_max", dpt_max)

print()
print("Difference between predicted and computed final momentum (rms), relative:")
print("dpx_rms/px_max", dpx_rms / px_max)
print("dpy_rms/py_max", dpy_rms / py_max)
print("dpt_rms/pt_max", dpt_rms / pt_max)

# Test maximum error:
atol = 5.1e-2
print(f"  tol={atol}")

assert np.allclose(
    [dpx_rms / px_max, dpy_rms / py_max],
    [0.0, 0.0],
    atol=atol,
)

# Longitudinal kick is sensitive to noise (relax tolerance):
atol = 0.5
print(f"  tol={atol}")

assert np.allclose(
    [dpt_rms / pt_max],
    [0.0],
    atol=atol,
)


print()
print("Difference between predicted and computed final momentum (max), relative:")
print("dpx_max/px_max", dpx_max / px_max)
print("dpy_max/py_max", dpy_max / py_max)
print("dpt_max/pt_max", dpt_max / pt_max)

# Test maximum error:
atol = 5.1e-2 if is_double else 1.2e-1
print(f"  tol={atol}")

assert np.allclose(
    [dpx_max / px_max, dpy_max / py_max],
    [0.0, 0.0],
    atol=atol,
)

# Longitudinal kick is sensitive to noise (relax tolerance):
atol = 1.1 if is_double else 3.0
print(f"  tol={atol}")

assert np.allclose(
    [dpt_max / pt_max],
    [0.0],
    atol=atol,
)
