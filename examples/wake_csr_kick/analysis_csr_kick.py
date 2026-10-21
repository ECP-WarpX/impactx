#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import numpy as np
import openpmd_api as io
import scipy.constants as sc
import scipy.special as sp

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
r_e = 1.0 / (4.0 * np.pi * eps0) * qe**2 / (me * clite**2)

# Basic beam parameters
bunch_charge_C = 1.0e-9
kin_energy_MeV = 5.0e3
gamma = 1.0e6 * kin_energy_MeV / me_eV + 1.0
beta_gamma_2 = gamma**2 - 1.0
beta = np.sqrt(beta_gamma_2) / gamma
Nb = bunch_charge_C / qe
p0_eVc = me_eV * beta * gamma

# Bend parameters
radius = 10.3462283686195526

# Beam distribution parameters
sigmact = 5.5555553e-8
sigmaz = beta * sigmact

# Simulation slice length
L = 0.500194828041958

# CSR intensity parameter (units of length, meters)
kappa_eV_m13 = 2.0 * r_e * me_eV / (3.0 ** (1.0 / 3.0) * radius ** (2.0 / 3.0))
Kscale_eV_per_m = Nb * kappa_eV_m13 / sigmaz ** (4.0 / 3.0)

# Initial particle data

ti = initial_sort["position_t"]
pti = initial_sort["momentum_t"]
zi = -beta * ti
zi_scaled = zi / sigmaz

# Predicted momentum kick, from J. Qiang, Phys. Rev. Accel. Beams 28, 114602 (2025), eqs. (31-32)
Gamma1 = sp.gamma(5.0 / 6.0)
Gamma2 = sp.gamma(4.0 / 3.0)
zeta = -(zi_scaled**2) / 2.0
Hyperg1 = sp.hyp1f1(2.0 / 3.0, 1.0 / 2.0, zeta)
Hyperg2 = sp.hyp1f1(7.0 / 6.0, 3.0 / 2.0, zeta)
wake = (
    Kscale_eV_per_m
    / (4.0 * np.sqrt(np.pi))
    * (
        -(2.0 ** (4.0 / 3.0)) * Gamma1 * Hyperg1
        + 2.0 ** (5.0 / 6.0) * zi_scaled * Gamma2 * Hyperg2
    )
)

dE_predicted = wake * L / p0_eVc
pt_predicted = -dE_predicted

# Maximum momentum kick
pt_max = pt_predicted.abs().max()

# Final particle data

tf = final_sort["position_t"]
ptf = final_sort["momentum_t"]
dEf = -ptf


# Difference between value and prediction

dpt_max = (ptf - pt_predicted).abs().max()
dpt_rms = np.sqrt(np.mean(np.square(ptf - pt_predicted)))

print()
print("Difference between predicted and computed final momentum, absolute rms:")
print("dpt_rms", dpt_rms)

print()
print("Difference between predicted and computed final momentum, absolute max:")
print("dpt_max", dpt_max)

print()
print("Difference between predicted and computed final momentum (rms), relative:")
print("dpt_rms/pt_max", dpt_rms / pt_max)

# Test maximum error:
# Longitudinal kick is sensitive to noise (relax tolerance):
atol = 0.02
print(f"  tol={atol}")

assert np.allclose(
    [dpt_rms / pt_max],
    [0.0],
    atol=atol,
)


print()
print("Difference between predicted and computed final momentum (max), relative:")
print("dpt_max/pt_max", dpt_max / pt_max)

# Test maximum error:
# Longitudinal kick is sensitive to noise (relax tolerance):
atol = 0.04 if is_double else 3.0
print(f"  tol={atol}")

assert np.allclose(
    [dpt_max / pt_max],
    [0.0],
    atol=atol,
)
