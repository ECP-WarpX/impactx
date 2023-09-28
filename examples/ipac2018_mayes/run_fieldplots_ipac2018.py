#!/usr/bin/env python3
#
# Copyright 2022 ImpactX contributors
# Authors: Axel Huebl, Christopher E. Mayes
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

import amrex.space3d as amr
from impactx import ImpactX, RefPart, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.n_cell = [128, 128, 256]  # 128^3 or higher for reasonable convergence
sim.particle_shape = 2  # B-spline order
sim.space_charge = True
sim.dynamic_size = True
sim.prob_relative = 3.0

# beam diagnostics
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = False

# domain decomposition & space charge mesh
sim.init_grids()

# load a cold 10 MeV electron beam
mass_MeV = 0.510998950  # electron mass in MeV/c^2
energy_MeV = 10.0  # reference energy (total)
bunch_charge_C = 1.0e-9  # charge in C
npart = int(1000000)  # number of macro particles

#   reference particle (Mayes Fig. 1 is for + charge)
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(mass_MeV).set_energy_MeV(energy_MeV)

#   particle bunch
r = 1.0  # aspect ratio = sigma_z / sigma_perp:  range 0.01 to 10
sigma_r = 1.0e-3  # fixed at 1 mm
sigma_z = r * sigma_r
gamma = energy_MeV / mass_MeV
beta = (1.0 - (1.0 / gamma) ** 2) ** 0.5
print(gamma, beta)
c0 = 2.99792458e8  # speed of light in m/s
sigma_t = sigma_z / (beta * gamma)  # recall t is implicitly scaled by c0
print(f"sigma_t={sigma_t}m")

distr = distribution.Gaussian(
    sigmaX=sigma_r,
    sigmaY=sigma_r,
    sigmaT=sigma_t,
    sigmaPx=0.0,
    sigmaPy=0.0,
    sigmaPt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
sim.lattice.extend([monitor])

# run simulation
sim.evolve()

# plotting

# theory data from eq. (1) in Mayes, sampled to text files in
txt = np.loadtxt("Ex_Mayes.dat")
x_theory, E_x_theory = txt[:, 0], txt[:, 3]

# simulation data
F_x = sim.space_charge_field(lev=0, comp="x")  # [V/m] in the lab frame

gm = sim.Geom(lev=0)
dr = gm.data().CellSize()
dV = np.prod(dr)

half_x, half_y, half_z = [n // 2 for n in sim.n_cell]  # order: x,y,z

# plot data slices
f = plt.figure()
ax = f.gca()
ng = F_x.nGrowVect
q_e = -1.602176634e-19
for mfi in F_x:
    bx = mfi.validbox()
    rbx = amr.RealBox(bx, dr, gm.ProbLo())

    arr_np = F_x.array(mfi).to_numpy(copy=True)  # indices: x, y, z, comp

    # shift box to zero-based local mfi index space
    half_y_local = half_y - bx.lo_vect[1]
    half_z_local = half_z - bx.lo_vect[2]
    bx.shift(bx.lo_vect * -1)

    # check if the current tile contains the half-y plane
    if half_y_local < 0 or half_y_local > arr_np.shape[1]:
        continue
    # check if the current tile contains the half-z plane
    if half_z_local < 0 or half_z_local > arr_np.shape[2]:
        continue

    comp = 0
    lineout = arr_np[ng[0] : -ng[0], half_y_local, half_z_local, comp]  # w/o guard, N
    E_x = lineout / gamma / 1.0e6  # E_x(x[m]) [MV/m]
    ax.plot(
        np.linspace(rbx.lo(0) / sigma_r, rbx.hi(0) / sigma_r, E_x.shape[0]),
        E_x,
        "ko",
        markersize=3,
    )
xmin = -6.0
xmax = 6.0
ymin = -5.0
ymax = 5.0
ax.set_xlim([xmin, xmax])
ax.set_ylim([ymin, ymax])
ax.plot(x_theory, E_x_theory, label="theory")  # E_x(sigma_r) [MV/m]
# cb = f.colorbar(im)
# cb.set_label(r"charge density  [C/m$^3$]")
ax.set_xlabel(r"$x$  [$\sigma_x$]")
ax.set_ylabel(r"$E_x$  [MV/m]")
ax.legend()
# plt.savefig("charge_deposition.png")
plt.show()

# clean shutdown
#   note: timers turned stop here
del sim
amr.finalize()
