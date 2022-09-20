#!/usr/bin/env python3
#
# Copyright 2022 ImpactX contributors
# Authors: Axel Huebl, Christopher E. Mayes
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np

import amrex
from impactx import ImpactX, RefPart, distribution, elements

pp_amr = amrex.ParmParse("amr")
# pp_amr.addarr("n_cell", [128, 128, 128])
pp_amr.addarr("n_cell", [40, 40, 32])

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = True
sim.dynamic_size = True
sim.prob_relative = 1.0

# beam diagnostics
sim.diagnostics = False
sim.slice_step_diagnostics = False

# domain decomposition & space charge mesh
sim.init_grids()

# load a 10 MeV electron beam with an initial
# unnormalized rms emittance of 2 nm
energy_MeV = 10  # reference energy
mass_MeV = 0.510998950  # electron mass in MeV/c^2
bunch_charge_C = 1.0e-9  # used with space charge
npart = int(10000)  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(mass_MeV).set_energy_MeV(energy_MeV)

#   particle bunch
r = 1.0  # aspect ratio = sigma_z / sigma_perp = 0.01 - 10
sigma_z = 1.0e-6  # mm
sigma_r = sigma_z / r
gamma = energy_MeV / mass_MeV
beta = (1.0 - (1.0 / gamma) ** 2) ** 0.5
c0 = 2.99792458e8  # speed of light in m/s
sigma_t = sigma_z / (c0 * beta)
print(f"sigma_t={sigma_t}s")

distr = distribution.Waterbag(
    sigmaX=sigma_r,
    sigmaY=sigma_r,
    sigmaT=sigma_t,
    sigmaPx=0.0,
    sigmaPy=0.0,
    sigmaPt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# design the accelerator lattice
sim.lattice.append(elements.Drift(ds=1.0, nslice=30))

# run simulation
sim.evolve()

# calculate phase space
#   hack:
import matplotlib.pyplot as plt

num_plots_per_row = 2
f, axs = plt.subplots(1, num_plots_per_row, figsize=(7, 2))
ax_x_px = axs[0]
ax_z_pz = axs[1]

pc = sim.particle_container()
lev = pc.GetParticles(0)
for tile_ind, pt in lev.items():
    # positions + id + cpuid
    aos = pt.GetArrayOfStructs()
    aos_arr = np.array(aos, copy=False)

    # momentum & particle weight
    real_arrays = pt.GetStructOfArrays().GetRealData()
    px = np.array(real_arrays[0], copy=False)
    pz = np.array(real_arrays[2], copy=False)

    print(f"tile_ind={tile_ind}, pt={pt}")
    print(f"aos_arr={aos_arr}, aos_arr.shape={aos_arr.shape}")
    print(f"px={px}, px.shape={px.shape}")

    ax_x_px.scatter(aos_arr[()]["x"], px)
    ax_z_pz.scatter(aos_arr[()]["z"], pz)
plt.show()

# MPI reduce phase space
#   TODO SUM up histograms via MPI

# plot phase space
# import matplotlib.pyplot as plt
# f = plt.figure()
# ax = plt.gca()
# ax.scatter()

# clean shutdown
#   note: timers turned stop here
del sim
amrex.finalize()
