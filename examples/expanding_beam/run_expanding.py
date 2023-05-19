#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import amrex.space3d as amr
from impactx import ImpactX, RefPart, distribution, elements

pp_amr = amr.ParmParse("amr")
pp_amr.add("max_level", 1)
pp_amr.addarr("n_cell", [56, 56, 48])

sim = ImpactX()

# set numerical parameters and IO control
#sim.n_cell = [56, 56, 48]
sim.particle_shape = 2  # B-spline order
sim.space_charge = True
sim.dynamic_size = True
sim.prob_relative = [3.0, 1.1]

# beam diagnostics
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = False

# domain decomposition & space charge mesh
sim.init_grids()

# load a 2 GeV electron beam with an initial
# unnormalized rms emittance of 2 nm
kin_energy_MeV = 250  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles (outside tests, use 1e5 or more)

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Kurth6D(
    sigmaX=4.472135955e-4,
    sigmaY=4.472135955e-4,
    sigmaT=9.12241869e-7,
    sigmaPx=0.0,
    sigmaPy=0.0,
    sigmaPt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
sim.lattice.extend([monitor, elements.Drift(ds=6.0, nslice=40), monitor])

# run simulation
sim.evolve()

# visualize
import numpy as np

rho = sim.rho(lev=0)
rs = rho.sum_unique(comp=0, local=False)

gm = sim.Geom(lev=0)
dr = gm.data().CellSize()
dV = np.prod(dr)

half_z = sim.n_cell[2] // 2  # order: x,y,z

import matplotlib.pyplot as plt

f = plt.figure()
ax = f.gca()
ng = rho.nGrowVect
for mfi in rho:
    bx = mfi.validbox()
    rbx = amr.RealBox(bx, dr, gm.ProbLo())

    arr = rho.array(mfi)
    arr_np = np.array(arr, copy=False)  # indices: comp, z, y, x

    # shift box to zero-based local mfi index space
    half_z_local = half_z - bx.lo_vect[2]
    bx.shift(bx.lo_vect * -1)
    # check if the current tile contains the half-z plane
    if half_z_local < 0 or half_z_local > arr_np.shape[2]:
        continue

    comp = 0
    mu = 1.0e6  # m->mu
    im = ax.imshow(
        # arr_np[comp, half_z, ...] * dV,  # including guard
        arr_np[comp, half_z_local, ng[1] : -ng[1], ng[0] : -ng[0]] * dV,  # w/o guard
        origin="lower",
        aspect="auto",
        extent=[rbx.lo(0) * mu, rbx.hi(0) * mu, rbx.lo(1) * mu, rbx.hi(1) * mu],
    )
    cb = f.colorbar(im)
    cb.set_label(r"charge density  [C/m$^3$]")
    ax.set_xlabel(r"$x$  [$\mu$m]")
    ax.set_ylabel(r"$y$  [$\mu$m]")
    save_png = False
    if save_png:
        plt.savefig("charge_deposition.png")
    else:
        plt.show()

# clean shutdown
del sim
amr.finalize()
