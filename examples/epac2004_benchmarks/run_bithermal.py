#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Marco Garten, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import amrex.space3d as amr
from impactx import ImpactX, RefPart, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.n_cell = [56, 56, 64]
sim.particle_shape = 2  # B-spline order
sim.space_charge = True
sim.dynamic_size = True
sim.prob_relative = 4.0

# beam diagnostics
sim.slice_step_diagnostics = False

# domain decomposition & space charge mesh
sim.init_grids()

# beam parameters
energy_MeV = 0.1  # reference energy
bunch_charge_C = 1.4285714285714285714e-10  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_energy_MeV(energy_MeV)

#   particle bunch
distr = distribution.Bithermal(
    k=6.283185307179586,
    kT=36.0e-6,
    kT_halo=900.0e-6,
    normalize=0.54226,
    normalize_halo=0.08195,
    w_halo=0.05,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
sim.lattice.append(monitor)

constf = elements.Constf(
    ds=10.0,
    kx=6.283185307179586,
    ky=6.283185307179586,
    kt=6.283185307179586,
    nslice=400,
)

sim.lattice.append(constf)
sim.lattice.append(monitor)

# run simulation
sim.evolve()

# clean shutdown
del sim
amr.finalize()
