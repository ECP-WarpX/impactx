#!/usr/bin/env python3
#
# Copyright 2022 ImpactX contributors
# Authors: Ryan Sandberg, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import amrex
from impactx import ImpactX, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.set_particle_shape(2)  # B-spline order
sim.set_space_charge(False)
#sim.set_diagnostics(False)  # benchmarking
sim.set_slice_step_diagnostics(True)

# domain decomposition & space charge mesh
sim.init_grids()

# load a 2.5 MeV proton beam
energy_MeV = 2.5  # reference energy
charge_C = 1.0e-9  # used with space charge
mass_MeV = 938.27208816
qm_qeeV = 1.0e-6/mass_MeV  # charge/mass
npart = 10000  # number of macro particles

distr = distribution.Waterbag(
    sigmaX = 2.0e-3,
    sigmaY = 2.0e-3,
    sigmaT = 1.0e-3,
    sigmaPx = 3.0e-4,
    sigmaPy = 3.0e-4,
    sigmaPt = 0.0,
    muxpx = 0.0,
    muypy = 0.0,
    mutpt = 0.0)
sim.add_particles(qm_qeeV, charge_C, distr, npart)

# set the energy in the reference particle
sim.particle_container().ref_particle() \
    .set_energy_MeV(energy_MeV, mass_MeV)

constEnd = elements.ConstF(
    ds = 0.0025,
    kx = 1.0,
    ky = 1.0,
    kt = 1.0e-12)
nllens = elements.NonlinearLens(
    knll = 2.0e-7,
    cnll = 0.01)
const = elements.ConstF(
    ds = 0.005,
    kx = 1.0,
    ky = 1.0,
    kt = 1.0e-12)
# design the accelerator lattice
num_lenses = 10
nllens_lattice = [constEnd] + [nllens, const] * (num_lenses-1) + [nllens, constEnd]

# assign a fodo segment
sim.lattice.extend(nllens_lattice)

# run simulation
sim.evolve()

# clean shutdown
del sim
amrex.finalize()
