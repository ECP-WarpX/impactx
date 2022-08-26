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

# load a 2 GeV electron beam with an initial
# unnormalized rms emittance of  nm
energy_MeV = 2.0e3  # reference energy
charge_C = 1.e-9  # used without space charge
mass_MeV = 0.510998950
qm_qeeV = -1.0e-6/mass_MeV  # charge/mass
npart = 10000  # number of macro particles

distr = distribution.Waterbag(
    sigmaX = 4.0e-3,
    sigmaY = 4.0e-3,
    sigmaT = 1.0e-3,
    sigmaPx = 3.0e-4,
    sigmaPy = 3.0e-4,
    sigmaPt = 2.0e-3)
sim.add_particles(qm_qeeV, charge_C, distr, npart)

# set the energy in the reference particle
sim.particle_container().ref_particle() \
    .set_energy_MeV(energy_MeV, mass_MeV)

# design the accelerator lattice
multipole = [
    elements.Multipole(multiple=2, K_normal=3.0, K_skew=0.0),
    elements.Multipole(multiple=3, K_normal=100.0, K_skew=-50.0),
    elements.Multipole(multiple=4, K_normal=65.0, K_skew=6.0)
]
# assign a fodo segment
sim.lattice.extend(multipole)

# run simulation
sim.evolve()

# clean shutdown
del sim
amrex.finalize()
