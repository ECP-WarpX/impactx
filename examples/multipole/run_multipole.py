#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Ryan Sandberg, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import amrex
from impactx import ImpactX, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 2 GeV electron beam with an initial
# unnormalized rms emittance of  nm
energy_MeV = 2.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used without space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_energy_MeV(energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    sigmaX=4.0e-3,
    sigmaY=4.0e-3,
    sigmaT=1.0e-3,
    sigmaPx=3.0e-4,
    sigmaPy=3.0e-4,
    sigmaPt=2.0e-3,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", "h5")

# design the accelerator lattice
multipole = [
    monitor,
    elements.Multipole(multiple=2, K_normal=3.0, K_skew=0.0),
    elements.Multipole(multiple=3, K_normal=100.0, K_skew=-50.0),
    elements.Multipole(multiple=4, K_normal=65.0, K_skew=6.0),
    monitor,
]
# assign a fodo segment
sim.lattice.extend(multipole)

# run simulation
sim.evolve()

# clean shutdown
del sim
amrex.finalize()
