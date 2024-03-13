#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Marco Garten, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from impactx import ImpactX, distribution

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 250 MeV proton beam with an initial
# horizontal rms emittance of 1 um and an
# initial vertical rms emittance of 2 um
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle().load_file("solenoid.madx")

#   particle bunch
distr = distribution.Waterbag(
    lambdaX=1.559531175539e-3,
    lambdaY=2.205510139392e-3,
    lambdaT=1.0e-3,
    lambdaPx=6.41218345413e-4,
    lambdaPy=9.06819680526e-4,
    lambdaPt=1.0e-3,
)
sim.add_particles(bunch_charge_C, distr, npart)

# design the accelerator lattice
sim.lattice.load_file("solenoid.madx", nslice=1)

# run simulation
sim.evolve()

# clean shutdown
sim.finalize()
