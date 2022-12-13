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
sim.particle_shape = 2  # B-spline order
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 2 GeV proton beam with an initial
# unnormalized rms emittance of 1 um in each
# coordinate plane
energy_MeV = 2.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_energy_MeV(energy_MeV)

#   particle bunch
distr = distribution.Kurth6D(
    sigmaX=1.0e-3,
    sigmaY=1.0e-3,
    sigmaT=3.369701494258956e-4,
    sigmaPx=1.0e-3,
    sigmaPy=1.0e-3,
    sigmaPt=2.9676219145931020e-3,
)
sim.add_particles(bunch_charge_C, distr, npart)

# design the accelerator lattice: here we just assign a single element
sim.lattice.append(elements.ConstF(ds=2.0, kx=1.0, ky=1.0, kt=1.0))

# run simulation
sim.evolve()

# clean shutdown
del sim
amrex.finalize()
