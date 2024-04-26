#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Ryan Sandberg, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

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
kin_energy_MeV = 2.0e3  # reference energy
bunch_charge_C = 1.0e-8  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Kurth6D(
    lambdaX=1.11e-3,
    lambdaY=1.11e-3,
    lambdaT=3.74036839224568e-4,
    lambdaPx=9.00900900901e-4,
    lambdaPy=9.00900900901e-4,
    lambdaPt=2.6735334467940146e-3,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
constf1 = elements.ConstF(ds=2.0, kx=0.7, ky=0.7, kt=0.7)
drift1 = elements.Drift(ds=1.0)
sim.lattice.extend([monitor, drift1, constf1, drift1, monitor])

# run simulation
sim.evolve()

# clean shutdown
sim.finalize()
