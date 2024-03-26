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

# load a 2.5 MeV proton beam
kin_energy_MeV = 2.5  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    lambdaX=2.0e-3,
    lambdaY=2.0e-3,
    lambdaT=1.0e-3,
    lambdaPx=3.0e-4,
    lambdaPy=3.0e-4,
    lambdaPt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")
monitor.nonlinear_lens_invariants = True
monitor.alpha = 0.0
monitor.beta = 1.0
monitor.tn = 0.4
monitor.cn = 0.01

# design the accelerator lattice
constEnd = elements.ConstF(ds=0.05, kx=1.0, ky=1.0, kt=1.0e-12)
nllens = elements.NonlinearLens(knll=4.0e-6, cnll=0.01)
const = elements.ConstF(ds=0.1, kx=1.0, ky=1.0, kt=1.0e-12)

num_lenses = 18
nllens_lattice = (
    [monitor, constEnd]
    + [nllens, const] * (num_lenses - 1)
    + [nllens, constEnd, monitor]
)

# add elements to the lattice segment
sim.lattice.extend(nllens_lattice)

# run simulation
sim.evolve()

# clean shutdown
sim.finalize()
