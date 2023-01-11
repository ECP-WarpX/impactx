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

# load a 2.5 MeV proton beam
energy_MeV = 2.5  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_energy_MeV(energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    sigmaX=2.0e-3,
    sigmaY=2.0e-3,
    sigmaT=1.0e-3,
    sigmaPx=3.0e-4,
    sigmaPy=3.0e-4,
    sigmaPt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
beam_monitor = elements.BeamMonitor("beam_monitor", "h5")

# design the accelerator lattice
constEnd = elements.ConstF(ds=0.0025, kx=1.0, ky=1.0, kt=1.0e-12)
nllens = elements.NonlinearLens(knll=2.0e-7, cnll=0.01)
const = elements.ConstF(ds=0.005, kx=1.0, ky=1.0, kt=1.0e-12)

num_lenses = 10
nllens_lattice = [constEnd] + [nllens, const] * (num_lenses - 1) + [nllens, constEnd]

# add elements to the lattice segment
sim.lattice.append(beam_monitor)
sim.lattice.extend(nllens_lattice)
sim.lattice.append(beam_monitor)

# run simulation
sim.evolve()

# clean shutdown
del sim
amrex.finalize()
