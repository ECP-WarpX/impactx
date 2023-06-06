#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import amrex.space3d as amr
from impactx import ImpactX, RefPart, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load initial beam
energy_MeV = 4.0e-3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_energy_MeV(energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    sigmaX=3.9984884770e-5,
    sigmaY=3.9984884770e-5,
    sigmaT=1.0e-3,
    sigmaPx=2.6623538760e-5,
    sigmaPy=2.6623538760e-5,
    sigmaPt=2.0e-5,
    muxpx=-0.846574929020762,
    muypy=0.846574929020762,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice)
ns = 25  # number of slices per ds in the element
period = [
    monitor,
    elements.ChrAcc(ds=0.038, ez=1.12188308693e-4, bz=1.0e-14, nslice=ns),
    elements.ExactSbend(ds=0.25, phi=180.0, B=1),
    elements.ChrAcc(ds=0.038, ez=1.12188308693e-4, bz=1.0e-14, nslice=ns),
    elements.ExactSbend(ds=0.25, phi=180.0, B=1),
    monitor,
]

sim.lattice.extend(period)

# number of periods through the lattice
sim.periods = 250

# run simulation
sim.evolve()

# clean shutdown
del sim
amr.finalize()
