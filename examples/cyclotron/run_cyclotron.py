#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
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

# load initial beam
kin_energy_MeV = 4.0e-3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    sigmaX=1.0e-3,
    sigmaY=1.0e-3,
    sigmaT=0.3,
    sigmaPx=2.0e-4,
    sigmaPy=2.0e-4,
    sigmaPt=2.0e-5,
    muxpx=-0.0,
    muypy=0.0,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice)
ns = 1  # number of slices per ds in the element
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
sim.periods = 150

# run simulation
sim.evolve()

# clean shutdown
sim.finalize()
