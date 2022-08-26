#!/usr/bin/env python3
#
# Copyright 2022 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import amrex
from impactx import ImpactX, RefPart, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.set_particle_shape(2)  # B-spline order
sim.set_space_charge(False)
# sim.set_diagnostics(False)  # benchmarking
sim.set_slice_step_diagnostics(True)

# domain decomposition & space charge mesh
sim.init_grids()

# load a 2 GeV electron beam with an initial
# unnormalized rms emittance of 2 nm
energy_MeV = 2.0e3  # reference energy
charge_C = 1.0e-9  # used with space charge
mass_MeV = 0.510998950
qm_qeeV = -1.0e-6 / mass_MeV  # charge/mass
npart = 10000  # number of macro particles

distr = distribution.Waterbag(
    sigmaX=3.9984884770e-5,
    sigmaY=3.9984884770e-5,
    sigmaT=1.0e-3,
    sigmaPx=2.6623538760e-5,
    sigmaPy=2.6623538760e-5,
    sigmaPt=2.0e-3,
    muxpx=-0.846574929020762,
    muypy=0.846574929020762,
    mutpt=0.0,
)
sim.add_particles(qm_qeeV, charge_C, distr, npart)

# set the energy in the reference particle
sim.particle_container().ref_particle().set_energy_MeV(energy_MeV, mass_MeV)

# design the accelerator lattice
ns = 25  # number of slices per ds in the element
fodo = [
    elements.Drift(ds=0.25, nslice=ns),
    elements.Quad(ds=1.0, k=1.0, nslice=ns),
    elements.Drift(ds=0.5, nslice=ns),
    elements.Quad(ds=1.0, k=-1.0, nslice=ns),
    elements.Drift(ds=0.25, nslice=ns),
]
# assign a fodo segment
sim.lattice.extend(fodo)

# run simulation
sim.evolve()

# clean shutdown
del sim
amrex.finalize()
