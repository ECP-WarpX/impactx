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
sim.set_diags_slice_step_diagnostics(True)

# domain decomposition & space charge mesh
sim.init_grids()
# access distributed particle beam storage
particles = sim.particle_container()

# load a 2 GeV electron beam with an initial
# unnormalized rms emittance of 2 nm
energy_MeV = 2.0e3  # reference energy
charge_C = 0.0  # assign zero weighting to particles
mass_MeV = 0.510998950
qm_qeeV = -1.0e-6/mass_MeV  # charge/mass
npart = 10000  # number of macro particles

distr = distribution.Waterbag(
    sigmaX = 3.9984884770e-5,
    sigmaY = 3.9984884770e-5,
    sigmaT = 1.0e-3,
    sigmaPx = 2.6623538760e-5,
    sigmaPy = 2.6623538760e-5,
    sigmaPt = 2.0e-3,
    muxpx = -0.846574929020762,
    muypy = 0.846574929020762,
    mutpt = 0.0)
sim.add_particles(qm_qeeV, charge_C, distr, npart)

# init reference particle
refPart = RefPart()
# make the next two lines a helper function?
refPart.pt = -energy_MeV / mass_MeV - 1.0
refPart.pz = (refPart.pt**2 - 1.0)**0.5
particles.set_ref_particle(refPart)

# design the accelerator lattice
ns = 25  # number of slices per ds in the element
fodo = [
    elements.Drift(ds=0.25, nslice=ns),
    elements.Quad(ds=1.0, k=1.0, nslice=ns),
    elements.Drift(ds=0.5, nslice=ns),
    elements.Quad(ds=1.0, k=-1.0, nslice=ns),
    elements.Drift(ds=0.25, nslice=ns)
]
# assign a fodo segment
sim.lattice.extend(fodo)

# run simulation
sim.evolve()

# clean shutdown
del sim
amrex.finalize()
