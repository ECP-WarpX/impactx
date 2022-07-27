#!/usr/bin/env python3
#
# Copyright 2022 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import amrex
from impactx import ImpactX, RefPart, distribution, elements

impactX = ImpactX()

impactX.set_particle_shape(2)
impactX.set_diags_slice_step_diagnostics(True)
impactX.init_grids()

# init particle beam
energy_MeV = 2.0e3
charge_C = 0.0  # assign zero weighting to particles
mass_MeV = 0.510998950
qm_qeeV = -1.0/0.510998950e6
npart = 10000

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
impactX.add_particles(qm_qeeV, charge_C, distr, npart)

# init reference particle
refPart = RefPart()
# make the next two lines a helper function?
refPart.pt = -energy_MeV / mass_MeV - 1.0
refPart.pz = (refPart.pt**2 - 1.0)**0.5
impactX.particle_container().set_ref_particle(refPart)

# init accelerator lattice
ns = 25  # number of slices per ds in the element
fodo = [
    elements.Drift(ds=0.25, nslice=ns),
    elements.Quad(ds=1.0, k=1.0, nslice=ns),
    elements.Drift(ds=0.5, nslice=ns),
    elements.Quad(ds=1.0, k=-1.0, nslice=ns),
    elements.Drift(ds=0.25, nslice=ns)
]
#  assign a fodo segment
impactX.lattice.extend(fodo)

# run simulation
impactX.evolve()

# clean shutdown
del impactX
amrex.finalize()
