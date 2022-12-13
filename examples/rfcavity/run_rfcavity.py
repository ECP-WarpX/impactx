#!/usr/bin/env python3
#
# Copyright 2022 ImpactX contributors
# Authors: Marco Garten, Axel Huebl, Chad Mitchell
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
sim.slice_step_diagnostics = False

# domain decomposition & space charge mesh
sim.init_grids()

# load a 230 MeV electron beam with an initial
# unnormalized rms emittance of 1 mm-mrad in all
# three phase planes
energy_MeV = 230.0  # reference energy
bunch_charge_C = 1.0e-10  # used with space charge
npart = 10000  # number of macro particles (outside tests, use 1e5 or more)

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_energy_MeV(energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    sigmaX=0.352498964601e-3,
    sigmaY=0.207443478142e-3,
    sigmaT=0.70399950746e-4,
    sigmaPx=5.161852770e-6,
    sigmaPy=9.163582894e-6,
    sigmaPt=0.260528852031e-3,
    muxpx=0.5712386101751441,
    muypy=-0.514495755427526,
    mutpt=-5.05773e-10,
)
sim.add_particles(bunch_charge_C, distr, npart)

# design the accelerator lattice

#   Drift elements
dr1 = elements.Drift(ds=0.4, nslice=1)
dr2 = elements.Drift(ds=0.032997, nslice=1)
#   RF cavity element
rf = elements.RFCavity(
    ds=1.31879807,
    escale=62.0,
    freq=1.3e9,
    phase=85.5,
    mapsteps=100,
    nslice=4,
)

sim.lattice.extend([dr1, dr2, rf, dr2, dr2, rf, dr2, dr2, rf, dr2, dr2, rf, dr2])

# run simulation
sim.evolve()

# clean shutdown
del sim
amrex.finalize()
