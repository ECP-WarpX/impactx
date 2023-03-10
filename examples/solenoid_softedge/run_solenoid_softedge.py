#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Chad Mitchell, Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import amrex
from impactx import ImpactX, RefPart, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = False

# domain decomposition & space charge mesh
sim.init_grids()

# load a 250 MeV proton beam with an initial
# horizontal rms emittance of 1 um and an
# initial vertical rms emittance of 2 um
energy_MeV = 250.0  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_energy_MeV(energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    sigmaX=1.559531175539e-3,
    sigmaY=2.205510139392e-3,
    sigmaT=1.0e-3,
    sigmaPx=6.41218345413e-4,
    sigmaPy=9.06819680526e-4,
    sigmaPt=1.0e-3,
)
sim.add_particles(bunch_charge_C, distr, npart)

# design the accelerator lattice
sol = elements.SoftSolenoid(
    ds=6.0,
    bscale=1.233482899483985,
    cos_coefficients=[
        0.350807812299706,
        0.323554693720069,
        0.260320578919415,
        0.182848575294969,
        0.106921016050403,
        4.409581845710694e-002,
        -9.416427163897508e-006,
        -2.459452716865687e-002,
        -3.272762575737291e-002,
        -2.936414401076162e-002,
        -1.995780078926890e-002,
        -9.102893342953847e-003,
        -2.456410658713271e-006,
        5.788233017324325e-003,
        8.040408292420691e-003,
        7.480064552867431e-003,
        5.230254569468851e-003,
        2.447614547094685e-003,
        -1.095525090532255e-006,
        -1.614586867387170e-003,
        -2.281365457438345e-003,
        -2.148709081338292e-003,
        -1.522541739363011e-003,
        -7.185505862719508e-004,
        -6.171194824600157e-007,
        4.842109305036943e-004,
        6.874508102002901e-004,
        6.535550288205728e-004,
        4.648795813759210e-004,
        2.216564722797528e-004,
        -4.100982995210341e-007,
        -1.499332112463395e-004,
        -2.151538438342482e-004,
        -2.044590946652016e-004,
        -1.468242784844341e-004,
    ],
    sin_coefficients=[
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    ],
    mapsteps=200,
    nslice=4,
)

sim.lattice.extend([sol])

# run simulation
sim.evolve()

# clean shutdown
del sim
amrex.finalize()
