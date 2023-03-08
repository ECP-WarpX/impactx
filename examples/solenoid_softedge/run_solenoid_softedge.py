#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Marco Garten, Axel Huebl, Chad Mitchell
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
sim.slice_step_diagnostics = True

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
sim.lattice.append(elements.SoftSolenoid(
    ds=6.0,
    bscale=1.233482899483985,
    cos_coefficients=[
             0.350807812299706,
             0.323554693720069,
             0.260320578919415,
             0.182848575294969,
             0.106921016050403,
             4.409581845710694E-002,
            -9.416427163897508E-006,
            -2.459452716865687E-002,
            -3.272762575737291E-002,
            -2.936414401076162E-002,
            -1.995780078926890E-002,
            -9.102893342953847E-003,
            -2.456410658713271E-006,
             5.788233017324325E-003,
             8.040408292420691E-003,
             7.480064552867431E-003,
             5.230254569468851E-003,
             2.447614547094685E-003,
            -1.095525090532255E-006,
            -1.614586867387170E-003,
            -2.281365457438345E-003,
            -2.148709081338292E-003,
            -1.522541739363011E-003,
            -7.185505862719508E-004,
            -6.171194824600157E-007,
             4.842109305036943E-004,
             6.874508102002901E-004,
             6.535550288205728E-004,
             4.648795813759210E-004,
             2.216564722797528E-004,
            -4.100982995210341E-007,
            -1.499332112463395E-004,
            -2.151538438342482E-004,
            -2.044590946652016E-004,
            -1.468242784844341E-004
    ],
    sin_coefficients=[
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0
    ],
    mapsteps=200,
    nslice=4)
)

# run simulation
sim.evolve()

# clean shutdown
del sim
amrex.finalize()
