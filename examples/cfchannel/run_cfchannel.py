#!/usr/bin/env python3
#
# Copyright 2022 ImpactX contributors
# Authors: Marco Garten, Axel Huebl, Chad Mitchell
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

# load a 2 GeV proton beam with an initial
# normalized transverse rms emittance of 1 um
energy_MeV = 2.0e3  # reference energy
charge_C = 1.0e-9  # used with space charge
mass_MeV = 938.27208816
qm_qeeV = -1.0e-6 / mass_MeV  # charge/mass
npart = 10000  # number of macro particles

distr = distribution.Waterbag(
    sigmaX=1.0e-3,
    sigmaY=1.0e-3,
    sigmaT=3.369701494258956e-4,
    sigmaPx=1.0e-3,
    sigmaPy=1.0e-3,
    sigmaPt=2.9676219145931020e-3,
    muxpx=0.0,
    muypy=0.0,
    mutpt=0.0,
)
sim.add_particles(qm_qeeV, charge_C, distr, npart)

# set the energy in the reference particle
sim.particle_container().ref_particle().set_energy_MeV(energy_MeV, mass_MeV)

# design the accelerator lattice
sim.lattice.append(elements.ConstF(
    ds=2.0, kx=1.0, ky=1.0, kt=1.0
))

# run simulation
sim.evolve()

# clean shutdown
del sim
amrex.finalize()
