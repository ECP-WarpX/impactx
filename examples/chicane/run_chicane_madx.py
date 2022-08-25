#!/usr/bin/env python3
#
# Copyright 2022 ImpactX contributors
# Authors: Marco Garten, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import os

import amrex
from impactx import (ImpactX, MADXParser, RefPart, distribution, elements,
                     madx2impactx)

sim = ImpactX()

# set numerical parameters and IO control
sim.set_particle_shape(2)  # B-spline order
sim.set_space_charge(False)
#sim.set_diagnostics(False)  # benchmarking
sim.set_slice_step_diagnostics(True)

# domain decomposition & space charge mesh
sim.init_grids()

# @TODO make this better ... shouldn't need that. But somehow the working directory is
dir_path = os.path.dirname(os.path.realpath(__file__))

try:
    madx = MADXParser()

    madx.parse('{}/input_chicane.madx'.format(dir_path))

except Exception as e:
    print(f"Unexpected {e = }, {type(e) = }")
    raise

# load a 40 MeV electron beam with an initial
# normalized transverse rms emittance of 1 um

# @TODO read CHARGE (single particle charge) from MADX as well
charge_C = 1.0e-9  # used with space charge
# @TODO get from dictionary b/c MADX knows P
mass_MeV = 0.510998950
qm_qeeV = -1.0e-6/mass_MeV  # charge/mass
npart = 10000  # number of macro particles

distr = distribution.Waterbag(
    sigmaX = 2.2951017632e-5,
    sigmaY = 1.3084093142e-5,
    sigmaT = 5.5555553e-8,
    sigmaPx = 1.598353425e-6,
    sigmaPy = 2.803697378e-6,
    sigmaPt = 2.000000000e-6,
    muxpx = 0.933345606203060,
    muypy = 0.933345606203060,
    mutpt = 0.999999961419755)
sim.add_particles(qm_qeeV, charge_C, distr, npart)

# design the accelerator lattice
ns = 25  # number of slices per ds in the element


# print summary
print(madx)

# MADX default energy is in GeV
energy_MeV = float(madx.getEtot()) * 1e3
beamline = madx.getBeamline()

# set the energy in the reference particle
sim.particle_container().ref_particle() \
    .set_energy_MeV(energy_MeV, mass_MeV)

chicane = madx2impactx(beamline)
# assign a fodo segment
sim.lattice.extend(chicane)

# run simulation
sim.evolve()

# clean shutdown
del sim
amrex.finalize()
