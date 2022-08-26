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
                     madx2impactx_lattice, madx2impactx_beam)

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


# print summary
print(madx)

beamline = madx.getBeamline()

ref_particle_dict = madx2impactx_beam(
    madx.getParticle(),  # if particle species is known, mass, charge, and potentially energy are set to default
    # TODO MADX parser needs to extract charge if it's given,
    # TODO MADX parser needs to extract mass if it's given,
    energy=float(madx.getEtot())  # MADX default energy is in GeV
)


# load a 40 MeV electron beam with an initial
# normalized transverse rms emittance of 1 um

# @TODO read CHARGE (single particle charge) from MADX as well
charge_C = 1.0e-9  # used with space charge

qm_qeeV = ref_particle_dict['mass']/ref_particle_dict['charge']  # charge/mass
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

# set the energy in the reference particle
sim.particle_container().ref_particle() \
    .set_energy_MeV(ref_particle_dict['energy'], ref_particle_dict['mass'])

chicane = madx2impactx_lattice(beamline)
# assign a fodo segment
sim.lattice.extend(chicane)

# run simulation
sim.evolve()

# clean shutdown
del sim
amrex.finalize()
