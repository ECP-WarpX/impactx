#!/usr/bin/env python3
#
# Copyright 2022 ImpactX contributors
# Authors: Marco Garten, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import os

import amrex
from impactx import (
    ImpactX,
    MADXParser,
    RefPart,
    distribution,
    elements,
    madx2impactx_beam,
    madx2impactx_lattice,
)

sim = ImpactX()

# set numerical parameters and IO control
sim.set_particle_shape(2)  # B-spline order
sim.set_space_charge(False)
# sim.set_diagnostics(False)  # benchmarking
sim.set_slice_step_diagnostics(True)

# domain decomposition & space charge mesh
sim.init_grids()

# @TODO make this better ... shouldn't need that
dir_path = os.path.dirname(os.path.realpath(__file__))

try:
    madx = MADXParser()

    madx.parse(f"{dir_path}/fodo.madx")

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
    energy=float(madx.getEtot()),  # MADX default energy is in GeV
)

# load a 2 GeV electron beam with an initial
# unnormalized rms emittance of 2 nm

# @TODO read CHARGE (single particle charge) from MADX as well
charge_C = 1.0e-9  # used with space charge

qm_qeeV = ref_particle_dict["mass"] / ref_particle_dict["charge"]  # charge/mass
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

# design the accelerator lattice
ns = 25  # number of slices per ds in the element

# set the energy in the reference particle
sim.particle_container().ref_particle().set_energy_MeV(
    ref_particle_dict["energy"], ref_particle_dict["mass"]
)

fodo = madx2impactx_lattice(beamline)
# assign a fodo segment
sim.lattice.extend(fodo)

# run simulation
sim.evolve()

# clean shutdown
del sim
amrex.finalize()
