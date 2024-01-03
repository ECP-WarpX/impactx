#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import amrex.space3d as amr
from impactx import ImpactX, RefPart, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 2 GeV electron beam with an initial
# unnormalized rms emittance of 2 nm
kin_energy_MeV = 2.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 100000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    sigmaX=1.16098260008648811e-3,
    sigmaY=1.16098260008648811e-3,
    sigmaT=1.0e-3,
    sigmaPx=0.580491300043e-3,
    sigmaPy=0.580491300043e-3,
    sigmaPt=2.0e-3,
    muxpx=0.0,
    muypy=0.0,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice)
ns = 1  # number of slices per ds in the element
lattice = [
    monitor,
    elements.Quad(ds=1.0, k=0.25, dx=0.003, dy=0.0, rotation=30.0, nslice=ns),
    monitor,
]
sim.lattice.extend(lattice)

# run simulation
sim.evolve()

# clean shutdown
del sim
amr.finalize()
