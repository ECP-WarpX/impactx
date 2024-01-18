#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Chad Mitchell, Axel Huebl
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

# load a 5 GeV electron beam with an initial
# normalized transverse rms emittance of 1 um
kin_energy_MeV = 2.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    sigmaX=5.0e-6,  # 5 um
    sigmaY=8.0e-6,  # 8 um
    sigmaT=0.0599584916,  # 200 ps
    sigmaPx=2.5543422003e-9,  # exn = 50 pm-rad
    sigmaPy=1.5964638752e-9,  # eyn = 50 pm-rad
    sigmaPt=9.0e-4,  # approximately dE/E
    muxpx=0.0,
    muypy=0.0,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
ns = 25  # number of slices per ds in the element

bend = [
    monitor,
    elements.CFbend(ds=0.5, rc=7.613657587094493, k=-7.057403, nslice=ns),
    monitor,
]

# assign a lattice segment
sim.lattice.extend(bend)

# run simulation
sim.evolve()

# clean shutdown
del sim
amr.finalize()
