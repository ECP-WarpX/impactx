#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Chad Mitchell, Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import math

import amrex.space3d as amr
from impactx import ImpactX, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
sim.slice_step_diagnostics = False

# domain decomposition & space charge mesh
sim.init_grids()

# load initial beam parameters
kin_energy_MeV = 0.8e3  # reference energy
bunch_charge_C = 0.0

#   reference particle
ref = sim.beam.ref
ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)
qm_eev = 1.0 / 938.27208816 / 1e6  # electron charge/mass in e / eV

#   particle bunch
beam = sim.beam

dx = [1, 0, 0, 0, 0, 0]
dpx = [0, 1, 0, 0, 0, 0]
dy = [0, 0, 1, 0, 0, 0]
dpy = [0, 0, 0, 1, 0, 0]
dt = [0, 0, 0, 0, 1, 0]
dpt = [0, 0, 0, 0, 0, 1]

if amr.ParallelDescriptor.IOProcessor():
    beam.add_n_particles(dx, dy, dt, dpx, dpy, dpt, qm_eev, bunch_charge_C)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice)
ns = 1  # number of slices per ds in the element

# parameters
edge_angle = math.pi / 8.0
dipedge0 = elements.DipEdge(
    name="dipedge0",
    psi=edge_angle,
    rc=10.0,
    g=1.0,
    K2=0.5,
    model="linear_slac75",
    location="entry",
)

sim.lattice.append(monitor)
sim.lattice.load_file("dipedge.madx", nslice=ns)  # to use the MAD-X parser
# sim.lattice.append(dipedge0)  # to use the equivalent parameters specified above
sim.lattice.append(monitor)

# return linear map
R_map = sim.lattice.transfer_map(ref).to_numpy()
print(R_map)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
