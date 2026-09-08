#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Marco Garten, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np

from impactx import ImpactX, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
sim.short_range_wakefields = True
sim.wake_bins = 4
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = False

# domain decomposition & space charge mesh
sim.init_grids()

# load a 230 MeV electron beam with an initial
# unnormalized rms emittance of 1 mm-mrad in all
# three phase planes
kin_energy_MeV = 230.0  # reference energy
bunch_charge_C = 1.0e-10  # used with space charge
npart = 10000  # number of macro particles (outside tests, use 1e5 or more)

#   reference particle
ref = sim.beam.ref
ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    lambdaX=0.352498964601e-3,
    lambdaY=0.207443478142e-3,
    lambdaT=0.70399950746e-4,
    lambdaPx=5.161852770e-6,
    lambdaPy=9.163582894e-6,
    lambdaPt=0.260528852031e-3,
    muxpx=0.5712386101751441,
    muypy=-0.514495755427526,
    mutpt=-5.05773e-10,
)
sim.add_particles(bunch_charge_C, distr, npart)

# design the accelerator lattice

# access RF cavity on-axis field data
data_in = np.loadtxt("onaxis_data.in")
z = data_in[:, 0]
ez_onaxis = data_in[:, 1]
ncoef = 25

# access RF cavity short-range structure wake data
wake_data_in = np.loadtxt("wake_data.in")
z_wake = wake_data_in[:, 0]
wake_x_data = wake_data_in[:, 1]
wake_y_data = wake_data_in[:, 2]
wake_z_data = wake_data_in[:, 3]

print("z_wake = ")
print(z_wake)

#   Drift elements
dr1 = elements.Drift(name="dr1", ds=0.4, nslice=1)
dr2 = elements.Drift(name="dr2", ds=0.032997, nslice=1)
#   RF cavity element
rf = elements.RFCavity(
    name="rf",
    ds=1.31879807,
    escale=62.0,
    z=z,
    field_on_axis=ez_onaxis,
    z_data=z_wake,
    wake_x_data=wake_x_data,
    wake_y_data=wake_y_data,
    wake_z_data=wake_z_data,
    ncoef=ncoef,
    freq=1.3e9,
    phase=85.5,
    mapsteps=100,
    nslice=4,
)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

sim.lattice.extend(
    [
        monitor,
        dr1,
        dr2,
        rf,
        dr2,
        dr2,
        rf,
        dr2,
        dr2,
        rf,
        dr2,
        dr2,
        rf,
        dr2,
        monitor,
    ]
)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
