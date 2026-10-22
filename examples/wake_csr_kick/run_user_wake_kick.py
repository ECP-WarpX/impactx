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
sim.wake_bins = 150
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = False

# domain decomposition & space charge mesh
sim.init_grids()

# load a 230 MeV electron beam with an initial
# unnormalized rms emittance of 1 mm-mrad in all
# three phase planes
kin_energy_MeV = 230.0  # reference energy
bunch_charge_C = 1.0e-10  # used with space charge
npart = 1000000  # number of macro particles (outside tests, use 1e5 or more)

#   reference particle
ref = sim.beam.ref
ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.KVdist(
# distr = distribution.Waterbag(
# distr = distribution.Triangle(
# distr = distribution.Gaussian(
    lambdaX=0.0,
    lambdaY=0.0,
    lambdaT=1.0e-3,
    lambdaPx=0.0,
    lambdaPy=0.0,
    lambdaPt=0.0,
    muxpx=0.0,
    muypy=0.0,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# design the accelerator lattice

# access RF cavity on-axis field data
data_in = np.loadtxt("onaxis_data.in")
z = data_in[:, 0]
ez_onaxis = data_in[:, 1]
ncoef = 25

# specify the RF cavity short-range structure wake data using an expression
zmin = (
    0.0  # lower value of longitudinal separation used in the wake function (in meters)
)
zmax = (
    0.01  # upper value of longitudinal separation used in the wake function (in meters)
)
nz = 5000  # number of longitudinal separation sampling points to be used
K0 = 4.30355e13  # first fitted parameter
s0 = 0.00175  # second fitted parameter (in meters)
z_wake = np.linspace(zmin, zmax, nz)
wake_x_data = np.zeros(nz)
wake_y_data = np.zeros(nz)
wake_z_data = K0 * np.exp(
    -np.sqrt(z_wake / s0)
)

#   RF cavity element
rf = elements.RFCavity(
    name="rf",
    ds=1.31879807,
    escale=0.0,
    z=z,
    field_on_axis=ez_onaxis,
    z_data=z_wake,
    wake_x_data=wake_x_data,
    wake_y_data=wake_y_data,
    wake_z_data=wake_z_data,
    ncoef=ncoef,
    freq=1.3e9,
    phase=85.5,
    mapsteps=400,
    nslice=1,
)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

sim.lattice.extend(
    [
        monitor,
        rf,
        monitor,
    ]
)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
