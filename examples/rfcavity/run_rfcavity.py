#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Marco Garten, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from impactx import ImpactX, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = False
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
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

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

#   Drift elements
dr1 = elements.Drift(name="dr1", ds=0.4, nslice=1)
dr2 = elements.Drift(name="dr2", ds=0.032997, nslice=1)
#   RF cavity element
rf = elements.RFCavity(
    name="rf",
    ds=1.31879807,
    escale=62.0,
    freq=1.3e9,
    phase=85.5,
    cos_coefficients=[
        0.1644024074311037,
        -0.1324009958969339,
        4.3443060026047219e-002,
        8.5602654094946495e-002,
        -0.2433578169042885,
        0.5297150596779437,
        0.7164884680963959,
        -5.2579522442877296e-003,
        -5.5025369142193678e-002,
        4.6845673335028933e-002,
        -2.3279346335638568e-002,
        4.0800777539657775e-003,
        4.1378326533752169e-003,
        -2.5040533340490805e-003,
        -4.0654981400000964e-003,
        9.6630592067498289e-003,
        -8.5275895985990214e-003,
        -5.8078747006425020e-002,
        -2.4044337836660403e-002,
        1.0968240064697212e-002,
        -3.4461179858301418e-003,
        -8.1201564869443749e-004,
        2.1438992904959380e-003,
        -1.4997753525697276e-003,
        1.8685171825676386e-004,
    ],
    sin_coefficients=[
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    ],
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
