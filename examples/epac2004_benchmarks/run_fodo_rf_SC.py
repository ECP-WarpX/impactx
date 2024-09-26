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
sim.n_cell = [56, 56, 64]
sim.particle_shape = 2  # B-spline order
sim.space_charge = True
sim.dynamic_size = True
sim.prob_relative = [4.0]

# beam diagnostics
sim.slice_step_diagnostics = False

# domain decomposition & space charge mesh
sim.init_grids()

# beam parameters
kin_energy_MeV = 250.0  # reference energy
bunch_charge_C = 1.42857142857142865e-10  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Kurth6D(
    lambdaX=9.84722273e-4,
    lambdaY=6.96967278e-4,
    lambdaT=4.486799242214e-03,
    lambdaPx=0.0,
    lambdaPy=0.0,
    lambdaPt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
sim.lattice.append(monitor)

#   Quad elements
fquad = elements.Quad(name="fquad", ds=0.15, k=2.4669749766168163, nslice=6)
dquad = elements.Quad(name="dquad", ds=0.3, k=-2.4669749766168163, nslice=12)

#   Drift element
dr = elements.Drift(name="dr", ds=0.1, nslice=4)

#   RF cavity elements
gapa1 = elements.RFCavity(
    name="gapa1",
    ds=1.0,
    escale=0.042631556991578,
    freq=7.0e8,
    phase=45.0,
    cos_coefficients=[
        0.120864178375839,
        -0.044057987631337,
        -0.209107290958498,
        -0.019831226655815,
        0.290428111491964,
        0.381974267375227,
        0.276801212694382,
        0.148265085353012,
        0.068569351192205,
        0.0290155855315696,
        0.011281649986680,
        0.004108501632832,
        0.0014277644197320,
        0.000474212125404,
        0.000151675768439,
        0.000047031436898,
        0.000014154595193,
        4.154741658e-6,
        1.191423909e-6,
        3.348293360e-7,
        9.203061700e-8,
        2.515007200e-8,
        6.478108000e-9,
        1.912531000e-9,
        2.925600000e-10,
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

gapb1 = elements.RFCavity(
    name="gapb1",
    ds=1.0,
    escale=0.042631556991578,
    freq=7.0e8,
    phase=-1.0,
    cos_coefficients=[
        0.120864178375839,
        -0.044057987631337,
        -0.209107290958498,
        -0.019831226655815,
        0.290428111491964,
        0.381974267375227,
        0.276801212694382,
        0.148265085353012,
        0.068569351192205,
        0.0290155855315696,
        0.011281649986680,
        0.004108501632832,
        0.0014277644197320,
        0.000474212125404,
        0.000151675768439,
        0.000047031436898,
        0.000014154595193,
        4.154741658e-6,
        1.191423909e-6,
        3.348293360e-7,
        9.203061700e-8,
        2.515007200e-8,
        6.478108000e-9,
        1.912531000e-9,
        2.925600000e-10,
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

lattice_no_drifts = [fquad, gapa1, dquad, gapb1, fquad]

#   set first lattice element
sim.lattice.append(lattice_no_drifts[0])

#   intersperse all remaining elements of the lattice with a drift element
for element in lattice_no_drifts[1:]:
    sim.lattice.extend([dr, element])

sim.lattice.append(monitor)

# run simulation
sim.evolve()

# clean shutdown
sim.finalize()
