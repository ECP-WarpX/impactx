#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from impactx import ImpactX, distribution, elements

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
kin_energy_MeV = 100.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 100000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Gaussian(
    lambdaX=1.288697604e-6,
    lambdaY=1.288697604e-6,
    lambdaT=1.0e-6,
    lambdaPx=3.965223396e-6,
    lambdaPy=3.965223396e-6,
    lambdaPt=0.01,  # 1% energy spread
    muxpx=0.0,
    muypy=0.0,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice)
ns = 25  # number of slices per ds in the element

# Drift elements
dr1 = elements.ChrDrift(ds=1.0, nslice=ns)
dr2 = elements.ChrDrift(ds=2.0, nslice=ns)

# Plasma lens elements
q1 = elements.ChrPlasmaLens(
    ds=0.331817852986604588, k=996.147787384348956, unit=1, nslice=ns
)
q2 = elements.ChrPlasmaLens(
    ds=0.176038957633108457, k=528.485181135649785, unit=1, nslice=ns
)
q3 = elements.ChrPlasmaLens(
    ds=1.041842576046930486, k=3127.707468391874166, unit=1, nslice=ns
)
q4 = elements.ChrPlasmaLens(
    ds=0.334367090894399520, k=501.900417308233112, unit=1, nslice=ns
)
q5 = elements.ChrPlasmaLens(
    ds=1.041842576046930486, k=3127.707468391874166, unit=1, nslice=ns
)
q6 = elements.ChrPlasmaLens(
    ds=0.176038957633108457, k=528.485181135649785, unit=1, nslice=ns
)
q7 = elements.ChrPlasmaLens(
    ds=0.331817852986604588, k=996.147787384348956, unit=1, nslice=ns
)

# q1 = elements.ChrPlasmaLens(ds=0.331817852986604588, k=2.98636067687944129, unit=0, nslice=ns)
# q2 = elements.ChrPlasmaLens(ds=0.176038957633108457, k=1.584350618697976110, unit=0, nslice=ns)
# q3 = elements.ChrPlasmaLens(ds=1.041842576046930486, k=9.37658318442237437, unit=0, nslice=ns)
# q4 = elements.ChrPlasmaLens(ds=0.334367090894399520, k=1.50465190902479784, unit=0, nslice=ns)
# q5 = elements.ChrPlasmaLens(ds=1.041842576046930486, k=9.37658318442237437, unit=0, nslice=ns)
# q6 = elements.ChrPlasmaLens(ds=0.176038957633108457, k=1.584350618697976110, unit=0, nslice=ns)
# q7 = elements.ChrPlasmaLens(ds=0.331817852986604588, k=2.98636067687944129, unit=0, nslice=ns)

lattice_line = [
    monitor,
    dr1,
    q1,
    dr2,
    q2,
    dr2,
    q3,
    dr2,
    q4,
    dr2,
    q5,
    dr2,
    q6,
    dr2,
    q7,
    dr1,
    monitor,
]

# define the lattice
sim.lattice.extend(lattice_line)

# run simulation
sim.evolve()

# clean shutdown
sim.finalize()
