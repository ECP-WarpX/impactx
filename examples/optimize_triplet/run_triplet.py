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
# unnormalized rms emittance of 5 nm
kin_energy_MeV = 2.0e3  # reference energy
bunch_charge_C = 100.0e-12  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    sigmaX=2.0e-4,
    sigmaY=2.0e-4,
    sigmaT=3.1622776602e-5,
    sigmaPx=1.1180339887e-5,
    sigmaPy=1.1180339887e-5,
    sigmaPt=3.1622776602e-5,
    muxpx=0.894427190999916,
    muypy=-0.894427190999916,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice)
ns = 10  # number of slices per ds in the element
period = [
    monitor,
    elements.Drift(ds=2.7, nslice=ns),
    elements.Quad(ds=0.1, k=-3.5, nslice=ns),
    elements.Drift(ds=1.4, nslice=ns),
    elements.Quad(ds=0.2, k=2.75, nslice=ns),
    elements.Drift(ds=1.4, nslice=ns),
    elements.Quad(ds=0.1, k=-3.5, nslice=ns),
    elements.Drift(ds=2.7, nslice=ns),
    monitor,
]

sim.lattice.extend(period)

# number of periods through the lattice
sim.periods = 1

# run simulation
sim.evolve()

# clean shutdown
sim.finalize()
