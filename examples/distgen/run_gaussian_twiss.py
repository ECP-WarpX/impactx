#!/usr/bin/env python3
#
# Copyright 2022-2024 ImpactX contributors
# Authors: Marco Garten, Axel Huebl, CHad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from impactx import ImpactX, distribution, elements, twiss

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
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Gaussian(
    **twiss(
        beta_x=2.8216194100262637,
        beta_y=2.8216194100262637,
        beta_t=0.5,
        emitt_x=2e-09,
        emitt_y=2e-09,
        emitt_t=2e-06,
        alpha_x=-1.5905003499999992,
        alpha_y=1.5905003499999992,
        alpha_t=0.0
    )
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice)
ns = 25  # number of slices per ds in the element
fodo = [
    monitor,
    elements.Drift(ds=0.25, nslice=ns),
    monitor,
    elements.Quad(ds=1.0, k=1.0, nslice=ns),
    monitor,
    elements.Drift(ds=0.5, nslice=ns),
    monitor,
    elements.Quad(ds=1.0, k=-1.0, nslice=ns),
    monitor,
    elements.Drift(ds=0.25, nslice=ns),
    monitor,
]
# assign a fodo segment
sim.lattice.extend(fodo)

# run simulation
sim.evolve()

# clean shutdown
sim.finalize()
