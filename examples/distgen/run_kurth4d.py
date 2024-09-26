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

# load a 2 GeV proton beam with an initial
# normalized transverse rms emittance of 1 um
kin_energy_MeV = 2.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Kurth4D(
    lambdaX=1.0e-3,
    lambdaY=1.0e-3,
    lambdaT=1.0e-3,
    lambdaPx=1.0e-3,
    lambdaPy=1.0e-3,
    lambdaPt=2.0e-3,
    muxpx=-0.0,
    muypy=0.0,
    mutpt=0.0,
)

sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice)
constf = [
    monitor,
    elements.ConstF(name="constf1", ds=2.0, kx=1.0, ky=1.0, kt=1.0e-4),
    monitor,
]

# assign a constf segment
sim.lattice.extend(constf)

# run simulation
sim.evolve()

# clean shutdown
sim.finalize()
