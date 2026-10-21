#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

# The kicker example, with more grids than MPI ranks: on two MPI ranks, each
# rank owns two of the four grids. Without space charge, the particles are not
# redistributed, so they stay on the grid that they were added to.

from impactx import ImpactX, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition: four grids
sim.n_cell = [16, 16, 64]
sim.max_grid_size = [16]
sim.init_grids()

# load a 2 GeV electron beam
kin_energy_MeV = 2.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used without space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.beam.ref
ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    lambdaX=4.0e-3,
    lambdaY=4.0e-3,
    lambdaT=1.0e-3,
    lambdaPx=3.0e-4,
    lambdaPy=3.0e-4,
    lambdaPt=2.0e-3,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("m1", backend="h5")

# design the accelerator lattice
kicklattice = [
    monitor,
    elements.Kicker(name="kick1", xkick=2.0e-3, ykick=0.0, unit="dimensionless"),
    elements.Kicker(name="kick2", xkick=0.0, ykick=3.0e-3, unit="dimensionless"),
    monitor,
]
# assign a lattice
sim.lattice.extend(kicklattice)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
