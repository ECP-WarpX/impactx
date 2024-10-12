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
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 250 MeV proton beam with an initial
# unnormalized rms emittance of 1 mm-mrad in all
# three phase planes
kin_energy_MeV = 250.0  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    lambdaX=3.131948925200e-3,
    lambdaY=1.148450209423e-3,
    lambdaT=2.159922887089e-3,
    lambdaPx=3.192900088357e-4,
    lambdaPy=8.707386631090e-4,
    lambdaPt=4.62979491526e-4,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
sim.lattice.append(monitor)
#   Quad elements
quad1 = elements.Quad(name="quad1", ds=0.15, k=2.5)
quad2 = elements.Quad(name="quad2", ds=0.3, k=-2.5)
#   Drift element
drift1 = elements.Drift(name="drift1", ds=1.0)
#   Short RF cavity element
shortrf1 = elements.Buncher(name="shortrf1", V=0.01, k=15.0)

lattice_no_drifts = [quad1, shortrf1, quad2, shortrf1, quad1]
#   set first lattice element
sim.lattice.append(lattice_no_drifts[0])
#   intersperse all remaining elements of the lattice with a drift element
for element in lattice_no_drifts[1:]:
    sim.lattice.extend([drift1, element])

sim.lattice.append(monitor)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
