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
kin_energy_MeV = 20.0  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    sigmaX=0.5e-3,
    sigmaY=0.5e-3,
    sigmaT=5.0e-3,
    sigmaPx=1.0e-5,
    sigmaPy=1.0e-5,
    sigmaPt=4.0e-6,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
sim.lattice.append(monitor)
#   Short RF cavity element
shortrf1 = elements.ShortRF(V=1000.0, freq=1.3e9, phase=-89.5)
#   Drift element
drift1 = elements.Drift(ds=1.7)

sim.lattice.extend([shortrf1, drift1])

sim.lattice.append(monitor)

# run simulation
sim.evolve()

# clean shutdown
sim.finalize()
