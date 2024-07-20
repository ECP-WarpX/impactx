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
sim.n_cell = [48, 48, 40]  # [72, 72, 64] for increased precision
sim.particle_shape = 2  # B-spline order
sim.space_charge = True
sim.prob_relative = [3.0]
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 2 GeV proton beam with an initial
# normalized transverse rms emittance of 1 um
kin_energy_MeV = 2.0e3  # reference energy
bunch_charge_C = 1.0e-8  # used with space charge
npart = 10000  # number of macro particles; use 1e5 for increased precision

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    lambdaX=1.2154443728379865788e-3,
    lambdaY=1.2154443728379865788e-3,
    lambdaT=4.0956844276541331005e-4,
    lambdaPx=8.2274435782286157175e-4,
    lambdaPy=8.2274435782286157175e-4,
    lambdaPt=2.4415943602685364584e-3,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
nslice = 50  # use 1e5 for increased precision

# design the accelerator lattice
sim.lattice.extend(
    [
        monitor,
        elements.ConstF(ds=2.0, kx=1.0, ky=1.0, kt=1.0, nslice=nslice),
        monitor,
    ]
)

# run simulation
sim.evolve()

# clean shutdown
sim.finalize()
