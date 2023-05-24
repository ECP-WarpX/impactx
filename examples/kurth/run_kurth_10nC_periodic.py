#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Ryan Sandberg, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import amrex.space3d as amr
from impactx import ImpactX, distribution, elements

pp_amr = amr.ParmParse("amr")
pp_amr.addarr("n_cell", [48, 48, 40])  # use [72, 72, 72] for increased precision

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = True
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 2 GeV proton beam with an initial
# unnormalized rms emittance of 1 um in each
# coordinate plane
energy_MeV = 2.0e3  # reference energy
bunch_charge_C = 1.0e-8  # used with space charge
npart = 10000  # number of macro particles; use 1e5 for increased precision

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_energy_MeV(energy_MeV)

#   particle bunch
distr = distribution.Kurth6D(
    sigmaX=1.46e-3,
    sigmaY=1.46e-3,
    sigmaT=4.9197638312420749e-4,
    sigmaPx=6.84931506849e-4,
    sigmaPy=6.84931506849e-4,
    sigmaPt=2.0326178944803812e-3,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
nslice = 20  # use 30 for increased precision
constf1 = elements.ConstF(ds=2.0, kx=0.7, ky=0.7, kt=0.7, nslice=nslice)
drift1 = elements.Drift(ds=1.0, nslice=nslice)
sim.lattice.extend([monitor, drift1, constf1, drift1, monitor])

# run simulation
sim.evolve()

# clean shutdown
del sim
amr.finalize()
