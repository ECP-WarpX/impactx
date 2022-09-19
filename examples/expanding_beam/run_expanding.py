#!/usr/bin/env python3
#
# Copyright 2022 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import amrex
from impactx import ImpactX, RefPart, distribution, elements

pp_amr = amrex.ParmParse("amr")
pp_amr.addarr("n_cell", [40, 40, 32])

sim = ImpactX()

# set numerical parameters and IO control
sim.set_particle_shape(2)  # B-spline order
sim.set_space_charge(True)
sim.dynamic_size = True
sim.prob_relative = 1.0

# beam diagnostics
# sim.set_diagnostics(False)  # benchmarking
sim.set_slice_step_diagnostics(False)

# domain decomposition & space charge mesh
sim.init_grids()

# load a 2 GeV electron beam with an initial
# unnormalized rms emittance of 2 nm
energy_MeV = 250  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_energy_MeV(energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    sigmaX=4.472135955e-4,
    sigmaY=4.472135955e-4,
    sigmaT=9.12241869e-7,
    sigmaPx=0.0,
    sigmaPy=0.0,
    sigmaPt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# design the accelerator lattice
sim.lattice.append(elements.Drift(ds=6.0, nslice=30))

# run simulation
sim.evolve()

# clean shutdown
del sim
amrex.finalize()
