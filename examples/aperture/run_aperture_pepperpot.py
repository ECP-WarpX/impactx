#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import amrex.space3d as amr
from impactx import ImpactX, distribution, elements

# work-around for https://github.com/ECP-WarpX/impactx/issues/499
pp_amrex = amr.ParmParse("amrex")
pp_amrex.add("the_arena_is_managed", 1)

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True
sim.particle_lost_diagnostics_backend = "h5"

# domain decomposition & space charge mesh
sim.init_grids()

# load a 250 MeV proton beam with an initial
# horizontal rms emittance of 1 um and an
# initial vertical rms emittance of 2 um
kin_energy_MeV = 250.0  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 1000000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    lambdaX=1.559531175539e-3,
    lambdaY=2.205510139392e-3,
    lambdaT=1.0e-3,
    lambdaPx=6.41218345413e-4,
    lambdaPy=9.06819680526e-4,
    lambdaPt=1.0e-3,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
sim.lattice.extend(
    [
        monitor,
        elements.Drift(name="drift", ds=0.123),
        elements.Aperture(
            name="pepperpot",
            xmax=1.5e-4,
            ymax=1.0e-4,
            repeat_x=1.0e-3,
            repeat_y=1.0e-3,
            shape="rectangular",
        ),
        monitor,
    ]
)

# run simulation
sim.evolve()

# clean shutdown
sim.finalize()
