#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Marco Garten, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from impactx import ImpactX, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = False
sim.csr = True
sim.csr_bins = 150
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 5 GeV electron beam with an initial
# normalized transverse rms emittance of 1 um
kin_energy_MeV = 5.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 1000000  # number of macro particles

#   reference particle
ref = sim.beam.ref
ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Gaussian(
    lambdaX=0.0,
    lambdaY=0.0,
    lambdaT=5.5555553e-8,
    lambdaPx=0.0,
    lambdaPy=0.0,
    lambdaPt=0.0,
    muxpx=0.0,
    muypy=0.0,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
ns = 1  # number of slices per ds in the element
rc = 10.3462283686195526  # bend radius (meters)
psi = 0.048345620280243  # pole face rotation angle (radians)
lb = 0.500194828041958  # bend arc length (meters)

# Bend elements
sbend1 = elements.Sbend(name="sbend1", ds=lb, rc=-rc, nslice=ns)
sbend2 = elements.Sbend(name="sbend2", ds=lb, rc=rc, nslice=ns)

# assign a segment with the first half of the lattice
sim.lattice.append(monitor)
sim.lattice.append(sbend1)
sim.lattice.append(monitor)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
