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
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Gaussian(
    lambdaX=2.2951017632e-5,
    lambdaY=1.3084093142e-5,
    lambdaT=5.5555553e-8,
    lambdaPx=1.598353425e-6,
    lambdaPy=2.803697378e-6,
    lambdaPt=2.000000000e-6,
    muxpx=0.933345606203060,
    muypy=0.933345606203060,
    mutpt=0.999999961419755,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
ns = 25  # number of slices per ds in the element
rc = 10.3462283686195526  # bend radius (meters)
psi = 0.048345620280243  # pole face rotation angle (radians)
lb = 0.500194828041958  # bend arc length (meters)

# Drift elements
dr1 = elements.Drift(name="dr1", ds=5.0058489435, nslice=ns)
dr2 = elements.Drift(name="dr2", ds=1.0, nslice=ns)
dr3 = elements.Drift(name="dr3", ds=2.0, nslice=ns)

# Bend elements
sbend1 = elements.Sbend(name="sbend1", ds=lb, rc=-rc, nslice=ns)
sbend2 = elements.Sbend(name="sbend2", ds=lb, rc=rc, nslice=ns)

# Dipole Edge Focusing elements
dipedge1 = elements.DipEdge(name="dipedge1", psi=-psi, rc=-rc, g=0.0, K2=0.0)
dipedge2 = elements.DipEdge(name="dipedge2", psi=psi, rc=rc, g=0.0, K2=0.0)

lattice_half = [sbend1, dipedge1, dr1, dipedge2, sbend2]
# assign a segment with the first half of the lattice
sim.lattice.append(monitor)
sim.lattice.extend(lattice_half)
sim.lattice.append(dr2)
lattice_half.reverse()
# extend the lattice by a reversed half
sim.lattice.extend(lattice_half)
sim.lattice.append(dr3)
sim.lattice.append(monitor)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
