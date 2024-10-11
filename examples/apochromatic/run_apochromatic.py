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

# load a 2 GeV electron beam with an initial
# unnormalized rms emittance of 2 nm
kin_energy_MeV = 100.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 100000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Gaussian(
    lambdaX=1.288697604e-6,
    lambdaY=1.288697604e-6,
    lambdaT=1.0e-6,
    lambdaPx=3.965223396e-6,
    lambdaPy=3.965223396e-6,
    lambdaPt=0.01,  # 1% energy spread
    muxpx=0.0,
    muypy=0.0,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice)
ns = 25  # number of slices per ds in the element

# Drift elements
dr1 = elements.ChrDrift(name="dr1", ds=1.0, nslice=ns)
dr2 = elements.ChrDrift(name="dr2", ds=10.0, nslice=ns)

# Quad elements
q1 = elements.ChrQuad(name="q1", ds=1.2258333333, k=0.5884, nslice=ns)
q2 = elements.ChrQuad(name="q2", ds=1.5677083333, k=-0.7525, nslice=ns)
q3 = elements.ChrQuad(name="q3", ds=1.205625, k=0.5787, nslice=ns)
q4 = elements.ChrQuad(name="q4", ds=1.2502083333, k=-0.6001, nslice=ns)
q5 = elements.ChrQuad(name="q5", ds=1.2502083333, k=0.6001, nslice=ns)
q6 = elements.ChrQuad(name="q6", ds=1.205625, k=-0.5787, nslice=ns)
q7 = elements.ChrQuad(name="q7", ds=1.5677083333, k=0.7525, nslice=ns)
q8 = elements.ChrQuad(name="q8", ds=1.2258333333, k=-0.5884, nslice=ns)

lattice_line = [monitor, dr1, q1, q2, q3, dr2, q4, q5, dr2, q6, q7, q8, dr1, monitor]

# define the lattice
sim.lattice.extend(lattice_line)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
