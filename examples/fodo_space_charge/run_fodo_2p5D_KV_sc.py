#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from impactx import ImpactX, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.max_level = 0
sim.n_cell = [32, 32, 1]
sim.blocking_factor_x = [16]
sim.blocking_factor_y = [16]
sim.blocking_factor_z = [1]

sim.particle_shape = 2  # B-spline order
sim.space_charge = "2p5D"
sim.poisson_solver = "fft"
sim.dynamic_size = True
sim.prob_relative = [1.1]

# special parameters for the 2.5D space charge solver
sim.space_charge_num_longitudinal_bins = 2
sim.space_charge_apply_longitudinal_kick = False

# beam diagnostics
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# reference energy
kin_energy_MeV = 6.7  # reference energy

#   reference particle
ref = sim.beam.ref
ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)

#  bunch charge
bunch_charge_C = 1.0e-12
npart = 10000  # number of macro particles (outside tests, use 1e5 or more)

#   particle bunch
distr = distribution.KVdist(
    lambdaX=3.22523868393e-4,
    lambdaY=3.22523868393e-4,
    lambdaT=1.73205080757e-4,
    lambdaPx=1.1641443538999e-3,
    lambdaPy=1.1641443538999e-3,
    lambdaPt=0.0,
    muxpx=0.926836840601929,
    muypy=-0.926836840601929,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice)
ns = 100  # number of slices per ds in the element
fodo = [
    monitor,
    elements.Drift(name="drift1", ds=7.44e-2, nslice=ns),
    elements.Quad(name="quad1", ds=6.10e-2, k=-103.12574100336, nslice=ns),
    elements.Drift(name="drift2", ds=14.88e-2, nslice=ns),
    elements.Quad(name="quad2", ds=6.10e-2, k=103.12574100336, nslice=ns),
    elements.Drift(name="drift3", ds=7.44e-2, nslice=ns),
    monitor,
]
# assign a fodo segment
sim.lattice.extend(fodo)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
