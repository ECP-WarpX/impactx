#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Ryan Sandberg, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np

import amrex
from impactx import (
    Config,
    ImpactX,
    ImpactXParIter,
    RefPart,
    TransformationDirection,
    coordinate_transformation,
    distribution,
    elements,
)

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 1 GeV electron beam with an initial
# unnormalized rms emittance of 2 nm
energy_MeV = 1e3  # reference energy
energy_gamma = energy_MeV / 0.510998950
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
pc = sim.particle_container()
ref = pc.ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_energy_MeV(energy_MeV)

#   particle bunch
distr = distribution.Gaussian(
    sigmaX=3e-6,
    sigmaY=3e-6,
    sigmaT=1e-2,
    sigmaPx=1.33 / energy_gamma,
    sigmaPy=1.33 / energy_gamma,
    sigmaPt=100 / energy_gamma,
    muxpx=-0.5,
    muypy=0.4,
    mutpt=0.8,
)
sim.add_particles(bunch_charge_C, distr, npart)
# number of slices per ds in each lattice element
ns = 1

rbc_s0 = pc.reduced_beam_characteristics()
coordinate_transformation(pc, TransformationDirection.to_fixed_t)
rbc_t = pc.reduced_beam_characteristics()
coordinate_transformation(pc, TransformationDirection.to_fixed_s)
rbc_s = pc.reduced_beam_characteristics()

# clean shutdown
del sim

# assert that forward-inverse transformation of the beam leaves beam unchanged
atol = 1e-14
rtol = 1e-10
for key, val in rbc_s0.items():
    if not np.isclose(val, rbc_s[key], rtol=rtol,atol=atol):
        print(f'initial[{key}]={val}, final[{key}]={rbc_s[key]} not equal')
    assert np.isclose(val, rbc_s[key], rtol=rtol,atol=atol)
# assert that the t-based beam is different, at least in the following keys:
large_st_diff_keys = ['beta_x', 'beta_y', 'emittance_y', 'emittance_x', 'sig_y', 'sig_x', 't_mean']
for key in large_st_diff_keys:
    rel_error = (rbc_s0[key] - rbc_t[key]) / rbc_s0[key]
    assert abs(rel_error) > 1
