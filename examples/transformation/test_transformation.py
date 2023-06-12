#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Ryan Sandberg, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np

import amrex
from impactx import (  # NOQA
    Config,
    ImpactX,
    ImpactXParIter,
    RefPart,
    TransformationDirection,
    coordinate_transformation,
    distribution,
    elements,
)


def get_particle_data(pc):
    for lvl in range(pc.finest_level + 1):
        for pti in ImpactXParIter(pc, level=lvl):
            aos = pti.aos()
            aos_arr = np.array(aos, copy=False)

            soa = pti.soa()
            real_arrays = soa.GetRealData()
            data_arr = np.vstack(
                [aos_arr["x"], aos_arr["y"], aos_arr["z"], real_arrays[:3]]
            ).T
    return data_arr


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

data_arr_si = get_particle_data(pc)
coordinate_transformation(pc, TransformationDirection.to_fixed_t)
data_arr_t = get_particle_data(pc)
coordinate_transformation(pc, TransformationDirection.to_fixed_s)
data_arr_s = get_particle_data(pc)

# clean shutdown
del sim
amrex.finalize()

# assert that the t-based beam is different
assert np.max(abs(data_arr_t - data_arr_si)) > 0.1
# assert that forward-inverse transformation of the beam leaves beam unchanged
assert np.max(abs(data_arr_s - data_arr_si)) < 1e-15
