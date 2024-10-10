#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Ryan Sandberg, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import math

from impactx import ImpactX, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 2.5 MeV proton beam
kin_energy_MeV = 2.5  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    lambdaX=1.397456296195e-003,
    lambdaY=1.397456296195e-003,
    lambdaT=1.0e-4,
    lambdaPx=1.256184325020e-003,
    lambdaPy=1.256184325020e-003,
    lambdaPt=0.0,
    muxpx=0.8090169943749474,
    muypy=0.8090169943749474,
    mutpt=0.0,
)

sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")
monitor.nonlinear_lens_invariants = True
monitor.alpha = 1.376381920471173
monitor.beta = 1.892632003628881
monitor.tn = 0.4
monitor.cn = 0.01

sim.lattice.append(monitor)

# defining parameters of the nonlinear lens
lens_length = 1.8
num_lenses = 18
tune_advance = 0.3
c_parameter = 0.01
t_strength = 0.4
ds = lens_length / num_lenses

# drift elements
ds_half = ds / 2.0
dr = elements.Drift(name="dr", ds=ds_half)

# define the nonlinear lens segments
for j in range(0, num_lenses):
    s = -lens_length / 2.0 + ds_half + j * ds
    beta_star = lens_length / 2.0 * 1.0 / math.tan(math.pi * tune_advance)
    beta = beta_star * (
        1.0 + (2.0 * s * math.tan(math.pi * tune_advance) / lens_length) ** 2
    )
    knll_s = t_strength * c_parameter**2 * ds / beta
    cnll_s = c_parameter * math.sqrt(beta)
    nllens = elements.NonlinearLens(name="nllens" + str(j), knll=knll_s, cnll=cnll_s)
    segments = [dr, nllens, dr]
    sim.lattice.extend(segments)

# focusing lens
const = elements.ConstF(
    name="const",
    ds=1.0e-8,
    kx=12060.113295833,
    ky=12060.113295833,
    kt=1.0e-12,
    nslice=1,
)
sim.lattice.append(const)
sim.lattice.append(monitor)

# number of periods
sim.periods = 1

# run simulation
sim.evolve()

# clean shutdown
sim.finalize()
