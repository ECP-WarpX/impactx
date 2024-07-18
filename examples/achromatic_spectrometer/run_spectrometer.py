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

# load a 1 GeV electron beam
kin_energy_MeV = 1.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    lambdaX=3.162277660e-6,
    lambdaY=3.162277660e-6,
    lambdaT=1.0e-3,
    lambdaPx=3.16227766017e-4,
    lambdaPy=3.16227766017e-4,
    lambdaPt=2.0e-2,
    muxpx=0.0,
    muypy=0.0,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
ns = 25  # number of slices per ds in the element

# specify thick tapered plasma lens element
lens_length = 0.02  # length in m
num_lenses = 10
focal_length = 0.5  # focal length in m
dtaper = 11.488289081903567  # 1/(horizontal dispersion in m)
ds = lens_length / num_lenses
dk = 1.0 / (focal_length * num_lenses)

# drifts appearing the drift-kick sequence
ds_half = ds / 2.0
dr = elements.Drift(ds=ds_half, nslice=ns)

# define the lens segments
thick_lens = []
for _ in range(0, num_lenses):
    pl = elements.TaperedPL(k=dk, taper=dtaper, unit=0)
    segment = [dr, pl, dr]
    thick_lens.extend(segment)

bend = elements.ExactSbend(ds=1.0, phi=10.0, B=0.0, nslice=ns)
drift = elements.Drift(ds=1.0, nslice=ns)

# specify the lattice sequence
sim.lattice.append(monitor)
sim.lattice.append(bend)
sim.lattice.extend(thick_lens)
sim.lattice.append(drift)
sim.lattice.append(monitor)

# run simulation
sim.evolve()

# clean shutdown
sim.finalize()
