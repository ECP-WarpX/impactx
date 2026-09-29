#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Ryan Sandberg, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
import transformation_utilities as pycoord

import amrex.space3d as amr
from impactx import ImpactX, elements

################

N_part = int(1e5)
beam_radius = 2e-3
sigr = 500e-6
sigpx = 10
sigpy = 10
# fixed seed for deterministic test results
rng = np.random.default_rng(seed=42)
px = rng.normal(0, sigpx, N_part)
py = rng.normal(0, sigpy, N_part)
theta = 2 * np.pi * rng.random(N_part)
r = np.abs(rng.normal(beam_radius, sigr, N_part))
x = r * np.cos(theta)
y = r * np.sin(theta)
z_mean = 0
pz_mean = 2e4
z_std = 1e-3
pz_std = 2e2
zpz_std = -0.18
zpz_cov_list = [[z_std**2, zpz_std], [zpz_std, pz_std**2]]
z, pz = rng.multivariate_normal([0, 0], zpz_cov_list, N_part).T
pz += pz_mean

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

energy_gamma = np.sqrt(1 + pz_mean**2)
energy_MeV = 0.510998950 * energy_gamma  # reference energy
bunch_charge_C = 10.0e-12  # used with space charge
q_e_C = 1.60217663e-19

#   reference particle
ref = sim.beam.ref
ref.set_species("electron").set_kin_energy_MeV(energy_MeV)
qm_eev = -1.0 / 0.510998950 / 1e6  # electron charge/mass in e / eV
ref.z = 0

beam = sim.beam

# Note for MPI-parallel simulations:
#   `beam.add_n_particles(...)` is local to the MPI rank, spatial
#   locality does not matter. Thus, you can add particles at any
#   MPI rank, e.g., equally chuncked up for perfect load balancing.
#
#   You do NOT want to add the same unique particle at multiple
#   MPI ranks.
#
#   When ImpactX needs to sort particles spatially, it will
#   redistribute them over MPI ranks automatically during tracking.
#
#   In the example here, we add all particles from one MPI rank.
#   This is simple but not scalable -- for many particles just
#   add 1/N unique particles per MPI rank.

# The particles are given at fixed t, so we add them inside a
# `with beam.at_fixed_t():` block: ImpactX transforms them to
# fixed s when the block ends, using the design energy of `ref`.
with beam.at_fixed_t():
    if amr.ParallelDescriptor.IOProcessor():
        # particle coordinates relative to the reference particle, at fixed t
        dx, dy, dz, dpx, dpy, dpz = pycoord.to_ref_part_t_from_global_t(
            ref, x, y, z, px, py, pz
        )

        # here we use equal particle weighting, but you can assign any weight to each particle
        w = np.ones_like(dx) * (bunch_charge_C / q_e_C / N_part)

        # This call has two options:
        # A) reassign equal weighting according to bunch_charge_C
        # B) use the particle weighting from the input array w
        beam.add_n_particles(
            dx, dy, dz, dpx, dpy, dpz, qm_eev, bunch_charge=bunch_charge_C
        )
        # ok, let's clear all particles and do option B
        beam.clear_particles()

        beam.add_n_particles(dx, dy, dz, dpx, dpy, dpz, qm_eev, w=w)

# build the accelerator lattice
monitor = elements.BeamMonitor("monitor", backend="h5")
sim.lattice.extend(
    [
        monitor,
        elements.Drift(name="drift", ds=0.01),
        monitor,
    ]
)

sim.track_particles()

# clean shutdown
sim.finalize()
