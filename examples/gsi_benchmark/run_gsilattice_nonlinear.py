#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Chad Mitchell, Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from distribution_gsi import get_reference_params

import amrex.space3d as amr
from impactx import ImpactX, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

#  set reference particle
ref = sim.beam.ref
kin_energy_MeV, bunch_charge_C, charge_qe = get_reference_params()
ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)
qm_eev = ref.charge_qe / (ref.mass_MeV * 1.0e6)  # electron charge/mass in e / eV


# set test particles
pc = sim.particle_container()

#  add test particles
if amr.ParallelDescriptor.IOProcessor():
    dx = np.linspace(0, 4.31e-2, 20)
    zero_arr = np.linspace(0, 0.0, 20)
    pc.add_n_particles(
        dx, zero_arr, zero_arr, zero_arr, zero_arr, zero_arr, qm_eev, bunch_charge=0.0
    )

# init accelerator lattice
ns = 1  # number of slices per ds in the element

# Drift elements
dr1 = elements.Drift(name="dr1", ds=0.6450000, nslice=ns)
dr2 = elements.Drift(name="dr2", ds=0.9700000, nslice=ns)
dr3 = elements.Drift(name="dr3", ds=6.8390117, nslice=ns)
dr4 = elements.Drift(name="dr4", ds=0.6000000, nslice=ns)
dr5 = elements.Drift(name="dr5", ds=0.7098000, nslice=ns)
dr6 = elements.Drift(name="dr6", ds=0.4998000, nslice=ns)

# Bend elements
rc = 10.0
sbend1 = elements.Sbend(name="sbend1", ds=2.6179938779914944, rc=rc)
e1 = elements.DipEdge(name="e1", psi=0.127409035395586, rc=rc, g=0.07, K2=0.5)
e2 = elements.DipEdge(name="e2", psi=0.127409035395586, rc=rc, g=0.07, K2=0.5)

# Quad elements
qs1f = elements.Quad(name="qs1f", ds=1.0400000, k=0.311872401, nslice=ns)
qs2d = elements.Quad(name="qs2d", ds=1.0400000, k=-0.496504354, nslice=ns)
qs3t = elements.Quad(name="qs3t", ds=0.4804000, k=0.62221964, nslice=ns)

# Sextupole elements
K2L = 0.2 * 3.0
# K2L = 0.2
sextupole = elements.Multipole(name="sextupole", multipole=3, K_normal=K2L, K_skew=0.0)

# Short RF element for bunching:
rf = elements.ShortRF(name="rf", V=2.3184782e-8, freq=1.0e4, phase=-90.0)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# Lines of interest
cell = [dr1, e1, sbend1, e2, dr2, e1, sbend1, e2, dr3, qs1f, dr4, qs2d, dr5, qs3t, dr6]
chain = 11 * cell

# Construct lattice
sim.lattice.append(monitor)
sim.lattice.append(rf)
sim.lattice.extend(cell)
sim.lattice.append(sextupole)
sim.lattice.extend(chain)
sim.lattice.append(rf)

# number of turns in the ring
sim.periods = 1000

sarr = []
test_data = []
mm_scale = 1.0e3


def hook_after_period(sim):
    s = sim.beam.ref.s
    sarr.append(s)
    beam = sim.beam.to_df()
    # Filter on particle weight (collect test particles only)
    for row in beam[beam["weighting"] == 0.0].itertuples():
        # collect test particle data
        test_data.append(
            [s, row.idcpu, row.position_x * mm_scale, row.momentum_x * mm_scale]
        )


sim.hook["after_period"] = hook_after_period


# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()

df = pd.DataFrame(test_data, columns=["s", "id", "x", "px"])
sorted_df = df.sort_values(by="id")

do_plot = True  # True to generate a plot of the test particle orbits
if do_plot:
    import matplotlib.pyplot as plt

    n = len(sarr)
    for i in range(0, len(df), n):
        subset = sorted_df.iloc[i : i + n]
        plt.scatter(subset["x"], subset["px"], s=5)

    plt.xlabel("x [mm]", fontsize=12)
    plt.ylabel("px [mrad]", fontsize=12)
    plt.title("Test Particles: Horizontal Coordinates")
    plt.show()
