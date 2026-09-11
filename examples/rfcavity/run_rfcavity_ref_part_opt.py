#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Marco Garten, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
from phase_opt import optimize

from impactx import ImpactX, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# reference kinetic energy (initial)
kin_energy_MeV = 230.0  # reference energy

#   reference particle
ref = sim.beam.ref
ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)

# design the accelerator lattice

# access RF cavity on-axis field data
data_in = np.loadtxt("onaxis_data.in")
z = data_in[:, 0]
ez_onaxis = data_in[:, 1]
ncoef = 25

#   Drift elements
dr1 = elements.Drift(name="dr1", ds=0.4, nslice=1)
dr2 = elements.Drift(name="dr2", ds=0.032997, nslice=1)

#   RF cavity elements: four physically distinct cavities, built from one template so
#   that the on-axis field is fitted once
rf_tpl = elements.RFCavity(
    name="rf1",
    ds=1.31879807,
    escale=20.0,
    z=z,
    field_on_axis=ez_onaxis,
    ncoef=ncoef,
    freq=1.3e9,
    phase=0.0,
    mapsteps=100,
    nslice=4,
)
rf = [rf_tpl.copy(name=f"rf{i}") for i in range(1, 5)]


# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

sim.lattice.extend(
    [
        monitor,
        dr1,
        dr2,
        rf[0],
        dr2,
        dr2,
        rf[1],
        dr2,
        dr2,
        rf[2],
        dr2,
        dr2,
        rf[3],
        dr2,
        monitor,
    ]
)


def hook_before_element(sim):
    element = sim.tracking_element
    if type(element) is elements.RFCavity:
        beam = sim.beam
        ref = beam.ref
        print(
            f"  Beam at s={ref.s:.2f}m, t={ref.t:.2f}s, gamma at entry={ref.gamma:.2f}",
            flush=True,
        )

        # Interpret input RF phase as measured relative to max accelerating phase:
        phase_shift = element.phase

        # Find RF phase that maximizes energy gain:
        phase_opt, e_gain = optimize(ref, element)

        # Reset input RF phase to appropriate value for tracking:
        element.phase = phase_opt + phase_shift
        print(
            f"  RF cavity updated (reset) values of phase={phase_opt:.2f}, energy gain (MeV) ={e_gain:.2f}",
            flush=True,
        )


sim.hook["before_element"] = hook_before_element

# run simulation
sim.track_reference(ref)

# clean shutdown
sim.finalize()
