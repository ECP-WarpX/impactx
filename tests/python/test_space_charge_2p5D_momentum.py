#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
import pytest

from impactx import Config, ImpactX, elements


@pytest.mark.parametrize("num_bins", [2, 33, 129])
@pytest.mark.parametrize("space_charge", ["Gauss2p5D", "2p5D"])
def test_space_charge_2p5D_momentum_conservation(space_charge, num_bins):
    """
    The longitudinal 2.5D space-charge kick conserves the total momentum.

    The longitudinal kick of a particle is the slope of the line charge density
    at its position times a function of its transverse position. Here, every
    sample of a skewed longitudinal profile carries the same transverse points,
    so the net kick is the sum of the line charge density slopes over all
    longitudinal samples times a transverse sum. The first sum vanishes for a
    momentum-conserving charge deposition, slope and gather.
    """
    sim = ImpactX()

    sim.particle_shape = 2
    sim.space_charge = space_charge
    sim.space_charge_apply_longitudinal_kick = True
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    if space_charge == "2p5D":
        sim.poisson_solver = "fft"
        sim.space_charge_num_longitudinal_bins = num_bins
        sim.max_level = 0
        sim.n_cell = [64, 64, 1]
        sim.blocking_factor_x = [16]
        sim.blocking_factor_y = [16]
        sim.blocking_factor_z = [1]
        sim.prob_relative = [1.1]
    else:
        sim.space_charge_gauss_nint = 101
        sim.space_charge_gauss_charge_z_bins = num_bins
        sim.space_charge_gauss_taylor_delta = 0.01
        sim.space_charge_gauss_long_scale = 6.0
    sim.init_grids()

    kin_energy_MeV = 100.0
    bunch_charge_C = 1.0e-9
    sim.beam.ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)
    qm_eev = -1.0 / 0.51099895000e6  # electron charge/mass in e / eV

    # the same transverse points at each longitudinal sample, all at rest
    rng = np.random.default_rng(seed=42)
    num_transverse = 32
    num_longitudinal = 2048
    x_m = rng.normal(0.0, 4.0e-5, num_transverse)
    y_m = rng.normal(0.0, 4.0e-5, num_transverse)
    g = rng.gamma(2.0, 1.0, num_longitudinal)
    t_j = (g - 2.0) / np.sqrt(2.0) * 1.0e-3  # skewed, rms 1 mm

    x = np.tile(x_m, num_longitudinal)
    y = np.tile(y_m, num_longitudinal)
    t = np.repeat(t_j, num_transverse)
    zeros = np.zeros_like(x)
    sim.beam.add_n_particles(x, y, t, zeros, zeros, zeros, qm_eev, bunch_charge_C)

    # one slice: one space-charge kick, then a drift that does not change pt
    sim.lattice.append(elements.Drift(name="drift", ds=1.0, nslice=1))
    sim.track_particles()

    # sums in double precision: only the round-off of the kick remains
    pt = sim.beam.to_df(local=True)["momentum_t"].to_numpy().astype(np.float64)
    sim.finalize()

    rms_pt = np.sqrt(np.mean(pt**2))
    net_pt = np.mean(pt)
    print(
        f"{space_charge} num_bins={num_bins}: mean(pt)/rms(pt) = {net_pt / rms_pt:.3e}"
    )

    # the beam got a longitudinal kick ...
    assert rms_pt > 0.0
    # ... with a zero sum up to round-off; in single precision, this is dominated by
    # the summation of the charge deposit, which is largest for few bins (up to ~1e-4)
    rtol = 1.0e-3 if Config.precision == "SINGLE" else 1.0e-10
    assert abs(net_pt) < rtol * rms_pt
