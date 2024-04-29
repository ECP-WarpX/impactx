#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Ryan Sandberg, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
import pytest
from impactx import CoordSystem, ImpactX, coordinate_transformation, distribution


def test_transformation():
    """
    This test ensures s->t and t->s transformations
    do round-trip.
    """
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
    kin_energy_MeV = 1e3  # reference energy
    energy_gamma = kin_energy_MeV / 0.510998950 + 1.0
    bunch_charge_C = 1.0e-9  # used with space charge
    npart = 10000  # number of macro particles

    #   reference particle
    pc = sim.particle_container()
    ref = pc.ref_particle()
    ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

    #   particle bunch
    distr = distribution.Gaussian(
        lambdaX=3e-6,
        lambdaY=3e-6,
        lambdaT=1e-2,
        lambdaPx=1.33 / energy_gamma,
        lambdaPy=1.33 / energy_gamma,
        lambdaPt=100 / energy_gamma,
        muxpx=-0.5,
        muypy=0.4,
        mutpt=0.8,
    )
    sim.add_particles(bunch_charge_C, distr, npart)
    rbc_s0 = pc.reduced_beam_characteristics()

    # this must fail: we cannot transform from s to s
    with pytest.raises(Exception):
        coordinate_transformation(pc, direction=CoordSystem.s)

    # transform to t
    coordinate_transformation(pc, direction=CoordSystem.t)
    rbc_t = pc.reduced_beam_characteristics()

    # this must fail: we cannot transform from t to t
    with pytest.raises(Exception):
        coordinate_transformation(pc, direction=CoordSystem.t)

    # transform back to s
    coordinate_transformation(pc, direction=CoordSystem.s)
    rbc_s = pc.reduced_beam_characteristics()

    # finalize simulation
    sim.finalize()
    del sim

    # assert that forward-inverse transformation of the beam leaves beam unchanged
    atol = 1e-14
    rtol = 1e-10
    for key, val in rbc_s0.items():
        if not np.isclose(val, rbc_s[key], rtol=rtol, atol=atol):
            print(f"initial[{key}]={val}, final[{key}]={rbc_s[key]} not equal")
        assert np.isclose(val, rbc_s[key], rtol=rtol, atol=atol)
    # assert that the t-based beam is different, at least in the following keys:
    large_st_diff_keys = [
        "beta_x",
        "beta_y",
        "emittance_y",
        "emittance_x",
        "sig_y",
        "sig_x",
        "t_mean",
    ]
    for key in large_st_diff_keys:
        rel_error = (rbc_s0[key] - rbc_t[key]) / rbc_s0[key]
        assert abs(rel_error) > 1
