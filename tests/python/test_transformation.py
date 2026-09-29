#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Ryan Sandberg, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
import pytest

from impactx import (
    Config,
    CoordSystem,
    ImpactX,
    coordinate_transformation,
    distribution,
)


def _assert_moments_close(expected, actual):
    """Compare two beam_moments() dicts to round-trip precision."""
    if Config.precision == "SINGLE":
        atol = 1e-6
        rtol = 5e-6
    else:
        atol = 1e-14
        rtol = 1e-10
    for key, val in expected.items():
        if not np.isclose(val, actual[key], rtol=rtol, atol=atol):
            print(f"initial[{key}]={val}, final[{key}]={actual[key]} not equal")
            assert False


def _make_beam():
    """Create a simulation with a long, correlated 1 GeV electron beam at fixed s."""
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
    beam = sim.beam
    ref = beam.ref
    ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

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

    return sim, beam


def test_transformation():
    """
    This test ensures s->t and t->s transformations
    do round-trip.
    """
    sim, beam = _make_beam()
    rbc_s0 = beam.beam_moments()

    # this must fail: we cannot transform from s to s
    with pytest.raises(Exception):
        coordinate_transformation(beam, direction=CoordSystem.s)

    # transform to t
    coordinate_transformation(beam, direction=CoordSystem.t)
    rbc_t = beam.beam_moments()

    # this must fail: we cannot transform from t to t
    with pytest.raises(Exception):
        coordinate_transformation(beam, direction=CoordSystem.t)

    # transform back to s
    coordinate_transformation(beam, direction=CoordSystem.s)
    rbc_s = beam.beam_moments()

    # finalize simulation
    sim.finalize()

    # assert that forward-inverse transformation of the beam leaves beam unchanged
    _assert_moments_close(rbc_s0, rbc_s)
    # assert that the t-based beam is different, at least in the following keys:
    large_st_diff_keys = [
        "beta_x",
        "beta_y",
        "emittance_y",
        "emittance_x",
        "sigma_y",
        "sigma_x",
        "mean_t",
    ]
    for key in large_st_diff_keys:
        rel_error = (rbc_s0[key] - rbc_t[key]) / rbc_s0[key]
        assert abs(rel_error) > 1


def test_at_fixed_t():
    """
    ``beam.at_fixed_t()`` matches the paired coordinate_transformation calls
    and always returns the beam to fixed s.
    """
    sim, beam = _make_beam()
    rbc_s0 = beam.beam_moments()

    # same result as the explicit transformation
    coordinate_transformation(beam, direction=CoordSystem.t)
    rbc_t_explicit = beam.beam_moments()
    coordinate_transformation(beam, direction=CoordSystem.s)

    assert beam.coord_system == CoordSystem.s
    with beam.at_fixed_t() as beam_t:
        assert beam_t is beam
        assert beam.coord_system == CoordSystem.t
        rbc_t = beam.beam_moments()
    assert beam.coord_system == CoordSystem.s
    _assert_moments_close(rbc_t_explicit, rbc_t)
    _assert_moments_close(rbc_s0, beam.beam_moments())

    # an exception inside the block is re-raised, with the beam back at fixed s
    class Oops(Exception):
        pass

    with pytest.raises(Oops):
        with beam.at_fixed_t():
            raise Oops()
    assert beam.coord_system == CoordSystem.s
    _assert_moments_close(rbc_s0, beam.beam_moments())

    # entering requires fixed s, so blocks cannot be nested
    with beam.at_fixed_t():
        with pytest.raises(RuntimeError, match="must be at fixed s"):
            with beam.at_fixed_t():
                pass
        assert beam.coord_system == CoordSystem.t
    assert beam.coord_system == CoordSystem.s

    # transforming back by hand inside the block is an error on exit
    with pytest.raises(RuntimeError, match="transformed out of fixed t"):
        with beam.at_fixed_t():
            coordinate_transformation(beam, direction=CoordSystem.s)
    assert beam.coord_system == CoordSystem.s

    # ... but does not hide an exception that is already in flight
    with pytest.raises(Oops):
        with beam.at_fixed_t():
            coordinate_transformation(beam, direction=CoordSystem.s)
            raise Oops()
    assert beam.coord_system == CoordSystem.s
    _assert_moments_close(rbc_s0, beam.beam_moments())

    sim.finalize()
