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
    ImpactXParIter,
    coordinate_transformation,
    distribution,
    elements,
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

    def raise_oops():
        # raising via a call keeps the code after pytest.raises reachable for linters
        raise Oops()

    with pytest.raises(Oops):
        with beam.at_fixed_t():
            raise_oops()
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
            raise_oops()
    assert beam.coord_system == CoordSystem.s
    _assert_moments_close(rbc_s0, beam.beam_moments())

    sim.finalize()


def test_names_follow_coord_system():
    """
    The longitudinal attributes are named after the coordinates they hold:
    ``position_t``/``momentum_t`` at fixed s, ``position_z``/``momentum_z`` at fixed t.
    """
    sim, beam = _make_beam()

    names_s = {"position_t", "momentum_t"}
    names_t = {"position_z", "momentum_z"}

    def assert_names(present, absent):
        for names in (
            set(beam.real_soa_names),
            set(beam.to_df(local=True).columns),
            *(set(pti.soa().to_xp().real) for pti in ImpactXParIter(beam, level=0)),
        ):
            assert present <= names, f"{present} not in {names}"
            assert not absent & names, f"{absent & names} in {names}"

    assert_names(names_s, names_t)
    df_s = beam.to_df(local=True)

    with beam.at_fixed_t():
        assert_names(names_t, names_s)
        df_t = beam.to_df(local=True)
    assert_names(names_s, names_t)

    # the explicit transformation renames as well
    coordinate_transformation(beam, direction=CoordSystem.t)
    assert_names(names_t, names_s)
    coordinate_transformation(beam, direction=CoordSystem.s)
    assert_names(names_s, names_t)

    # the transverse attributes keep their names; z and pz are different data than t and pt
    assert np.array_equal(df_s["weighting"], df_t["weighting"])
    assert not np.allclose(df_s["position_t"], df_t["position_z"])

    # consumers of fixed-s data refuse fixed-t data instead of mislabeling it
    monitor = elements.BeamMonitor("monitor_fixed_t", backend="h5")
    with beam.at_fixed_t():
        with pytest.raises(RuntimeError, match="must be at fixed s"):
            beam.plot_phasespace()
        with pytest.raises(RuntimeError, match="must be at fixed s"):
            monitor.push(beam)
    monitor.finalize()

    sim.finalize()
