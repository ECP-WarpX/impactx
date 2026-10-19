#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import pytest

from impactx import ImpactX, distribution, elements

NSLICE = 4
DS = 0.5  # drift length in m


class HookError(RuntimeError):
    """Raised on purpose by a hook of the Programmable element."""


def _run(raise_in):
    """Track through one Programmable drift whose hook ``raise_in`` fails once.

    :param raise_in: the hook that raises: "push", "ref_particle" or "beam_particles"
    :return: how often each hook was called until the error reached the caller
    """
    sim = ImpactX()

    sim.space_charge = False
    sim.slice_step_diagnostics = False
    sim.diagnostics = False

    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(2.0e3)

    distr = distribution.Waterbag(
        lambdaX=3.9984884770e-5,
        lambdaY=3.9984884770e-5,
        lambdaT=1.0e-3,
        lambdaPx=2.6623538760e-5,
        lambdaPy=2.6623538760e-5,
        lambdaPt=2.0e-3,
        muxpx=-0.846574929020762,
        muypy=0.846574929020762,
        mutpt=0.0,
    )
    sim.add_particles(1.0e-9, distr, 10000)

    element = elements.Programmable(name="d1")
    element.ds = DS
    element.nslice = NSLICE

    fail = {"now": True}
    calls = {"push": 0, "ref_particle": 0, "beam_particles": 0}

    def called(hook):
        calls[hook] += 1
        if fail["now"] and hook == raise_in:
            raise HookError(f"raised in {hook}")

    def beam_particles(pti, refpart):
        called("beam_particles")
        soa = pti.soa().to_xp()
        soa.real["position_x"][:] += (DS / NSLICE) * soa.real["momentum_x"][:]

    def ref_particle(refpart):
        called("ref_particle")
        refpart.z = refpart.z + DS / NSLICE
        refpart.s = refpart.s + DS / NSLICE

    element.beam_particles = beam_particles
    element.ref_particle = ref_particle
    if raise_in == "push":
        # replaces the two hooks above, which push the reference particle and the beam
        element.push = lambda pc, step, period: called("push")

    sim.lattice.append(element)

    with pytest.raises(HookError, match=f"raised in {raise_in}"):
        sim.track_particles()
    calls_until_error = dict(calls)

    # the failed run leaves nothing behind: once the hook works, tracking does, too
    assert sim.tracking_element is None
    fail["now"] = False
    sim.track_particles()

    sim.finalize()

    return calls_until_error


@pytest.mark.parametrize("raise_in", ["push", "ref_particle"])
def test_exception_in_a_serial_hook_reaches_the_caller(raise_in):
    """An exception raised in a hook that ImpactX calls serially arrives in Python."""
    calls = _run(raise_in)

    assert calls[raise_in] == 1
    # the beam is pushed after the reference particle: it was never reached
    assert calls["beam_particles"] == 0


def test_exception_in_beam_particles_reaches_the_caller():
    """An exception raised in ``beam_particles`` arrives in Python.

    ImpactX calls this hook once per particle tile, from inside an OpenMP parallel
    region. An exception that left that region ended the process instead of reaching
    the caller.
    """
    calls = _run("beam_particles")

    assert calls["ref_particle"] == 1
    # the remaining tiles are skipped: a hook that just failed is not called again
    assert calls["beam_particles"] == 1
