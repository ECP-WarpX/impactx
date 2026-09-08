#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import math

import pytest

from impactx import Config, ImpactX, RefPart, elements

KIN_ENERGY_MEV = 250.0

# The reference particle carries pt = -gamma, and a ShortRF advances it by
# ``pt -= V * cos(phase)`` (src/elements/ShortRF.H). The gamma gain across the
# cavity is therefore exactly V * cos(phase), independent of the incoming energy.
V = 0.5
PHASE_DEG = -30.0
GAMMA_GAIN = V * math.cos(math.radians(PHASE_DEG))

if Config.precision == "SINGLE":
    RTOL = 5.0e-5
else:
    RTOL = 1.0e-12


def _ref():
    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(KIN_ENERGY_MEV)
    return ref


def _lattice_with_short_rf():
    """QF -- cavity -- QF, reusing one element name either side of the cavity."""
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Quad(name="QF", ds=0.3, k=1.0),
            elements.Drift(name="d1", ds=0.5),
            elements.ShortRF(name="cav", V=V, freq=1.3e9, phase=PHASE_DEG),
            elements.Drift(name="d2", ds=0.5),
            elements.Quad(name="QF", ds=0.3, k=1.0),
        ]
    )
    return lattice


def _lattice_with_rf_cavity():
    """A thick RFCavity: ShortRF is not the only element that accelerates."""
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Drift(name="up", ds=0.1),
            elements.RFCavity(
                name="rf",
                ds=1.0,
                escale=0.04,
                freq=7.0e8,
                phase=45.0,
                cos_coefficients=[2.0],
                sin_coefficients=[0.0],
                mapsteps=100,
                nslice=4,
            ),
            elements.Drift(name="down", ds=0.1),
        ]
    )
    return lattice


def _track_reference_through(lattice):
    """The reference particle at the lattice exit, from actual reference tracking."""
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = "false"
    sim.slice_step_diagnostics = False
    sim.diagnostics = False
    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(KIN_ENERGY_MEV)
    sim.lattice.extend(lattice)

    try:
        sim.track_reference(ref)
        return ref.s, ref.gamma, ref.rigidity_Tm
    finally:
        sim.finalize()


def test_map_trace_carries_reference_particle():
    """Every map_trace entry reports the reference particle at that element's exit."""
    ref = _ref()
    lattice = _lattice_with_short_rf()

    trace = lattice.map_trace(ref)

    # one entry per element, plus the leading <start> entry
    assert len(trace) == len(lattice) + 1
    assert all("ref" in entry for entry in trace)

    # the leading entry carries the incoming reference particle unchanged
    assert trace[0]["type"] == "<start>"
    assert trace[0]["ref"].gamma == pytest.approx(ref.gamma, rel=RTOL)

    # s on the entry and on its reference particle agree
    for entry in trace:
        assert entry["ref"].s == pytest.approx(entry["s"], rel=RTOL, abs=1.0e-12)

    # the caller's reference particle is not modified in place
    assert ref.gamma == pytest.approx(_ref().gamma, rel=RTOL)
    assert ref.s == pytest.approx(0.0, abs=1.0e-12)


def test_reference_energy_changes_across_short_rf():
    """The gamma gain across a ShortRF matches the analytic V*cos(phase)."""
    lattice = _lattice_with_short_rf()
    trace = lattice.map_trace(_ref())
    by_name = {entry["name"]: entry for entry in trace if entry["name"]}

    gamma_before = by_name["d1"]["ref"].gamma
    gamma_after = by_name["cav"]["ref"].gamma

    assert gamma_after - gamma_before == pytest.approx(GAMMA_GAIN, rel=RTOL)

    # a drift does not change the reference energy
    assert by_name["d2"]["ref"].gamma == pytest.approx(gamma_after, rel=RTOL)


@pytest.mark.parametrize(
    "make_lattice, last_element",
    [
        (_lattice_with_short_rf, ("QF", 2)),
        (_lattice_with_rf_cavity, ("down", 1)),
    ],
    ids=["ShortRF", "RFCavity"],
)
def test_ref_at_matches_track_reference(make_lattice, last_element):
    """ref_at() at the lattice exit reproduces what tracking the reference produces.

    Parametrized over a thin and a thick accelerating element on purpose: anything
    that reconstructs the reference energy from ShortRF alone passes the first case
    and fails the second.
    """
    lattice = make_lattice()
    name, occurrence = last_element

    walked = lattice.ref_at(_ref(), name, occurrence=occurrence)
    s, gamma, rigidity_Tm = _track_reference_through(lattice)

    assert walked.s == pytest.approx(s, rel=RTOL)
    assert walked.gamma == pytest.approx(gamma, rel=RTOL)
    assert walked.rigidity_Tm == pytest.approx(rigidity_Tm, rel=RTOL)


def test_rigidity_at_differs_across_cavity():
    """rigidity_at() is the strength/field conversion constant, and a cavity moves it."""
    ref = _ref()
    lattice = _lattice_with_short_rf()

    brho_upstream = lattice.rigidity_at(ref, "QF")
    brho_downstream = lattice.rigidity_at(ref, "QF", occurrence=2)

    # Brho = m * beta * gamma * c / q carries the sign of the charge, so it is
    # negative for an electron; the cavity raises its magnitude.
    assert brho_upstream < 0.0
    assert abs(brho_downstream) > abs(brho_upstream)

    # Same magnet name, different rigidity: exactly the case that makes a single
    # lattice-entry rigidity the wrong conversion constant downstream of a cavity.
    assert brho_downstream != pytest.approx(brho_upstream, rel=1.0e-6)

    # Brho = m * beta * gamma * c / q, so it follows the reference particle.
    after = lattice.ref_at(ref, "QF", occurrence=2)
    assert brho_downstream == pytest.approx(after.rigidity_Tm, rel=RTOL)


def test_ref_at_name_errors():
    """Unknown names and out-of-range occurrences are errors, not silent picks."""
    ref = _ref()
    lattice = _lattice_with_short_rf()

    with pytest.raises(ValueError, match="No element named 'nope'"):
        lattice.ref_at(ref, "nope")

    with pytest.raises(ValueError, match="occurs 2 time"):
        lattice.ref_at(ref, "QF", occurrence=3)

    with pytest.raises(ValueError, match="occurrence must be >= 1"):
        lattice.ref_at(ref, "QF", occurrence=0)
