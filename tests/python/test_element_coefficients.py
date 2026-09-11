#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""Retuning the array-valued parameters of an element after it was built.

These were read-only while the lattice copied elements: a write through the user's handle
could not reach the element being tracked. The lattice shares elements now, so it can.
"""

import numpy as np
import pytest

from impactx import ImpactX, elements

# element factory -> (first array property, second array property, paired setter)
CASES = [
    (
        lambda: elements.SoftQuadrupole(
            ds=1.0, gscale=1.0, cos_coefficients=[1.0, 2.0], sin_coefficients=[0.0, 3.0]
        ),
        "cos_coefficients",
        "sin_coefficients",
        "set_coefficients",
    ),
    (
        lambda: elements.SoftSolenoid(
            ds=1.0, bscale=1.0, cos_coefficients=[1.0, 2.0], sin_coefficients=[0.0, 3.0]
        ),
        "cos_coefficients",
        "sin_coefficients",
        "set_coefficients",
    ),
    (
        lambda: elements.RFCavity(
            ds=1.0,
            escale=1.0,
            freq=1.0e9,
            phase=0.0,
            cos_coefficients=[1.0, 2.0],
            sin_coefficients=[0.0, 3.0],
        ),
        "cos_coefficients",
        "sin_coefficients",
        "set_coefficients",
    ),
    (
        lambda: elements.ExactMultipole(ds=1.0, k_normal=[1.0, 2.0], k_skew=[0.0, 3.0]),
        "k_normal",
        "k_skew",
        "set_coefficients",
    ),
    (
        lambda: elements.ExactCFbend(ds=1.0, k_normal=[1.0, 2.0], k_skew=[0.0, 3.0]),
        "k_normal",
        "k_skew",
        "set_coefficients",
    ),
]

IDS = ["SoftQuadrupole", "SoftSolenoid", "RFCavity", "ExactMultipole", "ExactCFbend"]


@pytest.mark.parametrize("make,first,second,setter", CASES, ids=IDS)
def test_arrays_are_readable(make, first, second, setter):
    el = make()
    assert getattr(el, first) == [1.0, 2.0]
    assert getattr(el, second) == [0.0, 3.0]


@pytest.mark.parametrize("make,first,second,setter", CASES, ids=IDS)
def test_setting_one_array_keeps_the_other(make, first, second, setter):
    el = make()
    setattr(el, first, [5.0, 6.0])

    assert getattr(el, first) == [5.0, 6.0]
    assert getattr(el, second) == [0.0, 3.0]


@pytest.mark.parametrize("make,first,second,setter", CASES, ids=IDS)
def test_paired_setter_can_change_the_length(make, first, second, setter):
    el = make()
    getattr(el, setter)([1.0, 2.0, 3.0], [0.0, 0.0, 0.0])

    assert getattr(el, first) == [1.0, 2.0, 3.0]
    assert getattr(el, second) == [0.0, 0.0, 0.0]


@pytest.mark.parametrize("make,first,second,setter", CASES, ids=IDS)
def test_mismatched_lengths_are_rejected_and_change_nothing(
    make, first, second, setter
):
    el = make()

    # through a single property, measured against the array it keeps
    with pytest.raises(ValueError, match="same length"):
        setattr(el, first, [1.0])

    # and through the paired setter
    with pytest.raises(ValueError, match="same length"):
        getattr(el, setter)([1.0, 2.0], [1.0])

    assert getattr(el, first) == [1.0, 2.0]
    assert getattr(el, second) == [0.0, 3.0]


def test_polygon_vertices():
    """The vertex arrays follow the same pattern; their validation is tested with the
    other array-valued parameters."""

    poly = elements.PolygonAperture(
        vertices_x=[0.0, 1.0, 1.0, 0.0, 0.0], vertices_y=[0.0, 0.0, 1.0, 1.0, 0.0]
    )
    assert poly.vertices_x == [0.0, 1.0, 1.0, 0.0, 0.0]

    poly.set_vertices([0.0, 2.0, 2.0, 0.0, 0.0], [0.0, 0.0, 2.0, 2.0, 0.0])
    assert poly.vertices_x == [0.0, 2.0, 2.0, 0.0, 0.0]
    assert poly.vertices_y == [0.0, 0.0, 2.0, 2.0, 0.0]


def test_retuning_reaches_the_element_in_the_lattice():
    """The reason these setters were held back until the lattice shared elements."""

    sq = elements.SoftQuadrupole(
        ds=1.0, gscale=1.0, cos_coefficients=[1.0, 2.0], sin_coefficients=[0.0, 3.0]
    )
    lattice = elements.KnownElementsList([sq])

    sq.set_coefficients([9.0, 9.0], [0.0, 0.0])

    assert lattice[0].cos_coefficients == [9.0, 9.0]


@pytest.fixture
def push_once():
    """Push a deterministic beam, resetting its particles and reference on every call."""
    sim = ImpactX()
    sim.particle_shape = 2
    sim.n_cell = [8, 8, 8]
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    def run(element):
        beam = sim.beam
        beam.clear_particles()
        beam.ref.reset()
        beam.ref.set_species("electron").set_kin_energy_MeV(100.0)
        # Identical local particles on each rank; no random sampling or host access to
        # device arrays. add_n_particles uploads the inputs on GPU builds.
        beam.add_n_particles(
            [0.0, 1.0e-3, 3.0e-3],
            [0.0, 2.0e-3, -1.0e-3],
            [0.0, 3.0e-3, -2.0e-3],
            [1.0e-4, 2.0e-4, -1.0e-4],
            [2.0e-4, -1.0e-4, 3.0e-4],
            [3.0e-4, 1.0e-4, -2.0e-4],
            -1.0 / 0.510998950e6,
            1.0e-12,
        )
        element.push(beam)
        # to_df copies device data to the host and synchronizes before comparison.
        phase_space = beam.to_df(local=True)[
            [
                "position_x",
                "position_y",
                "position_t",
                "momentum_x",
                "momentum_y",
                "momentum_t",
            ]
        ].to_numpy(copy=True)
        reference = np.array(
            [
                getattr(beam.ref, name)
                for name in ("x", "y", "z", "t", "px", "py", "pz", "pt", "s")
            ]
        )
        return (
            phase_space,
            reference,
            beam.total_number_of_particles(only_valid=True, only_local=True),
        )

    try:
        yield run
    finally:
        sim.finalize()


@pytest.mark.parametrize(
    "name,settings",
    [
        ("SoftQuadrupole", {"gscale": 1.0}),
        ("SoftSolenoid", {"bscale": 1.0}),
        ("RFCavity", {"escale": 1.0, "freq": 1.0e9, "phase": 10.0}),
        ("ExactMultipole", {}),
        ("ExactCFbend", {}),
    ],
    ids=IDS,
)
@pytest.mark.parametrize(
    "first,second",
    [
        ([2.0, 3.0], [0.0, 0.1]),
        ([2.0, 3.0, 0.4], [0.0, 0.1, 0.2]),
        ([2.0], [0.0]),
    ],
    ids=["same_length", "grow", "shrink"],
)
def test_coefficient_update_changes_the_particle_push(
    name, settings, first, second, push_once
):
    """A warmed-up element must push like a fresh element after its arrays change."""
    cls = getattr(elements, name)
    keys = (
        ("k_normal", "k_skew")
        if name in ("ExactMultipole", "ExactCFbend")
        else ("cos_coefficients", "sin_coefficients")
    )
    element = cls(ds=0.1, **settings, **{keys[0]: [1.0, 2.0], keys[1]: [0.0, 3.0]})

    # Populate the execution-space cache before changing values or reallocating arrays.
    before, _, _ = push_once(element)
    element.set_coefficients(first, second)
    after, ref_after, alive = push_once(element)

    fresh = cls(ds=0.1, **settings, **{keys[0]: first, keys[1]: second})
    expected, ref_expected, expected_alive = push_once(fresh)

    assert np.isfinite(after).all()
    assert np.isfinite(ref_after).all()
    assert not np.array_equal(before, after), (
        "The update must change the particle motion"
    )
    # Both pushes execute identical arithmetic on the same backend; no CPU/GPU
    # cross-comparison or statistical tolerance is needed.
    np.testing.assert_array_equal(after, expected)
    np.testing.assert_array_equal(ref_after, ref_expected)
    assert alive == expected_alive == 3


@pytest.mark.parametrize(
    "vertices_x,vertices_y,expected_alive",
    [
        (
            [-0.0015, 0.0015, 0.0015, -0.0015, -0.0015],
            [-0.004, -0.004, 0.004, 0.004, -0.004],
            2,
        ),
        (
            [-0.0015, 0.0, 0.0015, 0.0015, -0.0015, -0.0015],
            [-0.004, -0.004, -0.004, 0.004, 0.004, -0.004],
            2,
        ),
        (
            [-0.0015, 0.0015, 0.0, -0.0015],
            [-0.004, -0.004, 0.004, -0.004],
            1,
        ),
    ],
    ids=["same_length", "grow", "shrink"],
)
def test_vertex_update_changes_particle_losses(
    vertices_x, vertices_y, expected_alive, push_once
):
    """Updated polygon vertices must change which particles are transmitted."""
    polygon = elements.PolygonAperture(
        vertices_x=[-0.01, 0.01, 0.01, -0.01, -0.01],
        vertices_y=[-0.01, -0.01, 0.01, 0.01, -0.01],
    )
    _, _, before_alive = push_once(polygon)
    assert before_alive == 3

    polygon.set_vertices(vertices_x, vertices_y)
    _, _, after_alive = push_once(polygon)
    fresh = elements.PolygonAperture(vertices_x=vertices_x, vertices_y=vertices_y)
    _, _, fresh_alive = push_once(fresh)

    assert after_alive == fresh_alive == expected_alive
