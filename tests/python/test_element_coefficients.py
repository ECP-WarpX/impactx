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

    with pytest.raises(ValueError, match="same length"):
        setattr(el, first, [1.0])

    assert getattr(el, first) == [1.0, 2.0]
    assert getattr(el, second) == [0.0, 3.0]


def test_polygon_vertices():
    poly = elements.PolygonAperture(
        vertices_x=[0.0, 1.0, 1.0, 0.0, 0.0], vertices_y=[0.0, 0.0, 1.0, 1.0, 0.0]
    )
    assert poly.vertices_x == [0.0, 1.0, 1.0, 0.0, 0.0]

    poly.set_vertices([0.0, 2.0, 2.0, 0.0, 0.0], [0.0, 0.0, 2.0, 2.0, 0.0])
    assert poly.vertices_x == [0.0, 2.0, 2.0, 0.0, 0.0]

    # the polygon must stay closed
    with pytest.raises(ValueError, match="first and last vertex"):
        poly.set_vertices([0.0, 1.0, 1.0], [0.0, 0.0, 1.0])

    assert poly.vertices_x == [0.0, 2.0, 2.0, 0.0, 0.0]


def test_retuning_reaches_the_element_in_the_lattice():
    """The reason these setters were held back until the lattice shared elements."""

    sim = ImpactX()
    sim.particle_shape = 2
    sim.n_cell = [8, 8, 8]

    sq = elements.SoftQuadrupole(
        ds=1.0, gscale=1.0, cos_coefficients=[1.0, 2.0], sin_coefficients=[0.0, 3.0]
    )
    sim.lattice.append(sq)

    sq.set_coefficients([9.0, 9.0], [0.0, 0.0])

    assert sim.lattice[0].cos_coefficients == [9.0, 9.0]
