#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
"""Lattice operations stay linear in the number of elements.

A long beamline -- a ring, or anything that has been through
``insert_element_every_ds`` -- has many thousands of elements. An operation that moves
the tail of the lattice once per position is quadratic, which is unnoticeable in a unit
test and ruinous at scale. These tests compare the cost at two sizes rather than a wall
time, so they mean the same thing on any machine.

Timing on a shared build machine is noisy, so the two sizes are far apart: growing the
lattice fourfold costs a linear operation 4x and a quadratic one 16x. Judging against 8
therefore leaves a factor of two of headroom on either side. Comparing neighbouring
sizes would leave only the span from 2 to 4, which noise alone can cover.
"""

import time

import pytest

from impactx import elements

#: how much longer the larger lattice is
SIZE_RATIO = 4

#: linear costs SIZE_RATIO, quadratic costs SIZE_RATIO squared; judge halfway between
LINEAR_TOLERANCE = 8.0

SMALL = 2000


def lattice_of(n):
    lattice = elements.KnownElementsList()
    lattice.extend([elements.Drift(ds=0.1, name=f"d{i}") for i in range(n)])
    return lattice


def time_it(setup, operation, repeats=5):
    """Best of a few runs of ``operation``, with ``setup`` excluded.

    ``setup`` returns the arguments ``operation`` is called with. Everything the
    operation needs -- the lattice, and any elements written into it -- is built there,
    outside the timed region: constructing an element costs more than moving a handle,
    so timing the construction would hide a quadratic move behind a linear build.

    The fastest run is the one least disturbed by whatever else the machine is doing.
    """

    best = float("inf")
    for _ in range(repeats):
        arguments = setup()
        start = time.perf_counter()
        operation(*arguments)
        best = min(best, time.perf_counter() - start)
    return best


def growth(operation, make, small_n=SMALL):
    """How much the cost grows when the lattice grows by @see SIZE_RATIO.

    ``make(n)`` returns the arguments ``operation`` is timed with for ``n`` elements.
    """

    small = time_it(lambda: make(small_n), operation)
    large = time_it(lambda: make(SIZE_RATIO * small_n), operation)
    # a floor keeps a very fast operation from dividing two timer-resolution numbers
    return large / max(small, 1e-5)


@pytest.mark.parametrize(
    ("name", "key_of"),
    [
        ("delete_leading_half", lambda n: slice(0, n // 2)),
        ("delete_every_second", lambda n: slice(None, None, 2)),
        ("delete_all", lambda n: slice(None)),
        ("delete_reversed", lambda n: slice(None, None, -1)),
    ],
)
def test_slice_deletion_is_linear(name, key_of):
    measured = growth(
        lambda lattice, key: lattice.__delitem__(key),
        make=lambda n: (lattice_of(n), key_of(n)),
    )

    assert measured < LINEAR_TOLERANCE


@pytest.mark.parametrize(
    ("name", "key_of", "count_of"),
    [
        ("replace_all", lambda n: slice(None), lambda n: n),
        ("replace_leading_half", lambda n: slice(0, n // 2), lambda n: n // 2),
        ("prepend", lambda n: slice(0, 0), lambda n: n),
        (
            "replace_every_second",
            lambda n: slice(None, None, 2),
            lambda n: len(range(0, n, 2)),
        ),
    ],
)
def test_slice_assignment_is_linear(name, key_of, count_of):
    def make(n):
        replacements = [elements.Drift(ds=0.2) for _ in range(count_of(n))]
        return lattice_of(n), key_of(n), replacements

    measured = growth(
        lambda lattice, key, replacements: lattice.__setitem__(key, replacements),
        make=make,
    )

    assert measured < LINEAR_TOLERANCE


def test_building_a_lattice_is_linear():
    def make(n):
        return elements.KnownElementsList(), [
            elements.Drift(ds=0.1, name=f"d{i}") for i in range(n)
        ]

    measured = growth(lambda lattice, drifts: lattice.extend(drifts), make=make)

    assert measured < LINEAR_TOLERANCE


#: the filtered delete needs a longer lattice than the others before its cost separates:
#: the linear part of the work it used to do dominated until the lattice was very long
FILTERED_DELETE_SIZE = 16000


def test_filtered_delete_is_linear():
    """A scattered selection is the case that hides the cost.

    When every position is selected the removal is contiguous from the back and looks
    linear whatever the implementation; alternating kinds is what exposes it.
    """

    def alternating(n):
        lattice = elements.KnownElementsList()
        lattice.extend(
            [
                elements.Drift(ds=0.1, name=f"d{i}")
                if i % 2 == 0
                else elements.Quad(ds=0.1, k=1.0, name=f"q{i}")
                for i in range(n)
            ]
        )
        return lattice

    # the selection is taken outside the timed region: only the delete is measured
    measured = growth(
        lambda selection: selection.delete(),
        make=lambda n: (alternating(n).select(kind="Drift"),),
        small_n=FILTERED_DELETE_SIZE,
    )

    assert measured < LINEAR_TOLERANCE
