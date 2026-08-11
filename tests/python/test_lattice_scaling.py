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


def build(n):
    lattice = elements.KnownElementsList()
    lattice.extend([elements.Drift(ds=0.1, name=f"d{i}") for i in range(n)])
    return lattice


def time_it(setup, operation, repeats=5):
    """Best of a few runs of ``operation``, with ``setup`` excluded.

    The fastest run is the one least disturbed by whatever else the machine is doing.
    """

    best = float("inf")
    for _ in range(repeats):
        subject = setup()
        start = time.perf_counter()
        operation(subject)
        best = min(best, time.perf_counter() - start)
    return best


def growth(operation, make=build, small_n=SMALL):
    """How much the cost grows when the lattice grows by @see SIZE_RATIO."""

    small = time_it(lambda: make(small_n), operation)
    large = time_it(lambda: make(SIZE_RATIO * small_n), operation)
    # a floor keeps a very fast operation from dividing two timer-resolution numbers
    return large / max(small, 1e-5)


@pytest.mark.parametrize(
    ("name", "operation"),
    [
        ("delete_leading_half", lambda lat: lat.__delitem__(slice(0, len(lat) // 2))),
        ("delete_every_second", lambda lat: lat.__delitem__(slice(None, None, 2))),
        ("delete_all", lambda lat: lat.__delitem__(slice(None))),
        ("delete_reversed", lambda lat: lat.__delitem__(slice(None, None, -1))),
    ],
)
def test_slice_deletion_is_linear(name, operation):
    measured = growth(operation)

    assert measured < LINEAR_TOLERANCE


@pytest.mark.parametrize(
    ("name", "operation"),
    [
        (
            "replace_all",
            lambda lat: lat.__setitem__(
                slice(None), [elements.Drift(ds=0.2) for _ in range(len(lat))]
            ),
        ),
        (
            "replace_leading_half",
            lambda lat: lat.__setitem__(
                slice(0, len(lat) // 2),
                [elements.Drift(ds=0.2) for _ in range(len(lat) // 2)],
            ),
        ),
        (
            "prepend",
            lambda lat: lat.__setitem__(
                slice(0, 0), [elements.Drift(ds=0.2) for _ in range(len(lat))]
            ),
        ),
        (
            "replace_every_second",
            lambda lat: lat.__setitem__(
                slice(None, None, 2),
                [elements.Drift(ds=0.2) for _ in range(len(range(0, len(lat), 2)))],
            ),
        ),
    ],
)
def test_slice_assignment_is_linear(name, operation):
    measured = growth(operation)

    assert measured < LINEAR_TOLERANCE


def test_building_a_lattice_is_linear():
    measured = growth(build, make=lambda n: n)

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

    measured = growth(
        lambda lattice: lattice.select(kind="Drift").delete(),
        make=alternating,
        small_n=FILTERED_DELETE_SIZE,
    )

    assert measured < LINEAR_TOLERANCE


# --- the operations above must also still be correct -------------------------------


def names_of(lattice):
    return [element.name for element in lattice]


@pytest.mark.parametrize(
    ("count", "key"),
    [
        (8, slice(0, 4)),
        (8, slice(None, None, 2)),
        (8, slice(None)),
        (8, slice(None, None, -1)),
        (8, slice(3, 7)),
        (9, slice(None, None, 3)),
    ],
)
def test_slice_deletion_matches_a_list(count, key):
    lattice = build(count)
    reference = [f"d{i}" for i in range(count)]

    del lattice[key]
    del reference[key]

    assert names_of(lattice) == reference


@pytest.mark.parametrize(
    ("count", "key", "replacements"),
    [
        (6, slice(0, 3), 3),
        (6, slice(0, 3), 5),
        (6, slice(0, 3), 1),
        (6, slice(None), 2),
        (6, slice(0, 0), 2),
        (6, slice(6, 6), 2),
        (6, slice(2, 4), 0),
    ],
)
def test_contiguous_slice_assignment_matches_a_list(count, key, replacements):
    lattice = build(count)
    reference = [f"d{i}" for i in range(count)]
    new = [elements.Drift(ds=0.2, name=f"n{i}") for i in range(replacements)]

    lattice[key] = new
    reference[key] = [element.name for element in new]

    assert names_of(lattice) == reference


def test_a_rebuild_keeps_the_exact_objects():
    """The elements that survive a slice edit are the same objects, not copies."""

    class Tagged(elements.Drift):
        def __init__(self, **kwargs):
            super().__init__(**kwargs)
            self.tag = "kept"

    lattice = elements.KnownElementsList()
    keep_first = Tagged(ds=0.1, name="a")
    keep_last = Tagged(ds=0.2, name="b")
    lattice.extend([keep_first, elements.Drift(ds=0.3, name="gone"), keep_last])

    del lattice[1:2]

    assert lattice[0] is keep_first
    assert lattice[1] is keep_last
    assert lattice[1].tag == "kept"
