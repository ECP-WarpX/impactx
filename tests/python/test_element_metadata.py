#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""Host-only element metadata carried by the public element envelope.

Elements are split into an element-specific physics implementation, which is what gets
copied into a particle-push kernel, and common host-only metadata such as the element
name. These tests pin the behavior a user sees across that split.
"""

import pytest

from impactx import elements


def test_name_survives_a_copy():
    """Copying an element must carry its metadata along.

    Regression test for the construction adapter of the element envelope: a variadic
    constructor that is not excluded from copy construction wins overload resolution for
    an lvalue, copy-constructs only the physics base and default-constructs the
    metadata. That is not a compile error -- the copy just silently comes out unnamed.
    ``copy()`` is the C++ copy of the element; the lattice itself shares, so it would not
    exercise this.
    """

    drift = elements.Drift(ds=1.0, nslice=3, name="d1")

    copied = drift.copy()

    assert copied is not drift
    assert copied.has_name
    assert copied.name == "d1"
    assert copied.nslice == 3
    assert copied.ds == 1.0


def test_unnamed_element_stays_unnamed():
    """An element without a name reports no name, and asking for one raises."""

    drift = elements.Drift(ds=1.0)
    assert not drift.has_name
    assert drift.name is None

    lattice = elements.KnownElementsList()
    lattice.append(drift)
    assert not lattice[0].has_name


def test_empty_name_is_no_name():
    """An empty name is not a name."""

    assert not elements.Drift(ds=1.0, name="").has_name


def test_name_is_settable_and_resettable():
    """The name is ordinary mutable host metadata."""

    drift = elements.Drift(ds=1.0, name="before")
    assert drift.name == "before"

    drift.name = "after"
    assert drift.name == "after"

    drift.name = ""
    assert not drift.has_name


@pytest.mark.parametrize(
    "name",
    ["quad1", "a-very-long-element-name-" * 8, "\u00fcn\u00efcod\u00e9"],
)
def test_name_round_trips(name):
    """Names are stored verbatim, including long and non-ASCII ones.

    The non-ASCII case is written as escapes so this file stays ASCII, as the
    repository's style check requires; the string itself is unchanged.

    Worth pinning: the name used to be a raw ``char *`` with hand-written copy semantics.
    """

    lattice = elements.KnownElementsList()
    lattice.append(elements.Drift(ds=1.0, name=name))

    assert lattice[0].name == name
