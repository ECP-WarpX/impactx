#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""The lattice holds elements, it does not copy them.

``sim.lattice`` follows Python list rules: appending stores the object given, so the
caller keeps a handle to the very element that will be tracked. Distinct physical
elements come from constructing one per position.
"""

import gc

import pytest

from impactx import elements


def test_append_keeps_the_object():
    q = elements.Quad(ds=0.3, k=2.0, name="q1")

    lattice = elements.KnownElementsList()
    lattice.append(q)

    assert lattice[0] is q


def test_changes_are_visible_through_the_lattice():
    q = elements.Quad(ds=0.3, k=2.0, name="q1")
    lattice = elements.KnownElementsList()
    lattice.append(q)

    q.k = 3.0
    assert lattice[0].k == 3.0

    lattice[0].k = 4.0
    assert q.k == 4.0


def test_repeated_insertion_is_one_element_twice():
    """``[q, q]`` is two occurrences of one element, not two elements."""

    q = elements.Quad(ds=0.3, k=2.0, name="q1")
    lattice = elements.KnownElementsList()
    lattice.extend([q, q])

    assert len(lattice) == 2
    assert lattice[0] is q
    assert lattice[1] is q

    # one element, so retuning it retunes both occurrences
    q.k = 7.0
    assert lattice[0].k == 7.0
    assert lattice[1].k == 7.0


def test_fresh_construction_gives_independent_elements():
    """The documented way to get distinct physical elements."""

    lattice = elements.KnownElementsList()
    lattice.append(elements.Quad(ds=0.3, k=2.0))
    lattice.append(elements.Quad(ds=0.3, k=2.0))

    assert lattice[0] is not lattice[1]

    lattice[0].k = 9.0
    assert lattice[1].k == 2.0


def test_element_survives_losing_every_external_reference():
    """The lattice keeps the element alive, and keeps it the same object."""

    lattice = elements.KnownElementsList()
    lattice.append(elements.Quad(ds=0.3, k=2.0, name="only-in-lattice"))
    gc.collect()

    assert lattice[0].name == "only-in-lattice"
    assert lattice[0] is lattice[0]


def test_python_subclass_and_attributes_survive():
    """A Python subclass is not sliced away to its C++ base.

    The subclass, its ``__dict__`` and its callbacks live in the Python wrapper, so the
    lattice has to keep that wrapper -- sharing the C++ element is not enough.
    """

    class Tagged(elements.Programmable):
        def __init__(self, tag):
            super().__init__(ds=0.0)
            self.tag = tag

    lattice = elements.KnownElementsList()
    lattice.append(Tagged("kept"))
    gc.collect()

    assert type(lattice[0]).__name__ == "Tagged"
    assert lattice[0].tag == "kept"


def test_constructor_shares_too():
    """Building a lattice from a list keeps those objects, like ``extend`` does."""

    q = elements.Quad(ds=0.3, k=2.0)
    d = elements.Drift(ds=1.0)

    lattice = elements.KnownElementsList([q, d])

    assert lattice[0] is q
    assert lattice[1] is d


def test_setitem_replaces_with_the_object_given():
    q = elements.Quad(ds=0.3, k=2.0)
    d = elements.Drift(ds=1.0)

    lattice = elements.KnownElementsList([q])
    lattice[0] = d

    assert lattice[0] is d


def test_negative_indexing():
    q = elements.Quad(ds=0.3, k=2.0)
    d = elements.Drift(ds=1.0)
    lattice = elements.KnownElementsList([q, d])

    assert lattice[-1] is d
    assert lattice[-2] is q

    with pytest.raises(IndexError):
        _ = lattice[2]
    with pytest.raises(IndexError):
        _ = lattice[-3]


def test_iteration_yields_the_same_objects():
    q = elements.Quad(ds=0.3, k=2.0)
    d = elements.Drift(ds=1.0)
    lattice = elements.KnownElementsList([q, d, q])

    seen = list(lattice)
    assert seen[0] is q
    assert seen[1] is d
    assert seen[2] is q


def test_pop_back_returns_the_object():
    q = elements.Quad(ds=0.3, k=2.0)
    d = elements.Drift(ds=1.0)
    lattice = elements.KnownElementsList([q, d])

    assert lattice.pop_back() is d
    assert len(lattice) == 1
    assert lattice[0] is q


def test_appending_a_non_element_is_rejected():
    lattice = elements.KnownElementsList()

    with pytest.raises(TypeError):
        lattice.append("not an element")


def test_extend_is_all_or_nothing():
    """A bad entry anywhere leaves the lattice untouched."""

    q = elements.Quad(ds=0.3, k=2.0)
    lattice = elements.KnownElementsList([q])

    with pytest.raises(TypeError):
        lattice.extend([elements.Drift(ds=1.0), "not an element"])

    assert len(lattice) == 1
    assert lattice[0] is q


def test_two_lattices_can_share_one_element():
    q = elements.Quad(ds=0.3, k=2.0)

    a = elements.KnownElementsList([q])
    b = elements.KnownElementsList([q])

    assert a[0] is b[0] is q

    del q
    gc.collect()

    # still one element, still shared
    assert a[0] is b[0]
    a[0].k = 5.0
    assert b[0].k == 5.0


def test_del_removes_one_occurrence():
    """Removing a position does not remove the element from other positions."""

    q = elements.Quad(ds=0.3, k=2.0)
    d = elements.Drift(ds=1.0)
    lattice = elements.KnownElementsList([q, d, q])

    del lattice[0]

    assert len(lattice) == 2
    assert lattice[0] is d
    assert lattice[1] is q


def test_insert():
    q = elements.Quad(ds=0.3, k=2.0)
    d = elements.Drift(ds=1.0)
    lattice = elements.KnownElementsList([q])

    lattice.insert(0, d)
    assert lattice[0] is d
    assert lattice[1] is q

    # out-of-range clamps, as for a Python list
    m = elements.Marker("m")
    lattice.insert(99, m)
    assert lattice[-1] is m


def test_membership_and_counting_are_identity_based():
    """Two elements with equal parameters are still two different elements."""

    q = elements.Quad(ds=0.3, k=2.0)
    twin = elements.Quad(ds=0.3, k=2.0)  # equal by value, a different element
    lattice = elements.KnownElementsList([q, q])

    assert q in lattice
    assert twin not in lattice

    assert lattice.count(q) == 2
    assert lattice.count(twin) == 0
    assert lattice.index(q) == 0

    with pytest.raises(ValueError):
        lattice.index(twin)


def test_remove_takes_the_first_occurrence():
    q = elements.Quad(ds=0.3, k=2.0)
    d = elements.Drift(ds=1.0)
    lattice = elements.KnownElementsList([q, d, q])

    lattice.remove(q)

    assert len(lattice) == 2
    assert lattice[0] is d
    assert lattice[1] is q

    with pytest.raises(ValueError):
        lattice.remove(elements.Quad(ds=0.3, k=2.0))


def test_reversed_iteration():
    q = elements.Quad(ds=0.3, k=2.0)
    d = elements.Drift(ds=1.0)
    lattice = elements.KnownElementsList([q, d])

    back = list(reversed(lattice))
    assert back[0] is d
    assert back[1] is q


def _numbered(n):
    """A lattice of n quads with k = 0..n-1, plus the elements themselves."""
    els = [elements.Quad(ds=0.1 * (i + 1), k=float(i)) for i in range(n)]
    return elements.KnownElementsList(els), els


def test_slice_read_shares_elements():
    lattice, els = _numbered(5)

    part = lattice[1:4]
    assert [e.k for e in part] == [1.0, 2.0, 3.0]
    # a slice is a new lattice over the same elements, as for a Python list
    assert part[0] is els[1]


def test_slice_read_with_step_and_negative_bounds():
    lattice, _ = _numbered(5)

    assert [e.k for e in lattice[::2]] == [0.0, 2.0, 4.0]
    assert [e.k for e in lattice[-2:]] == [3.0, 4.0]
    assert len(lattice[10:]) == 0


def test_slice_delete():
    lattice, _ = _numbered(5)

    del lattice[1:3]
    assert [e.k for e in lattice] == [0.0, 3.0, 4.0]


def test_contiguous_slice_assignment_may_resize():
    lattice, _ = _numbered(5)
    d = elements.Drift(ds=1.0, name="new")

    lattice[1:3] = [d]

    assert len(lattice) == 4
    assert lattice[1] is d
    assert [type(e).__name__ for e in lattice] == ["Quad", "Drift", "Quad", "Quad"]


def test_extended_slice_assignment_requires_matching_length():
    lattice, _ = _numbered(4)

    lattice[::2] = [elements.Drift(ds=1.0), elements.Drift(ds=2.0)]
    assert [type(e).__name__ for e in lattice] == ["Drift", "Quad", "Drift", "Quad"]

    with pytest.raises(ValueError):
        lattice[::2] = [elements.Drift(ds=1.0)]


def test_slice_assignment_is_all_or_nothing():
    lattice, els = _numbered(3)

    with pytest.raises(TypeError):
        lattice[0:2] = [elements.Drift(ds=1.0), "not an element"]

    assert len(lattice) == 3
    assert lattice[0] is els[0]


def test_slicing_a_subclass_gives_a_plain_lattice():
    """As slicing a ``list`` subclass gives a plain ``list``.

    A subclass may take constructor arguments there is nothing to pass, and may mean
    something of its own by ``extend``; neither is the slice's business.
    """

    class Beamline(elements.KnownElementsList):
        def __init__(self, label):
            super().__init__()
            self.label = label

        def extend(self, added):
            raise AssertionError("a slice must not go through the subclass")

    beamline = Beamline("ring")
    elements.KnownElementsList.extend(
        beamline, [elements.Drift(ds=0.1, name=f"d{i}") for i in range(4)]
    )

    sliced = beamline[0:2]

    assert type(sliced) is elements.KnownElementsList
    assert [element.name for element in sliced] == ["d0", "d1"]
    assert sliced[0] is beamline[0]
    assert beamline.label == "ring"
