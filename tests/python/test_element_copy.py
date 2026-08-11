#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""``element.copy()`` gives a distinct element with the same configuration.

This is how a user asks for an independent element when they already have one:
``lattice.append(q)`` adds another occurrence of ``q``, ``lattice.append(q.copy())``
adds a second element.
"""

import pytest

from impactx import elements


def test_copy_is_a_distinct_element():
    q = elements.Quad(ds=0.3, k=2.0, name="q1")
    c = q.copy()

    assert c is not q
    assert c.ds == q.ds and c.k == q.k and c.name == q.name

    c.k = 9.0
    assert q.k == 2.0


def test_copy_in_a_lattice_is_a_second_element():
    q = elements.Quad(ds=0.3, k=2.0)
    lattice = elements.KnownElementsList([q, q.copy()])

    assert lattice[0] is q
    assert lattice[1] is not q

    lattice[1].k = 5.0
    assert q.k == 2.0


def test_every_element_type_can_be_copied():
    without = [
        name
        for name in dir(elements)
        if name[0].isupper()
        and name not in ("KnownElementsList", "FilteredElementsList")
        and not hasattr(getattr(elements, name), "copy")
    ]
    assert without == []


def test_dynamic_element_copy_owns_its_arrays():
    """The coefficient arrays must not be shared between two elements."""

    sq = elements.SoftQuadrupole(
        ds=1.0, gscale=1.0, cos_coefficients=[1.0, 2.0], sin_coefficients=[0.0, 3.0]
    )
    c = sq.copy()

    assert c is not sq
    assert c.to_dict()["cos_coefficients"] == [1.0, 2.0]
    assert c.to_dict()["sin_coefficients"] == [0.0, 3.0]


def test_beam_monitor_copy_does_not_share_the_output_series():
    """Two monitors sharing one open series would interleave into the same iterations."""

    monitor = elements.BeamMonitor("mon", backend="h5")
    c = monitor.copy()

    assert c is not monitor
    assert c.name == monitor.name

    monitor.finalize()
    c.finalize()


def test_python_subclass_must_say_what_a_copy_means():
    """Refusing beats silently returning a plain base element."""

    class Tagged(elements.Programmable):
        def __init__(self):
            super().__init__(ds=0.0)
            self.tag = "x"

    with pytest.raises(TypeError, match="copy"):
        Tagged().copy()


def test_subclass_can_define_its_own_copy():
    class Tagged(elements.Programmable):
        def __init__(self, tag):
            super().__init__(ds=0.0)
            self.tag = tag

        def copy(self):
            return Tagged(self.tag)

    original = Tagged("kept")
    c = original.copy()

    assert c is not original
    assert c.tag == "kept"


def test_filtered_delete_leaves_other_elements_untouched():
    """Unselected elements keep their identity, subclass and attributes."""

    class Tagged(elements.Programmable):
        def __init__(self, tag):
            super().__init__(ds=0.0)
            self.tag = tag

    keep = Tagged("survivor")
    lattice = elements.KnownElementsList(
        [keep, elements.Quad(ds=0.3, k=1.0, name="drop"), elements.Drift(ds=1.0)]
    )

    lattice.select(name="drop").delete()

    assert lattice[0] is keep
    assert type(lattice[0]).__name__ == "Tagged"
    assert lattice[0].tag == "survivor"


def test_replace_each_uses_one_copy_per_position():
    lattice = elements.KnownElementsList(
        [
            elements.Quad(ds=0.3, k=1.0, name="a"),
            elements.Quad(ds=0.3, k=2.0, name="b"),
        ]
    )

    lattice.select(kind="Quad").replace_each(elements.Drift(ds=1.0))

    assert [type(e).__name__ for e in lattice] == ["Drift", "Drift"]
    assert lattice[0] is not lattice[1]
    assert [e.name for e in lattice] == ["a", "b"]


class TestCopyWithOverrides:
    """``copy()`` takes the differences that make the copy a different element."""

    def test_an_override_applies_to_the_copy_only(self):
        template = elements.Quad(ds=1.0, k=1.0, name="q")

        derived = template.copy(k=2.0, name="q2")

        assert derived.k == 2.0
        assert derived.name == "q2"
        assert template.k == 1.0
        assert template.name == "q"

    def test_a_template_serves_a_scan(self):
        template = elements.Quad(ds=1.0, k=1.0, name="q")

        scan = [template.copy(k=k) for k in (0.8, 0.9, 1.0)]

        assert [element.k for element in scan] == [0.8, 0.9, 1.0]
        assert len({id(element) for element in scan}) == 3
        assert template.k == 1.0

    def test_copy_without_overrides_is_unchanged(self):
        template = elements.Quad(ds=1.0, k=1.0, name="q")

        assert template.copy().k == 1.0

    def test_a_parameter_the_element_does_not_have_is_reported(self):
        """A mistyped name must not be silently ignored."""

        template = elements.Quad(ds=1.0, k=1.0)

        with pytest.raises(AttributeError):
            template.copy(kk=2.0)

    def test_an_element_carrying_arrays_keeps_them(self):
        element = elements.SoftQuadrupole(
            ds=0.1, gscale=1.0, cos_coefficients=[2.0], sin_coefficients=[0.0]
        )

        derived = element.copy(gscale=3.0)

        assert derived.gscale == 3.0
        assert list(derived.cos_coefficients) == list(element.cos_coefficients)

    def test_a_subclass_is_still_refused(self):
        class MyQuad(elements.Quad):
            pass

        with pytest.raises(TypeError, match="copy"):
            MyQuad(ds=1.0, k=1.0).copy(k=2.0)
