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

from impactx import Config, ImpactX, distribution, elements


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
    """The coefficient arrays must not be shared between two elements.

    Writing through either side must leave the other as it was; equal contents alone
    would also hold for one shared array.
    """

    sq = elements.SoftQuadrupole(
        ds=1.0, gscale=1.0, cos_coefficients=[1.0, 2.0], sin_coefficients=[0.0, 3.0]
    )
    c = sq.copy()
    assert c is not sq

    c.set_coefficients([9.0, 9.0], [8.0, 8.0])
    assert sq.cos_coefficients == [1.0, 2.0]
    assert sq.sin_coefficients == [0.0, 3.0]

    sq.cos_coefficients = [5.0, 6.0]
    assert c.cos_coefficients == [9.0, 9.0]


def _track_monitors(name, tail, tmp_path, monkeypatch):
    """Track a few particles through ``[head, drift, tail]`` and return the
    reference-particle positions written to the series ``name``.

    ``tail`` is built from the head monitor: the monitor itself, or a copy of it.
    """

    io = pytest.importorskip("openpmd_api")
    from pathlib import Path

    monkeypatch.chdir(tmp_path)

    sim = ImpactX()
    sim.particle_shape = 2
    sim.slice_step_diagnostics = False
    sim.init_grids()
    sim.beam.ref.set_species("electron").set_kin_energy_MeV(2.0e3)
    sim.add_particles(
        1.0e-9,
        distribution.Waterbag(
            lambdaX=4.0e-5,
            lambdaY=4.0e-5,
            lambdaT=1.0e-3,
            lambdaPx=2.7e-5,
            lambdaPy=2.7e-5,
            lambdaPt=2.0e-3,
        ),
        16,
    )

    head = elements.BeamMonitor(name, backend="h5")
    sim.lattice.extend([head, elements.Drift(ds=0.5), tail(head)])
    try:
        sim.track_particles()
    finally:
        sim.finalize()

    (path,) = sorted(Path("diags/openPMD").glob(f"{name}.*"))
    series = io.Series(str(path), io.Access.read_linear)
    return [
        iteration.particles["beam"].get_attribute("s_ref")
        for iteration in series.read_iterations()
    ]


@pytest.mark.skipif(not Config.have_openpmd, reason="built without openPMD")
def test_beam_monitor_copy_is_a_second_monitor(tmp_path, monkeypatch):
    """A copy of a monitor records the beam like the monitor itself would.

    Monitors of one name write one series, so a monitor at the head and the tail of a
    lattice -- the same object, or a copy of it -- records both passes. The copy starts
    without the open series and the per-pass state of the original, and opens its own on
    first use; reusing them would have the tail write into the head's pass.
    """

    aliased = _track_monitors("aliased", lambda head: head, tmp_path, monkeypatch)
    copied = _track_monitors("copied", lambda head: head.copy(), tmp_path, monkeypatch)

    assert aliased == pytest.approx([0.0, 0.5])
    assert copied == aliased


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

    def test_a_mistyped_name_is_reported_on_a_dynamic_attribute_element(self):
        """`Programmable` takes attributes of its own, which hid mistyped names.

        Regression test: the overrides were applied with `setattr`, and an element that
        accepts dynamic attributes took a mistyped name as a new attribute instead of
        rejecting it -- leaving the parameter it was meant for unchanged, where every
        other element type raises.
        """

        template = elements.Programmable(ds=1.0, nslice=3)

        with pytest.raises(AttributeError, match="nslcie"):
            template.copy(nslcie=5)

        # the name it was meant for still works
        assert template.copy(nslice=5).nslice == 5

    def test_a_setting_the_copy_cannot_own_is_refused(self):
        """A BeamMonitor's Twiss settings are keyed by its name, which a copy shares.

        Regression test: `monitor.copy(beta=5.0)` applied the override to state the
        original reads too, so it silently retuned the original's diagnostics. Its name
        has no setter, so the copy cannot be given values of its own either.
        """

        monitor = elements.BeamMonitor("mon_shared")
        monitor.beta = 2.0

        with pytest.raises(ValueError, match="shares with its copy"):
            monitor.copy(beta=5.0)

        assert monitor.beta == 2.0

        # copying without those overrides is still fine
        assert type(monitor.copy()).__name__ == "BeamMonitor"

    def test_a_name_is_checked_without_running_its_getter(self):
        """Checking an override's name must not call the property's getter.

        Regression test: the name was looked up on the instance, which runs the getter,
        and a fresh `BeamMonitor`'s getters refuse until something is configured -- so a
        perfectly good name was reported as no such attribute. Asking the type finds the
        property without calling it, which is why a fresh monitor reaches the rule below
        rather than an `AttributeError`.
        """

        fresh = elements.BeamMonitor("mon_fresh")

        with pytest.raises(ValueError, match="shares with its copy"):
            fresh.copy(beta=5.0)

    def test_an_element_carrying_arrays_keeps_them(self):
        element = elements.SoftQuadrupole(
            ds=0.1, gscale=1.0, cos_coefficients=[2.0], sin_coefficients=[0.0]
        )

        derived = element.copy(gscale=3.0)

        assert derived.gscale == 3.0
        assert list(derived.cos_coefficients) == list(element.cos_coefficients)
