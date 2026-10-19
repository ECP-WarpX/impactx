#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""``element.copy()`` gives a distinct element with the same configuration.

``lattice.append(q)`` adds another occurrence of ``q``; ``lattice.append(q.copy())``
adds a second element. Copying every element type is covered by
``test_element_serialization.py``.
"""

import pytest

from impactx import Config, ImpactX, distribution, elements


def test_copy_is_a_distinct_element():
    q = elements.Quad(ds=0.3, k=2.0, nslice=3, name="q1")
    c = q.copy()

    assert c is not q
    assert (c.ds, c.k, c.nslice, c.name) == (q.ds, q.k, q.nslice, q.name)
    assert (c.ds, c.k, c.nslice, c.name) == (pytest.approx(0.3), 2.0, 3, "q1")

    c.k = 9.0
    assert q.k == 2.0


def test_copy_owns_its_arrays():
    sq = elements.SoftQuadrupole(
        ds=1.0, gscale=1.0, cos_coefficients=[1.0, 2.0], sin_coefficients=[0.0, 3.0]
    )
    c = sq.copy()

    c.set_coefficients([9.0, 9.0], [8.0, 8.0])
    assert sq.cos_coefficients == [1.0, 2.0]
    assert sq.sin_coefficients == [0.0, 3.0]

    sq.cos_coefficients = [5.0, 6.0]
    assert c.cos_coefficients == [9.0, 9.0]


def test_overrides_apply_to_the_copy_only():
    template = elements.Quad(ds=1.0, k=1.0, name="q")

    derived = template.copy(k=2.0, name="q2")
    assert (derived.k, derived.name) == (2.0, "q2")
    assert (template.k, template.name) == (1.0, "q")

    scan = [template.copy(k=k) for k in (0.8, 0.9, 1.0)]
    assert [element.k for element in scan] == pytest.approx([0.8, 0.9, 1.0])

    # paired arrays are given together
    rf = elements.RFCavity(
        ds=1.0,
        escale=1.0,
        freq=1.0e9,
        phase=0.0,
        cos_coefficients=[2.0],
        sin_coefficients=[0.0],
    )
    longer = rf.copy(cos_coefficients=[2.0, 0.1], sin_coefficients=[0.0, 0.2])
    assert longer.cos_coefficients == pytest.approx([2.0, 0.1])
    assert rf.cos_coefficients == [2.0]


@pytest.mark.parametrize(
    "template",
    [
        pytest.param(lambda: elements.Quad(ds=1.0, k=1.0), id="Quad"),
        # accepts arbitrary attributes, so a typo must not become a new attribute
        pytest.param(
            lambda: elements.Programmable(ds=1.0, nslice=3), id="Programmable"
        ),
    ],
)
def test_an_unknown_override_raises(template):
    with pytest.raises(AttributeError, match="nslcie"):
        template().copy(nslcie=5)


def test_a_beam_monitor_copy_cannot_change_shared_settings():
    """A monitor's Twiss settings are keyed by its name, which the copy shares."""

    monitor = elements.BeamMonitor("mon_shared")
    monitor.beta = 2.0

    with pytest.raises(ValueError, match="shares with its copy"):
        monitor.copy(beta=5.0)
    assert monitor.beta == 2.0


def test_a_python_subclass_defines_its_own_copy():
    class Tagged(elements.Drift):
        def __init__(self, tag):
            super().__init__(ds=1.0)
            self.tag = tag

    with pytest.raises(TypeError, match="copy"):
        Tagged("x").copy()

    class Copyable(Tagged):
        def copy(self):
            return Copyable(self.tag)

    original = Copyable("kept")
    c = original.copy()
    assert c is not original
    assert c.tag == "kept"


def test_names():
    drift = elements.Drift(ds=1.0)
    assert not drift.has_name
    assert drift.name is None
    assert not drift.copy().has_name

    assert not elements.Drift(ds=1.0, name="").has_name

    drift.name = "after"
    assert drift.name == "after"
    drift.name = ""
    assert not drift.has_name

    unicode_name = "\u00fcn\u00efcod\u00e9-" * 8
    assert elements.Drift(ds=1.0, name=unicode_name).copy().name == unicode_name


def _monitor_positions(name, tail, tmp_path, monkeypatch):
    """Reference positions written by a monitor at the head, and ``tail(head)`` at the
    end, of ``[head, drift, tail]``."""

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
    positions = []
    for iteration in series.read_iterations():
        positions.append(iteration.particles["beam"].get_attribute("s_ref"))
    series.close()
    return positions


@pytest.mark.skipif(not Config.have_openpmd, reason="built without openPMD")
def test_a_monitor_at_head_and_tail_records_both_passes(tmp_path, monkeypatch):
    """The same monitor twice, or a monitor and its copy, write one series."""

    aliased = _monitor_positions("aliased", lambda head: head, tmp_path, monkeypatch)
    copied = _monitor_positions(
        "copied", lambda head: head.copy(), tmp_path, monkeypatch
    )

    assert aliased == pytest.approx([0.0, 0.5])
    assert copied == aliased
