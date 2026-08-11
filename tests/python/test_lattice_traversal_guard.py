#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""What a tracking hook may do to the lattice.

Retuning an element is supported and is the point of the hooks. Changing which elements
the lattice contains, while tracking is walking it, is not.
"""

import pytest

from impactx import ImpactX, elements


@pytest.fixture()
def sim():
    s = ImpactX()
    s.particle_shape = 2
    s.n_cell = [8, 8, 8]
    s.slice_step_diagnostics = False
    s.init_grids()
    s.beam.ref.set_species("electron").set_kin_energy_MeV(100.0)
    yield s
    s.finalize()


def _track(sim, hook):
    sim.hook["before_element"] = hook
    sim.track_reference(sim.beam.ref)


def test_retuning_an_element_is_allowed(sim):
    sim.lattice.append(elements.Drift(ds=0.5, name="d"))

    seen = {}

    def hook(s):
        s.tracking_element.ds = 0.25
        seen["ds"] = s.tracking_element.ds

    _track(sim, hook)
    assert seen["ds"] == 0.25


@pytest.mark.parametrize(
    "edit",
    [
        pytest.param(lambda lat: lat.append(elements.Drift(ds=0.1)), id="append"),
        pytest.param(lambda lat: lat.extend([elements.Drift(ds=0.1)]), id="extend"),
        pytest.param(lambda lat: lat.insert(0, elements.Drift(ds=0.1)), id="insert"),
        pytest.param(lambda lat: lat.clear(), id="clear"),
        pytest.param(lambda lat: lat.pop_back(), id="pop_back"),
        pytest.param(lambda lat: lat.__delitem__(0), id="del"),
        pytest.param(
            lambda lat: lat.__setitem__(0, elements.Drift(ds=0.1)), id="setitem"
        ),
        pytest.param(
            lambda lat: lat.__setitem__(slice(0, 1), [elements.Drift(ds=0.1)]),
            id="setitem_slice",
        ),
    ],
)
def test_changing_the_sequence_is_rejected(sim, edit):
    sim.lattice.append(elements.Drift(ds=0.5, name="d"))

    seen = {}

    def hook(s):
        try:
            edit(s.lattice)
            seen["result"] = "allowed"
        except RuntimeError as e:
            seen["result"] = str(e)

    _track(sim, hook)

    assert "while tracking" in seen["result"]
    assert len(sim.lattice) == 1


def test_the_lattice_is_editable_again_after_tracking(sim):
    sim.lattice.append(elements.Drift(ds=0.5, name="d"))
    _track(sim, lambda s: None)

    sim.lattice.append(elements.Drift(ds=0.5, name="d2"))
    assert len(sim.lattice) == 2


def test_the_diagnostic_traversals_are_guarded(sim):
    """map_trace() and transfer_map() walk the lattice too, and a reference push can
    run a user callback that reaches it.

    Regression test: the diagnostic walk took no traversal guard, so appending from a
    ``Programmable.ref_particle`` callback reallocated the storage under the iterator
    that was walking it. That surfaced as "visit with an empty element handle" -- a read
    of freed memory that happened to hit a null check.
    """

    lattice = sim.lattice

    def grow(refpart):
        for _ in range(64):
            lattice.append(elements.Drift(ds=0.1))

    grower = elements.Programmable(ds=0.0)
    grower.ref_particle = grow
    lattice.extend([elements.Drift(ds=0.5), grower, elements.Drift(ds=0.5)])

    with pytest.raises(RuntimeError, match="while tracking through it"):
        lattice.transfer_map(sim.beam.ref, fallback_identity_map=True)


def test_a_diagnostic_traversal_leaves_the_lattice_editable(sim):
    """The guard is released when the walk ends, including when it ends in an error."""

    sim.lattice.extend(
        [elements.Drift(ds=0.5, name="d1"), elements.Drift(ds=0.5, name="d2")]
    )
    sim.lattice.transfer_map(sim.beam.ref, fallback_identity_map=True)

    sim.lattice.append(elements.Drift(ds=0.5, name="d3"))
    assert len(sim.lattice) == 3
