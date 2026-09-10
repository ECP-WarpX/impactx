#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
"""Edges of the lattice API: empty lattices, odd slices, and views that outlive things."""

import gc

import pytest

from impactx import ImpactX, elements


def drifts(*names):
    lattice = elements.KnownElementsList()
    for i, name in enumerate(names):
        lattice.append(elements.Drift(ds=0.1 * (i + 1), name=name))
    return lattice


def names_of(lattice):
    return [element.name for element in lattice]


@pytest.mark.parametrize(
    ("count", "key", "expected"),
    [
        (4, slice(None, None, -1), []),
        (5, slice(None, None, -2), ["d1", "d3"]),
        (5, slice(None, None, 2), ["d1", "d3"]),
        (4, slice(None, None, -2), ["d0", "d2"]),
        (6, slice(4, 1, -1), ["d0", "d1", "d5"]),
    ],
)
def test_delete_slice_matches_a_list(count, key, expected):
    """Deleting through a slice removes what the same slice removes from a list."""

    lattice = drifts(*[f"d{i}" for i in range(count)])
    reference = [f"d{i}" for i in range(count)]

    del lattice[key]
    del reference[key]

    assert names_of(lattice) == expected
    assert names_of(lattice) == reference


def test_delete_slice_keeps_the_owner_list_in_step():
    """A slice delete leaves as many wrappers as elements, so indexing stays correct."""

    lattice = drifts("d0", "d1", "d2", "d3")
    del lattice[::-1]

    assert len(lattice) == 0
    lattice.append(elements.Drift(ds=0.9, name="fresh"))
    assert names_of(lattice) == ["fresh"]


@pytest.mark.parametrize("method", ["erase", "index"])
def test_out_of_range_position_raises(method):
    """A position past the end is reported, not read."""

    lattice = drifts("d0")
    with pytest.raises(IndexError):
        if method == "erase":
            del lattice[5]
        else:
            _ = lattice[5]


def test_pop_on_an_empty_lattice_raises():
    lattice = elements.KnownElementsList()
    with pytest.raises(IndexError):
        lattice.pop_back()


def test_tracking_an_empty_lattice_raises():
    """An empty lattice is reported rather than walked off the front."""

    from impactx import distribution

    sim = ImpactX()
    sim.particle_shape = 2
    sim.slice_step_diagnostics = False
    sim.diagnostics = False
    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(2.0e3)
    distr = distribution.Waterbag(
        lambdaX=4.0e-5,
        lambdaY=4.0e-5,
        lambdaT=1.0e-3,
        lambdaPx=2.7e-5,
        lambdaPy=2.7e-5,
        lambdaPt=2.0e-3,
    )
    sim.init_envelope(ref, distr)

    assert len(sim.lattice) == 0
    with pytest.raises(RuntimeError, match="zero elements"):
        sim.track_envelope()

    sim.finalize()


class TestViewOutlivingItsSimulation:
    """A lattice reached through ``sim.lattice`` says so rather than reading freed memory."""

    @staticmethod
    def orphaned_view():
        sim = ImpactX()
        view = sim.lattice
        view.append(elements.Drift(ds=0.2, name="d"))
        del sim
        gc.collect()
        return view

    @pytest.mark.parametrize(
        "use",
        [
            pytest.param(len, id="len"),
            pytest.param(lambda v: v.size(), id="size"),
            pytest.param(lambda v: v.is_empty(), id="is_empty"),
            pytest.param(lambda v: v[0], id="getitem"),
            pytest.param(list, id="iterate"),
            pytest.param(lambda v: v.append(elements.Drift(ds=0.1)), id="append"),
        ],
    )
    def test_reports_instead_of_reading_freed_memory(self, use):
        view = self.orphaned_view()

        with pytest.raises(RuntimeError, match="no longer exists"):
            use(view)


def test_a_standalone_lattice_is_not_affected():
    """Only a view of a simulation's lattice has a parent to outlive."""

    lattice = drifts("d0", "d1")
    assert len(lattice) == 2
    assert names_of(lattice) == ["d0", "d1"]


def test_slice_assignment_from_a_generator_that_appends():
    """A replacement iterable may change the lattice while it is being consumed.

    Regression test: the slice bounds were computed before the iterable was consumed, so
    a generator that appended left the assignment writing at stale positions. Starting
    from ["initial"], this produced ["replacement", "added"] where a list gives
    ["initial", "replacement"]. ``list.__setitem__`` materializes first; so do we.
    """

    lattice = drifts("initial")
    reference = ["initial"]

    def grow_lattice():
        lattice.append(elements.Drift(ds=1.0, name="added"))
        yield elements.Drift(ds=1.0, name="replacement")

    def grow_list():
        reference.append("added")
        yield "replacement"

    lattice[-1:] = grow_lattice()
    reference[-1:] = grow_list()

    assert names_of(lattice) == ["initial", "replacement"]
    assert names_of(lattice) == reference


def test_slice_assignment_from_a_generator_that_clears():
    """The stale bounds could also index past the end after the generator shortened it.

    Regression test: with the bounds taken first, clearing inside the generator left the
    head loop reading positions that were no longer there, raising IndexError.
    """

    lattice = drifts("a", "b", "c")
    reference = ["a", "b", "c"]

    def clear_lattice():
        lattice.clear()
        yield elements.Drift(ds=1.0, name="new")

    def clear_list():
        reference.clear()
        yield "new"

    lattice[-1:] = clear_lattice()
    reference[-1:] = clear_list()

    assert names_of(lattice) == ["new"]
    assert names_of(lattice) == reference


def test_extended_slice_survives_a_finalizer_that_reads_the_lattice():
    """A displaced element's ``__del__`` may run while the replacement is still going.

    Regression test: assigning through an extended slice released each displaced element
    as it went -- both its handle and, through the owner list, the last reference to its
    Python wrapper. A subclass finalizer that read the lattice from there saw a generation
    that had already moved and rebuilt the owner list, so the remaining writes landed in a
    list that was no longer the lattice's. The second replacement lost its wrapper and came
    back as a plain ``Drift`` without its attributes.
    """

    lattice = elements.KnownElementsList()
    reads = []

    class Observed(elements.Drift):
        def __del__(self):
            reads.append(lattice[1].ds)

    class Tagged(elements.Drift):
        def __init__(self, ds):
            super().__init__(ds=ds)
            self.tag = "retained"

    lattice.extend([Observed(ds=1.0), elements.Drift(ds=2.0), elements.Drift(ds=3.0)])
    lattice[::2] = [Tagged(10.0), Tagged(30.0)]

    assert type(lattice[0]) is Tagged
    assert type(lattice[2]) is Tagged
    assert lattice[0].tag == "retained"
    assert lattice[2].tag == "retained"
    assert reads == [2.0]
