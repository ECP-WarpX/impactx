#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
"""A lattice view belongs to one simulation, and a selection to one arrangement of it."""

import gc

import pytest

from impactx import ImpactX, elements


def names_of(lattice):
    return [element.name for element in lattice]


def test_a_new_simulation_does_not_inherit_the_previous_lattice():
    """Two simulations never share a lattice, whatever the allocator does.

    A simulation that is gone frees an address that the next one may be given, and the
    Python wrapper is kept per address. The elements of the first simulation must not
    turn up in the second.
    """

    first = ImpactX()
    first.lattice.append(elements.Drift(ds=1.0, name="from_the_first"))
    del first
    gc.collect()

    for _ in range(8):
        second = ImpactX()

        assert len(second.lattice) == 0
        assert names_of(second.lattice) == []

        second.lattice.append(elements.Drift(ds=2.0, name="from_the_second"))
        assert names_of(second.lattice) == ["from_the_second"]

        second.finalize()


def test_a_view_of_a_gone_simulation_never_becomes_another_simulation():
    """A kept view must stop working, not start working on someone else's lattice.

    Returning ``sim.lattice`` from a helper is the ordinary way to hit this: the
    simulation is gone by the time the caller looks, and the next one may be allocated
    exactly where it was.
    """

    def first_run():
        simulation = ImpactX()
        simulation.lattice.append(elements.Drift(ds=1.0, name="first"))
        return simulation.lattice

    kept = first_run()
    gc.collect()

    second = ImpactX()
    second.lattice.append(elements.Drift(ds=2.0, name="second"))

    assert kept is not second.lattice
    with pytest.raises(RuntimeError, match="no longer exists"):
        len(kept)
    with pytest.raises(RuntimeError, match="no longer exists"):
        kept.append(elements.Drift(ds=99.0, name="injected"))

    assert names_of(second.lattice) == ["second"]

    second.finalize()


def test_each_simulation_gets_its_own_lattice_contents():
    """Two simulations alive at once keep their lattices apart."""

    first = ImpactX()
    second = ImpactX()

    first.lattice.append(elements.Drift(ds=1.0, name="first"))
    second.lattice.append(elements.Quad(ds=0.3, k=1.0, name="second"))

    assert names_of(first.lattice) == ["first"]
    assert names_of(second.lattice) == ["second"]

    first.finalize()
    second.finalize()


class TestSelectionsGoStaleWithTheLattice:
    """A selection is a list of positions, so any edit that moves elements voids it."""

    @staticmethod
    def lattice_of_six():
        lattice = elements.KnownElementsList()
        lattice.extend(
            [
                elements.Quad(ds=0.1, k=1.0, name="q0"),
                elements.Quad(ds=0.1, k=1.0, name="q1"),
                elements.Quad(ds=0.1, k=1.0, name="q2"),
                elements.Drift(ds=0.1, name="d0"),
                elements.Drift(ds=0.1, name="d1"),
                elements.Drift(ds=0.1, name="d2"),
            ]
        )
        return lattice

    @pytest.mark.parametrize(
        "edit",
        [
            pytest.param(lambda lat: lat.append(elements.Drift(ds=0.1)), id="append"),
            pytest.param(
                lambda lat: lat.insert(0, elements.Drift(ds=0.1)), id="insert"
            ),
            pytest.param(lambda lat: lat.__delitem__(0), id="del"),
            pytest.param(lambda lat: lat.__delitem__(slice(0, 2)), id="del_slice"),
            pytest.param(
                lambda lat: lat.__setitem__(slice(0, 2), [elements.Drift(ds=0.1)]),
                id="setitem_slice",
            ),
            pytest.param(lambda lat: lat.clear(), id="clear"),
            pytest.param(lambda lat: lat.pop_back(), id="pop_back"),
            pytest.param(lambda lat: lat.extend([elements.Drift(ds=0.1)]), id="extend"),
        ],
    )
    def test_reading_a_stale_selection_raises(self, edit):
        lattice = self.lattice_of_six()
        selection = lattice.select(kind="Quad")
        assert len(selection) == 3

        edit(lattice)

        with pytest.raises(RuntimeError, match="no longer valid"):
            len(selection)
        with pytest.raises(RuntimeError, match="no longer valid"):
            _ = selection[0]

    def test_writing_through_a_stale_selection_raises(self):
        """The dangerous case: the positions now name different elements."""

        lattice = self.lattice_of_six()
        selection = lattice.select(kind="Quad")

        lattice.insert(0, elements.Drift(ds=0.1, name="new_head"))

        with pytest.raises(RuntimeError, match="no longer valid"):
            selection.replace_with_drifts()

        # nothing was rewritten
        assert names_of(lattice) == ["new_head", "q0", "q1", "q2", "d0", "d1", "d2"]

    def test_retuning_an_element_does_not_void_a_selection(self):
        """Only moving elements invalidates; changing one in place does not."""

        lattice = self.lattice_of_six()
        selection = lattice.select(kind="Quad")

        lattice[0].k = 3.0

        assert len(selection) == 3
        assert selection[0].k == 3.0

    def test_a_fresh_selection_after_an_edit_is_usable(self):
        lattice = self.lattice_of_six()
        lattice.append(elements.Quad(ds=0.1, k=1.0, name="q3"))

        selection = lattice.select(kind="Quad")

        assert len(selection) == 4
        selection.replace_with_drifts()
        assert [type(element).__name__ for element in lattice] == ["Drift"] * 7


@pytest.mark.parametrize(
    "no_op",
    [
        pytest.param(lambda lat: lat.__delitem__(slice(3, 3)), id="del_empty_slice"),
        pytest.param(lambda lat: lat.__setitem__(slice(2, 2), []), id="assign_nothing"),
        pytest.param(lambda lat: lat.extend([]), id="extend_nothing"),
    ],
)
def test_an_edit_that_changes_nothing_keeps_selections_usable(no_op):
    """Only an edit that moves elements makes a selection describe something else."""

    lattice = elements.KnownElementsList()
    lattice.extend(
        [elements.Quad(ds=0.1, k=1.0, name=f"q{i}") for i in range(3)]
        + [elements.Drift(ds=0.1, name=f"d{i}") for i in range(3)]
    )
    selection = lattice.select(kind="Quad")
    before = lattice.generation

    no_op(lattice)

    assert lattice.generation == before
    assert len(selection) == 3


def test_an_empty_selection_leaves_other_selections_alone():
    lattice = elements.KnownElementsList()
    lattice.extend([elements.Quad(ds=0.1, k=1.0, name="q")])

    quads = lattice.select(kind="Quad")
    nothing = lattice.select(kind="Sbend")
    assert len(nothing) == 0

    nothing.delete()

    assert len(quads) == 1


def test_generation_counts_structural_edits_only():
    lattice = elements.KnownElementsList()
    start = lattice.generation

    lattice.append(elements.Drift(ds=0.1, name="d"))
    after_append = lattice.generation
    assert after_append != start

    lattice[0].ds = 0.5
    assert lattice.generation == after_append


class TestFilteredEditsAreAllOrNothing:
    """A rejected edit leaves the lattice as it was, as `Lattice` promises."""

    @staticmethod
    def three_quads():
        lattice = elements.KnownElementsList()
        lattice.extend([elements.Quad(ds=0.1, k=1.0, name=f"q{i}") for i in range(3)])
        return lattice

    def test_a_template_that_cannot_be_named_changes_nothing(self):
        """A `BeamMonitor` has no settable name, so the second position fails."""

        lattice = self.three_quads()

        with pytest.raises(AttributeError):
            lattice.select(kind="Quad").replace_each(elements.BeamMonitor("mon"))

        assert [type(element).__name__ for element in lattice] == ["Quad"] * 3

    def test_a_template_that_cannot_be_copied_changes_nothing(self):
        class MyDrift(elements.Drift):
            pass

        lattice = self.three_quads()

        with pytest.raises(TypeError, match="copy"):
            lattice.select(kind="Quad").replace_each(MyDrift(ds=0.5))

        assert [type(element).__name__ for element in lattice] == ["Quad"] * 3


class TestInsertEveryDsKeepsElements:
    @staticmethod
    def tagged_quad():
        class MyQuad(elements.Quad):
            def __init__(self, **kwargs):
                super().__init__(**kwargs)
                self.tag = "mine"

        return MyQuad(ds=0.5, k=1.0)

    def test_unsplit_elements_survive_the_source_lattice(self):
        """The result holds the elements, not wrappers minted from them later."""

        def build():
            source = elements.KnownElementsList(
                [elements.Drift(ds=1.0), self.tagged_quad()]
            )
            return elements.transformation.insert_element_every_ds(
                source, 1.0, elements.Marker("m")
            )

        result = build()
        gc.collect()

        assert type(result[2]).__name__ == "MyQuad"
        assert result[2].tag == "mine"

    def test_an_element_at_several_positions_is_matched_to_all_of_them(self):
        shared = elements.Quad(ds=0.2, k=1.0, name="shared")
        source = elements.KnownElementsList([shared, elements.Drift(ds=1.0), shared])

        result = elements.transformation.insert_element_every_ds(
            source, 1.0, elements.Marker("m")
        )

        held = [element for element in result if element is shared]
        assert len(held) == 2
