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
        iterator = iter(selection)
        assert next(iterator) is lattice[0]

        edit(lattice)

        with pytest.raises(RuntimeError, match="no longer valid"):
            len(selection)
        with pytest.raises(RuntimeError, match="no longer valid"):
            _ = selection[0]
        with pytest.raises(RuntimeError, match="no longer valid"):
            next(iterator)

    @pytest.mark.parametrize("bound", ["start", "stop", "step"])
    @pytest.mark.parametrize("edit", ["insert", "clear"])
    def test_slice_bounds_cannot_revalidate_a_stale_selection(self, bound, edit):
        lattice = self.lattice_of_six()
        selection = lattice.select(kind="Quad")

        class EditingBound:
            def __index__(self):
                if edit == "insert":
                    lattice.insert(0, elements.Drift(ds=0.1, name="new_head"))
                else:
                    lattice.clear()
                return 1

        bounds = {"start": None, "stop": None, "step": None}
        bounds[bound] = EditingBound()
        key = slice(bounds["start"], bounds["stop"], bounds["step"])
        with pytest.raises(RuntimeError, match="no longer valid"):
            selection[key]

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
        """Only moving elements invalidates; changing one in place does not.

        The lattice counts structural edits, and a retune is not one of them.
        """

        lattice = self.lattice_of_six()
        selection = lattice.select(kind="Quad")
        generation = lattice.generation
        iterator = iter(selection)
        assert next(iterator) is lattice[0]

        lattice[0].k = 3.0

        assert lattice.generation == generation
        assert len(selection) == 3
        assert selection[0].k == 3.0
        assert next(iterator) is lattice[1]
        assert list(iterator) == [lattice[2]]

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

    def test_a_template_that_copies_once_changes_nothing(self):
        """The failure can arrive after a good copy, not on the first one.

        Regression test: the replacements were checked only as they were installed, so a
        `copy()` that returned an element once and something else next left the earlier
        positions already replaced.
        """

        class CopiesOnce(elements.Drift):
            calls = 0

            def copy(self, **overrides):
                CopiesOnce.calls += 1
                return elements.Drift(ds=0.5) if CopiesOnce.calls == 1 else None

        # unnamed, so `keep_name` does not touch the replacement first and the bad
        # value reaches the point where the type is checked
        lattice = elements.KnownElementsList()
        lattice.extend([elements.Quad(ds=0.1, k=1.0) for _ in range(3)])

        with pytest.raises(TypeError, match="expected a lattice element"):
            lattice.select(kind="Quad").replace_each(CopiesOnce(ds=0.5))

        assert [type(element).__name__ for element in lattice] == ["Quad"] * 3


def test_matching_drifts_follow_the_element_kind_of_a_subclass():
    """`model="match"` picks the drift for the element's kind, subclass or not.

    Regression test: the tier came from the exact class name, so a subclass of
    `ExactSbend` matched nothing and was replaced with a plain `Drift`.
    """

    class MyBend(elements.ExactSbend):
        pass

    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.ExactSbend(ds=1.0, phi=30.0, B=0.0),
            MyBend(ds=1.0, phi=30.0, B=0.0),
        ]
    )

    lattice.select(kind="ExactSbend").replace_with_drifts(model="match")

    assert [type(element).__name__ for element in lattice] == ["ExactDrift"] * 2


def test_a_copy_that_edits_the_lattice_is_refused():
    """Building the replacements runs user code, which may move the positions.

    Regression test: the positions were taken before `copy()` ran and used afterwards
    without rechecking, so a template whose `copy()` inserted an element wrote the
    replacements over whatever had moved into those positions.
    """

    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Quad(ds=1.0, k=1.0, name="q0"),
            elements.Drift(ds=1.0, name="d0"),
            elements.Quad(ds=1.0, k=2.0, name="q1"),
        ]
    )

    class Meddler(elements.Drift):
        edited = False

        def copy(self, **overrides):
            if not Meddler.edited:
                Meddler.edited = True
                lattice.insert(0, elements.Drift(ds=0.1, name="inserted"))
            return elements.Drift(ds=0.5)

    selection = lattice.select(kind="Quad")

    with pytest.raises(RuntimeError, match="no longer valid"):
        selection.replace_each(Meddler(ds=0.5))

    # only what the callback itself did; nothing written at the stale positions
    assert [element.name for element in lattice] == ["inserted", "q0", "d0", "q1"]


def test_dropping_the_simulation_inside_a_callback_does_not_abort():
    """A callback may release the last reference to the simulation it is tracking in.

    Regression test: the lattice bindings checked that the simulation was still there but
    did not hold it, so a callback that dropped the last reference had the simulation
    destroyed mid-traversal. Its destructor finalizes, finalizing refuses while a
    traversal is in progress, and an exception out of a destructor ends the process --
    this aborted with SIGABRT.
    """

    import subprocess
    import sys

    program = """
from impactx import ImpactX, RefPart, elements

holder = [ImpactX()]
lattice = holder[0].lattice

ref = RefPart()
ref.set_species("electron").set_kin_energy_MeV(100)

dropper = elements.Programmable()
dropper.ref_particle = lambda refpart: holder.clear()
lattice.append(dropper)

lattice.transfer_map(ref, fallback_identity_map=True)
print("survived")
"""

    finished = subprocess.run(
        [sys.executable, "-c", program], capture_output=True, text=True, timeout=120
    )

    assert finished.returncode == 0, finished.stderr[-2000:]
    assert "survived" in finished.stdout


class TestFinalizersDuringFilteredEdits:
    """A displaced element's `__del__` runs while the edit is still going.

    Regression tests: the filtered edits released each displaced element as they went, so
    a subclass finalizer saw a half-finished lattice and could move the positions the rest
    of the edit was about to write at.
    """

    def test_delete_never_shows_an_empty_lattice(self):
        lattice = elements.KnownElementsList()
        seen = []

        class Watcher(elements.Drift):
            def __del__(self):
                seen.append([element.name for element in lattice])

        lattice.extend(
            [Watcher(ds=1.0, name="remove"), elements.Drift(ds=1.0, name="keep")]
        )
        lattice.select(name="remove").delete()

        assert [element.name for element in lattice] == ["keep"]
        assert seen == [["keep"]]

    def test_replace_each_is_not_derailed_by_a_finalizer(self):
        lattice = elements.KnownElementsList()

        class Watcher(elements.Drift):
            def __del__(self):
                lattice.insert(0, elements.Drift(ds=1.0, name="added"))

        lattice.extend(
            [
                Watcher(ds=1.0, name="a"),
                elements.Drift(ds=1.0, name="keep"),
                elements.Drift(ds=1.0, name="b"),
            ]
        )
        lattice.select(name=["a", "b"]).replace_each(elements.Quad(ds=1.0, k=1.0))

        assert [(e.name, type(e).__name__) for e in lattice] == [
            ("added", "Drift"),
            ("a", "Quad"),
            ("keep", "Drift"),
            ("b", "Quad"),
        ]

    def test_a_property_getter_that_edits_the_lattice_is_refused(self):
        """`replace_with_drifts` reads the elements it replaces to pick their drift.

        Regression test: on a Python subclass a property getter is user code and can edit
        the lattice, and the positions were taken before it ran. `replace_each` rechecked
        afterwards; this one did not, and wrote the replacements at stale positions.
        """

        lattice = elements.KnownElementsList()

        class MutatingQuad(elements.Quad):
            edited = False

            @property
            def ds(self):
                if not MutatingQuad.edited:
                    MutatingQuad.edited = True
                    lattice.insert(0, elements.Drift(ds=0.25, name="inserted"))
                return super().ds

        lattice.extend(
            [
                MutatingQuad(ds=1.0, k=1.0, name="a"),
                elements.Drift(ds=1.0, name="keep"),
                elements.Quad(ds=1.0, k=1.0, name="b"),
            ]
        )

        with pytest.raises(RuntimeError, match="no longer valid"):
            lattice.select(kind="Quad").replace_with_drifts()

        assert [element.name for element in lattice] == ["inserted", "a", "keep", "b"]

    def test_the_returned_selection_goes_stale_if_a_finalizer_edits(self):
        """The selection handed back is stamped before the displaced elements are let go.

        Regression test: it was built after the release, so a finalizer that edited the
        lattice left it carrying positions from before the edit and the generation from
        after -- it looked valid and named the wrong elements. Built beforehand, that edit
        invalidates it, as it invalidates any other selection.
        """

        lattice = elements.KnownElementsList()

        class Watcher(elements.Quad):
            def __del__(self):
                lattice.insert(0, elements.Quad(ds=1.0, k=2.0, name="added"))

        lattice.extend(
            [
                Watcher(ds=1.0, k=1.0, name="a"),
                elements.Drift(ds=1.0, name="keep"),
                elements.Quad(ds=1.0, k=1.0, name="b"),
            ]
        )

        replaced = lattice.select(kind="Quad").replace_with_drifts()

        with pytest.raises(RuntimeError, match="no longer valid"):
            list(replaced)

    def test_the_returned_selection_is_usable_when_nothing_edits(self):
        """The ordinary case still hands back a selection over the replacements."""

        lattice = elements.KnownElementsList()
        lattice.extend(
            [
                elements.Quad(ds=1.0, k=1.0, name="a"),
                elements.Drift(ds=1.0, name="keep"),
                elements.Quad(ds=1.0, k=1.0, name="b"),
            ]
        )

        replaced = lattice.select(kind="Quad").replace_with_drifts()

        assert [element.name for element in replaced] == ["a", "b"]

    def test_replace_with_drifts_is_not_derailed_by_a_finalizer(self):
        lattice = elements.KnownElementsList()

        class Watcher(elements.Quad):
            def __del__(self):
                lattice.insert(0, elements.Drift(ds=1.0, name="added"))

        lattice.extend(
            [
                Watcher(ds=1.0, k=1.0, name="a"),
                elements.Drift(ds=1.0, name="keep"),
                elements.Quad(ds=1.0, k=1.0, name="b"),
            ]
        )
        lattice.select(kind="Quad").replace_with_drifts()

        assert [(e.name, type(e).__name__) for e in lattice] == [
            ("added", "Drift"),
            ("a", "Drift"),
            ("keep", "Drift"),
            ("b", "Drift"),
        ]


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


@pytest.mark.parametrize("chained", [False, True], ids=["direct", "chained"])
@pytest.mark.parametrize("edit", ["insert", "clear"])
@pytest.mark.parametrize("position", [0, 1], ids=["first", "last"])
def test_selection_rejects_edits_during_matching(chained, edit, position):
    """Matching a subclass property must not stamp old positions as current."""
    lattice = elements.KnownElementsList()

    class MutatingQuad(elements.Quad):
        edited = False

        @property
        def name(self):
            if not self.edited:
                self.edited = True
                if edit == "insert":
                    lattice.insert(0, elements.Drift(ds=0.5, name="inserted"))
                else:
                    lattice.clear()
            return super().name

    row = [elements.Quad(ds=1.0, k=1.0, name="keep")]
    row.insert(position, MutatingQuad(ds=1.0, k=1.0, name="selected"))
    lattice.extend(row)
    source = lattice.select(kind="Quad") if chained else lattice

    with pytest.raises(RuntimeError, match="no longer valid"):
        source.select(name="selected").delete()

    # Only the getter's edit took effect; no stale selection deleted a different element.
    expected = (
        ["inserted"] + [element.name for element in row] if edit == "insert" else []
    )
    assert names_of(lattice) == expected


def test_owner_release_callback_does_not_keep_a_simulation_cycle_alive():
    """The C++ cleanup callback must not hide a strong Python ownership cycle."""
    import weakref

    sim = ImpactX()
    element = elements.Programmable()
    element.sim = sim
    sim.lattice.append(element)
    observed_sim = weakref.ref(sim)
    observed_element = weakref.ref(element)

    del sim, element
    gc.collect()

    assert observed_sim() is None
    assert observed_element() is None
