#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""``sim.lattice``: the lattice of one simulation, tracked, edited and selected from."""

import gc
import math
import weakref

import pytest

from impactx import ImpactX, distribution, elements


def names_of(lattice):
    return [element.name for element in lattice]


def waterbag():
    return distribution.Waterbag(
        lambdaX=1.0e-4,
        lambdaY=1.0e-4,
        lambdaT=1.0e-3,
        lambdaPx=1.0e-5,
        lambdaPy=1.0e-5,
        lambdaPt=1.0e-3,
    )


@pytest.fixture
def sim():
    simulation = ImpactX()
    simulation.particle_shape = 2
    simulation.n_cell = [8, 8, 8]
    simulation.slice_step_diagnostics = False
    simulation.diagnostics = False
    simulation.init_grids()
    simulation.beam.ref.set_species("electron").set_kin_energy_MeV(100.0)
    yield simulation
    simulation.finalize()


def track_reference(s):
    s.track_reference(s.beam.ref)


def track_envelope(s):
    s.init_envelope(s.beam.ref, waterbag())
    s.track_envelope()


def track_particles(s):
    s.add_particles(1.0e-9, waterbag(), 16)
    s.track_particles()


TRACKERS = [
    pytest.param(track_reference, id="reference"),
    pytest.param(track_envelope, id="envelope"),
    pytest.param(track_particles, id="particles"),
]


# --- one lattice per simulation ---------------------------------------------------------


def test_each_simulation_has_its_own_lattice():
    first = ImpactX()
    first.lattice.append(elements.Drift(ds=1.0, name="first"))
    second = ImpactX()
    second.lattice.append(elements.Drift(ds=2.0, name="second"))

    assert names_of(first.lattice) == ["first"]
    assert names_of(second.lattice) == ["second"]
    first.finalize()
    second.finalize()

    # a new simulation starts empty, even when allocated where an old one was
    for _ in range(4):
        fresh = ImpactX()
        assert len(fresh.lattice) == 0
        fresh.finalize()


def test_a_lattice_outliving_its_simulation_raises():
    def build():
        simulation = ImpactX()
        simulation.lattice.append(elements.Drift(ds=1.0, name="first"))
        return simulation.lattice

    kept = build()
    gc.collect()

    for use in (len, list, lambda v: v.append(elements.Drift(ds=1.0))):
        with pytest.raises(RuntimeError, match="no longer exists"):
            use(kept)


def test_an_element_referring_to_its_simulation_is_collected():
    """A reference cycle through the lattice must not keep the simulation alive."""

    simulation = ImpactX()
    element = elements.Programmable()
    element.sim = simulation
    simulation.lattice.append(element)
    observed_sim = weakref.ref(simulation)
    observed_element = weakref.ref(element)

    del simulation, element
    gc.collect()

    assert observed_sim() is None
    assert observed_element() is None


def test_a_hook_referring_to_its_simulation_is_collected():
    """finalize() releases the hooks, and with them a cycle through the simulation."""

    def run():
        # the closure keeps the cell of `simulation`; `del` would clear it and hide a cycle
        simulation = ImpactX()
        simulation.hook["before_element"] = lambda _: simulation.beam
        simulation.finalize()
        assert len(simulation.hook) == 0
        return weakref.ref(simulation)

    observed_sim = run()
    gc.collect()

    assert observed_sim() is None


def test_assigning_the_lattice(sim):
    sim.lattice = [elements.Drift(ds=0.1, name="d0"), elements.Drift(ds=0.2, name="d1")]
    assert names_of(sim.lattice) == ["d0", "d1"]

    sim.lattice = sim.lattice
    assert names_of(sim.lattice) == ["d0", "d1"]

    # a rejected assignment leaves the previous lattice in place
    with pytest.raises(TypeError):
        sim.lattice = [elements.Quad(ds=0.3, k=1.0), "not an element"]
    assert names_of(sim.lattice) == ["d0", "d1"]


def test_a_lattice_from_an_inputs_file_replaces_the_python_one(sim, tmp_path):
    """Parsing replaces the lattice in C++; the Python view must see the new elements."""

    inputs = tmp_path / "one_drift.in"
    inputs.write_text(
        "lattice.elements = from_inputs\n"
        "from_inputs.type = drift\n"
        "from_inputs.ds = 1.23\n"
    )

    sim.lattice.append(elements.Quad(ds=0.3, k=2.0, name="from_python"))

    sim.load_inputs_file(str(inputs))
    sim.init_lattice_elements_from_inputs()

    assert names_of(sim.lattice) == ["from_inputs"]
    assert sim.lattice[0].ds == pytest.approx(1.23)


def test_finalize_empties_the_lattice():
    simulation = ImpactX()
    simulation.lattice.append(elements.Drift(ds=1.0))

    simulation.finalize()

    assert len(simulation.lattice) == 0


# --- tracking ---------------------------------------------------------------------------


def test_tracking_an_empty_lattice_raises(sim):
    sim.init_envelope(sim.beam.ref, waterbag())

    with pytest.raises(RuntimeError, match="zero elements"):
        sim.track_envelope()


def test_a_hook_can_retune_the_element_being_tracked(sim):
    drift = elements.Drift(ds=0.5)
    sim.lattice.append(drift)

    def hook(s):
        assert s.tracking_element is drift
        s.tracking_element.ds = 0.25

    sim.hook["before_element"] = hook
    track_reference(sim)

    assert sim.beam.ref.s == pytest.approx(0.25)
    assert sim.tracking_element is None


@pytest.mark.parametrize("track", TRACKERS)
def test_a_hook_cannot_change_the_sequence(sim, track):
    sim.lattice.append(elements.Drift(ds=0.5, name="d"))
    seen = []

    def hook(s):
        with pytest.raises(RuntimeError, match="while tracking"):
            s.lattice.append(elements.Drift(ds=0.1))
        seen.append(len(s.lattice))

    sim.hook["before_element"] = hook
    track(sim)

    assert seen == [1]
    # editable again once tracking is done
    sim.lattice.append(elements.Drift(ds=0.5, name="d2"))
    assert names_of(sim.lattice) == ["d", "d2"]


def _envelope_through(second_cavity):
    """Envelope moments after ``[rf, drift, second_cavity(rf), drift]``."""

    simulation = ImpactX()
    simulation.particle_shape = 2
    simulation.slice_step_diagnostics = False
    simulation.diagnostics = False
    simulation.init_grids()

    ref = simulation.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(100.0)
    simulation.init_envelope(ref, waterbag())

    rf = elements.RFCavity(
        ds=1.0,
        escale=20.0,
        freq=1.3e9,
        phase=-89.5,
        cos_coefficients=[2.0],
        sin_coefficients=[0.0],
        mapsteps=10,
    )
    drift = elements.Drift(ds=0.5)
    simulation.lattice.extend([rf, drift, second_cavity(rf), drift])
    try:
        simulation.track_envelope()
        # per-particle extrema are NaN for an envelope, and NaN never compares equal
        return {
            name: value
            for name, value in simulation.envelope.beam_moments(ref).items()
            if not math.isnan(value)
        }
    finally:
        simulation.finalize()


def test_an_element_at_two_positions_tracks_like_two_copies():
    """Per-slice state of an element (the cavity's linear map) must not leak between
    the positions it occupies."""

    shared = _envelope_through(lambda rf: rf)
    copied = _envelope_through(lambda rf: rf.copy())

    assert shared == copied


# --- selections -------------------------------------------------------------------------


def test_a_selection_goes_stale_when_elements_move():
    lattice = elements.KnownElementsList(
        [elements.Quad(ds=0.1, k=1.0, name=f"q{i}") for i in range(3)]
        + [elements.Drift(ds=0.1, name="d0")]
    )
    quads = lattice.select(kind="Quad")

    # retuning moves nothing, so the selection stays usable
    lattice[0].k = 3.0
    assert len(quads) == 3
    assert quads[0].k == 3.0

    lattice.insert(0, elements.Drift(ds=0.1, name="new_head"))

    for use in (len, list, lambda s: s[0], lambda s: s.replace_with_drifts()):
        with pytest.raises(RuntimeError, match="no longer valid"):
            use(quads)
    # nothing was written through the stale positions
    assert names_of(lattice) == ["new_head", "q0", "q1", "q2", "d0"]


def test_filtered_edits_keep_subclassed_elements():
    class Tagged(elements.Drift):
        def __init__(self, tag):
            super().__init__(ds=1.0)
            self.tag = tag

    class MyBend(elements.ExactSbend):
        pass

    keep = Tagged("survivor")
    lattice = elements.KnownElementsList(
        [
            keep,
            elements.Quad(ds=0.3, k=1.0, name="drop"),
            MyBend(ds=1.0, phi=30.0, B=0.0, name="bend"),
        ]
    )

    lattice.select(name="drop").delete()
    assert lattice[0] is keep and lattice[0].tag == "survivor"

    # a subclass is still of its element kind: selected by it, and given its drift
    lattice.select(kind="ExactSbend").replace_with_drifts(model="match")
    assert lattice[0] is keep
    assert type(lattice[1]) is elements.ExactDrift
    assert lattice[1].name == "bend"


def test_replace_each_needs_a_template_it_can_copy():
    class MyDrift(elements.Drift):
        pass

    lattice = elements.KnownElementsList(
        [elements.Quad(ds=0.1, k=1.0, name=f"q{i}") for i in range(3)]
    )

    with pytest.raises(TypeError, match="copy"):
        lattice.select(kind="Quad").replace_each(MyDrift(ds=0.5))
    assert [type(element) for element in lattice] == [elements.Quad] * 3

    # a replacement per position, not one element repeated
    lattice.select(kind="Quad").replace_each(elements.Drift(ds=0.5))
    assert names_of(lattice) == ["q0", "q1", "q2"]
    assert lattice[0] is not lattice[1]


def test_insert_element_every_ds_keeps_the_elements_it_does_not_split():
    class MyQuad(elements.Quad):
        tag = "mine"

    shared = MyQuad(ds=0.2, k=1.0, name="shared")
    source = elements.KnownElementsList([shared, elements.Drift(ds=1.0), shared])

    result = elements.transformation.insert_element_every_ds(
        source, 1.0, elements.Marker("m")
    )
    del source
    gc.collect()

    held = [element for element in result if element is shared]
    assert len(held) == 2
    assert held[0].tag == "mine"
