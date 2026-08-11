#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
"""``sim.lattice`` reports the elements that are actually there, and edits are all-or-nothing."""

import pytest

from impactx import ImpactX, elements


@pytest.fixture
def sim():
    simulation = ImpactX()
    simulation.particle_shape = 2
    simulation.slice_step_diagnostics = False
    simulation.diagnostics = False
    yield simulation
    simulation.finalize()


def names_of(lattice):
    return [element.name for element in lattice]


def test_lattice_follows_a_replacement_of_the_same_length(sim, tmp_path, monkeypatch):
    """A lattice built in C++ replaces one built in Python, even at the same length.

    Parsing an input file clears the lattice and fills it again without Python being
    told. The count is one before and one after, so the number of elements alone cannot
    reveal that the element changed.
    """

    inputs = tmp_path / "one_drift.in"
    inputs.write_text(
        "\n".join(
            [
                "beam.npart = 100",
                "beam.units = static",
                "beam.kin_energy = 2.0e3",
                "beam.charge = 1.0e-9",
                "beam.particle = electron",
                "beam.distribution = waterbag",
                "beam.lambdaX = 1.0e-4",
                "beam.lambdaY = 1.0e-4",
                "beam.lambdaT = 1.0e-3",
                "beam.lambdaPx = 1.0e-5",
                "beam.lambdaPy = 1.0e-5",
                "beam.lambdaPt = 1.0e-3",
                "lattice.elements = from_cxx",
                "from_cxx.type = drift",
                "from_cxx.ds = 1.23",
                "",
            ]
        )
    )

    monkeypatch.chdir(tmp_path)
    sim.init_grids()
    sim.lattice.append(elements.Quad(ds=0.3, k=2.0, name="from_python"))
    assert names_of(sim.lattice) == ["from_python"]

    sim.load_inputs_file(str(inputs))
    sim.init_lattice_elements_from_inputs()

    assert len(sim.lattice) == 1
    assert names_of(sim.lattice) == ["from_cxx"]
    assert sim.lattice[0].ds == pytest.approx(1.23)


def test_assigning_a_bad_element_leaves_the_lattice_alone(sim):
    """A rejected assignment is not allowed to destroy what was there."""

    sim.init_grids()
    sim.lattice.extend(
        [
            elements.Drift(ds=0.1, name="d0"),
            elements.Drift(ds=0.2, name="d1"),
        ]
    )

    with pytest.raises(TypeError):
        sim.lattice = [elements.Quad(ds=0.3, k=1.0, name="q"), "not an element"]

    assert names_of(sim.lattice) == ["d0", "d1"]


def test_assigning_the_lattice_to_itself_keeps_it(sim):
    """``sim.lattice = sim.lattice`` is a no-op, not a way to empty the lattice."""

    sim.init_grids()
    sim.lattice.extend(
        [
            elements.Drift(ds=0.1, name="d0"),
            elements.Drift(ds=0.2, name="d1"),
        ]
    )

    sim.lattice = sim.lattice

    assert names_of(sim.lattice) == ["d0", "d1"]


def test_finalize_releases_elements_without_init_grids(sim):
    """``finalize()`` empties the lattice whether or not grids were ever initialized.

    The elements have been finalized by then, so leaving them in place would keep
    elements that are done with -- a beam monitor with its series closed -- trackable.
    """

    sim.lattice.append(elements.Drift(ds=1.0, name="d"))
    assert len(sim.lattice) == 1

    sim.finalize()

    assert len(sim.lattice) == 0
