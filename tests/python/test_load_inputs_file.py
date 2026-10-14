#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""``sim.load_inputs_file``: the parameters of an inputs file, as for the impactx executable."""

import subprocess
import sys

import pytest

import amrex.space3d as amr
from impactx import ImpactX


def write_inputs(path, amrex_verbose=1):
    """An inputs file with one drift, which also sets the AMReX verbosity."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"amrex.verbose = {amrex_verbose}\n"
        "algo.particle_shape = 1\n"
        "diag.enable = false\n"
        "diag.slice_step_diagnostics = false\n"
        "lattice.elements = d1\n"
        "d1.type = drift\n"
        "d1.ds = 1.23\n"
    )
    return path


@pytest.fixture
def sim():
    """A simulation that is finalized at the end, also of a failing test."""
    simulation = ImpactX()
    yield simulation
    simulation.finalize()


@pytest.mark.manages_amrex
@pytest.mark.parametrize("amrex_verbose", [0, 1])
def test_load_inputs_file_initializes_amrex(sim, tmp_path, capfd, amrex_verbose):
    """Loading a file without an initialized AMReX initializes AMReX with it.

    AMReX reads amrex.verbose while it initializes: its start-up message shows
    whether the file came in time. A path with "=" and blanks, as in a parameter
    scan, is still read as a file.
    """
    inputs = write_inputs(tmp_path / "a0=2 scan" / "inputs", amrex_verbose)
    assert not amr.initialized()

    # a setting that AMReX reads while it initializes
    sim.abort_on_unused_inputs = 0
    sim.particle_shape = 3

    sim.load_inputs_file(inputs)

    assert amr.initialized()
    assert ("Initializing AMReX" in capfd.readouterr().out) == bool(amrex_verbose)
    # ImpactX defaults for AMReX parameters, where the file sets none
    assert sim.omp_threads == "nosmt"

    # the file comes after properties set before, and before those set after
    assert sim.particle_shape == 1
    sim.particle_shape = 2
    assert sim.particle_shape == 2

    sim.init_grids()
    sim.init_lattice_elements_from_inputs()
    assert len(sim.lattice) == 1
    assert sim.lattice[0].ds == pytest.approx(1.23)

    sim.finalize()
    assert not amr.initialized()


def test_warning_threshold_before_amrex():
    """The warning logger, which needs AMReX, takes its settings before AMReX is initialized.

    The warning logger keeps its threshold for the rest of the process, so this
    runs in a separate one.
    """
    subprocess.run(
        [
            sys.executable,
            "-c",
            "from impactx import ImpactX; ImpactX().abort_on_warning_threshold = 'high'",
        ],
        check=True,
    )


@pytest.mark.manages_amrex
def test_finalize_without_init_grids(sim, tmp_path):
    """finalize() shuts down the AMReX that loading a file initialized."""
    sim.load_inputs_file(write_inputs(tmp_path / "inputs"))
    assert amr.initialized()

    sim.finalize()
    assert not amr.initialized()


@pytest.mark.parametrize(
    "amrex_initialized",
    [pytest.param(False, marks=pytest.mark.manages_amrex), True],
    ids=["amrex_not_initialized", "amrex_initialized"],
)
@pytest.mark.parametrize(
    "problem, exception",
    [
        ("missing", FileNotFoundError),
        ("directory", IsADirectoryError),
        ("malformed", RuntimeError),
    ],
)
def test_unreadable_inputs_file(sim, tmp_path, amrex_initialized, problem, exception):
    """A file that cannot be read raises, and leaves AMReX as it was."""
    inputs = tmp_path / "inputs"
    if problem == "directory":
        inputs.mkdir()
    elif problem == "malformed":
        inputs.write_text("algo.particle_shape =\n")

    with pytest.raises(exception) as error:
        sim.load_inputs_file(inputs)
    if issubclass(exception, OSError):
        assert error.value.filename == str(inputs)
    assert amr.initialized() == amrex_initialized

    sim.load_inputs_file(write_inputs(tmp_path / "valid" / "inputs"))
    assert sim.particle_shape == 1
