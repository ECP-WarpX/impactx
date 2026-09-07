#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""
Test the lifetime of a simulation: that inputs can be set before the simulation
is initialized, and that multiple simulations can run in the same process, e.g.,
for parameter scans and optimization loops that create one simulation per
iteration.

Each scenario runs in a fresh subprocess that does not import mpi4py, so that
ImpactX controls the MPI lifetime itself. In the pytest process, ``conftest.py``
already initializes MPI (via mpi4py) and AMReX, which would hide regressions in
the ImpactX-owned lifetime.
"""

import subprocess
import sys

import pytest

# a small, quick simulation: tracked in a helper function so we can run it twice
SIM_HELPER = """
from impactx import ImpactX, distribution, elements

distr = distribution.Waterbag(
    lambdaX=3.9984884770e-5,
    lambdaY=3.9984884770e-5,
    lambdaT=1.0e-3,
    lambdaPx=2.6623538760e-5,
    lambdaPy=2.6623538760e-5,
    lambdaPt=2.0e-3,
    muxpx=-0.846574929020762,
    muypy=0.846574929020762,
    mutpt=0.0,
)


def run_sim():
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.tiny_profiler = False
    sim.init_grids()

    sim.beam.ref.set_species("electron").set_kin_energy_MeV(2.0e3)
    sim.add_particles(1.0e-9, distr, 1000)
    sim.lattice.extend(
        [
            elements.Drift(name="d1", ds=0.25, nslice=1),
            elements.Quad(name="q1", ds=1.0, k=1.0, nslice=1),
            elements.Drift(name="d2", ds=0.5, nslice=1),
        ]
    )
    sim.track_particles()

    rbc = sim.beam.reduced_beam_characteristics()
    sim.finalize()

    return rbc
"""


def run_snippet(snippet):
    """Run a Python snippet in a subprocess.

    Returns the completed process and its combined output.
    """
    proc = subprocess.run(
        [sys.executable, "-c", snippet],
        capture_output=True,
        text=True,
        timeout=600,
    )
    return proc, proc.stdout + " " + proc.stderr


@pytest.mark.manages_amrex
def test_two_simulations_in_one_process():
    """Two simulations, one after the other, in a single process."""
    snippet = (
        SIM_HELPER
        + """
first = run_sim()
second = run_sim()

# same input, same random seed: identical results
assert first["sig_x"] == second["sig_x"], (first["sig_x"], second["sig_x"])
print("TWO_RUNS_OK")
"""
    )
    proc, output = run_snippet(snippet)
    assert proc.returncode == 0, output
    assert "TWO_RUNS_OK" in output, output


@pytest.mark.manages_amrex
def test_inputs_do_not_leak_between_simulations():
    """A finalized simulation does not pass its inputs on to the next one."""
    snippet = """
from impactx import ImpactX

sim = ImpactX()
sim.particle_shape = 2
sim.space_charge = True
sim.dynamic_size = True
sim.prob_relative = [3.0]
sim.n_cell = [16, 16, 16]
sim.init_grids()
assert sim.space_charge == "3D"
sim.finalize()
del sim

# a new simulation does not inherit the inputs of the previous one
sim = ImpactX()
for prop in ["space_charge", "prob_relative", "n_cell"]:
    try:
        value = getattr(sim, prop)
    except RuntimeError:
        pass  # not set: as expected for a fresh simulation
    else:
        raise AssertionError(f"{prop} leaked from the previous simulation: {value}")
sim.finalize()
print("NO_LEAK_OK")
"""
    proc, output = run_snippet(snippet)
    assert proc.returncode == 0, output
    assert "NO_LEAK_OK" in output, output


@pytest.mark.manages_amrex
def test_externally_initialized_amrex_is_not_finalized():
    """ImpactX only finalizes AMReX if it initialized AMReX itself."""
    snippet = """
import amrex.space3d as amr
from impactx import ImpactX

amr.initialize(["amrex.verbose=0", "particles.do_tiling=1"])

sim = ImpactX()
sim.particle_shape = 2
sim.init_grids()
sim.finalize()
del sim

assert amr.initialized(), "AMReX was finalized by ImpactX, but not initialized by it"

# a second simulation in the same AMReX context
sim = ImpactX()
sim.particle_shape = 2
sim.init_grids()
sim.finalize()
del sim

assert amr.initialized()
amr.finalize()
print("EXTERNAL_AMREX_OK")
"""
    proc, output = run_snippet(snippet)
    assert proc.returncode == 0, output
    assert "EXTERNAL_AMREX_OK" in output, output


@pytest.mark.manages_amrex
def test_warning_inputs_before_init_grids():
    """Warning logger inputs can be set before the simulation is initialized."""
    snippet = """
from impactx import ImpactX

sim = ImpactX()
sim.abort_on_warning_threshold = "high"
sim.abort_on_unused_inputs = 0
sim.always_warn_immediately = 1
sim.particle_shape = 2
sim.init_grids()

assert sim.abort_on_warning_threshold == "high"
assert sim.abort_on_unused_inputs == 0

sim.finalize()
print("WARNING_INPUTS_OK")
"""
    proc, output = run_snippet(snippet)
    assert proc.returncode == 0, output
    assert "WARNING_INPUTS_OK" in output, output
