#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell, Eric Stern
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from pathlib import Path

import pytest

from impactx import Config, ImpactX, distribution, elements

io = pytest.importorskip("openpmd_api")

if not Config.have_openpmd:
    pytest.skip("ImpactX was compiled without openPMD support", allow_module_level=True)


def build_sim(lattice_elements, npart):
    """Build a small simulation with the given lattice elements."""
    sim = ImpactX()

    sim.particle_shape = 2
    sim.slice_step_diagnostics = False
    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(2.0e3)

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
    sim.add_particles(1.0e-9, distr, npart)

    sim.lattice.extend(lattice_elements)

    return sim


def track_and_finalize(sim):
    """Track particles and always finalize the simulation.

    Finalizing in all paths is essential: if tracking raises, pytest keeps
    the exception traceback (and thus ``sim``) alive beyond this test and
    beyond this test's AMReX lifetime. A later, delayed ``~ImpactX`` would
    then call ``amrex::Finalize()`` in the middle of an unrelated test.
    """
    try:
        sim.track_particles()
    finally:
        sim.finalize()


def check_series(name, particles, beam_moments, npart):
    """Validate the series written by a monitor with the given flags.

    Reads with linear access, which supports all iteration encodings.
    """
    files = sorted(Path("diags/openPMD").glob(f"{name}.*"))
    assert files, f"no openPMD series found for monitor '{name}'"

    series = io.Series(str(files[0]), io.Access.read_linear)
    n_iterations = 0
    for it in series.read_iterations():
        n_iterations += 1
        beam = it.particles["beam"]

        # reference particle attributes: always present
        for attr in ("s_ref", "beta_ref", "gamma_ref", "mass_ref", "charge_ref"):
            assert attr in beam.attributes

        # beam moments attributes: present only if enabled
        # this includes the current period (turn)
        for attr in ("sig_x", "sig_y", "sig_t", "emittance_x", "charge_C", "period"):
            assert (attr in beam.attributes) == beam_moments
        if beam_moments:
            assert beam.get_attribute("period") == 0  # single pass

        # per-particle records: present only if enabled
        if particles:
            assert "id" in beam
            assert "positionOffset" in beam
            assert beam["id"][io.Record_Component.SCALAR].shape == [npart]
        else:
            # zero-extent constant record for openPMD-api 0.17+ readers
            assert sorted(beam) == ["empty"]
            assert beam["empty"][io.Record_Component.SCALAR].shape == [0]
    assert n_iterations == 2
    series.close()


@pytest.mark.parametrize(
    "particles,beam_moments",
    [(True, True), (True, False), (False, True), (False, False)],
)
def test_beam_monitor_flags(particles, beam_moments):
    """
    This tests that particle output and beam-moments output of the
    BeamMonitor element can be independently turned on and off.
    """
    npart = 512
    name = f"mon_p{int(particles)}_m{int(beam_moments)}"

    monitor = elements.BeamMonitor(
        name=name, particles=particles, beam_moments=beam_moments
    )
    assert monitor.particles == particles
    assert monitor.beam_moments == beam_moments
    d = monitor.to_dict()
    assert d["particles"] == particles
    assert d["beam_moments"] == beam_moments

    sim = build_sim([monitor, elements.Drift(name="d1", ds=0.25), monitor], npart)
    track_and_finalize(sim)

    check_series(name, particles, beam_moments, npart)


def test_beam_monitor_moments_only_bp4():
    """
    The ADIOS2 BP4 engine drops attribute-only iterations written through
    I/O steps; this tests the random access fall-back used for bp4.
    """
    # skip based on the writer library: the reading openpmd_api module
    # may support ADIOS2 (e.g., PyPI wheels) even if ImpactX does not
    if not Config.openpmd_backends.get("adios2", False):
        pytest.skip("ImpactX was compiled without openPMD ADIOS2 support")

    npart = 512
    name = "mon_bp4"
    monitor = elements.BeamMonitor(
        name=name, backend="bp4", encoding="g", particles=False
    )

    sim = build_sim([monitor, elements.Drift(name="d1", ds=0.25), monitor], npart)
    track_and_finalize(sim)

    check_series(name, particles=False, beam_moments=True, npart=npart)


def test_beam_monitor_invariants_require_particles():
    """
    This tests that requesting nonlinear lens invariants (per-particle
    columns) without particle output fails loudly.
    """
    monitor = elements.BeamMonitor(name="mon_err", particles=False)
    monitor.nonlinear_lens_invariants = True

    sim = build_sim([monitor], npart=512)
    with pytest.raises(RuntimeError, match="nonlinear_lens_invariants"):
        track_and_finalize(sim)
