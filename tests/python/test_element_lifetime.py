#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""
Tests for element lifetime management with dynamic GPU data.

These tests verify that elements with dynamic GPU data (like SoftSolenoid,
RFCavity, etc.) can be:
1. Created before AMReX/ImpactX is initialized
2. Cleaned up properly even if never added to a lattice
3. Reused across multiple simulations (with AMReX finalize/init cycles)
"""

import pytest

from impactx import ImpactX, distribution, elements


def run_minimal_simulation(sim, lattice_elements):
    """
    Run a minimal simulation with given lattice elements.

    Parameters
    ----------
    sim : ImpactX
        The simulation object (must be initialized with init_grids)
    lattice_elements : list
        List of beamline elements to add to the lattice
    """
    kin_energy_MeV = 250.0
    bunch_charge_C = 1.0e-9
    npart = 100

    ref = sim.beam.ref
    ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_kin_energy_MeV(kin_energy_MeV)

    distr = distribution.Waterbag(
        lambdaX=1.0e-3,
        lambdaY=1.0e-3,
        lambdaT=1.0e-3,
        lambdaPx=1.0e-4,
        lambdaPy=1.0e-4,
        lambdaPt=1.0e-3,
    )
    sim.add_particles(bunch_charge_C, distr, npart)

    sim.lattice.extend(lattice_elements)
    sim.track_particles()


def create_soft_solenoid(name="sol"):
    """Create a SoftSolenoid element with minimal coefficients."""
    return elements.SoftSolenoid(
        name=name,
        ds=1.0,
        bscale=0.1,
        cos_coefficients=[0.5, 0.3, 0.1],
        sin_coefficients=[0.0, 0.0, 0.0],
        mapsteps=10,
        nslice=1,
    )


def create_soft_quadrupole(name="quad"):
    """Create a SoftQuadrupole element with minimal coefficients."""
    return elements.SoftQuadrupole(
        name=name,
        ds=0.5,
        gscale=1.0,
        cos_coefficients=[1.0, 0.5],
        sin_coefficients=[0.0, 0.0],
        mapsteps=10,
        nslice=1,
    )


@pytest.mark.manages_amrex
def test_element_creation_before_amrex_init():
    """
    Test that elements with dynamic GPU data can be created before AMReX init.

    This verifies that host-side data storage works without AMReX being
    initialized, and that the element can be used once AMReX is started.
    """
    # Create element BEFORE AMReX/ImpactX is initialized
    sol = create_soft_solenoid("sol_pre_init")
    quad = create_soft_quadrupole("quad_pre_init")

    # Now initialize ImpactX (which initializes AMReX) and run simulation
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    run_minimal_simulation(sim, [sol, quad])

    sim.finalize()


@pytest.mark.manages_amrex
def test_unused_element_cleanup():
    """
    Test that elements never added to lattice are cleaned up without error.

    This verifies that sim.finalize() handles elements that were created
    but never used in a simulation.
    """
    # Create elements before ImpactX init
    unused_sol = create_soft_solenoid("unused_sol")
    unused_quad = create_soft_quadrupole("unused_quad")

    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    # Create more elements after ImpactX init
    another_unused = create_soft_solenoid("another_unused")

    # Run simulation with DIFFERENT element, not the ones created above
    used_sol = create_soft_solenoid("used_sol")
    run_minimal_simulation(sim, [used_sol])

    # finalize() should NOT error even though unused_sol, unused_quad,
    # and another_unused were never added to lattice
    sim.finalize()

    # Keep references alive to ensure they weren't prematurely destroyed
    assert unused_sol is not None
    assert unused_quad is not None
    assert another_unused is not None


@pytest.mark.manages_amrex
def test_element_reuse_across_simulations():
    """
    Test that elements can be reused across multiple simulations.

    This verifies that release_gpu_data() preserves host data, allowing
    elements to be uploaded to GPU again in a subsequent simulation.
    """
    # Create element before first simulation
    reusable_sol = create_soft_solenoid("reusable_sol")

    # First simulation
    sim1 = ImpactX()
    sim1.particle_shape = 2
    sim1.space_charge = False
    sim1.diagnostics = False
    sim1.slice_step_diagnostics = False
    sim1.init_grids()

    run_minimal_simulation(sim1, [reusable_sol])

    sim1.finalize()

    # Second simulation - reuse the same element
    sim2 = ImpactX()
    sim2.particle_shape = 2
    sim2.space_charge = False
    sim2.diagnostics = False
    sim2.slice_step_diagnostics = False
    sim2.init_grids()

    # Use the SAME element again - its host data should still be valid
    run_minimal_simulation(sim2, [reusable_sol])

    sim2.finalize()


@pytest.mark.manages_amrex
def test_a_failed_run_does_not_leave_an_element_behind():
    """A run that ends in an exception must let go of the element it was tracking.

    ``sim.tracking_element`` holds a share of the element being tracked. A hook that raises
    used to leave that share in place: the release sat after the element loop, so the
    exception carried straight past it, and finalize() never let go either. A finalized
    simulation then still named an element it had already finished with.
    """

    sim = ImpactX()
    sim.particle_shape = 2
    sim.slice_step_diagnostics = False
    sim.init_grids()
    sim.beam.ref.set_species("electron").set_kin_energy_MeV(100.0)

    stopper = elements.Programmable()

    def stop(refpart):
        raise RuntimeError("stop tracking")

    stopper.ref_particle = stop
    sim.lattice.append(stopper)

    with pytest.raises(RuntimeError, match="stop tracking"):
        sim.track_reference(sim.beam.ref)

    # the traversal releases its share on the way out, exception or not
    assert sim.tracking_element is None

    # and finalize() leaves nothing naming an element either
    sim.finalize()
    assert sim.tracking_element is None
