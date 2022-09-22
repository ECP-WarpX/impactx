#!/usr/bin/env python3
#
# Copyright 2022 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import pytest

import amrex
from impactx import ImpactX, RefPart, distribution, elements


def test_impactx_fodo_file():
    """
    This tests an equivalent to main.cpp in C++
    """
    sim = ImpactX()

    sim.load_inputs_file("examples/fodo/input_fodo.in")

    sim.init_grids()
    sim.init_beam_distribution_from_inputs()
    sim.init_lattice_elements_from_inputs()

    sim.evolve()


def test_impactx_nofile():
    """
    This tests using ImpactX without an inputs file
    """
    sim = ImpactX()

    sim.particle_shape = 2
    sim.slice_step_diagnostics = True
    sim.init_grids()

    # init particle beam
    energy_MeV = 2.0e3
    bunch_charge_C = 1.0e-9
    npart = 10000

    #   reference particle
    ref = sim.particle_container().ref_particle()
    ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_energy_MeV(energy_MeV)

    #   particle bunch
    distr = distribution.Waterbag(
        sigmaX=3.9984884770e-5,
        sigmaY=3.9984884770e-5,
        sigmaT=1.0e-3,
        sigmaPx=2.6623538760e-5,
        sigmaPy=2.6623538760e-5,
        sigmaPt=2.0e-3,
        muxpx=-0.846574929020762,
        muypy=0.846574929020762,
        mutpt=0.0,
    )
    sim.add_particles(bunch_charge_C, distr, npart)

    assert sim.particle_container().TotalNumberOfParticles() == npart

    # init accelerator lattice
    fodo = [
        elements.Drift(0.25),
        elements.Quad(1.0, 1.0),
        elements.Drift(0.5),
        elements.Quad(1.0, -1.0),
        elements.Drift(0.25),
    ]
    #  assign a fodo segment
    # sim.lattice = fodo

    #  add 4 more FODO segments
    for i in range(4):
        sim.lattice.extend(fodo)

    # add 2 more drifts
    for i in range(4):
        sim.lattice.append(elements.Drift(0.25))

    print(sim.lattice)
    print(len(sim.lattice))
    assert len(sim.lattice) > 5

    sim.evolve()


def test_impactx_noparticles():
    """
    This tests using ImpactX without particles:
    must throw a user-friendly runtime error
    """
    sim = ImpactX()

    sim.particle_shape = 2
    sim.init_grids()

    # init particle beam
    energy_MeV = 2.0e3

    #   reference particle
    ref = sim.particle_container().ref_particle()
    ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_energy_MeV(energy_MeV)
    #   particle bunch: init intentionally missing

    # init accelerator lattice
    sim.lattice.append(elements.Drift(0.5))

    with pytest.raises(
        RuntimeError, match="No particles found. Cannot run evolve without a beam."
    ):
        sim.evolve()


def test_impactx_noshape():
    """
    This tests using ImpactX without particle shape:
    must throw a user-friendly runtime error
    """
    sim = ImpactX()

    # "sim.particle_shape = order" is intentionally missing

    with pytest.raises(
        RuntimeError,
        match="particle_shape is not set, cannot initialize grids with guard cells.",
    ):
        sim.init_grids()

    with pytest.raises(
        RuntimeError,
        match="particle_shape is not set yet",
    ):
        print(sim.particle_shape)

    # correct the mistake and keep going
    sim.particle_shape = 2
    sim.init_grids()


def test_impactx_resting_refparticle():
    """
    This tests using ImpactX with a resting reference particle:
    must throw a user-friendly runtime error
    """
    sim = ImpactX()

    sim.particle_shape = 2
    sim.init_grids()

    # init particle beam
    #   reference particle: init intentionally missing
    #   particle bunch
    gaussian = distribution.Gaussian(
        sigmaX=4.0e-5,
        sigmaY=5.0e-5,
        sigmaT=1.0e-3,
        sigmaPx=1.0e-5,
        sigmaPy=3.0e-5,
        sigmaPt=2.0e-3,
    )
    with pytest.raises(
        RuntimeError,
        match="add_particles: Reference particle charge not yet set!",
    ):
        sim.add_particles(bunch_charge=0.0, distr=gaussian, npart=10)

    sim.lattice.append(elements.Drift(0.25))

    with pytest.raises(
        RuntimeError,
        match="The reference particle energy is zero. Not yet initialized?",
    ):
        sim.evolve()


def test_impactx_no_elements():
    """
    This tests using ImpactX without a beamline lattice:
    must throw a user-friendly runtime error
    """
    sim = ImpactX()

    sim.load_inputs_file("examples/fodo/input_fodo.in")

    sim.init_grids()
    sim.init_beam_distribution_from_inputs()
    # intentionally skipped: no lattice initialized

    with pytest.raises(
        RuntimeError,
        match="Beamline lattice has zero elements. Not yet initialized?",
    ):
        sim.evolve()


def test_impactx_change_resolution():
    """
    This test checks we can change the grid resolution.
    This is currently a work-around because we cannot yet change the cells
    after the simulation object as been created.
    """
    pp_amr = amrex.ParmParse("amr")
    pp_amr.addarr("n_cell", [16, 24, 32])

    sim = ImpactX()

    # Future:
    # sim.ncell = [16, 24, 32]
    # sim.domain = amrex.RealBox([1., 2., 3.], [4., 5., 6.])

    sim.particle_shape = 2
    sim.slice_step_diagnostics = False
    sim.diagnostics = False
    sim.init_grids()

    assert sim.n_cell == [16, 24, 32]

    rho = sim.rho(lev=0)
    assert rho.nComp == 1
    assert rho.size == 1
    assert rho.num_comp == 1
    # assert rho.nGrowVect == [2, 2, 2]
    print(f"rho.nGrowVect={rho.nGrowVect}")
    assert iter(rho).length > 0
    assert not rho.is_all_cell_centered
    assert rho.is_all_nodal
