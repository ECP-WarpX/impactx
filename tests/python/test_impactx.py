#!/usr/bin/env python3
#
# Copyright 2022-2023 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
import pytest
from conftest import basepath

from impactx import ImpactX, distribution, elements

# FIXME in AMReX via https://github.com/AMReX-Codes/amrex/pull/3727
# def test_impactx_module():
#    """
#    Tests the basic modules we provide.
#    """
#    print(f"version={impactx.__version__}")
#    assert impactx.__version__  # version must not be empty


def test_impactx_fodo_file():
    """
    This tests an equivalent to main.cpp in C++
    """
    sim = ImpactX()

    sim.load_inputs_file(basepath + "/examples/fodo/input_fodo.in")

    sim.init_grids()
    sim.init_beam_distribution_from_inputs()
    sim.init_lattice_elements_from_inputs()

    sim.evolve()

    # validate the results
    beam = sim.particle_container()
    num_particles = beam.total_number_of_particles()
    assert num_particles == 10000
    atol = 0.0  # ignored
    rtol = 2.2 * num_particles**-0.5  # from random sampling of a smooth distribution

    # in situ calculate the reduced beam characteristics
    rbc = beam.reduced_beam_characteristics()

    # see examples/fodo/analysis_fodo.py
    print("charge=", rbc["charge_C"])
    assert np.allclose(
        [
            rbc["sig_x"],
            rbc["sig_y"],
            rbc["sig_t"],
            rbc["emittance_x"],
            rbc["emittance_y"],
            rbc["emittance_t"],
            rbc["charge_C"],
        ],
        [
            7.5451170454175073e-005,
            7.5441588239210947e-005,
            9.9775878164077539e-004,
            1.9959540393751392e-009,
            2.0175015289132990e-009,
            2.0013820193294972e-006,
            -1.0e-9,
        ],
        rtol=rtol,
        atol=atol,
    )

    # finalize simulation
    sim.finalize()


def test_impactx_nofile():
    """
    This tests using ImpactX without an inputs file
    """
    sim = ImpactX()

    sim.particle_shape = 2
    sim.slice_step_diagnostics = True
    sim.init_grids()

    # init particle beam
    kin_energy_MeV = 2.0e3
    bunch_charge_C = 1.0e-9
    npart = 10000

    #   reference particle
    ref = sim.particle_container().ref_particle()
    ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

    #   particle bunch
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
    sim.add_particles(bunch_charge_C, distr, npart)

    assert sim.particle_container().total_number_of_particles() == npart

    # init accelerator lattice
    fodo = [
        elements.Drift(name="d1", ds=0.25),
        elements.Quad(name="q1", ds=1.0, k=1.0),
        elements.Drift(name="d2", ds=0.5),
        elements.Quad(name="q2", ds=1.0, k=-1.0),
        elements.Drift(name="d3", ds=0.25),
    ]
    #  assign a fodo segment
    # sim.lattice = fodo

    #  add 4 more FODO segments
    for i in range(4):
        sim.lattice.extend(fodo)

    # add 2 more drifts
    for i in range(4):
        sim.lattice.append(elements.Drift(name="d" + str(4 + i), ds=0.25))

    print(sim.lattice)
    print(len(sim.lattice))
    assert len(sim.lattice) > 5

    sim.evolve()

    # finalize simulation
    sim.finalize()


def test_impactx_noparticles():
    """
    This tests using ImpactX without particles:
    must throw a user-friendly runtime error
    """
    sim = ImpactX()

    sim.particle_shape = 2
    sim.init_grids()

    # init particle beam
    kin_energy_MeV = 2.0e3

    #   reference particle
    ref = sim.particle_container().ref_particle()
    ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)
    #   particle bunch: init intentionally missing

    # init accelerator lattice
    sim.lattice.append(elements.Drift(ds=0.5))

    with pytest.raises(
        RuntimeError, match="No particles found. Cannot run evolve without a beam."
    ):
        sim.evolve()

    # finalize simulation
    sim.finalize()


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
        match="algo.particle_shape is not set yet",
    ):
        print(sim.particle_shape)

    # correct the mistake and keep going
    sim.particle_shape = 2
    sim.init_grids()

    # finalize simulation
    sim.finalize()


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
        lambdaX=4.0e-5,
        lambdaY=5.0e-5,
        lambdaT=1.0e-3,
        lambdaPx=1.0e-5,
        lambdaPy=3.0e-5,
        lambdaPt=2.0e-3,
    )
    with pytest.raises(
        RuntimeError,
        match="add_particles: Reference particle charge not yet set!",
    ):
        sim.add_particles(bunch_charge=0.0, distr=gaussian, npart=10)

    sim.lattice.append(elements.Drift(ds=0.25))

    with pytest.raises(
        RuntimeError,
        match="The reference particle energy is zero. Not yet initialized?",
    ):
        sim.evolve()

    # finalize simulation
    sim.finalize()


def test_impactx_no_elements():
    """
    This tests using ImpactX without a beamline lattice:
    must throw a user-friendly runtime error
    """
    sim = ImpactX()

    sim.load_inputs_file(basepath + "/examples/fodo/input_fodo.in")

    sim.init_grids()
    sim.init_beam_distribution_from_inputs()
    # intentionally skipped: no lattice initialized

    with pytest.raises(
        RuntimeError,
        match="Beamline lattice has zero elements. Not yet initialized?",
    ):
        sim.evolve()

    # finalize simulation
    sim.finalize()


def test_impactx_change_resolution():
    """
    This test checks we can change the grid resolution.
    This is currently a work-around because we cannot yet change the cells
    after the simulation object as been created.
    """
    sim = ImpactX()

    sim.n_cell = [16, 24, 32]
    sim.particle_shape = 2
    sim.slice_step_diagnostics = False
    sim.diagnostics = False
    sim.init_grids()

    assert sim.n_cell == [16, 24, 32]

    rho = sim.rho(lev=0)
    assert rho.nComp == 1
    assert rho.size == 1
    assert rho.num_comp == 1
    # assert rho.n_grow_vect == [2, 2, 2]
    print(f"rho.n_grow_vect={rho.n_grow_vect}")
    assert iter(rho).length > 0
    assert not rho.is_all_cell_centered
    assert rho.is_all_nodal

    # finalize simulation
    sim.finalize()
