#!/usr/bin/env python3
#
# Copyright 2022-2023 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from impactx import ImpactX, distribution, elements, push


def test_element_push():
    """
    This tests using ImpactX without a lattice.
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

    pc = sim.particle_container()
    assert pc.total_number_of_particles() == npart

    # init accelerator lattice
    fodo = [
        elements.Drift(0.25),
    ]
    sim.lattice.extend(fodo)

    sim.evolve()

    # Push manually through a few elements
    elements.Quad(1.0, 1.0).push(pc)
    elements.Drift(0.5).push(pc)
    elements.Quad(1.0, -1.0).push(pc)

    # alternative formulation
    push(pc, elements.Drift(0.25))

    # finalize simulation
    sim.finalize()
