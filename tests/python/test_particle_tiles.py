#!/usr/bin/env python3
#
# Copyright 2022-2024 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from impactx import ImpactX, ImpactXParIter, amr, distribution, elements


def test_particle_tiles():
    """
    This tests shows low-level, AMReX particle iteration.
    See other tests for high-level particle access.
    """
    sim = ImpactX()

    sim.particle_shape = 2
    sim.space_charge = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    # init particle beam
    kin_energy_MeV = 2.0e3
    bunch_charge_C = 1.0e-9
    npart = 10000

    #   reference particle
    pc = sim.particle_container()
    ref = pc.ref_particle()
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

    assert pc.total_number_of_particles() == npart

    # init accelerator lattice
    fodo = [
        elements.Drift(name="d1", ds=0.25),
        elements.Quad(name="q1", ds=1.0, k=1.0),
        elements.Drift(name="d2", ds=0.5),
        elements.Quad(name="q2", ds=1.0, k=-1.0),
        elements.Drift(name="d3", ds=0.25),
    ]
    sim.lattice.extend(fodo)

    # simulate
    sim.evolve()

    # access local particles
    for lvl in range(pc.finest_level + 1):
        for pti in ImpactXParIter(pc, level=lvl):
            soa = pti.soa()
            real_arrays = soa.get_real_data()
            print(real_arrays)

    # finalize simulation
    sim.finalize()


if __name__ == "__main__":
    test_particle_tiles()

    # clean simulation shutdown
    if amr.initialized():
        amr.finalize()
