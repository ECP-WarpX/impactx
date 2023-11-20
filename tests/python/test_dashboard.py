#!/usr/bin/env python3
#
# Copyright 2022-2023 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import asyncio

from impactx import ImpactX, distribution, elements, push


def test_dashboard():
    """
    This tests using ImpactX with an interactive dashboard.
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
    assert pc.TotalNumberOfParticles() == npart

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

    async def run():
        # start interactive dashboard
        dashboard = sim.dashboard()

        # run simulation
        sim.evolve()
        # TODO / Idea:
        # - add callbacks to a python "async def" function that calls
        #     await asyncio.sleep(0)
        #   to yield control to the event loop
        # - await sim.evolve_async()
        # - OR make evolve return an Awaitable object that has a
        #   __await__ method for each iteration
        # - add an option of the sorts of "pause" to evolve,
        #   fall into event loop yield busy loop and listen to step
        #   or evolve events / state changes
        # - await sim.evolve_async(pause=True)
        #   before calling the blocking dashboard is an option, too
        #   -> but, using the dashboard blocking might not be ideal for
        #      Jupyter usage...?

        await dashboard

    asyncio.run(run())


if __name__ == "__main__":
    test_dashboard()
