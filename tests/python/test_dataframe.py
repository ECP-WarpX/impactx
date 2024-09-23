#!/usr/bin/env python3
#
# Copyright 2022-2023 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import importlib

import matplotlib.pyplot as plt
import pytest

from impactx import ImpactX, amr, distribution, elements


@pytest.mark.skipif(
    importlib.util.find_spec("pandas") is None, reason="pandas is not available"
)
def test_df_pandas(save_png=True):
    """
    This tests using ImpactX and Pandas Dataframes
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
        elements.Drift("d1", 0.25),
        elements.Quad("q1", 1.0, 1.0),
        elements.Drift("d2", 0.5),
        elements.Quad("q2", 1.0, -1.0),
        elements.Drift("d3", 0.25),
    ]
    sim.lattice.extend(fodo)

    # simulate
    sim.evolve()

    # check local particles
    df = pc.to_df(local=True)
    print(df)

    # ensure the column heads are correctly labeled
    assert df.columns.tolist() == [
        "idcpu",
        "position_x",
        "position_y",
        "position_t",
        "momentum_x",
        "momentum_y",
        "momentum_t",
        "qm",
        "weighting",
    ]

    # compare number of global particles
    # FIXME
    # df = pc.to_df(local=False)
    # if df is not None:
    #    assert npart == len(df)
    #    assert df.columns.tolist() == ['idcpu', 'position_x', 'position_y', 'position_t', 'momentum_x', 'momentum_y', 'momentum_t', 'qm', 'weighting']

    # plot
    fig = pc.plot_phasespace()

    #   note: figure data available on MPI rank zero
    if fig is not None:
        fig.savefig("phase_space.png")
        if save_png:
            fig.savefig("phase_space.png")
        else:
            plt.show()

    # finalize simulation
    sim.finalize()


if __name__ == "__main__":
    test_df_pandas(save_png=False)

    # clean simulation shutdown
    amr.finalize()
