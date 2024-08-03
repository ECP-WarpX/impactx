#!/usr/bin/env python3
#
# Copyright 2022-2023 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from trame.app import get_server

server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller
import importlib

from impactx import ImpactX

from ...Input.distributionParametersCard.distributionMain import (
    save_distribution_parameters_to_file,
)
from ...Input.latticeConfigurationCard.latticeMain import save_latticeElements_to_file
from ..plot_phase_space.phaseSpaceSettings import adjusted_settings_plot

def run_simulation(save_png=True):
    """
    This tests using ImpactX and Pandas Dataframes
    """
    sim = ImpactX()

    sim.particle_shape = state.particle_shape
    sim.space_charge = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    # init particle beam
    kin_energy_MeV = state.kin_energy_MeV
    bunch_charge_C = state.bunch_charge_C
    npart = state.npart

    #   reference particle
    pc = sim.particle_container()
    ref = pc.ref_particle()
    ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

    #   particle bunch
    # distr = distribution.Waterbag(
    #     lambdaX=3.9984884770e-5,
    #     lambdaY=3.9984884770e-5,
    #     lambdaT=1.0e-3,
    #     lambdaPx=2.6623538760e-5,
    #     lambdaPy=2.6623538760e-5,
    #     lambdaPt=2.0e-3,
    #     muxpx=-0.846574929020762,
    #     muypy=0.846574929020762,
    #     mutpt=0.0,
    # )
    distr = save_distribution_parameters_to_file()

    sim.add_particles(bunch_charge_C, distr, npart)

    assert pc.total_number_of_particles() == npart

    # init accelerator lattice
    # fodo = [
    #     elements.Drift(0.25),
    #     elements.Quad(1.0, 1.0),
    #     elements.Drift(0.5),
    #     elements.Quad(1.0, -1.0),
    #     elements.Drift(0.25),
    # ]
    fodo = save_latticeElements_to_file()

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

    # fig = pc.plot_phasespace()
    fig = adjusted_settings_plot(pc)

    return fig
