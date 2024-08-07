#!/usr/bin/env python3
#
# Copyright 2022-2023 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------

from ...trame_setup import setup_server

server, state, ctrl = setup_server()

import base64
import io


from impactx import ImpactX, Config

from ...Input.distributionParameters.distributionMain import (
    save_distribution_parameters,
)
from ...Input.latticeConfiguration.latticeMain import save_lattice_elements
from ..plot_PhaseSpaceProjections.phaseSpaceSettings import adjusted_settings_plot

# Call MPI_Init and MPI_Finalize only once:
if Config.have_mpi:
    from mpi4py import MPI  # noqa


def fig_to_base64(fig):
    """
    Puts png in trame-compatible form
    """
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


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
    distr = save_distribution_parameters()

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
    fodo = save_lattice_elements()

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
    fig_original = pc.plot_phasespace()

    if fig_original is not None:
        image_base64 = fig_to_base64(fig_original)
        state.image_data = f"data:image/png;base64, {image_base64}"

    sim.finalize()

    return fig
