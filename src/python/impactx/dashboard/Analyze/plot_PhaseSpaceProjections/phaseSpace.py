"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from ...trame_setup import setup_server

server, state, ctrl = setup_server()

import base64
import io

from impactx import Config, ImpactX

from ...Input.distributionParameters.distributionMain import distribution_parameters
from ...Input.latticeConfiguration.latticeMain import lattice_elements
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


def run_simulation():
    """
    This tests using ImpactX and Pandas Dataframes
    """
    sim = ImpactX()

    # space charge selections
    sim.n_cell = state.n_cell
    sim.particle_shape = state.particle_shape
    sim.space_charge = state.space_charge
    sim.dynamic_size = state.dynamic_size

    sim.slice_step_diagnostics = True
    sim.init_grids()

    # init particle beam
    kin_energy_MeV = state.kin_energy_MeV
    bunch_charge_C = state.bunch_charge_C
    npart = state.npart

    #   reference particle
    pc = sim.particle_container()
    ref = pc.ref_particle()
    ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

    distribution = distribution_parameters()
    sim.add_particles(bunch_charge_C, distribution, npart)

    lattice_configuration = lattice_elements()

    sim.lattice.extend(lattice_configuration)

    # simulate
    sim.evolve()

    fig = adjusted_settings_plot(pc)
    fig_original = pc.plot_phasespace()

    if fig_original is not None:
        image_base64 = fig_to_base64(fig_original)
        state.image_data = f"data:image/png;base64, {image_base64}"

    sim.finalize()

    return fig
