"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from trame.widgets import vuetify

from ..trame_setup import setup_server

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Common toolbar elements
# -----------------------------------------------------------------------------


class ToolbarElements:
    """
    Helper functions to create
    Vuetify UI elements for toolbar.
    """

    @staticmethod
    def plot_options():
        vuetify.VSelect(
            v_model=("active_plot", "1D plots over s"),
            items=("plot_options",),
            label="Select plot to view",
            id="select_plot",
            hide_details=True,
            dense=True,
            style="max-width: 250px",
            disabled=("disableRunSimulationButton", True),
        )

    @staticmethod
    def run_simulation_button():
        vuetify.VBtn(
            "Run Simulation",
            id="run_simulation_button",
            style="background-color: #00313C; color: white; margin: 0 20px;",
            click=ctrl.run_simulation,
            disabled=("disableRunSimulationButton", True),
        )

    @staticmethod
    def show_simulation_complete():
        vuetify.VAlert(
                "Simulation Complete",
                v_model=("simulation_complete", False),
                id="simulation_complete",
                type="success",
                dense=True,
                classes="mt-4",
            )


# -----------------------------------------------------------------------------
# Content
# -----------------------------------------------------------------------------


class Toolbars:
    """
    Builds section toolbars for various pages.
    """

    @staticmethod
    def run_toolbar():
        """
        Builds toolbar for the 'Run' page.
        """

        vuetify.VSpacer(),
        ToolbarElements.show_simulation_complete(),
        ToolbarElements.run_simulation_button(),

    @staticmethod
    def analyze_toolbar():
        """
        Builds toolbar for the 'Analyze' page.
        """

        vuetify.VSpacer()
        ToolbarElements.plot_options()
