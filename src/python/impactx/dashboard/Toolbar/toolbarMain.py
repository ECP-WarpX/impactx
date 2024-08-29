"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from trame.widgets import vuetify

from ..trame_setup import setup_server

server, state, ctrl = setup_server()

state.show_dashboard_alert = True

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
            hide_details=True,
            dense=True,
            style="max-width: 250px",
            disabled=("disableRunSimulationButton", True),
        )

    @staticmethod
    def run_simulation_button():
        vuetify.VBtn(
            "Run Simulation",
            style="background-color: #00313C; color: white; margin: 0 20px;",
            click=ctrl.run_simulation,
            disabled=("disableRunSimulationButton", True),
        )

    @staticmethod
    def dashboard_info():
        """
        Creates an alert box with dashboard information.
        """

        vuetify.VAlert(
            "ImpactX Dashboard is provided as a preview and continues to be developed. "
            "Thus, it may not yet include all the features available in ImpactX.",
            type="info",
            dense=True,
            dismissible=True,
            v_model=("show_dashboard_alert", True),
            classes="mt-4",
        )


class Toolbars:
    """
    Builds section toolbars for various pages.
    """

    @staticmethod
    def input_toolbar():
        """
        Builds toolbar for the 'Input' page.
        """

        (ToolbarElements.dashboard_info(),)

    @staticmethod
    def run_toolbar():
        """
        Builds toolbar for the 'Run' page.
        """

        (ToolbarElements.dashboard_info(),)
        (vuetify.VSpacer(),)
        (ToolbarElements.run_simulation_button(),)

    @staticmethod
    def analyze_toolbar():
        """
        Builds toolbar for the 'Analyze' page.
        """

        (ToolbarElements.dashboard_info(),)
        vuetify.VSpacer()
        ToolbarElements.plot_options()
