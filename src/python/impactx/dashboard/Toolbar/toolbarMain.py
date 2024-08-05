from trame.widgets import vuetify

# from ..Analyze.plot_PhaseSpaceProjections.phaseSpace import outputTerminal
from ..Analyze.analyzeFunctions import AnalyzeFunctions
from ..Input.trameFunctions import TrameFunctions
from ..trame_setup import setup_server
from .exportTemplate import retrieve_state_content

# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Trigger
# -----------------------------------------------------------------------------


@ctrl.trigger("export")
def on_export_click():
    return retrieve_state_content()


# -----------------------------------------------------------------------------
# Common toolbar elements
# -----------------------------------------------------------------------------

state.selectedVisualization = None

TERMINAL_BUTTON_STYLES = {
    "background-color": "#2E86C1",
    "color": "white",
    "margin": "0 20px",
}


class ToolbarElements:
    """
    Helper functions to create
    Vuetify UI elements for toolbar.
    """

    @staticmethod
    def select_visualization():
        vuetify.VCombobox(
            placeholder="Select Visualization",
            v_model=("selectedVisualization",),
            items=(["Twiss Phase Space Ellipses"],),
            clearable=True,
            dense=True,
            hide_details=True,
            style="max-width: 250px",
        )

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
    def export_input_data():
        vuetify.VIcon(
            "mdi-download",
            style="color: #00313C; margin: 0 10px;",
            click="utils.download('input.in', trigger('export'), 'text/plain')",
            disabled=("disableRunSimulationButton", True),
        )

    @staticmethod
    def switch_theme():
        vuetify.VSwitch(
            v_model="$vuetify.theme.dark",
            hide_details=True,
        )

    @staticmethod
    def file_upload():
        """
        Allows users to upload file, but nothing more than that.
        """

        vuetify.VFileInput(
            label="Upload Input File",
            clearable=True,
            chips=True,
            show_size=True,
            dense=True,
            hide_details=True,
            style="max-width: 175px;",
        )

    def run_simulation():
        # ctrl.terminal_println("Running simulation...")
        # ctrl.terminal_println("Simulation complete.")
        AnalyzeFunctions.outputTerminal()

    @staticmethod
    def kill_button():
        return TrameFunctions.create_button("Kill")

    @staticmethod
    def stop_button():
        return TrameFunctions.create_button("Stop")

    @staticmethod
    def start_button():
        vuetify.VBtn(
            "START",
            style=TERMINAL_BUTTON_STYLES,
            classes="mx-1",
            click=ToolbarElements.run_simulation,
        )

    @staticmethod
    def checkbox_2d():
        vuetify.VCheckbox(
            label="2D",
            hide_details=True,
        )

    @staticmethod
    def checkbox_3d():
        vuetify.VCheckbox(label="3D", classes="px-2", hide_details=True)


# -----------------------------------------------------------------------------
# Content
# -----------------------------------------------------------------------------


class Toolbars:
    """
    Builds section toolbars for various pages.
    """

    @staticmethod
    def input_toolbar():
        """
        Builds toolbar for the 'Input' page.
        """

        # ToolbarElements.file_upload()
        vuetify.VSpacer()
        # ToolbarElements.select_visualization()
        ToolbarElements.run_simulation_button()
        # ToolbarElements.export_input_data()
        # ToolbarElements.switch_theme()

    @staticmethod
    def run_toolbar():
        """
        Builds toolbar for the 'Run' page.
        """

        ToolbarElements.stop_button(),
        ToolbarElements.start_button(),
        ToolbarElements.kill_button(),
        vuetify.VSpacer(),
        ToolbarElements.run_simulation_button(),
        # ToolbarElements.export_input_data(),
        # ToolbarElements.switch_theme(),

    @staticmethod
    def analyze_toolbar():
        """
        Builds toolbar for the 'Analyze' page.
        """

        vuetify.VSpacer()
        # ToolbarElements.checkbox_2d()
        # ToolbarElements.checkbox_3d()
        ToolbarElements.plot_options()
        ToolbarElements.run_simulation_button()
        # ToolbarElements.export_input_data()
        # ToolbarElements.switch_theme()
