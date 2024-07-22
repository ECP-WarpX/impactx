import os
from trame.app import get_server
from trame.widgets import vuetify, plotly, matplotlib

from Analyze.analyzeFunctions import analyzeFunctions
from Analyze.plot_phase_space.phaseSpace import run_simulation
from Analyze.plot_over_s.overS import line_plot_1d
from Run.optimize_triplet.run_triplet import run_optimize_triplet

# -----------------------------------------------------------------------------
# Start server
# -----------------------------------------------------------------------------

server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller

# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

# Call plot_over_s
def plot_over_s():
    fig = line_plot_1d(state.selected_headers, state.filtered_data)
    ctrl.plotly_figure_update(fig)

PLOTS = {
    "Plot Over S": plot_over_s,
    "Phase Space Plots": None,
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def validPlotOptions(simulationClicked):
    if simulationClicked:
        return list(PLOTS.keys())
    else:
        return ["Run Simulation To See Options"]

def load_dataTable_data():
    combined_files= analyzeFunctions.combine_files(REDUCED_BEAM_DATA, REF_PARTICLE_DATA)
    combined_files_data_converted_to_dictionary_format = analyzeFunctions.convert_to_dict(combined_files)
    data, headers = combined_files_data_converted_to_dictionary_format
    state.all_data = data
    state.all_headers = headers

# -----------------------------------------------------------------------------
# Defaults
# -----------------------------------------------------------------------------

BASE_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REDUCED_BEAM_DATA = os.path.join(BASE_PATH, 'diags', 'reduced_beam_characteristics.0.0')
REF_PARTICLE_DATA= os.path.join(BASE_PATH, 'diags', 'ref_particle.0.0')
DEFAULT_HEADERS = ["step", "s", "alpha_x", "alpha_y", "alpha_t"]

state.selected_headers = DEFAULT_HEADERS
state.plot_options = validPlotOptions(simulationClicked=False)
state.show_table = False
state.active_plot = None
state.filtered_data = []
state.all_data = []
state.all_headers = []

# -----------------------------------------------------------------------------
# Functions to update table/plot
# -----------------------------------------------------------------------------

def update_data_table():
    """
    Combines reducedBeam and refParticle files
    and updates data table upon column selection by user
    """
    load_dataTable_data()
    state.filtered_data = analyzeFunctions.filter_data(state.all_data, state.selected_headers)
    state.filtered_headers = analyzeFunctions.filter_headers(state.all_headers, state.selected_headers)

def update_plot():
    """
    Performs actions to display correct information,
    based on the plot optin selected by the user
    """
    if state.active_plot == "Plot Over S":
        update_data_table()
        ctrl.plotly_figure_update(plot_over_s())
        state.show_table = True
    elif state.active_plot == "Phase Space Plots":
        state.show_table = False
        # ctrl.matplotlib_figure_update(state.simulation_data) 

# -----------------------------------------------------------------------------
# State changes
# -----------------------------------------------------------------------------

@state.change("selected_headers")
def on_header_selection_change(selected_headers, **kwargs):
    state.filtered_headers = analyzeFunctions.filter_headers(state.all_headers, selected_headers)
    state.filtered_data = analyzeFunctions.filter_data(state.all_data, selected_headers)
 
@state.change("filtered_data", "active_plot")
def on_filtered_data_change(**kwargs):
    update_plot() 

@ctrl.add("run_simulation")
def run_simulation_and_store():
    workflow = state.selectedWorkflow
    state.plot_options = validPlotOptions(simulationClicked=True)
    if workflow == "DataFrameTest":
        state.simulation_data = run_simulation()
        update_plot()
    elif workflow == "Optimize Triplet":
        run_optimize_triplet()
    load_dataTable_data()
    
# -----------------------------------------------------------------------------
# GUI
# -----------------------------------------------------------------------------

class AnalyzeSimulation:
        
    def card():
        with vuetify.VContainer():
            with vuetify.VCard(v_if=("show_table")):
                with vuetify.VCol(style="width: 500px;"):
                    vuetify.VSelect(
                        v_model=("selected_headers",),
                        items=("all_headers",),
                        label="Select data to view",
                        multiple=True,
                    )
                    vuetify.VDivider()
                    vuetify.VDataTable(
                        headers=("filtered_headers",),
                        items=("filtered_data",),
                        header_class="centered-header",
                        dense=True,
                        height="250px",
                    )

    def plot():
        with vuetify.VContainer(v_if="active_plot === 'Plot Over S'", style="height: 90vh; width: 100vh;"):
            plotly_figure = plotly.Figure(display_mode_bar="true", config={"responsive": True})
            ctrl.plotly_figure_update = plotly_figure.update
        # with vuetify.VLayout(v_if="active_plot === 'Phase Space Plots'"):
            # matplotlib_figure = matplotlib.Figure(style="position: absolute")
            # ctrl.matplotlib_figure_update = matplotlib_figure.update
