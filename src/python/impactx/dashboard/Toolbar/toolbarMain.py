
from trame.app import get_server
from trame.widgets import vuetify
from Toolbar.exportTemplate import retrieve_state_content

from Input.trameFunctions import trameFunctions
# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------

server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller

state.selectedWorkflow = "Optimize Triplet"
state.isSelectedWorkflow = None
state.selectedVisualization = None

# -----------------------------------------------------------------------------
# Trigger
# -----------------------------------------------------------------------------

@ctrl.trigger("export")
def on_export_click():
    return retrieve_state_content()

@state.change("selectedWorkflow")
def on_selectedWorkflow_change(selectedWorkflow, **kwargs):
    # print(f"Selected workflow is {selectedWorkflow}")
    if selectedWorkflow == None:
        state.isSelectedWorkflow = "Please select a workflow"

# -----------------------------------------------------------------------------
# Common toolbar elements
# -----------------------------------------------------------------------------

class toolbarElements:
    
    def select_workflow():
        vuetify.VCombobox(
        placeholder="Select Workflow",
        v_model=("selectedWorkflow",),
        items=(["DataFrameTest", "Optimize Triplet"],),
        clearable=True,
        error_messages=("isSelectedWorkflow",),
        dense=True,
        hide_details=True,
        style="max-width: 175px",
        classes="mr-2",
    )
    
    def select_visualization():
        vuetify.VCombobox(
            placeholder="Select Visualization",
            v_model=("selectedVisualization",),
            items=(["Twiss Phase Space Ellipses", "Lattice Visualization"],),
            clearable=True,
            dense=True,
            hide_details=True,
            style="max-width: 250px",
        )

    def plot_options():
        vuetify.VSelect(
            v_model=("active_plot", "1D plots over s"),
            items=("plot_options",),
            label="Select plot to view",
            hide_details=True,
            dense=True,
            style="max-width: 250px",
            disabled=("disableRunSimulationButton",True),
        )

    def run_simulation_button():
        vuetify.VBtn(
            "Run Simulation",
            style="background-color: #00313C; color: white; margin: 0 20px;",
            click=ctrl.run_simulation,
            disabled=("disableRunSimulationButton",True),
        )
    
    def export_input_data():
        vuetify.VIcon(
            "mdi-download",
            style="color: #00313C; margin: 0 10px;",
            click="utils.download('input.in', trigger('export'), 'text/plain')",
            disabled=("disableRunSimulationButton",True),
        )

    def switch_theme():
        vuetify.VSwitch(
            v_model="$vuetify.theme.dark",
            hide_details=True,
        )

    def file_upload():
        vuetify.VFileInput(
            #Allows users to upload file, but nothing more than that.
            label="Upload Input File",
            clearable=True,
            chips=True,
            show_size=True,
            dense=True,
            hide_details=True,
            style="max-width: 175px;",
        )

    def kill_button():
        return trameFunctions.create_button("Kill")
            
    def stop_button():
        return trameFunctions.create_button("Stop")

    def start_button():
        return trameFunctions.create_button("Start")

    def checkbox_2d():
        vuetify.VCheckbox(
            label="2D",
            hide_details=True,
        )
    def checkbox_3d():
        vuetify.VCheckbox(
            label="3D", 
            classes="px-2",
            hide_details=True)
    
# -----------------------------------------------------------------------------
# Content
# -----------------------------------------------------------------------------

class toolbars:

    def runToolbar():

        toolbarElements.stop_button(),
        toolbarElements.start_button(),
        toolbarElements.kill_button(),
        vuetify.VSpacer(),
        toolbarElements.select_workflow(),
        toolbarElements.run_simulation_button(),
        toolbarElements.export_input_data(),
        toolbarElements.switch_theme(),
    
    def analyzeToolbar():
        vuetify.VSpacer()
        toolbarElements.checkbox_2d()
        toolbarElements.checkbox_3d()
        toolbarElements.select_workflow()
        toolbarElements.plot_options()
        toolbarElements.run_simulation_button()
        toolbarElements.export_input_data()
        toolbarElements.switch_theme()
        
    def latticeToolbar():
        toolbarElements.file_upload()
        vuetify.VSpacer()
        toolbarElements.select_workflow()
        toolbarElements.select_visualization()
        toolbarElements.run_simulation_button()
        toolbarElements.export_input_data()
        toolbarElements.switch_theme()
