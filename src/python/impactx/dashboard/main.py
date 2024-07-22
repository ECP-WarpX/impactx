from trame.app import get_server
from trame.ui.vuetify import SinglePageWithDrawerLayout
from trame.widgets import vuetify, router, xterm
from trame.ui.router import RouterViewLayout

from Input.trameFunctions import trameFunctions

from Toolbar.toolbarMain import toolbars

from Input.inputParametersCard.inputMain import inputParameters
from Input.distributionParametersCard.distributionMain import distributionParameters
from Input.latticeConfigurationCard.latticeMain import latticeConfiguration
from Analyze.plotsMain import AnalyzeSimulation
from Input.Visualiztion.twiss_phase_space_ellipse.x_px import visualizeTwiss
from Optimize.optimizeMain import Optimize

# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------

server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller

# -----------------------------------------------------------------------------
# Router Views
# -----------------------------------------------------------------------------

inputParameters = inputParameters()

with RouterViewLayout(server, "/Input"):
    with vuetify.VContainer(fluid=True):
        with vuetify.VRow():
            with vuetify.VCol(cols="auto", classes="pa-2"):
                with vuetify.VRow(no_gutters=True):
                    with vuetify.VCol(cols="auto", classes="pa-2"):
                        inputParameters.card()
                    with vuetify.VCol(cols="auto", classes="pa-2"):
                        distributionParameters.card()
                with vuetify.VRow(no_gutters=True):
                    with vuetify.VCol(cols="auto", classes="pa-2"):
                        latticeConfiguration.card()
            with vuetify.VCol(cols="auto", classes="pa-2", v_if="selectedVisualization == 'Twiss Phase Space Ellipses'"):
                with vuetify.VRow(no_gutters=True):
                    with vuetify.VCol(cols="auto", classes="pa-2"):
                        visualizeTwiss.card_x()
                    with vuetify.VCol(cols="auto", classes="pa-2"):
                        visualizeTwiss.card_t()
                with vuetify.VRow(no_gutters=True):
                    with vuetify.VCol(cols="auto", classes="pa-2"):
                        visualizeTwiss.card_y()

with RouterViewLayout(server, "/Analyze"):
        with vuetify.VContainer(fluid=True):
            with vuetify.VRow(no_gutters=True, classes="fill-height"):
                with vuetify.VCol(cols="auto", classes="pa-2 fill-height"):
                    AnalyzeSimulation.card()
                with vuetify.VCol(classes="pa-2 d-flex align-center justify-center fill-height"):
                    AnalyzeSimulation.plot()

with RouterViewLayout(server, "/Optimize"):
        with vuetify.VContainer(fluid=True):
            with vuetify.VRow(no_gutters=True, classes="fill-height"):
                with vuetify.VCol(cols="auto", classes="pa-2 fill-height"):
                    Optimize.card()

# -----------------------------------------------------------------------------
# GUI
# -----------------------------------------------------------------------------

def application():
    with SinglePageWithDrawerLayout(server) as layout:
        layout.title.hide()
        with layout.toolbar:
            with vuetify.Template(v_if="$route.path == '/Input'"):
                toolbars.latticeToolbar()
            with vuetify.Template(v_if="$route.path == '/Analyze'"):
                toolbars.analyzeToolbar()
            with vuetify.Template(v_if="$route.path == '/Run'"):
                toolbars.runToolbar()
            
        with layout.drawer as drawer:
            drawer.width = 200
            with vuetify.VList():
                vuetify.VSubheader("Simulation")
            trameFunctions.create_route("Input","mdi-file-edit")
            trameFunctions.create_route("Optimize","mdi-trending-up")
            trameFunctions.create_route("Run", "mdi-play")
            trameFunctions.create_route("Analyze", "mdi-chart-box-multiple")

        with layout.content:
            router.RouterView()
            with xterm.XTerm(shell=["/bin/bash"], v_if="$route.path == '/Run'") as term:
                ctrl.clear = term.clear
    return layout

class JupyterMainApplication:
    def __init__(self):
        self.ui = self.generate_ui()

    def generate_ui(self):
        return application()

application()
# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    server.start()

'''
from main import JupyterMainApplication

# Create new application instance
app = JupyterMainApplication()

# Let's start the server by waiting for its UI to be ready
await app.ui.ready

# Put the UI into the resulting cell
app.ui
'''
