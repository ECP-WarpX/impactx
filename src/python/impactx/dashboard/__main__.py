import sys

from trame.app import get_server
from trame.ui.router import RouterViewLayout
from trame.ui.vuetify import SinglePageWithDrawerLayout
from trame.widgets import router, vuetify, xterm

from .Analyze.plotsMain import AnalyzeSimulation
from .Input.distributionParameters.distributionMain import DistributionParameters
from .Input.inputParameters.inputMain import InputParameters
from .Input.latticeConfiguration.latticeMain import LatticeConfiguration
from .Input.trameFunctions import trameFunctions
from .Input.Visualiztion.twiss_phase_space_ellipse.x_px import VisualizeTwiss
from .Optimize.optimizeMain import Optimize
from .Toolbar.toolbarMain import Toolbars
from .start import main
# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------

server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller

# -----------------------------------------------------------------------------
# Router Views
# -----------------------------------------------------------------------------

inputParameters = InputParameters()

with RouterViewLayout(server, "/Input"):
    with vuetify.VContainer(fluid=True):
        with vuetify.VRow():
            with vuetify.VCol(cols="auto", classes="pa-2"):
                with vuetify.VRow(no_gutters=True):
                    with vuetify.VCol(cols="auto", classes="pa-2"):
                        inputParameters.card()
                    with vuetify.VCol(cols="auto", classes="pa-2"):
                        DistributionParameters.card()
                with vuetify.VRow(no_gutters=True):
                    with vuetify.VCol(cols="auto", classes="pa-2"):
                        LatticeConfiguration.card()
            with vuetify.VCol(
                cols="auto",
                classes="pa-2",
                v_if="selectedVisualization == 'Twiss Phase Space Ellipses'",
            ):
                with vuetify.VRow(no_gutters=True):
                    with vuetify.VCol(cols="auto", classes="pa-2"):
                        VisualizeTwiss.card_x()
                    with vuetify.VCol(cols="auto", classes="pa-2"):
                        VisualizeTwiss.card_t()
                with vuetify.VRow(no_gutters=True):
                    with vuetify.VCol(cols="auto", classes="pa-2"):
                        VisualizeTwiss.card_y()

with RouterViewLayout(server, "/Analyze"):
    with vuetify.VContainer(fluid=True):
        with vuetify.VRow(no_gutters=True, classes="fill-height"):
            with vuetify.VCol(cols="auto", classes="pa-2 fill-height"):
                AnalyzeSimulation.card()
            with vuetify.VCol(
                classes="pa-2 d-flex align-center justify-center fill-height"
            ):
                AnalyzeSimulation.plot()

with RouterViewLayout(server, "/Optimize"):
    with vuetify.VContainer(fluid=True):
        with vuetify.VRow(no_gutters=True, classes="fill-height"):
            with vuetify.VCol(cols="auto", classes="pa-2 fill-height"):
                Optimize.card()

# -----------------------------------------------------------------------------
# GUI
# -----------------------------------------------------------------------------


def application ():
    with SinglePageWithDrawerLayout(server) as layout:
        layout.title.hide()
        with layout.toolbar:
            with vuetify.Template(v_if="$route.path == '/Input'"):
                Toolbars.latticeToolbar()
            with vuetify.Template(v_if="$route.path == '/Analyze'"):
                Toolbars.analyzeToolbar()
            with vuetify.Template(v_if="$route.path == '/Run'"):
                Toolbars.runToolbar()

        with layout.drawer as drawer:
            drawer.width = 200
            with vuetify.VList():
                vuetify.VSubheader("Simulation")
            trameFunctions.create_route("Input", "mdi-file-edit")
            trameFunctions.create_route("Optimize", "mdi-trending-up")
            trameFunctions.create_route("Run", "mdi-play")
            trameFunctions.create_route("Analyze", "mdi-chart-box-multiple")

        with layout.content:
            router.RouterView()
            with xterm.XTerm(v_if="$route.path == '/Run'") as term:
                ctrl.clear = term.clear
                ctrl.terminal_println = term.writeln
    return layout

application()
# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    sys.exit(main())