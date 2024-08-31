from trame.widgets import vuetify

from ..generalFunctions import generalFunctions
from ...trame_setup import setup_server

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------


class csrConfiguration:

    @staticmethod
    def card():
        """
        Creates UI content for CSR configuration
        """

        with vuetify.VCard(v_show="CSR"):
            with vuetify.VCardTitle("CSR Configuration"):
                vuetify.VSpacer()
                vuetify.VIcon(
                    "mdi-information",
                    classes="ml-2",
                    click=lambda: generalFunctions.documentation("csrConfiguration"),
                    style="color: #00313C;",
                )
            vuetify.VDivider()
            with vuetify.VCardText():
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(cols=12):
                        vuetify.VTextField(
                            label="CSR Bins",
                            v_model=("CSR_bins",150),
                            type="number",
                            dense=True,
                        )
