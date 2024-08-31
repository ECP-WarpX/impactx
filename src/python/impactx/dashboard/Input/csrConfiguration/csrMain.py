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
        Creates UI content for CSR.
        """

        with vuetify.VCard(v_show="CSR", style="width: 170px;"):
            with vuetify.VCardTitle("CSR"):
                vuetify.VSpacer()
                vuetify.VIcon(
                    "mdi-information",
                    classes="ml-2",
                    click=lambda: generalFunctions.documentation("CSR"),
                    style="color: #00313C;",
                )
            vuetify.VDivider()
            with vuetify.VCardText():
                with vuetify.VRow(classes="py-2"):
                    with vuetify.VCol():
                        vuetify.VTextField(
                            label="CSR Bins",
                            v_model=("csr_bins", 150),
                            type="number",
                            dense=True,
                        )
