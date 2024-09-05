from trame.widgets import vuetify

from ...trame_setup import setup_server
from ..generalFunctions import generalFunctions

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
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(classes="py-0"):
                        vuetify.VTextField(
                            label="CSR Bins",
                            v_model=("csr_bins", 150),
                            type="number",
                            dense=True,
                        )
