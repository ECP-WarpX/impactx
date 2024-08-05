from trame.widgets import vuetify

# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------

from ..trame_setup import setup_server
server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Helpful
# -----------------------------------------------------------------------------


class Optimize:
    """
    User-Input section for the Optimize workflow.
    This is a base implementation without functionality
    thus far. (8/4/24)
    """

    @staticmethod
    def card():
        """
        User-Input section for Optimize workflow.
        """

        with vuetify.VCard(style="width: 680px; height: 300px"):
            with vuetify.VCardTitle("Optimizer Settings"):
                vuetify.VSpacer()
                vuetify.VIcon(
                    "mdi-information",
                    style="color: #00313C;",
                )
            vuetify.VDivider()
            with vuetify.VCardText():
                with vuetify.VRow():
                    with vuetify.VCol():
                        vuetify.VTextField(
                            label="Objective",
                            style="width: 500px",
                            dense=True,
                        )
                with vuetify.VRow():
                    with vuetify.VCol():
                        vuetify.VTextField(
                            label="Adjustment1",
                            style="width: 200px",
                            dense=True,
                        )
                    with vuetify.VCol():
                        vuetify.VTextField(
                            label="Adjustment2",
                            style="width: 200px",
                            dense=True,
                        )
                    with vuetify.VCol():
                        vuetify.VTextField(
                            label="Adjustment3",
                            style="width: 200px",
                            dense=True,
                        )
