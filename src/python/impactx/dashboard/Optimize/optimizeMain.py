from trame.app import get_server
from trame.widgets  import vuetify

# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------

server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller

# -----------------------------------------------------------------------------
# Helpful
# -----------------------------------------------------------------------------

class Optimize:
    
    def card():
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