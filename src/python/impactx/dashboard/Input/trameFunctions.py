from trame.app import get_server
from trame.widgets import vuetify

# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------

server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller

# -----------------------------------------------------------------------------
# Code
# -----------------------------------------------------------------------------

terminal_button_styles = {
    "background-color": "#2E86C1",
    "color": "white",
    "margin": "0 20px",
}


class trameFunctions:

    @staticmethod
    def create_route(route_title, mdi_icon):
        """
        Creates a route with a specified title and icon.
        :param route_title: The title of the route.
        :param mdi_icon: The icon to be used for the route.
        """

        state[route_title] = False  # Does not display route by default

        to = f"/{route_title}"
        click = f"{route_title} = true"

        with vuetify.VListItem(to=to, click=click):
            with vuetify.VListItemIcon():
                vuetify.VIcon(mdi_icon)
            with vuetify.VListItemContent():
                vuetify.VListItemTitle(route_title)

    @staticmethod
    def create_section(title, content, expand_section_index):
        with vuetify.VExpansionPanels(v_model=(expand_section_index,), accordion=True):
            with vuetify.VExpansionPanel():
                with vuetify.VExpansionPanelHeader():
                    vuetify.VCardText(title)
                with vuetify.VExpansionPanelContent():
                    vuetify.VCardText(content)

    @staticmethod
    def create_button(label):
        """
        Creates a Vuetify button component with the specified label and styles.
        :param label: The name of the button.
        """
                
        return vuetify.VBtn(
            label,
            style=terminal_button_styles,
            classes="mx-1",
        )
