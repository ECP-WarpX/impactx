from trame.widgets import vuetify

from ..trame_setup import setup_server

# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Code
# -----------------------------------------------------------------------------

TERMINAL_BUTTON_STYLES = {
    "background-color": "#2E86C1",
    "color": "white",
    "margin": "0 20px",
}


class TrameFunctions:
    """
    Contains functions containing Vuetify
    components.
    """

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
    def create_button(label):
        """
        Creates a Vuetify button component with the specified label and styles.
        :param label: The name of the button.
        """

        return vuetify.VBtn(
            label,
            style=TERMINAL_BUTTON_STYLES,
            classes="mx-1",
        )
