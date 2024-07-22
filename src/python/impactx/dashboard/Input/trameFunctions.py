from trame.app import get_server
from trame.widgets  import vuetify

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

class trameFunctions():

    def create_route(route_title, mdi_icon):
        state[route_title] = False # Does not display route by default

        to = f"/{route_title}"
        click =  f"{route_title} = true"

        with vuetify.VListItem(to=to, click=click):
            with vuetify.VListItemIcon():
                vuetify.VIcon(mdi_icon)
            with vuetify.VListItemContent():
                vuetify.VListItemTitle(route_title)
                
    def create_section(title, content, expand_section_index):
        with vuetify.VExpansionPanels(v_model=(expand_section_index,), accordion=True):
            with vuetify.VExpansionPanel():
                with vuetify.VExpansionPanelHeader():
                    vuetify.VCardText(title)
                with vuetify.VExpansionPanelContent():
                    vuetify.VCardText(content)
                        
    def create_button(label):
        return vuetify.VBtn(
            label,
            style=terminal_button_styles,
            classes="mx-1",
        )