from trame.app import get_server

# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------

server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def main():
    """
    Launches Trame application server
    """
    server.start()
    return 0
