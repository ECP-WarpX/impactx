from trame.app import get_server


def setup_server(client_type="vue2"):
    server = get_server(client_type=client_type)
    state, ctrl = server.state, server.controller
    return server, state, ctrl
