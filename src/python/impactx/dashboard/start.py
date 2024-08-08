from .trame_setup import setup_server

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def main():
    """
    Launches Trame application server
    """
    server.start()
    return 0
