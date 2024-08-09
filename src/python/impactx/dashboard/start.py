"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

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
