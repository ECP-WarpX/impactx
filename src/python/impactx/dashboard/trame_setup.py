"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from trame.app import get_server


def setup_server(client_type="vue2"):
    server = get_server(client_type=client_type)
    state, ctrl = server.state, server.controller
    return server, state, ctrl
