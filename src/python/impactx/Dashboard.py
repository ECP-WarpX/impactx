"""
This file is part of ImpactX

Copyright 2023 ImpactX contributors
Authors: Axel Huebl
License: BSD-3-Clause-LBNL
"""

import asyncio

import matplotlib.pyplot as plt
import numpy as np


def figure_size():
    if state.figure_size is None:
        return {}

    dpi = state.figure_size.get("dpi")
    rect = state.figure_size.get("size")
    w_inch = rect.get("width") / dpi
    h_inch = rect.get("height") / dpi

    return {
        "figsize": (w_inch, h_inch),
        "dpi": dpi,
    }


def FirstDemo():
    plt.close("all")
    fig, ax = plt.subplots(**figure_size())
    np.random.seed(0)
    ax.plot(
        np.random.normal(size=100), np.random.normal(size=100), "or", ms=10, alpha=0.3
    )
    ax.plot(
        np.random.normal(size=100), np.random.normal(size=100), "ob", ms=20, alpha=0.1
    )

    ax.set_xlabel("this is x")
    ax.set_ylabel("this is y")
    ax.set_title("Matplotlib Plot Rendered in D3!", size=14)
    ax.grid(color="lightgray", alpha=0.7)

    return fig


# Function to update histogram in a given Matplotlib widget
def update_histogram(widget, data):
    widget.figure.clear()
    ax = widget.figure.add_subplot(111)
    ax.hist(data.flatten(), bins=50)
    widget.draw_idle()
    # widget.update(fig1)


def create_app(server):
    from trame.ui.vuetify import SinglePageLayout
    from trame.widgets import matplotlib
    from trame.widgets.vuetify import VContainer, VSelect, VSpacer

    layout = SinglePageLayout(server)
    layout.title.set_text("Hello trame")
    # plot = matplotlib.Figure()

    return layout

    from trame.widgets.trame import SizeObserver

    ctrl = server.controller

    # create GUI layout
    with SinglePageLayout(server) as layout:
        layout.title.set_text("trame ❤️ matplotlib")

        with layout.toolbar:
            VSpacer()
            VSelect(
                v_model=("active_figure", "FirstDemo"),
                items=(
                    "figures",
                    [
                        {"text": "First Demo", "value": "FirstDemo"},
                        #                        {"text": "Multi Lines", "value": "MultiLines"},
                        #                        {"text": "Dots and Points", "value": "DotsandPoints"},
                        #                        {"text": "Moving Window Average", "value": "MovingWindowAverage"},
                        #                        {"text": "Subplots", "value": "Subplots"},
                    ],
                ),
                hide_details=True,
                dense=True,
            )

        with layout.content:
            with VContainer(fluid=True, classes="fill-height pa-0 ma-0"):
                with SizeObserver("figure_size"):
                    html_figure = matplotlib.Figure(style="position: absolute")
                    ctrl.update_figure = html_figure.update


def update_dashboard(sim):
    update_histogram(widget, data)


# -----------------------------------------------------------------------------
# Life Cycle events
# -----------------------------------------------------------------------------


def server_ready(**state):
    import json

    print("on_server_ready")
    print("  => current state:")
    print(json.dumps(state, indent=2))
    print("-" * 60)


def client_connected():
    print("on_client_connected")


def client_unmounted():
    print("on_client_unmounted")


def client_exited():
    print("on_client_exited")


def server_exited(**state):
    import json

    print("on_server_exited")
    print("  => current state:")
    print(json.dumps(state, indent=2))
    print("-" * 60)


async def start(server):
    await server.start(
        exec_mode="coroutine",
        open_browser=True,
    )


async def init_dashboard(sim):
    from trame.app import get_server

    sim.trame_server = get_server()
    sim.trame_server.client_type = "vue2"  # Until Matplotlib is ported to vue3

    # -----------------------------------------------------------------------------
    # Life Cycle registration
    # -----------------------------------------------------------------------------
    ctrl = sim.trame_server.controller

    ctrl.on_server_ready.add(server_ready)
    ctrl.on_client_connected.add(client_connected)
    ctrl.on_client_unmounted.add(client_unmounted)
    ctrl.on_client_exited.add(client_exited)
    ctrl.on_server_exited.add(server_exited)

    layout = create_app(sim.trame_server)

    await asyncio.create_task(start(sim.trame_server))


def register_dashboard(sim):
    """Simulation helper methods for the Dashboard"""

    # register member functions for ImpactX simulation class
    sim.dashboard = init_dashboard
    sim.update_dashboard = update_dashboard
