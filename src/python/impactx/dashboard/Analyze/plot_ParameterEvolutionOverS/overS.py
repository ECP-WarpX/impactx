"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

import plotly.graph_objects as go


def line_plot_1d(selected_headers, filtered_data):
    """
    Generates a 1D line plot using Plotly based on selected headers and filtered data.
    """

    x_axis = selected_headers[0] if len(selected_headers) > 1 else None
    y_axis = selected_headers[1:] if len(selected_headers) > 2 else None

    x = [row[x_axis] for row in filtered_data] if x_axis else []

    figure_data = []
    if y_axis:
        for column in y_axis:
            y = [row[column] for row in filtered_data]
            trace = go.Scatter(
                x=x,
                y=y,
                mode="lines+markers",
                name=column,
                line=dict(width=2),
                marker=dict(size=8),
            )
            figure_data.append(trace)

    return go.Figure(
        data=figure_data,
        layout=go.Layout(
            title="Plot Over S",
            xaxis=dict(title="s"),
            yaxis=dict(title=""),
            margin=dict(l=20, r=20, t=25, b=30),
        ),
    )
