"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""


def adjusted_settings_plot(self, num_bins=50, root_rank=0):
    """
    Plot the longitudinal and transverse phase space projections with matplotlib.

    Parameters
    ----------
    self : ImpactXParticleContainer_*
        The particle container class in ImpactX
    num_bins : int, default=50
        The number of bins for spatial and momentum directions per plot axis.
    root_rank : int, default=0
        MPI root rank to reduce to in parallel runs.

    Returns
    -------
    A matplotlib figure with containing the plot.
    For MPI-parallel ranks, the figure is only created on the root_rank.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats import gaussian_kde

    # Beam Characteristics
    # rbc = self.reduced_beam_characteristics()
    # update for plot unit system
    m2mm = 1.0e3
    rad2mrad = 1.0e3

    # Data Extraction
    df = self.to_df(local=True)

    # Matplotlib canvas: figure and plottable axes areas
    fig, axes = plt.subplots(1, 3, figsize=(12, 3))
    (ax_xpx, ax_ypy, ax_tpt) = axes

    # Plotting data points if df is not None
    if df is not None:
        # update for plot unit system
        df["position_x"] = df["position_x"].multiply(m2mm)
        df["position_y"] = df["position_y"].multiply(m2mm)
        df["position_t"] = df["position_t"].multiply(m2mm)
        df["momentum_x"] = df["momentum_x"].multiply(rad2mrad)
        df["momentum_y"] = df["momentum_y"].multiply(rad2mrad)
        df["momentum_t"] = df["momentum_t"].multiply(rad2mrad)

        def scatter_density_plot(ax, x, y, xlabel, ylabel, title):
            xy = np.vstack([x, y])
            z = gaussian_kde(xy)(xy)
            scatter = ax.scatter(x, y, c=z, cmap="viridis", alpha=0.5)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(title)
            ax.axis("equal")
            return scatter

        scatter_xpx = scatter_density_plot(
            ax_xpx,
            df["position_x"],
            df["momentum_x"],
            "Δ x [mm]",
            "Δ p_x [mrad]",
            "",
        )
        scatter_ypy = scatter_density_plot(
            ax_ypy,
            df["position_y"],
            df["momentum_y"],
            "Δ y [mm]",
            "Δ p_y [mrad]",
            "",
        )
        scatter_tpt = scatter_density_plot(
            ax_tpt,
            df["position_t"],
            df["momentum_t"],
            "Δ ct [mm]",
            "Delta p_t [p_0 . c]",
            "",
        )

        fig.colorbar(scatter_xpx, ax=ax_xpx, fraction=0.046, pad=0.04)
        fig.colorbar(scatter_ypy, ax=ax_ypy, fraction=0.046, pad=0.04)
        fig.colorbar(scatter_tpt, ax=ax_tpt, fraction=0.046, pad=0.04)

        fig.tight_layout()

    else:
        ax_xpx.text(
            0.5,
            0.5,
            "No data available",
            horizontalalignment="center",
            verticalalignment="center",
        )
        ax_ypy.text(
            0.5,
            0.5,
            "No data available",
            horizontalalignment="center",
            verticalalignment="center",
        )
        ax_tpt.text(
            0.5,
            0.5,
            "No data available",
            horizontalalignment="center",
            verticalalignment="center",
        )

    fig.canvas.manager.set_window_title("Phase Space")

    return fig
