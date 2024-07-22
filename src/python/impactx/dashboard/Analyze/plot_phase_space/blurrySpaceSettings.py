"""
This file is part of ImpactX

Copyright 2023 ImpactX contributors
Authors: Axel Huebl
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
    from quantiphy import Quantity

    # Beam Characteristics
    rbc = self.reduced_beam_characteristics()

    # update for plot unit system
    m2mm = 1.0e3
    rad2mrad = 1.0e3

    # Data Extraction
    df = self.to_df(local=True)

    # Matplotlib canvas: figure and plottable axes areas
    fig, axes = plt.subplots(1, 3, figsize=(16, 4), constrained_layout=True)
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

        xpx, x_edges, px_edges = np.histogram2d(
            df["position_x"],
            df["momentum_x"],
            bins=num_bins,
            range=[
                [rbc["x_min"] * m2mm, rbc["x_max"] * m2mm],
                [rbc["px_min"] * rad2mrad, rbc["px_max"] * rad2mrad],
            ],
        )

        ypy, y_edges, py_edges = np.histogram2d(
            df["position_y"],
            df["momentum_y"],
            bins=num_bins,
            range=[
                [rbc["y_min"] * m2mm, rbc["y_max"] * m2mm],
                [rbc["py_min"] * rad2mrad, rbc["py_max"] * rad2mrad],
            ],
        )

        tpt, t_edges, pt_edges = np.histogram2d(
            df["position_t"],
            df["momentum_t"],
            bins=num_bins,
            range=[
                [rbc["t_min"] * m2mm, rbc["t_max"] * m2mm],
                [rbc["pt_min"] * rad2mrad, rbc["pt_max"] * rad2mrad],
            ],
        )

        def plot_2d(hist, r_edges, p_edges, ax_rp):
            hist = np.ma.masked_where(hist == 0, hist)
            im = ax_rp.imshow(
                hist.T,
                origin="lower",
                aspect="auto",
                extent=[r_edges[0], r_edges[-1], p_edges[0], p_edges[-1]],
            )
            return fig.colorbar(im, ax=ax_rp)

        cbar_xpx = plot_2d(xpx, x_edges, px_edges, ax_xpx)
        cbar_ypy = plot_2d(ypy, y_edges, py_edges, ax_ypy)
        cbar_tpt = plot_2d(tpt, t_edges, pt_edges, ax_tpt)
    else:
        ax_xpx.text(0.5, 0.5, 'No data available', horizontalalignment='center', verticalalignment='center')
        ax_ypy.text(0.5, 0.5, 'No data available', horizontalalignment='center', verticalalignment='center')
        ax_tpt.text(0.5, 0.5, 'No data available', horizontalalignment='center', verticalalignment='center')

    # Annotations
    fig.canvas.manager.set_window_title("Phase Space")
    ax_xpx.set_xlabel("Δ x [mm]")
    ax_xpx.set_ylabel("Δ p_x [mrad]")
    ax_xpx.set_title("Longitudinal Phase Space")
    cbar_xpx.set_label("Charge [C/bin]")

    ax_ypy.set_xlabel("Δ y [mm]")
    ax_ypy.set_ylabel("Δ p_y [mrad]")
    ax_ypy.set_title("Transverse Phase Space (y)")
    cbar_ypy.set_label("Charge [C/bin]")

    ax_tpt.set_xlabel("Δ ct [mm]")
    ax_tpt.set_ylabel("Δ p_t [p_0 · c]")
    ax_tpt.set_title("Transverse Phase Space (t)")
    cbar_tpt.set_label("Charge [C/bin]")

    # Adding legends
    leg = ax_xpx.legend(
        title=(
            "ε_n,x = {}\n"
            "σ_x = {}\n"
            "β_x = {}\n"
            "α_x = {:.3g}"
        ).format(
            Quantity(rbc['emittance_x'], 'm').render(prec=3),
            Quantity(rbc['sig_x'], 'm').render(prec=3),
            Quantity(rbc['beta_x'], 'm').render(prec=3),
            rbc['alpha_x']
        ),
        loc="upper right",
        framealpha=0.8,
        handles=[]
    )
    leg._legend_box.sep = 0

    leg = ax_ypy.legend(
        title=(
            "ε_n,y = {}\n"
            "σ_y = {}\n"
            "β_y = {}\n"
            "α_y = {:.3g}"
        ).format(
            Quantity(rbc['emittance_y'], 'm').render(prec=3),
            Quantity(rbc['sig_y'], 'm').render(prec=3),
            Quantity(rbc['beta_y'], 'm').render(prec=3),
            rbc['alpha_y']
        ),
        loc="upper right",
        framealpha=0.8,
        handles=[]
    )
    leg._legend_box.sep = 0

    leg = ax_tpt.legend(
        title=(
            "ε_n,t = {}\n"
            "σ_ct = {}\n"
            "σ_pt = {:.3g}"
        ).format(
            Quantity(rbc['emittance_t'], 'm').render(prec=3),
            Quantity(rbc['sig_t'], 'm').render(prec=3),
            rbc['sig_pt']
        ),
        loc="upper right",
        framealpha=0.8,
        handles=[]
    )
    leg._legend_box.sep = 0

    return fig
