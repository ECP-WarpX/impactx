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
    from scipy.stats import gaussian_kde

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

        def scatter_density_plot(ax, x, y, xlabel, ylabel, title):
            xy = np.vstack([x, y])
            z = gaussian_kde(xy)(xy)
            scatter = ax.scatter(x, y, c=z, cmap='viridis', alpha=0.5)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(title)
            return scatter

        scatter_xpx = scatter_density_plot(
            ax_xpx, df["position_x"], df["momentum_x"],
            "Δ x [mm]", "Δ p_x [mrad]", "Longitudinal Phase Space"
        )
        scatter_ypy = scatter_density_plot(
            ax_ypy, df["position_y"], df["momentum_y"],
            "Δ y [mm]", "Δ p_y [mrad]", "Transverse Phase Space (y)"
        )
        scatter_tpt = scatter_density_plot(
            ax_tpt, df["position_t"], df["momentum_t"],
            "Δ ct [mm]", "Δ p_t [p_0 · c]", "Transverse Phase Space (t)"
        )

        fig.colorbar(scatter_xpx, ax=ax_xpx, label="Particle Density")
        fig.colorbar(scatter_ypy, ax=ax_ypy, label="Particle Density")
        fig.colorbar(scatter_tpt, ax=ax_tpt, label="Particle Density")

    else:
        ax_xpx.text(0.5, 0.5, 'No data available', horizontalalignment='center', verticalalignment='center')
        ax_ypy.text(0.5, 0.5, 'No data available', horizontalalignment='center', verticalalignment='center')
        ax_tpt.text(0.5, 0.5, 'No data available', horizontalalignment='center', verticalalignment='center')

    # Annotations
    fig.canvas.manager.set_window_title("Phase Space")

    # # Adding legends
    # def add_legend(ax, title):
    #     leg = ax.legend(
    #         title=title,
    #         loc="upper right",
    #         framealpha=0.8,
    #         handles=[]
    #     )
    #     leg._legend_box.sep = 0

    # add_legend(
    #     ax_xpx, (
    #         "ε_n,x = {}\n"
    #         "σ_x = {}\n"
    #         "β_x = {}\n"
    #         "α_x = {:.3g}"
    #     ).format(
    #         Quantity(rbc['emittance_x'], 'm').render(prec=3),
    #         Quantity(rbc['sig_x'], 'm').render(prec=3),
    #         Quantity(rbc['beta_x'], 'm').render(prec=3),
    #         rbc['alpha_x']
    #     )
    # )
    # add_legend(
    #     ax_ypy, (
    #         "ε_n,y = {}\n"
    #         "σ_y = {}\n"
    #         "β_y = {}\n"
    #         "α_y = {:.3g}"
    #     ).format(
    #         Quantity(rbc['emittance_y'], 'm').render(prec=3),
    #         Quantity(rbc['sig_y'], 'm').render(prec=3),
    #         Quantity(rbc['beta_y'], 'm').render(prec=3),
    #         rbc['alpha_y']
    #     )
    # )
    # add_legend(
    #     ax_tpt, (
    #         "ε_n,t = {}\n"
    #         "σ_ct = {}\n"
    #         "σ_pt = {:.3g}"
    #     ).format(
    #         Quantity(rbc['emittance_t'], 'm').render(prec=3),
    #         Quantity(rbc['sig_t'], 'm').render(prec=3),
    #         rbc['sig_pt']
    #     )
    # )

    return fig
