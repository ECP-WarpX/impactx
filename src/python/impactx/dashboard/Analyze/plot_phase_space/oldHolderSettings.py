### Only the single plot shows deltaX_delxaP_X

def adjusted_settings_plot(pc, num_bins=50, root_rank=0):
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats import gaussian_kde

    # Data Histogramming
    df = pc.to_df(local=True)

    if df is None:
        print("DataFrame is None")
        return None

    # Update for plot unit system
    m2mm = 1.0e3
    rad2mrad = 1.0e3

    # Convert positions to mm and momenta to mrad
    df.position_x = df.position_x.multiply(m2mm)
    df.momentum_x = df.momentum_x.multiply(rad2mrad)

    # Calculate the density of points using a 2D histogram
    xy = np.vstack([df.position_x, df.momentum_x])
    z = gaussian_kde(xy)(xy)

    # Simple scatter plot of position_x and momentum_x with density coloring
    fig, ax = plt.subplots(figsize=(8, 4), constrained_layout=True)
    scatter = ax.scatter(df.position_x, df.momentum_x, c=z, cmap='viridis', alpha=0.5)
    ax.set_xlabel("Delta x [mm]")
    ax.set_ylabel("Delta p_x [mrad]")

    # Adding a color bar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("Particle Density")

    # Adding a legend for the scatter plot (optional if needed for more information)
    # We can create custom legend if needed, but color bar usually serves the purpose for density.

    """
    # Beam Characteristics
    rbc = pc.reduced_beam_characteristics()

    # calculate local histograms
    xpx, x_edges, px_edges = np.histogram2d(
        df["position_x"],
        df["momentum_x"],
        bins=num_bins,
        range=[
            [rbc["x_min"] * m2mm, rbc["x_max"] * m2mm],
            [rbc["px_min"] * rad2mrad, rbc["px_max"] * rad2mrad],
        ],
    )

    # histograms per axis
    x = np.sum(xpx, axis=1)
    px = np.sum(xpx, axis=0)

    # Check histogram contents
    print(f"xpx:\n{xpx}\nx:\n{x}\npx:\n{px}\n")

    # Matplotlib canvas: figure and plottable axes areas
    fig, ax_xpx = plt.subplots(figsize=(8, 4), constrained_layout=True)

    #   projected axes
    ax_x, ax_px = ax_xpx.twinx(), ax_xpx.twiny()

    def plot_2d(hist, r, p, r_edges, p_edges, ax_r, ax_p, ax_rp):
        hist = np.ma.masked_where(hist == 0, hist)
        im = ax_rp.imshow(
            hist.T,
            origin="lower",
            aspect="auto",
            extent=[r_edges[0], r_edges[-1], p_edges[0], p_edges[-1]],
        )
        cbar = fig.colorbar(im, ax=ax_rp)

        r_mids = (r_edges[:-1] + r_edges[1:]) / 2
        p_mids = (p_edges[:-1] + p_edges[1:]) / 2
        ax_r.plot(r_mids, r, c="w", lw=0.8, alpha=0.7)
        ax_r.plot(r_mids, r, c="k", lw=0.5, alpha=0.7)
        ax_r.fill_between(r_mids, r, facecolor="k", alpha=0.2)
        ax_p.plot(p, p_mids, c="w", lw=0.8, alpha=0.7)
        ax_p.plot(p, p_mids, c="k", lw=0.5, alpha=0.7)
        ax_p.fill_betweenx(p_mids, p, facecolor="k", alpha=0.2)

        return cbar

    cbar_xpx = plot_2d(xpx, x, px, x_edges, px_edges, ax_x, ax_px, ax_xpx)

    def set_limits(r, p, r_edges, p_edges, ax_r, ax_p, ax_rp):
        pad = 0.1
        len_r = r_edges[-1] - r_edges[0]
        len_p = p_edges[-1] - p_edges[0]
        ax_rp.set_xlim(r_edges[0] - len_r * pad, r_edges[-1] + len_r * pad)
        ax_rp.set_ylim(p_edges[0] - len_p * pad, p_edges[-1] + len_p * pad)

        def on_xlims_change(axes):
            if not axes.xlim_reset_in_progress:
                pad = 6.0
                axes.xlim_reset_in_progress = True
                axes.set_xlim(0, np.max(p) * pad)
                axes.xlim_reset_in_progress = False

        ax_p.xlim_reset_in_progress = False
        ax_p.callbacks.connect("xlim_changed", on_xlims_change)
        on_xlims_change(ax_p)

        def on_ylims_change(axes):
            if not axes.ylim_reset_in_progress:
                pad = 6.0
                axes.ylim_reset_in_progress = True
                axes.set_ylim(0, np.max(r) * pad)
                axes.ylim_reset_in_progress = False

        ax_r.ylim_reset_in_progress = False
        ax_r.callbacks.connect("ylim_changed", on_ylims_change)
        on_ylims_change(ax_r)

    set_limits(x, px, x_edges, px_edges, ax_x, ax_px, ax_xpx)

    ax_xpx.set_xlabel("Delta x [mm]")
    ax_xpx.set_ylabel("Delta p_x [mrad]")
    cbar_xpx.set_label("Q [C/bin]")
    ax_x.set_yticks([])
    ax_px.set_xticks([])

    return fig
    """

    return fig
