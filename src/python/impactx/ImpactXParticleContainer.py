"""
This file is part of ImpactX

Copyright 2023 ImpactX contributors
Authors: Axel Huebl
License: BSD-3-Clause-LBNL
"""


def ix_pc_plot_mpl_phasespace(self, num_bins=50, root_rank=0):
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

    # Data Histogramming
    df = self.to_df(local=True)

    # calculate local histograms
    if df is None:
        xpx = np.zeros(
            (
                num_bins,
                num_bins,
            )
        )
        x_edges = np.linspace(rbc["x_min"] * m2mm, rbc["x_max"] * m2mm, num_bins + 1)
        px_edges = np.linspace(
            rbc["px_min"] * rad2mrad, rbc["px_max"] * rad2mrad, num_bins + 1
        )

        ypy = np.zeros(
            (
                num_bins,
                num_bins,
            )
        )
        y_edges = np.linspace(rbc["y_min"] * m2mm, rbc["y_max"] * m2mm, num_bins + 1)
        py_edges = np.linspace(
            rbc["py_min"] * rad2mrad, rbc["py_max"] * rad2mrad, num_bins + 1
        )

        tpt = np.zeros(
            (
                num_bins,
                num_bins,
            )
        )
        t_edges = np.linspace(rbc["t_min"] * m2mm, rbc["t_max"] * m2mm, num_bins + 1)
        pt_edges = np.linspace(
            rbc["pt_min"] * rad2mrad, rbc["pt_max"] * rad2mrad, num_bins + 1
        )
    else:
        # update for plot unit system
        # TODO: normalize to t/z to um and mc depending on s or t
        df.position_x = df.position_x.multiply(m2mm)
        df.position_y = df.position_y.multiply(m2mm)
        df.position_t = df.position_t.multiply(m2mm)

        df.momentum_x = df.momentum_x.multiply(rad2mrad)
        df.momentum_y = df.momentum_y.multiply(rad2mrad)
        df.momentum_t = df.momentum_t.multiply(rad2mrad)

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

    # MPI reduce
    #   nothing to do for non-MPI runs
    from inspect import getmodule

    ix = getmodule(self)
    if ix.Config.have_mpi:
        from mpi4py import MPI

        comm = MPI.COMM_WORLD  # TODO: get currently used ImpactX communicator here
        rank = comm.Get_rank()

        # MPI_Reduce the node-local histogram data
        combined_data = np.concatenate([xpx, ypy, tpt])
        summed_data = comm.reduce(
            combined_data,
            op=MPI.SUM,
            root=root_rank,
        )

        if rank != root_rank:
            return None

        [xpx, ypy, tpt] = np.split(
            summed_data,
            [
                len(xpx),
                len(xpx) + len(ypy),
            ],
        )

    # histograms per axis
    x = np.sum(xpx, axis=1)
    px = np.sum(xpx, axis=0)
    y = np.sum(ypy, axis=1)
    py = np.sum(ypy, axis=0)
    t = np.sum(tpt, axis=1)
    pt = np.sum(tpt, axis=0)

    # Matplotlib canvas: figure and plottable axes areas
    fig, axes = plt.subplots(1, 3, figsize=(16, 4), constrained_layout=True)
    (ax_xpx, ax_ypy, ax_tpt) = axes

    #   projected axes
    ax_x, ax_px = ax_xpx.twinx(), ax_xpx.twiny()
    ax_y, ax_py = ax_ypy.twinx(), ax_ypy.twiny()
    ax_t, ax_pt = ax_tpt.twinx(), ax_tpt.twiny()

    # Plotting
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
    cbar_ypy = plot_2d(ypy, y, py, y_edges, py_edges, ax_y, ax_py, ax_ypy)
    cbar_tpt = plot_2d(tpt, t, pt, t_edges, pt_edges, ax_t, ax_pt, ax_tpt)

    # Limits
    def set_limits(r, p, r_edges, p_edges, ax_r, ax_p, ax_rp):
        pad = 0.1
        len_r = r_edges[-1] - r_edges[0]
        len_p = p_edges[-1] - p_edges[0]
        ax_rp.set_xlim(r_edges[0] - len_r * pad, r_edges[-1] + len_r * pad)
        ax_rp.set_ylim(p_edges[0] - len_p * pad, p_edges[-1] + len_p * pad)

        # ensure zoom does not change value axis for projections
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
    set_limits(y, py, y_edges, py_edges, ax_y, ax_py, ax_ypy)
    set_limits(t, pt, t_edges, pt_edges, ax_t, ax_pt, ax_tpt)

    # Annotations
    fig.canvas.manager.set_window_title("Phase Space")
    ax_xpx.set_xlabel(r"$\Delta x$ [mm]")
    ax_xpx.set_ylabel(r"$\Delta p_x$ [mrad]")
    cbar_xpx.set_label(r"$Q$ [C/bin]")
    # ax_x.patch.set_alpha(0)
    ax_x.set_yticks([])
    ax_px.set_xticks([])

    ax_ypy.set_xlabel(r"$\Delta y$ [mm]")
    ax_ypy.set_ylabel(r"$\Delta p_y$ [mrad]")
    cbar_ypy.set_label(r"$Q$ [C/bin]")
    ax_y.set_yticks([])
    ax_py.set_xticks([])

    # TODO: update depending on s or t
    ax_tpt.set_xlabel(r"$\Delta ct$ [mm]")
    ax_tpt.set_ylabel(r"$\Delta p_t$ [$p_0\cdot c$]")
    cbar_tpt.set_label(r"$Q$ [C/bin]")
    ax_t.set_yticks([])
    ax_pt.set_xticks([])

    leg = ax_xpx.legend(
        title=r"$\epsilon_{n,x}=$"
        f"{Quantity(rbc['emittance_x'], 'm'):.3}\n"
        rf"$\sigma_x=${Quantity(rbc['sig_x'], 'm'):.3}"
        "\n"
        rf"$\beta_x=${Quantity(rbc['beta_x'], 'm'):.3}"
        "\n"
        rf"$\alpha_x=${rbc['alpha_x']:.3g}",
        loc="upper right",
        framealpha=0.8,
        handles=[],
    )
    leg._legend_box.sep = 0
    leg = ax_ypy.legend(
        title=r"$\epsilon_{n,y}=$"
        f"{Quantity(rbc['emittance_y'], 'm'):.3}\n"
        rf"$\sigma_y=${Quantity(rbc['sig_y'], 'm'):.3}"
        "\n"
        rf"$\beta_y=${Quantity(rbc['beta_y'], 'm'):.3}"
        "\n"
        rf"$\alpha_y=${rbc['alpha_y']:.3g}",
        loc="upper right",
        framealpha=0.8,
        handles=[],
    )
    leg._legend_box.sep = 0
    leg = ax_tpt.legend(
        title=r"$\epsilon_{n,t}=$"
        f"{Quantity(rbc['emittance_t'], 'm'):.3}\n"
        r"$\sigma_{ct}=$"
        f"{Quantity(rbc['sig_t'], 'm'):.3}\n"
        r"$\sigma_{pt}=$"
        f"{rbc['sig_pt']:.3g}",
        # TODO: I_peak, t_FWHM, ...
        loc="upper right",
        framealpha=0.8,
        handles=[],
    )
    leg._legend_box.sep = 0

    return fig


def register_ImpactXParticleContainer_extension(ixpc):
    """ImpactXParticleContainer helper methods"""
    # register member functions for ImpactXParticleContainer
    ixpc.plot_phasespace = ix_pc_plot_mpl_phasespace
