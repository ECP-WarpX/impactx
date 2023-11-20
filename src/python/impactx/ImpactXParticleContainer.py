"""
This file is part of ImpactX

Copyright 2023 ImpactX contributors
Authors: Axel Huebl
License: BSD-3-Clause-LBNL
"""


def ix_pc_to_df(self, local=True, comm=None, root_rank=0):
    """
    Copy all particles into a pandas.DataFrame

    Parameters
    ----------
    self : amrex.ParticleContainer_*
        A ParticleContainer class in pyAMReX
    local : bool
        MPI-local particles
    comm : MPI Communicator
        if local is False, this defaults to mpi4py.MPI.COMM_WORLD
    root_rank : MPI root rank to gather to
        if local is False, this defaults to 0

    Returns
    -------
    A concatenated pandas.DataFrame with particles from all levels.

    Returns None if no particles were found.
    If local=False, then all ranks but the root_rank will return None.
    """
    df = super(type(self), self).to_df(local=local, comm=comm, root_rank=root_rank)

    # rename columns according to our attribute names
    if df is not None:
        # todo: check if currently in fixed s or fixed t and pick name accordingly

        names = []
        for n in self.RealAoS.names_s:
            names.append(n)
        names.append("cpuid")
        for n in self.RealSoA.names_s:
            names.append(n)

        df.columns.values[0 : len(names)] = names

        # todo: also rename runtime attributes (e.g., "s_lost")
        # https://github.com/ECP-WarpX/impactx/pull/398

    return df


def ix_pc_plot_mpl_phasespace(self, num_bins=50):
    """
    ...
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import numpy as np

    # Matplotlib canvas: figure and plottable axes areas
    fig, axes = plt.subplots(1, 3, figsize=(16, 4), constrained_layout=True)
    (ax_xpx, ax_ypy, ax_tpt) = axes

    #   projected axes
    ax_x, ax_px = ax_xpx.twinx(), ax_xpx.twiny()
    ax_y, ax_py = ax_ypy.twinx(), ax_ypy.twiny()
    ax_t, ax_pt = ax_tpt.twinx(), ax_tpt.twiny()

    # Beam Characteristics
    rbc = self.reduced_beam_characteristics()

    # Data Histogramming
    # TODO: it will be actually scalable if we first histogram locally,
    #       then reduce the histrograms with MPI
    # TODO: normalize to um and mc depending on s or t
    df = self.to_df(local=False)

    # update for plot unit system
    df.position_x = df.position_x.multiply(1e6)
    df.position_y = df.position_y.multiply(1e6)

    xpx, x_edges, px_edges = np.histogram2d(
        df["position_x"],
        df["momentum_x"],
        bins=num_bins,
    )
    x = np.sum(xpx, axis=1)
    px = np.sum(xpx, axis=0)

    ypy, y_edges, py_edges = np.histogram2d(
        df["position_y"],
        df["momentum_y"],
        bins=num_bins,
    )
    y = np.sum(ypy, axis=1)
    py = np.sum(ypy, axis=0)

    tpt, t_edges, pt_edges = np.histogram2d(
        df["position_t"],
        df["momentum_t"],
        bins=num_bins,
    )
    t = np.sum(tpt, axis=1)
    pt = np.sum(tpt, axis=0)

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
    ax_xpx.set_xlabel(r"$\Delta x$ [$\mu$m]")
    ax_xpx.set_ylabel(r"$\Delta p_x$ [mc]")
    cbar_xpx.set_label(r"$Q$ [C/bin]")
    # ax_x.patch.set_alpha(0)
    ax_x.set_yticks([])
    ax_px.set_xticks([])

    ax_ypy.set_xlabel(r"$\Delta y$ [$\mu$m]")
    ax_ypy.set_ylabel(r"$\Delta p_y$ [mc]")
    cbar_ypy.set_label(r"$Q$ [C/bin]")
    ax_y.set_yticks([])
    ax_py.set_xticks([])

    # TODO: update depending on s or t
    ax_tpt.set_xlabel(r"$\Delta t$ [...]")
    ax_tpt.set_ylabel(r"$\Delta p_t$ [...]")
    cbar_tpt.set_label(r"$Q$ [C/bin]")
    ax_t.set_yticks([])
    ax_pt.set_xticks([])

    # TODO: write an auto-formatter that picks m, mu, n, p, f, k, M, G
    #       automatically
    ax_xpx.legend(
        title=r"$\epsilon_{n,x}=$"
        f"{rbc['emittance_x']*1e6:.3f} µm"
        "\n"
        rf"$\sigma_x=${rbc['sig_x']*1e6:.3f} µm"
        "\n"
        rf"$\beta_x=${rbc['beta_x']*1e3:.3f} mm"
        "\n"
        rf"$\alpha_x=${rbc['alpha_x']:.3f}",
        loc="upper right",
        framealpha=0.8,
        handles=[],
    )
    ax_ypy.legend(
        title=r"$\epsilon_{n,y}=$"
        f"{rbc['emittance_y']*1e6:.3f} µm"
        "\n"
        rf"$\sigma_x=${rbc['sig_y']*1e6:.3f} µm"
        "\n"
        rf"$\beta_x=${rbc['beta_y']*1e3:.3f} mm"
        "\n"
        rf"$\alpha_x=${rbc['alpha_y']:.3f}",
        loc="upper right",
        framealpha=0.8,
        handles=[],
    )
    ax_tpt.legend(
        title=r"$\epsilon_{n,t}=$"
        f"{rbc['emittance_t']*1e6:.3f}"
        r" MeV$\cdot$s"
        "\n"
        rf"$\sigma_t=${rbc['sig_t']*1e15:.3f} fs",
        # TODO: sigma_pz, I_peak, t_FWHM, ...
        loc="upper right",
        framealpha=0.8,
        handles=[],
    )

    return fig


def register_ImpactXParticleContainer_extension(ixpc):
    """ImpactXParticleContainer helper methods"""

    # register member functions for ImpactXParticleContainer
    ixpc.to_df = ix_pc_to_df
    ixpc.plot_phasespace = ix_pc_plot_mpl_phasespace
