"""
This file is part of ImpactX

Copyright 2023 ImpactX contributors
Authors: Axel Huebl
License: BSD-3-Clause-LBNL
"""

from contextlib import contextmanager

from ..impactx_pybind import CoordSystem, coordinate_transformation


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
    if self.coord_system != CoordSystem.s:
        raise RuntimeError(
            "plot_phasespace: the particles must be at fixed s, "
            f"but are at {self.coord_system}."
        )

    import matplotlib.pyplot as plt
    import numpy as np
    from quantiphy import Quantity

    # Beam Characteristics
    rbc = self.beam_moments()

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
        x_edges = np.linspace(rbc["min_x"] * m2mm, rbc["max_x"] * m2mm, num_bins + 1)
        px_edges = np.linspace(
            rbc["min_px"] * rad2mrad, rbc["max_px"] * rad2mrad, num_bins + 1
        )

        ypy = np.zeros(
            (
                num_bins,
                num_bins,
            )
        )
        y_edges = np.linspace(rbc["min_y"] * m2mm, rbc["max_y"] * m2mm, num_bins + 1)
        py_edges = np.linspace(
            rbc["min_py"] * rad2mrad, rbc["max_py"] * rad2mrad, num_bins + 1
        )

        tpt = np.zeros(
            (
                num_bins,
                num_bins,
            )
        )
        t_edges = np.linspace(rbc["min_t"] * m2mm, rbc["max_t"] * m2mm, num_bins + 1)
        pt_edges = np.linspace(
            rbc["min_pt"] * rad2mrad, rbc["max_pt"] * rad2mrad, num_bins + 1
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
                [rbc["min_x"] * m2mm, rbc["max_x"] * m2mm],
                [rbc["min_px"] * rad2mrad, rbc["max_px"] * rad2mrad],
            ],
        )

        ypy, y_edges, py_edges = np.histogram2d(
            df["position_y"],
            df["momentum_y"],
            bins=num_bins,
            range=[
                [rbc["min_y"] * m2mm, rbc["max_y"] * m2mm],
                [rbc["min_py"] * rad2mrad, rbc["max_py"] * rad2mrad],
            ],
        )

        tpt, t_edges, pt_edges = np.histogram2d(
            df["position_t"],
            df["momentum_t"],
            bins=num_bins,
            range=[
                [rbc["min_t"] * m2mm, rbc["max_t"] * m2mm],
                [rbc["min_pt"] * rad2mrad, rbc["max_pt"] * rad2mrad],
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


def ix_beam_moments_history(self):
    """
    Return the history of the beam as calculated by the reduced beam characteristics on every step.
    """
    import pandas as pd

    return pd.DataFrame(self.beam_moments_history_list())


def _as_real_device_vector(arr):
    """Convert an input into the AMReX real ``DeviceVector`` that
    ``AddNParticles`` expects, copying the data as needed.

    Accepts ``None`` (passed through), a NumPy or CuPy array (or array-like),
    or any pyAMReX ``PODVector``. A ``DeviceVector_real`` is passed through
    unchanged; any other allocator is copied via ``to_device`` (which picks
    the AMReX copy direction across memory spaces).
    """
    if arr is None:
        return None

    import amrex.space3d as amr

    # already the exact Gpu::DeviceVector type that AddNParticles expects
    if isinstance(arr, amr.DeviceVector_real):
        return arr
    # any other PODVector allocator: copy into a DeviceVector (CuPy-free)
    if type(arr).__name__.startswith("PODVector_"):
        return arr.to_device()
    # CuPy (and other device) arrays expose the CUDA array interface
    if hasattr(arr, "__cuda_array_interface__"):
        return amr.DeviceVector_real.from_cupy(arr)
    # NumPy arrays and array-likes (lists, tuples, ...)
    return amr.DeviceVector_real.from_numpy(arr)


def ix_pc_add_n_particles(
    self, x, y, t, px, py, pt, qm, bunch_charge=None, w=None, sx=None, sy=None, sz=None
):
    """
    Add new particles to the container for fixed s.

    The coordinate and weight arguments accept NumPy or CuPy arrays (or
    array-likes), as well as pyAMReX ``PODVector`` objects. Inputs are copied
    into device-compatible PODVectors as needed.

    Either the total charge (``bunch_charge``) or the weight of each particle
    (``w``) must be provided.

    Note: This can only be used *after* the grids have been created, i.e. after
    ``ImpactX.init_grids`` has been called.

    Parameters
    ----------
    x, y, t, px, py, pt : array_like
        Particle positions (x, y, time-of-flight c*t) and momenta.
    qm : float
        Charge over mass in 1/eV.
    bunch_charge : float, optional
        Total charge within a bunch in C.
    w : array_like, optional
        Weight of each particle: how many real particles to represent.
    sx, sy, sz : array_like, optional
        Spin components in x, y, z.
    """
    return self._add_n_particles(
        _as_real_device_vector(x),
        _as_real_device_vector(y),
        _as_real_device_vector(t),
        _as_real_device_vector(px),
        _as_real_device_vector(py),
        _as_real_device_vector(pt),
        qm,
        bunch_charge,
        _as_real_device_vector(w),
        _as_real_device_vector(sx),
        _as_real_device_vector(sy),
        _as_real_device_vector(sz),
    )


@contextmanager
def ix_pc_at_fixed_t(self):
    """
    Temporarily represent the beam at fixed t, e.g., to exchange particles with a code
    that uses time as the independent variable.

    On entry, the particle coordinates are transformed from fixed s to fixed t.
    When the ``with`` block is left, normally or through an exception, they are
    transformed back to fixed s.

    Inside the block, the particle arrays hold ``x, y, z, px, py, pz`` instead of
    ``x, y, t, px, py, pt``: ``z`` is the longitudinal position relative to the
    reference particle, in meters, ``px, py`` are the transverse momenta and ``pz``
    is the deviation from the reference momentum, all normalized by the reference
    momentum.
    The arrays are named accordingly: ``position_z`` and ``momentum_z`` replace
    ``position_t`` and ``momentum_t``.
    The arguments of ``add_n_particles`` keep their names ``t`` and ``pt``; inside the
    block, they take ``z`` and ``pz``.

    Both transformations take the design energy from ``self.ref.pt`` at the moment
    they run: the one on entry uses the reference particle as it is when the block
    starts, the one on exit uses the reference particle as it is when the block ends.
    If the block changes the reference energy, it must also write the particle
    coordinates relative to the new reference particle before the block ends.
    Only ``ref.pt`` is read, so update it (e.g., with ``ref.set_kin_energy_MeV``)
    rather than ``ref.pz`` alone.

    Parameters
    ----------
    self : ImpactXParticleContainer
        The particle container; its particles must be at fixed s.

    Yields
    ------
    ImpactXParticleContainer
        The same particle container, now at fixed t.

    Raises
    ------
    RuntimeError
        If the particles are not at fixed s on entry, or no longer at fixed t
        when the block ends without an exception.

    Examples
    --------
    >>> with sim.beam.at_fixed_t() as beam:
    ...     beam.add_n_particles(x, y, z, px, py, pz, qm, bunch_charge=charge_C)
    """
    if self.coord_system != CoordSystem.s:
        raise RuntimeError(
            "at_fixed_t: the particles must be at fixed s when entering the block, "
            f"but are at {self.coord_system}."
        )

    coordinate_transformation(self, CoordSystem.t)
    try:
        yield self
    except BaseException:
        # keep the original exception; only restore fixed s if it is still possible
        if self.coord_system == CoordSystem.t:
            coordinate_transformation(self, CoordSystem.s)
        raise

    if self.coord_system != CoordSystem.t:
        raise RuntimeError(
            "at_fixed_t: the particles were transformed out of fixed t inside the "
            "block. Do not call coordinate_transformation inside the block: leaving "
            "the block transforms the particles back to fixed s."
        )
    coordinate_transformation(self, CoordSystem.s)


def register_ImpactXParticleContainer_extension(ixpc):
    """ImpactXParticleContainer helper methods"""
    # register member functions for ImpactXParticleContainer
    ixpc.plot_phasespace = ix_pc_plot_mpl_phasespace
    ixpc.beam_moments_history = ix_beam_moments_history
    ixpc.add_n_particles = ix_pc_add_n_particles
    ixpc.at_fixed_t = ix_pc_at_fixed_t
