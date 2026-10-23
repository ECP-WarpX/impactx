"""
This file is part of ImpactX

Copyright 2023 ImpactX contributors
Authors: Axel Huebl
License: BSD-3-Clause-LBNL
"""

from __future__ import annotations

from contextlib import contextmanager

from impactx.impactx_pybind import CoordSystem, coordinate_transformation

__all__: list[str] = [
    "CoordSystem",
    "contextmanager",
    "coordinate_transformation",
    "ix_beam_moments_history",
    "ix_pc_add_n_particles",
    "ix_pc_at_fixed_t",
    "ix_pc_plot_mpl_phasespace",
    "register_ImpactXParticleContainer_extension",
]

def _as_real_device_vector(arr):
    """
    Convert an input into the AMReX real ``DeviceVector`` that
    ``AddNParticles`` expects, copying the data as needed.

    Accepts ``None`` (passed through), a NumPy or CuPy array (or array-like),
    or any pyAMReX ``PODVector``. A ``DeviceVector_real`` is passed through
    unchanged; any other allocator is copied via ``to_device`` (which picks
    the AMReX copy direction across memory spaces).
    """

def ix_beam_moments_history(self):
    """
    Return the history of the beam as calculated by the reduced beam characteristics on every step.
    """

def ix_pc_add_n_particles(
    self, x, y, t, px, py, pt, qm, bunch_charge=None, w=None, sx=None, sy=None, sz=None
):
    """
    Add new particles to the container.

    The coordinate and weight arguments accept NumPy or CuPy arrays (or
    array-likes), as well as pyAMReX ``PODVector`` objects. Inputs are copied
    into device-compatible PODVectors as needed.

    Either the total charge (``bunch_charge``) or the weight of each particle
    (``w``) must be provided.
    With ``bunch_charge``, each particle represents the same number of real
    particles with the charge of the reference particle, which must be set up
    before.

    Note: This can only be used *after* the grids have been created, i.e. after
    ``ImpactX.init_grids`` has been called.

    Parameters
    ----------
    x, y, t, px, py, pt : array_like
        Particle coordinates, relative to the reference particle.
        ``x`` and ``y`` are positions in m, ``px`` and ``py`` are momenta
        normalized by the reference momentum (unitless).
        At fixed s (the default), ``t`` is the time-of-flight c*t in m and
        ``pt`` is the energy deviation p_t = -Delta(gamma)/(beta_0 gamma_0)
        (unitless).
        Inside a ``beam.at_fixed_t()`` block, ``t`` takes the position z in m
        and ``pt`` takes the momentum deviation
        p_z = Delta(beta_z gamma)/(beta_0 gamma_0) (unitless).
    qm : float
        Charge over mass in 1/eV.
    bunch_charge : float, optional
        Magnitude of the total charge of the particles added in this call in C.
    w : array_like, optional
        Weight of each particle: how many real particles to represent.
    sx, sy, sz : array_like, optional
        Spin components in x, y, z.
    """

def ix_pc_at_fixed_t(*args, **kwds):
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

def register_ImpactXParticleContainer_extension(ixpc):
    """
    ImpactXParticleContainer helper methods
    """
