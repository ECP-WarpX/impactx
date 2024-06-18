"""

This file is part of ImpactX

Copyright 2023 ImpactX contributors
Authors: Axel Huebl
License: BSD-3-Clause-LBNL
"""

from __future__ import annotations

__all__ = ["ix_pc_plot_mpl_phasespace", "register_ImpactXParticleContainer_extension"]

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
