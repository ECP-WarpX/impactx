"""

This file is part of ImpactX

Copyright 2023 ImpactX contributors
Authors: Axel Huebl
License: BSD-3-Clause-LBNL
"""

from __future__ import annotations

__all__ = [
    "register_ImpactXParIter_extension",
    "soa",
    "soa_int_comps",
    "soa_real_comps",
]

def register_ImpactXParIter_extension(impactx_pybind):
    """
    ImpactXParIter helper methods
    """

def soa(self):
    """
    Get the StructOfArrays on the current tile

        Parameters
        ----------
        self : ImpactXParIter or ImpactXParConstIter
          used to query particle container component names

    """

def soa_int_comps(pti, num_comps):
    """

    Name the ImpactX int components in SoA.

    Parameters
    ----------
    pti : ImpactXParIter or ImpactXParConstIter
      used to query particle container component names
    num_comps : int
      number of components to generate names for.

    Returns
    -------
    A list of length num_comps with values.

    """

def soa_real_comps(pti, num_comps):
    """

    Name the ImpactX ParticleReal components in SoA.

    Parameters
    ----------
    pti : ImpactXParIter or ImpactXParConstIter
      used to query particle container component names
    num_comps : int
      number of components to generate names for.

    Returns
    -------
    A list of length num_comps with values.

    """
