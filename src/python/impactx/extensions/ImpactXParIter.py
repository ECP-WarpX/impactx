"""
This file is part of ImpactX

Copyright 2023 ImpactX contributors
Authors: Axel Huebl
License: BSD-3-Clause-LBNL
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
    names = pti.pc().RealSoA_names

    if len(names) != num_comps:
        raise RuntimeError(
            "num_comps ({num_comps}) is not equal to length of RealSoA_names ({len(names)})"
        )

    return names


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
    names = pti.pc().intSoA_names

    if len(names) != num_comps:
        raise RuntimeError(
            f"num_comps ({num_comps}) is not equal to length of intSoA_names ({len(names)})"
        )

    return names


def soa(self):
    """Get the StructOfArrays on the current tile

    Parameters
    ----------
    self : ImpactXParIter or ImpactXParConstIter
      used to query particle container component names
    """
    soa = super(type(self), self).soa()

    # overwrite name providers
    soa.soa_real_comps = lambda num_comps: soa_real_comps(self, num_comps)
    soa.soa_int_comps = lambda num_comps: soa_int_comps(self, num_comps)

    return soa


def register_ImpactXParIter_extension(impactx_pybind):
    """ImpactXParIter helper methods"""

    impactx_pybind.ImpactXParIter.soa = soa
    impactx_pybind.ImpactXParConstIter.soa = soa
