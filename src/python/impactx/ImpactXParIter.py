"""
This file is part of ImpactX

Copyright 2023 ImpactX contributors
Authors: Axel Huebl
License: BSD-3-Clause-LBNL
"""


def soa_real_comps(self, num_comps):
    """
    Name the ImpactX ParticleReal components in SoA.

    Parameters
    ----------
    self : impactx_pybind module
      used to query freestanding functions
    num_comps : int
      number of components to generate names for.

    Returns
    -------
    A list of length num_comps with values.
    """
    names = self.get_RealSoA_names(num_comps)

    if len(names) != num_comps:
        raise RuntimeError(
            "num_comps ({num_comps}) is not equal to length of RealSoA_names ({len(names)})"
        )

    return names


def soa_int_comps(self, num_comps):
    """
    Name the ImpactX int components in SoA.

    Parameters
    ----------
    self : impactx_pybind module
      used to query freestanding functions
    num_comps : int
      number of components to generate names for.

    Returns
    -------
    A list of length num_comps with values.
    """
    names = self.get_intSoA_names(num_comps)

    if len(names) != num_comps:
        raise RuntimeError(
            f"num_comps ({num_comps}) is not equal to length of intSoA_names ({len(names)})"
        )

    return names


def soa(self, impactx_pybind):
    """Get the StructOfArrays on the current tile"""
    soa = super(type(self), self).soa()

    # overwrite name providers
    soa.soa_real_comps = lambda num_comps: soa_real_comps(impactx_pybind, num_comps)
    soa.soa_int_comps = lambda num_comps: soa_int_comps(impactx_pybind, num_comps)

    return soa


def register_ImpactXParIter_extension(impactx_pybind):
    """ImpactXParIter helper methods"""

    impactx_pybind.ImpactXParIter.soa = lambda self: soa(self, impactx_pybind)
    impactx_pybind.ImpactXParConstIter.soa = lambda self: soa(self, impactx_pybind)
