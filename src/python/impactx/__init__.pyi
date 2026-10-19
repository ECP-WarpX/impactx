"""
impactx_pybind
--------------
.. currentmodule:: impactx_pybind

.. autosummary::
   :toctree: _generate
   ImpactX
   distribution
   elements

"""

from __future__ import annotations

import os as os

from amrex import space3d as amr
from amrex.space3d.amrex_3d_pybind import SmallMatrix_3x1_F_SI1_double as Vector3
from amrex.space3d.amrex_3d_pybind import SmallMatrix_3x6_F_SI1_double as Map3x6
from amrex.space3d.amrex_3d_pybind import SmallMatrix_6x6_F_SI1_double as Map6x6
from impactx.distribution_input_helpers import twiss
from impactx.extensions.Elements import register_elements_value_semantics
from impactx.extensions.ImpactXParticleContainer import (
    register_ImpactXParticleContainer_extension,
)
from impactx.extensions.KnownElementsList import register_KnownElementsList_extension
from impactx.extensions.RFCavity import register_RFCavity_extension
from impactx.extensions.SoftQuadrupole import register_SoftQuadrupole_extension
from impactx.extensions.SoftSolenoid import register_SoftSolenoid_extension
from impactx.fourier import fourier_coefficients
from impactx.impactx_pybind import (
    Config,
    CoordSystem,
    Envelope,
    ImpactX,
    ImpactXParConstIter,
    ImpactXParIter,
    ImpactXParticleContainer,
    ItemsView,
    KeysView,
    RefPart,
    UnorderedMap,
    ValuesView,
    coordinate_transformation,
    create_envelope,
    distribution,
    elements,
    flatten_charge_to_2D,
    push,
    reverse,
    wakeconvolution,
)
from impactx.impactx_pybind.elements import FilteredElementsList
from impactx.madx_to_impactx import read_beam

from . import (
    MADXParser,
    distribution_input_helpers,
    element_models,
    extensions,
    fourier,
    impactx_pybind,
    madx_to_impactx,
    plot,
)

__all__: list[str] = [
    "Config",
    "CoordSystem",
    "Envelope",
    "FilteredElementsList",
    "ImpactX",
    "ImpactXParConstIter",
    "ImpactXParIter",
    "ImpactXParticleContainer",
    "ItemsView",
    "KeysView",
    "MADXParser",
    "Map3x6",
    "Map6x6",
    "RefPart",
    "UnorderedMap",
    "ValuesView",
    "Vector3",
    "amr",
    "coordinate_transformation",
    "create_envelope",
    "distribution",
    "distribution_input_helpers",
    "element_models",
    "elements",
    "extensions",
    "flatten_charge_to_2D",
    "fourier",
    "fourier_coefficients",
    "impactx_pybind",
    "madx_to_impactx",
    "os",
    "plot",
    "push",
    "read_beam",
    "register_ImpactXParticleContainer_extension",
    "register_KnownElementsList_extension",
    "register_RFCavity_extension",
    "register_SoftQuadrupole_extension",
    "register_SoftSolenoid_extension",
    "register_elements_value_semantics",
    "reverse",
    "s",
    "t",
    "twiss",
    "wakeconvolution",
]
__author__: str = (
    "Axel Huebl, Chad Mitchell, Ryan Sandberg, Marco Garten, Ji Qiang, et al."
)
__license__: str = "BSD-3-Clause-LBNL"
__version__: str = "26.10"
s: impactx_pybind.CoordSystem
t: impactx_pybind.CoordSystem
