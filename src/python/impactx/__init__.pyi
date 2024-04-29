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
from impactx.impactx_pybind import (
    Config,
    CoordSystem,
    ImpactX,
    ImpactXParConstIter,
    RefPart,
    coordinate_transformation,
    distribution,
    elements,
    push,
)
from impactx.ImpactXParIter import register_ImpactXParIter_extension
from impactx.ImpactXParticleContainer import register_ImpactXParticleContainer_extension
from impactx.madx_to_impactx import read_beam, read_lattice

from . import (
    ImpactXParIter,
    ImpactXParticleContainer,
    MADXParser,
    impactx_pybind,
    madx_to_impactx,
)

__all__ = [
    "Config",
    "CoordSystem",
    "ImpactX",
    "ImpactXParConstIter",
    "ImpactXParIter",
    "ImpactXParticleContainer",
    "MADXParser",
    "RefPart",
    "amr",
    "coordinate_transformation",
    "cxx",
    "distribution",
    "elements",
    "impactx_pybind",
    "madx_to_impactx",
    "os",
    "push",
    "read_beam",
    "read_lattice",
    "register_ImpactXParIter_extension",
    "register_ImpactXParticleContainer_extension",
    "s",
    "t",
]
__author__: str = (
    "Axel Huebl, Chad Mitchell, Ryan Sandberg, Marco Garten, Ji Qiang, et al."
)
__license__: str = "BSD-3-Clause-LBNL"
__version__: str = "24.04"
s: impactx_pybind.CoordSystem  # value = <CoordSystem.s: 0>
t: impactx_pybind.CoordSystem  # value = <CoordSystem.t: 1>
cxx = impactx_pybind
