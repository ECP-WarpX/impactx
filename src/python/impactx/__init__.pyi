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
from amrex import space3d as amr
from impactx.extensions.ImpactXParIter import register_ImpactXParIter_extension
from impactx.extensions.ImpactXParticleContainer import register_ImpactXParticleContainer_extension
from impactx.impactx_pybind import Config
from impactx.impactx_pybind import CoordSystem
from impactx.impactx_pybind import ImpactX
from impactx.impactx_pybind import ImpactXParConstIter
from impactx.impactx_pybind import ImpactXParIter
from impactx.impactx_pybind import ImpactXParticleContainer
from impactx.impactx_pybind import RefPart
from impactx.impactx_pybind import coordinate_transformation
from impactx.impactx_pybind import distribution
from impactx.impactx_pybind import elements
from impactx.impactx_pybind import push
from impactx.madx_to_impactx import read_beam
from impactx.madx_to_impactx import read_lattice
import os as os
from . import MADXParser
from . import extensions
from . import impactx_pybind
from . import madx_to_impactx
__all__ = ['Config', 'CoordSystem', 'ImpactX', 'ImpactXParConstIter', 'ImpactXParIter', 'ImpactXParticleContainer', 'MADXParser', 'RefPart', 'amr', 'coordinate_transformation', 'cxx', 'distribution', 'elements', 'extensions', 'impactx_pybind', 'madx_to_impactx', 'os', 'push', 'read_beam', 'read_lattice', 'register_ImpactXParIter_extension', 'register_ImpactXParticleContainer_extension', 's', 't']
__author__: str = 'Axel Huebl, Chad Mitchell, Ryan Sandberg, Marco Garten, Ji Qiang, et al.'
__license__: str = 'BSD-3-Clause-LBNL'
__version__: str = '24.07'
s: impactx_pybind.CoordSystem  # value = <CoordSystem.s: 0>
t: impactx_pybind.CoordSystem  # value = <CoordSystem.t: 1>
cxx = impactx_pybind
