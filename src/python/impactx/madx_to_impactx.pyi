from __future__ import annotations

import typing
import warnings as warnings

import impactx.impactx_pybind
from impactx.impactx_pybind import RefPart, elements
from impactx.MADXParser import MADXParser

__all__ = [
    "MADXParser",
    "RefPart",
    "beam",
    "elements",
    "lattice",
    "read_beam",
    "read_lattice",
    "sc",
    "warnings",
]

class sc:
    """

    This class is used in lieu of scipy.constants
    to avoid a direct dependency on it.
    At the time of writing, this file was the only
    one requiring scipy in the ImpactX source outside
    of examples.

    """

    c: typing.ClassVar[float] = 299792458.0
    electron_volt: typing.ClassVar[float] = 1.602176634e-19
    m_e: typing.ClassVar[float] = 9.1093837015e-31
    m_p: typing.ClassVar[float] = 1.67262192369e-27
    m_u: typing.ClassVar[float] = 1.6605390666e-27
    physical_constants: typing.ClassVar[dict] = {
        "electron-muon mass ratio": (0.00483633169, "", 1.1e-10)
    }

def beam(particle, charge=None, mass=None, energy=None):
    """

    Function that converts a list of beam parameter dictionaries in the MADXParser format into ImpactX format

    Rules following https://mad.web.cern.ch/mad/releases/5.02.08/madxuguide.pdf pages 55f.

    :param str particle: reference particle name
    :param float charge: particle charge (proton charge units)
    :param float mass: particle mass (electron masses)
    :param float energy: total particle energy (GeV)
        - MAD-X default: 1 GeV
    :return dict: dictionary containing particle and beam attributes in ImpactX units

    """

def lattice(parsed_beamline, nslice=1):
    """

    Function that converts a list of elements in the MADXParser format into ImpactX format
    :param parsed_beamline: list of dictionaries
    :param nslice: number of ds slices per element
    :return: list of translated dictionaries

    """

def read_beam(ref: impactx.impactx_pybind.RefPart, madx_file):
    """

    Function that reads elements from a MAD-X file into a list of ImpactX.KnownElements
    :param RefPart ref: ImpactX reference particle (passed by reference)
    :param madx_file: file name to MAD-X file with beamline elements
    :return: list of ImpactX.KnownElements

    """

def read_lattice(madx_file, nslice=1):
    """

    Function that reads elements from a MAD-X file into a list of ImpactX.KnownElements
    :param madx_file: file name to MAD-X file with beamline elements
    :param nslice: number of ds slices per element
    :return: list of ImpactX.KnownElements

    """
