#!/usr/bin/env python3
#
# Copyright 2022 ImpactX contributors
# Authors: Marco Garten
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import warnings

import scipy.constants as sc

from impactx import elements


def madx2impactx_lattice(parsed_beamline,nslice=25):
    """
    Function that converts a list of elements in the MADXParser format into ImpactX format
    :param parsed_beamline: list of dictionaries
    :return: list of translated dictionaries
    """

    madx_to_impactx_dict = {
        "MARKER": "None",
        "DRIFT": "Drift",
        "SBEND": "Sbend",  # Sector Bending Magnet
        "QUAD": "Quad",  # Quadrupole
        "DIPEDGE": "DipEdge",
        "MULTIPOLE": "Multipole",
        "NLLENS": "NonlinearLens",
        # TODO Figure out how to identify these
        "ShortRF": "ShortRF",
        "ConstF": "ConstF",
    }

    impactx_beamline = []

    for d in parsed_beamline:

        if d['type'] in [k.casefold() for k in list(madx_to_impactx_dict.keys())]:
            if d['name'] == "drift":
                impactx_beamline.append(
                    elements.Drift(ds=d['l'], nslice=nslice)
                )
            elif d['name'] == "quadrupole":
                impactx_beamline.append(
                    elements.Quad(ds=d['l'], k=d['k1'], nslice=nslice)
                )
        else:
            raise NotImplementedError(
                "The beamline element named ",
                d['name'],
                "of type ",
                d['type'],
                "is not implemented in impactx.elements.",
                "Available elements are:",
                list(madx_to_impactx_dict.keys())
            )
    return impactx_beamline

def madx2impactx_beam(particle, charge=None, mass=None, energy=None):
    """
    Function that converts a list of beam parameter dictionaries in the MADXParser format into ImpactX format

    Rules following https://mad.web.cern.ch/mad/releases/5.02.08/madxuguide.pdf pages 55f.

    :param str particle: reference particle name
    :param float charge: particle charge (proton charge units)
    :param float mass: particle mass (electron masses)
    :param float energy: particle energy (GeV)
        - MAD-X default: 1 GeV
    :return dict: dictionary containing particle and beam attributes in ImpactX units
    """

    GeV2MeV = 1e-3
    kg2MeV = sc.c ** 2 / sc.electron_volt * 1e-6
    muon_mass = sc.physical_constants['electron-muon mass ratio'][0] / sc.m_e
    if energy is None:
        energy_MeV = 1e3  # MAD-X default is 1 GeV particle energy
    else:
        energy_MeV = energy * GeV2MeV

    impactx_beam = {
        'positron': {'mass': sc.m_e * kg2MeV, 'charge': 1},
        'electron': {'mass': sc.m_e * kg2MeV, 'charge': -1},
        'proton': {'mass': sc.m_p * kg2MeV, 'charge': 1},
        'antiproton': {'mass': sc.m_p * kg2MeV, 'charge': -1},
        'posmuon': {'mass': muon_mass * kg2MeV, 'charge': 1},  # positively charged muon (anti-muon)
        'negmuon': {'mass': muon_mass * kg2MeV, 'charge': -1},  # negatively charged muon
        'ion': {'mass': sc.m_u * kg2MeV, 'charge': 1},
        'generic': {'mass': mass, 'charge': charge}
    }


    if particle not in impactx_beam.keys():
        warnings.warn('Particle species name "' + particle + '" not in ' + impactx_beam.keys(), UserWarning)
        print('Choosing generic particle species, using provided `charge`, `mass` and `energy`.')
        _particle = 'generic'
    else:
        _particle = particle
        print('Choosing provided particle species "', particle, '", ignoring potentially provided `charge` and `mass` and setting defaults.')

    reference_particle = impactx_beam[_particle]
    reference_particle['energy'] = energy_MeV

    return reference_particle
