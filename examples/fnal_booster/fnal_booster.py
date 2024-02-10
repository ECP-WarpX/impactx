#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

# The following script is intended to parse the file 'fermi-booster-madx-sxf'
# corresponding to a thin-kick reduced model of the Fermilab Booster
# lattice (from F. Schmidt), for space charge studies in ImpactX.
# The input file is provided in Standard eXchange Format (SXF),
# described here:  https://doi.org/10.2172/1119545
# The present version is a draft for preliminary testing only.
# Authors:  C. Mitchell, A. Huebl


# Import the required Python packages

import math
import re

import amrex.space3d as amr
from impactx import ImpactX, distribution, elements

# Read the input file as a single string, since elements span multiple lines

with open("fermi-booster-madx-sxf_corr_new") as f:
    text = f.read()

# cut off the ring sequence wrapping lines (top/bottom of the file)
sequence = text[36:-45]

# Create list of element strings, which are separated by semicolon ;

element_list = sequence.split(";")

# print(element_list[0])
# print(element_list[1])


# Define some regular expressions (list to be expanded later)

# https://regex101.com
rx_dict = {
    "name": re.compile(r"^\s*(?P<name>[\w\.]+)\s+"),  # element name
    "type": re.compile(r"^\s*[\w\.]+\s+(?P<type>[\w]+)\s*\{"),  # element name
    "tag": re.compile(r"tag = (?P<tag>[\w\.]+)"),  # element tag (often same as name)
    "zloc": re.compile(r"at\s*=\s*(?P<zloc>\d*\.?\d*)[\s\n\r\}]+"),  # element location
    "body": re.compile(r"body\s*=\s*\{(?P<body>.+)\s*\}.*\}"),  # ...
    # multipole
    "lrad": re.compile(
        r"lrad\s*=\s*(?P<lrad>\d*\.?\d*)[\s\n\r\}]+"
    ),  # fictitious length
    "kl": re.compile(
        r"kl\s*=\s*\[(?P<kl>.*)\]"
    ),  # coefficients for the multipole strength
    # dipedge
    "e1": re.compile(r"e1\s*=\s*(?P<e1>[+-]?\d*\.?\d*)[\s\n|r\}]+"),  # ...
    "e2": re.compile(r"e2\s*=\s*(?P<e2>[+-]?\d*\.?\d*)[\s\n|r\}]+"),  # ...
    "h": re.compile(r"h\s*=\s*(?P<h>[+-]?\d*\.?\d*)[\s?(\n|\r)\}]?+"),  # ...
    "tilt": re.compile(r"tilt\s*=\s*(?P<tilt>[+-]?\d*\.?\d*)[\s\n|r\}]?+"),  # ...
    #   more in the MAD-X file? fint, hgap, h - always zero
    # rfcavity
    "volt": re.compile(r"volt\s*=\s*(?P<volt>[\w\.]+)\s+"),  # ...
    "harmon": re.compile(r"harmon\s*=\s*(?P<harmon>[\w\.]+)\s*"),  # ...
}


# Basic function for parsing a single element


def parse_one_group(element, key):
    """Returns regex match or None"""
    rx = rx_dict[key]
    element = element.replace("\n", " ")
    match = rx.search(element)
    if match:
        return match.group(key).strip()
    else:
        None


def parse_element(sim, element, zprev):
    """
    Input: SXF Element String
    Outpu: ImpactX Element
    """
    name = parse_one_group(element, "name")
    etyp = parse_one_group(element, "type")
    tag = parse_one_group(element, "tag")
    zloc = parse_one_group(element, "zloc")
    print(f"name={name}")
    print(f"etyp={etyp}")
    print(f"tag={tag}")

    # convert z-location to float
    if zloc is None:
        print("... empty zloc - ignored")
        zloc = zprev
    else:
        zloc = float(zloc)
        dz = zloc - zprev

        print(f"elements.Drift(ds={dz})")
        sim.lattice.append(elements.Drift(ds=dz))

    # tag

    # vkicker & hkicker: they all seem to be turned off in our file
    ignored_types = ["beambeam", "marker", "hkicker", "vkicker", None]

    if etyp in ignored_types:
        print("... ignored")

    elif etyp == "multipole":
        body = parse_one_group(element, "body")

        if body is not None:
            lrad = parse_one_group(body, "lrad")
            kl = parse_one_group(body, "kl")

            if kl is None:
                print("... empty kl - ignored")
            else:
                # simplify spaces to one
                kl = re.sub(r"\s+", " ", kl)
                # to list
                kl = kl.split(" ")
                # convert strings to floats
                kl = list(map(float, kl))
                print(f"lrad={lrad}")
                print(f"kl={kl}")

                for i in range(len(kl)):
                    if i == 0:
                        angledeg = kl[i] * 180.0 / math.pi
                        print(f"elements.ThinDipole(theta={angledeg}, rc=1.0)")
                        sim.lattice.append(elements.ThinDipole(theta=kl[i], rc=1.0))
                    else:
                        print(
                            f"elements.Multipole(multiple={i}, K_normal={kl[i]}, K_skew=0.0)"
                        )
                        sim.lattice.append(
                            elements.Multipole(multiple=i, K_normal=kl[i], K_skew=0.0)
                        )
        else:
            print("... empty body - ignored")

    elif etyp == "dipedge":
        body = parse_one_group(element, "body")
        # print(f"body={body}")

        if body is not None:
            e1 = parse_one_group(body, "e1")
            # e2 = parse_one_group(body, "e2")
            h = parse_one_group(body, "h")
            tilt = parse_one_group(body, "tilt")

            if e1 is None:
                print("... empty e1 - ignored")
            elif h is None:
                print("... empty h - ignored")
            else:
                e1 = float(e1)
                print(f"e1={e1}")
                h = float(h)
                print(f"h={h}")
                if tilt is None:
                    print("... empty tilt - ignored")
                    angledeg = 0.0
                else:
                    tilt = float(tilt)
                    print(f"tilt={tilt}")
                    angledeg = tilt * 180.0 / math.pi
                print(
                    f"elements.DipEdge(psi={e1}, rc=1.0, g={h}, K2=0.125, rotation={angledeg})"
                )
                sim.lattice.append(
                    elements.DipEdge(psi=e1, rc=1.0, g=h, K2=0.125, rotation=angledeg)
                )

        else:
            print("... empty body - ignored")

    elif etyp == "rfcavity":
        body = parse_one_group(element, "body")
        # print(f"body={body}")

        if body is not None:
            volt = parse_one_group(body, "volt")
            harmon = parse_one_group(body, "harmon")

            if volt is None or harmon is None:
                print("... empty volt or harmon - ignored")
            else:
                volt = float(volt)
                harmon = int(harmon)
                print(f"volt={volt}")
                print(f"harmon={harmon}")

                f0 = 0.45079517081715459e6
                Erest = 938.27208816
                print(
                    f"elements.ShortRF(V={volt/Erest}, freq={harmon*f0}, phase=-90.0)"
                )
                sim.lattice.append(
                    elements.ShortRF(V=volt / Erest, freq=harmon * f0, phase=-90.0)
                )

        else:
            print("... empty body - ignored")
    else:
        print(f"Unknown type: {etyp}")

    return zloc


# Set up simulation input

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# intial proton beam parameters
tot_energy_MeV = 1.3382720460000002e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
proton_mass_MeV = 938.27204600000003
kin_energy_MeV = tot_energy_MeV - proton_mass_MeV
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(proton_mass_MeV).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    sigmaX=1.288697604e-6,
    sigmaY=1.288697604e-6,
    sigmaT=1.0e-6,
    sigmaPx=3.965223396e-6,
    sigmaPy=3.965223396e-6,
    sigmaPt=0.01,  # 1% energy spread
    muxpx=0.0,
    muypy=0.0,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# Initial beam diagnostics
sim.lattice.append(elements.BeamMonitor("monitor", backend="h5"))

# Loop over elements

zprev = 0.0
for element in element_list:  # elements[0:12]:
    zloc = parse_element(sim, element, zprev)
    zprev = zloc

# Final beam diagnostics
sim.lattice.append(elements.BeamMonitor("monitor", backend="h5"))

# run simulation
sim.evolve()

# clean shutdown
del sim
amr.finalize()
