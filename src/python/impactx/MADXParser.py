#!/usr/bin/env python3
#
# Copyright 2022 ImpactX contributors
# Authors: Matthias Frey, Andreas Adelmann, Marco Garten
# License: BSD-3-Clause
#
# Changelog:
#   Oct 13th, 2017: original version for pyAcceLEGOrator elements
#   Aug 10th, 2022: adapted for standalone use
#
# -*- coding: utf-8 -*-
import os
import re
import warnings


class MADXParserError(Exception):
    pass


class MADXInputError(MADXParserError):
    def __init__(self, args, with_traceback):
        self.args = args
        self.with_traceback = with_traceback


class MADXInputWarning(UserWarning):
    pass


class MADXParser:
    """
    Simple MADX parser.
    It expects a single line per element.
    """

    def __init__(self):

        self.__drift = {"name": "", "l": 0.0, "type": "drift"}

        self.__drift_pattern = r"(.*):drift,(.*)=(.*);"

        self.__quadrupole = {"name": "", "l": 0.0, "k1": 0.0, "type": "quadrupole"}

        # don't count name and type --> len - 2
        self.__nQuad = 2 * (len(self.__quadrupole) - 2)

        self.quad_pattern = r"(.*):quadrupole,(.*)=(.*),(.*)=(.*);"

        self.__sbend = {
            "name": "",
            "l": 0.0,
            "angle": 0.0,
            "k1": 0.0,
            "e1": 0.0,
            "e2": 0.0,
            "type": "sbend",
        }

        self.__sbend_pattern = (
            r"(.*):sbend,(.*)=(.*),(.*)=(.*),(.*)=(.*),(.*)=(.*),(.*)=(.*);"
        )

        # don't count name and type --> len - 2
        # TODO add rbend
        self.__nDipole = 2 * (len(self.__sbend) - 2)

        self.__dipedge = {
            "name": "",
            "h": 0.0,
            "e1": 0.0,
            "fint": 0.0,
            "hgap": 0.0,
            "tilt": 0.0,
            "type": "dipedge",
        }

        self.__dipedge_pattern = (
            r"(.*):dipedge,(.*)=(.*),(.*)=(.*),(.*)=(.*),(.*)=(.*),(.*)=(.*);"
        )

        # don't count name and type --> len - 2
        self.__nDipedge = 2 * (len(self.__dipedge) - 2)

        self.beam = {
            "energy": 0.0,
            # TODO extend by 'PC'
            "particle": "",
        }

        self.__nBeam = 2 * len(self.beam)

        self.beam_pattern = r"beam,(.*)=(.*),(.*)=(.*);"

        self.__line = {
            "name": "",
            "elem": [],
        }

        self.__line_pattern = r"(.*):line=\(+(.*)\);"

        self.sequence = {"name": ""}

        self.seq_pattern = r"use,sequence=(.*);"

        self.__elements = []
        self.__lines = []

        self.__lattice = []

    # nonblank_lines inspired by
    # https://stackoverflow.com/questions/4842057/easiest-way-to-ignore-blank-lines-when-reading-a-file-in-python
    def nonblank_lines_to_lowercase(self, f):
        for l in f:
            line = l.rstrip().casefold()
            if line:
                yield line

    def parse(self, fn):
        """
        fn (str)    filename

        """

        if not os.path.isfile(fn):
            raise FileNotFoundError(f"File '{fn}' not found!")

        nLine = 0

        with open(fn, "r") as f:
            for line in self.nonblank_lines_to_lowercase(f):
                nLine += 1
                line = self._noWhitespace(line)

                if line[0] == "!":
                    # this is a comment
                    pass
                elif "drift" in line:

                    obj = re.match(self.__drift_pattern, line)

                    # first tag is name
                    self.__drift["name"] = obj.group(1)

                    if obj.group(2) in self.__drift:
                        self.__drift[obj.group(2)] = float(obj.group(3))
                    else:
                        raise MADXInputError(
                            "Drift",
                            "Line "
                            + str(nLine)
                            + ": Parameter "
                            + "'"
                            + obj.group(2)
                            + "'"
                            + " does not exist for drift.",
                        )

                    self.__elements.append(self.__drift.copy())

                elif "quadrupole" in line:

                    obj = re.match(self.quad_pattern, line)

                    # first tag is name
                    self.__quadrupole["name"] = obj.group(1)

                    for i in range(2, self.__nQuad + 2, 2):
                        if obj.group(i) in self.__quadrupole:
                            self.__quadrupole[obj.group(i)] = float(obj.group(i + 1))
                        else:
                            raise MADXInputError(
                                "Quadrupole",
                                "Line "
                                + str(nLine)
                                + ": Parameter "
                                + "'"
                                + obj.group(i)
                                + "'"
                                + " does not exist for quadrupole.",
                            )

                    self.__elements.append(self.__quadrupole.copy())

                elif "sbend" in line:

                    obj = re.match(self.__sbend_pattern, line)

                    # first tag is name
                    self.__sbend["name"] = obj.group(1)

                    for i in range(2, self.__nDipole + 2, 2):

                        if obj.group(i) in self.__sbend:
                            self.__sbend[obj.group(i)] = float(obj.group(i + 1))
                        else:
                            raise MADXInputError(
                                "Dipole",
                                "Line "
                                + str(nLine)
                                + ": Parameter "
                                + "'"
                                + obj.group(i)
                                + "'"
                                + " does not exist for dipole.",
                            )

                    self.__elements.append(self.__sbend.copy())

                elif "dipedge" in line:

                    obj = re.match(self.__dipedge_pattern, line)

                    # first tag is name
                    self.__dipedge["name"] = obj.group(1)

                    for i in range(2, self.__nDipedge + 2, 2):

                        if obj.group(i) in self.__dipedge:
                            self.__dipedge[obj.group(i)] = float(obj.group(i + 1))
                        else:
                            raise MADXInputError(
                                "DipEdge",
                                "Line "
                                + str(nLine)
                                + ": Parameter "
                                + "'"
                                + obj.group(i)
                                + "'"
                                + " does not exist for dipole edge.",
                            )

                    self.__elements.append(self.__dipedge.copy())

                elif "marker" in line:
                    pass

                elif "beam" in line:

                    obj = re.match(self.beam_pattern, line)

                    for i in range(1, self.__nBeam, 2):
                        if obj.group(i) in self.beam:
                            if obj.group(i + 1).isdigit():
                                self.beam[obj.group(i)] = float(obj.group(i + 1))
                            else:
                                self.beam[obj.group(i)] = obj.group(i + 1)

                        else:
                            raise MADXInputError(
                                "Beam",
                                "Line "
                                + str(nLine)
                                + ": Parameter "
                                + "'"
                                + obj.group(i)
                                + "'"
                                + " does not exist for beam.",
                            )

                elif "line" in line:

                    obj = re.match(self.__line_pattern, line)

                    self.__line["name"] = obj.group(1)

                    lines = obj.group(2).split(",")
                    newlines = []

                    # check multiplication and insert that many lines
                    for l in lines:
                        if "*" in l:
                            tmp = l.split("*")

                            n = 0
                            ll = ""

                            if tmp[0].isdigit():
                                n = int(tmp[0])
                                ll = tmp[1]
                            else:
                                n = int(tmp[1])
                                ll = tmp[0]

                            newlines += [ll] * n

                        else:
                            newlines.append(l)

                    self.__line["elem"] = newlines

                    self.__lines.append(self.__line.copy())

                elif "use" in line and "sequence" in line:

                    obj = re.match(self.seq_pattern, line)

                    self.sequence["name"] = obj.group(1)
                else:
                    raise MADXInputError("", "Error at line " + str(nLine))

        # 14. Oct. 2017,
        # https://stackoverflow.com/questions/7900882/extract-item-from-list-of-dictionaries
        start = [l for l in self.__lines if l["name"] == self.sequence["name"]][0]

        self._flatten(start)

        # we need to add start line
        self.__lattice = [start] + self.__lattice

        self.__lattice = self._combine(self.__lattice)

    def _flatten(self, line):
        """
        Find sublines.

        """

        name = line["name"]

        for l in self.__lines:
            if name in l["name"]:
                for ll in l["elem"]:
                    for lll in self.__lines:
                        if lll["name"] == ll:
                            self.__lattice.append(lll)
                            self._flatten(lll)

    def _combine(self, lattice):
        """
        Combine to one list of all basic
        elements.

        return a list of of element dictionaries
        """

        l1 = self.__lattice[0]

        for i in range(1, len(self.__lattice)):
            l2 = self.__lattice[i]

            for e in l1["elem"]:
                if l2["name"] == e:
                    idx = l1["elem"].index(e)
                    l1["elem"].remove(e)
                    l1["elem"] = l1["elem"][0:idx] + l2["elem"] + l1["elem"][idx:]

        return l1

    def _noWhitespace(self, string):
        """
        Remove white space from a string.

        14. Oct. 2017,
        https://stackoverflow.com/questions/3739909/how-to-strip-all-whitespace-from-string

        """
        return string.replace(" ", "")

    def __str__(self):

        if self.__lattice:
            length = 0.0

            # drift, dipole, quadrupole, dipedge
            nTypes = [0, 0, 0, 0]

            for elem in self.__lattice["elem"]:
                for e in self.__elements:
                    if elem == e["name"]:
                        length += e["l"]

                        if e["type"] == "drift":
                            nTypes[0] += 1
                        elif e["type"] == "sbend":
                            nTypes[1] += 1
                        elif e["type"] == "quadrupole":
                            nTypes[2] += 1
                        elif e["type"] == "dipedge":
                            nTypes[3] += 1
                        break

            sign = "*" * 70
            info = (
                sign
                + "\n"
                + "MADX-Parser information:\n"
                + "         length:\t"
                + str(length)
                + " [m]\n"
                + "      #elements:\t"
                + str(len(self.__lattice["elem"]))
                + "\n"
                + "            *       #drift:\t"
                + str(nTypes[0])
                + "\n"
                + "            *      #dipole:\t"
                + str(nTypes[1])
                + "\n"
                + "            *  #quadrupole:\t"
                + str(nTypes[2])
                + "\n"
                + "            *     #dipedge:\t"
                + str(nTypes[3])
                + "\n"
                + "           beam:\t\n"
                + "            *     particle:\t"
                + self.beam["particle"]
                + "\n"
                + "            * total energy:\t"
                + str(self.beam["energy"])
                + " [GeV]\n"
                + sign
            )

            return info
        else:
            return "No information available."

    def getBeamline(self):
        if self.__lattice:

            beamline = []

            for elem in self.__lattice["elem"]:
                for e in self.__elements:
                    if elem == e["name"]:
                        if e["type"] == "drift":
                            print("Drift L= " + str(e["l"]))
                            beamline.append(e)
                        elif e["type"] == "sbend":
                            print("Sbend L= ", e["l"], " angle = ", e["angle"])
                            beamline.append(e)
                        elif e["type"] == "quadrupole":
                            print("Quadrupole L= ", e["l"], " k1 = ", e["k1"])
                            beamline.append(e)
                        elif e["type"] == "dipedge":
                            print("Dipedge H= ", e["h"], " E1 = ", e["e1"])
                            beamline.append(e)
                        else:
                            print("Skipping element type " + "'" + e["type"] + "'")
            return beamline

        else:
            return []

    def getParticle(self):
        particle = self.beam["particle"]
        known_particles = [
            "positron",
            "electron",
            "proton",
            "antiproton",
            "posmuon",
            "negmuon",
            "ion",
        ]

        if particle not in known_particles:
            warning_message = (
                "No particle type " + "'" + self.beam["particle"] + "' available."
            )
            warnings.warn(warning_message, MADXInputWarning)
        else:
            pass

        return particle

    def getEtot(self):
        return self.beam["energy"]
