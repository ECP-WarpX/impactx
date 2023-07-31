# The following script is intended to parse the file 'fermi-booster-madx-sxf'
# corresponding to a thin-kick reduced model of the Fermilab Booster
# lattice (from F. Schmidt), for space charge studies in ImpactX.
# The input file is provided in Standard eXchange Format (SXF),
# described here:  https://doi.org/10.2172/1119545
# The present version is a draft for preliminary testing only.
# Authors:  C. Mitchell, A. Huebl


# Import the required Python packages

import re

import pandas as pd

# Read the input file as a single string, since elements span multiple lines

with open("fermi-booster-madx-sxf") as f:
    text = f.read()


# Create list of element strings, which are separated by semicolon ;

elements = text.split(";")


# Define some regular expressions (list to be expanded later)

rx_dict = {
    "tag": re.compile(r"tag = (?P<tag>.*)"),  # element name
    "zloc": re.compile(r"at = (?P<zloc>\d*\.?\d*)\n"),  # element location
    "lrad": re.compile(r"lrad =  (?P<lrad>\d*\.?\d*) "),  # fictitious length
}


# Basic function for parsing a single element


def _parse_line(element):
    for key, rx in rx_dict.items():
        match = rx.search(element)
        if match:
            return key, match
    return None, None


# Loop over elements and store data (test keys separately for now)

tagdata = []
zlocdata = []
lraddata = []

for i in range(len(elements)):
    key, match = _parse_line(elements[i])
    if key == "tag":
        tag = match.group("tag")
        tagdata.append(tag)
    if key == "zloc":
        zloc = match.group("zloc")
        zlocdata.append(zloc)
    if key == "lrad":
        lrad = match.group("lrad")
        lraddata.append(lrad)

print(len(tagdata))
print(tagdata)

print(len(zlocdata))
print(zlocdata)

print(len(lraddata))
print(lraddata)


# Analyze data using Pandas
# To follow
# data = pd.DataFrame(data)
