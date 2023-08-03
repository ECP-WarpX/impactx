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

# cut off the ring sequence wrapping lines (top/bottom of the file)
sequence = text[36:-45]

# Create list of element strings, which are separated by semicolon ;

elements = sequence.split(";")

#print(elements[0])
#print(elements[1])


# Define some regular expressions (list to be expanded later)

# https://regex101.com
rx_dict = {
    "name": re.compile(r"^\s*(?P<name>[\w\.]+)\s+"),  # element name
    "type": re.compile(r"^\s*[\w\.]+\s+(?P<type>[\w]+)\s*\{"),  # element name
    "tag": re.compile(r"tag = (?P<tag>[\w\.]+)"),  # element tag (often same as name)
    "zloc": re.compile(r"at = (?P<zloc>.*)[\s\n\r}]+"),  # element location

    "body": re.compile(r"body\s*=\s*\{(?P<body>.+)\s*\}.*\}"),  # ...
    # multipole
    "lrad": re.compile(r"lrad\s*=\s*(?P<lrad>\d*\.?\d*)[\s\n\r\}]+"),  # fictitious length
    "kl": re.compile(r"kl\s*=\s*\[(?P<kl>.*)\]"),  # coefficients for the multipole strength
    # dipedge
    "e1": re.compile(r"e1\s*=\s*(?P<e1>.*)\s*"),  # ...
    "e2": re.compile(r"e2\s*=\s*(?P<e2>.*)\s*"),  # ...
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

def parse_element(element):
    """
    Input: SXF Element String
    Outpu: ImpactX Element
    """
    name = parse_one_group(element, "name")
    etyp = parse_one_group(element, "type")
    tag = parse_one_group(element, "tag")
    print(f"name={name}")
    print(f"etyp={etyp}")
    print(f"tag={tag}")
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
                kl = re.sub(r'\s+', ' ', kl)
                # to list
                kl = kl.split(" ")
                # convert strings to floats
                kl = list(map(float, kl))
                print(f"lrad={lrad}")
                print(f"kl={kl}")
            
                # TODO: create ImpactX.Multiple and return
        else:
            print("... empty body - ignored")

    elif etyp == "dipedge":
        body = parse_one_group(element, "body")
        #print(f"body={body}")

        if body is not None:
            e1 = parse_one_group(body, "e1")
            #e2 = parse_one_group(body, "e2")

            if e1 is None:
                print("... empty e1 - ignored")
            else:
                e1 = float(e1)
                print(f"e1={e1}")
                # TODO: create ImpactX.Dipedge and return
        else:
            print("... empty body - ignored")

    elif etyp == "rfcavity":
        body = parse_one_group(element, "body")
        #print(f"body={body}")
        
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
                # TODO: create ImpactX.ThinRF? and return
        else:
            print("... empty body - ignored")
    else:
        print(f"Unknown type: {etyp}")

# Loop over elements and store data (test keys separately for now)

tagdata = []
zlocdata = []
lraddata = []

for element in elements: #elements[0:12]:
    parse_element(element)
    # TODO: append returned element to a beamline (as a list)
    print()

#print(len(tagdata))
#print(tagdata)

#print(len(zlocdata))
#print(zlocdata)

#print(len(lraddata))
#print(lraddata)


# Analyze data using Pandas
# To follow
# data = pd.DataFrame(data)
