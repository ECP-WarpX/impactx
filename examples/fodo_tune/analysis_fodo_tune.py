#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import numpy as np
import openpmd_api as io
import PyNAFF as pnf

# Collect beam data series
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)

# Specify time series for particle j
j = 5
print(f"output for particle index = {j}")

# Create array of TBT data values
x = []
px = []
n = 0
for k_i, i in series.iterations.items():
    beam = i.particles["beam"]
    turn = beam.to_df()
    x.append(turn["position_x"][j])
    px.append(turn["momentum_x"][j])
    n = n + 1

# Output number of periods in data series
nturns = len(x)
print(f"number of periods = {nturns}")
print()

# Approximate the tune and closed orbit using the 4-turn formula:

# from x data only
argument = (x[0] - x[1] + x[2] - x[3])/(2.0*(x[1] - x[2]))
tune1 = np.arccos(argument)/(2.0*np.pi)
print(f"tune output from 4-turn formula, using x data = {tune1}")

# from px data only
argument = (px[0] - px[1] + px[2] - px[3])/(2.0*(px[1] - px[2]))
tune2 = np.arccos(argument)/(2.0*np.pi)
print(f"tune output from 4-turn formula, using px data = {tune2}")

# orbit offset
xco = (x[1]**2 - x[2]**2 + x[1]*x[3] - x[2]*x[0])/(3.0*(x[1] - x[2]) + x[3] - x[0])
pxco = (px[1]**2 - px[2]**2 + px[1]*px[3] - px[2]*px[0])/(3.0*(px[1] - px[2]) + px[3] - px[0])
print(f"orbit offset from 4-turn formula: (x[m], px[rad]) = ({xco},{pxco})")

print()

# Approximate the tune by using NAFF on the entire data series:

output = pnf.naff(x, turns=nturns, nterms=4, skipTurns=0, getFullSpectrum=True, window=1)
print(f"tune output from NAFF, using x data:")
print(output[0,1], output[1,1])

output = pnf.naff(px, turns=nturns, nterms=4, skipTurns=0, getFullSpectrum=True, window=1)
print(f"tune output from NAFF, using px data:")
print(output[0,1], output[1,1])

# Normalize TBT data using linear CS parameters
betax = 2.82161941
alphax = -1.59050035
z = []
for n in range(0,nturns):
    xn = x[n]/np.sqrt(betax)
    pxn = px[n]*np.sqrt(betax) + x[n]*alphax/np.sqrt(betax)
    z.append(complex(xn,-pxn))

output = pnf.naff(z, turns=nturns, nterms=4, skipTurns=0, getFullSpectrum=True, window=1)
tune3 = output[0,1]
print(f"tune output from NAFF, using normalized x-ipx data: {tune3}")


rtol = 1.0e-3
print(f"  rtol={rtol}")

assert np.allclose(
    [tune1,tune2,tune3],
    [ 
        0.1883,
        0.1883,
        0.1883,
    ],
    rtol=rtol,
)
