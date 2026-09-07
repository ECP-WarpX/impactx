#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import matplotlib.pyplot as plt
import numpy as np
import openpmd_api as io
import PyNAFF as pnf

# Collect beam data series
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)


def get_horizontal_tune(particle_id):

    # Create array of TBT data values
    x = []
    px = []
    n = 0
    for k_i, i in series.iterations.items():
        beam = i.particles["beam"]
        turn = beam.to_df()
        x.append(turn["position_x"][particle_id])
        px.append(turn["momentum_x"][particle_id])
        n = n + 1

    # Number of periods in data series
    nturns = len(x)

    # matched Twiss functions (from input)
    alpha_x = 1.29174698
    beta_x = 12.79711091

    # Normalize TBT data using Twiss functions:
    z = []
    for n in range(0, nturns):
        xn = x[n] / np.sqrt(beta_x)
        pxn = px[n] * np.sqrt(beta_x) + x[n] * alpha_x / np.sqrt(beta_x)
        z.append(complex(xn, -pxn))

    # Approximate the tune by using NAFF on the entire data series:

    # Option to use raw x-values only:
    # output = pnf.naff(
    #    x, turns=nturns, nterms=4, skipTurns=0, getFullSpectrum=True, window=1
    # )

    output = pnf.naff(
        z, turns=nturns, nterms=4, skipTurns=0, getFullSpectrum=True, window=1
    )
    tune = output[0, 1]

    return tune


# Accesss initial x-coordinates:
beam_initial = series.iterations[1].particles["beam"].to_df()

mm_scale = 1.0e3
x = beam_initial["position_x"] * mm_scale
print("length")
print(len(x))

tune_x = []
for j in range(len(x)):
    tune = get_horizontal_tune(j)
    tune_x.append(tune)

plt.figure(figsize=(8, 6))
plt.scatter(x, tune_x, c="red", s=3)
plt.xlabel("x [mm]", fontsize=12)
plt.ylabel("tune", fontsize=12)
plt.title("Tune vs. Horizontal Position")
plt.show()
