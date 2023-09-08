#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import argparse
import glob
import re

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import openpmd_api as io
import pandas as pd

# collect final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
final_beam = series.iterations[last_step].particles["beam"].to_df()

# scaling to figure units
beta = 0.9986935469557160
gamma = 19.569511835591836
ErMeV = 0.510998950
scalex = 1.0e3 # x [m] -> x [mm]
scalet = -beta*1.0e3 # ct [m] -> z [m]
scalepx = beta*gamma*ErMeV # px/p0 -> px [MeV/c]
scalept = -gamma*ErMeV # pt/p0 -> pz [MeV/c]

# beam phase space scatter plots
num_plots_per_row = 2
fig, axs = plt.subplots(
    1, num_plots_per_row, figsize=(10, 4.0)
)

# plot final x-px distribution
ax = axs[(0)]
ax.scatter(
    final_beam.position_x.multiply(scalex),
    final_beam.momentum_x.multiply(scalepx),
    s = 0.5,
)
ax.tick_params(labelsize=12)
ax.set_xlabel(r"$x$ [mm]", fontsize=18)
ax.set_ylabel(r"$p_x$ [MeV/c]", fontsize=18)

# plot final z-pz distribution
ax = axs[(1)]
ax.scatter(
    final_beam.position_t.multiply(scalet),
    final_beam.momentum_t.multiply(scalept),
    s = 0.5,
)
ax.tick_params(labelsize=12)
ax.set_xlabel(r"$z$ [mm]", fontsize=18)
ax.set_ylabel(r"$p_z$ [MeV/c]", fontsize=18)

# store plots
fig.tight_layout()
plt.savefig("expanding_scatter.png")
