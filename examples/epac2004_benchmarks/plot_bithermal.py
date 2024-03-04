#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import argparse
from math import pi

import matplotlib.pyplot as plt
import numpy as np
import openpmd_api as io

# options to run this script
parser = argparse.ArgumentParser(description="Plot the Bithermal benchmark.")
parser.add_argument(
    "--save-png", action="store_true", help="non-interactive run: save to PNGs"
)
args = parser.parse_args()


# initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial_beam = series.iterations[1].particles["beam"].to_df()
final_beam = series.iterations[last_step].particles["beam"].to_df()

# Constants
w1 = 0.95
w2 = 0.05
bg = 0.0146003
Min = 0.0
Max = 0.025
Np = 100000001
n = 300


# Function for radius calculation
def r(x, y, z):
    return np.sqrt(x**2 + y**2 + z**2)


# Calculate radius and bin data
initial_radii = r(
    bg * initial_beam["position_t"],
    initial_beam["position_x"],
    initial_beam["position_y"],
)
initial_hist, bin_edges = np.histogram(initial_radii, bins=n, range=(Min, Max))
initial_bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

final_radii = r(
    bg * final_beam["position_t"], final_beam["position_x"], final_beam["position_y"]
)
final_hist, _ = np.histogram(final_radii, bins=n, range=(Min, Max))

# dr (m)
initial_r = initial_hist / (Np * (bin_edges[1] - bin_edges[0]))
final_r = final_hist / (Np * (bin_edges[1] - bin_edges[0]))

# Plotting
plt.figure(figsize=(10, 6))
plt.xscale("linear")
plt.yscale("log")
plt.xlim([Min, Max])
plt.ylim([0.1, 1.6e6])
plt.xlabel("r (m)", fontsize=30)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.grid(True)

# Plot the data
plt.plot(
    initial_bin_centers,
    initial_r / (4.0 * pi * (initial_bin_centers) ** 2),
    label="Initial beam",
    linewidth=2,
)
plt.plot(
    initial_bin_centers,
    final_r / (4.0 * pi * (initial_bin_centers) ** 2),
    label="Final beam",
    linewidth=2,
    linestyle="dotted",
)

# Show plot
plt.legend(fontsize=20)

plt.tight_layout()
if args.save_png:
    plt.savefig("bithermal.png")
else:
    plt.show()
