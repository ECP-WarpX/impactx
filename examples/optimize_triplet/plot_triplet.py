#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np

# options to run this script
parser = argparse.ArgumentParser(description="Plot the quadrupole triplet benchmark.")
parser.add_argument(
    "--save-png", action="store_true", help="non-interactive run: save to PNGs"
)
args = parser.parse_args()

# read reduced diagnostics
rdc_name = "diags/reduced_beam_characteristics.0"
if os.path.exists(rdc_name):
    data = np.loadtxt(rdc_name, skiprows=1)
else:  # OpenMP
    data = np.loadtxt(rdc_name + ".0", skiprows=1)
s = data[:, 1]
beta_x = data[:, 33]
beta_y = data[:, 34]

xMin = 0.0
xMax = 8.6
yMin = 0.0
yMax = 65.0

# Plotting
plt.figure(figsize=(10, 6))
plt.xscale("linear")
plt.yscale("linear")
plt.xlim([xMin, xMax])
# plt.ylim([yMin, yMax])
plt.xlabel("s (m)", fontsize=30)
plt.ylabel("CS Twiss beta (m)", fontsize=30)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.grid(True)

# Plot the data
plt.plot(s, beta_x, "b", label="Horizontal", linewidth=2, linestyle="solid")

plt.plot(s, beta_y, "r", label="Vertical", linewidth=2, linestyle="solid")

# Show plot
plt.legend(fontsize=20)

plt.tight_layout()
if args.save_png:
    plt.savefig("triplet.png")
else:
    plt.show()
