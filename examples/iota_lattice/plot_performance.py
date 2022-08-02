#!/usr/bin/env python3
#
# Copyright 2022 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import argparse

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy.stats import moment

matplotlib.rcParams.update({
    'font.size': 14
})


# options to run this script
parser = argparse.ArgumentParser(description='Plot the FODO benchmark.')
parser.add_argument('--save-png', action="store_true",
    help='non-interactive run: save to PNGs')
args = parser.parse_args()


# measured performance data
d = {
    'measurement': [
        "1 core",
        "8 cores",
        "16 cores",
        "32 cores",
        "1 A100",
        "2 A100"
    ],
    # in seconds
    'timing': [
        621.5,
        180.0,
        147.1,
        130.3,
        40.51,
        23.31
    ]
}
df = pd.DataFrame(data=d)
#print(df)

# calc speedup
ref_loc = df.loc[df['measurement'] == "16 cores"]
print(ref_loc["timing"].values[0])
df['speedup'] = ref_loc["timing"].values[0]/df["timing"]
df['label_padding'] = df["timing"].multiply(0.)

print(df)

# print beam position over s
f = plt.figure(figsize=(7, 2))
ax = f.gca()

p1 = ax.bar(df.measurement, df.speedup)
ax.bar_label(
    p1,
    labels=[f'{x:,.1f}x' for x in df.speedup],
#    padding = df.label_padding,
#    label_type = 'center'
)

# 1x reference line
ax.hlines(
    y=1.0,
    xmin=-1, xmax=len(df)-0.5,
    colors='grey', linestyles=':', lw=2,
    label='1x'
)

#ax.set_xlabel(r"$z$ [m]")
ax.set_ylabel("speedup")
ax.set_xlim([-0.5, len(df)-0.5])
ax.set_ylim([0, 8])
ax.set_yticks([0, 2, 4, 6, 8])
#ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.tight_layout()
if args.save_png:
    plt.savefig("iotalattice_performance.eps")
else:
    plt.show()
