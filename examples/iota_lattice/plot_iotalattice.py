#!/usr/bin/env python3
#
# Copyright 2022 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import argparse
import glob
import re

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy.stats import moment

matplotlib.rcParams.update({
    'font.size': 14
})


def get_moments(beam):
    """Calculate standard deviations of beam position & momenta
    and emittance values

    Returns
    -------
    sigx, sigy, sigt, emittance_x, emittance_y, emittance_t
    """
    sigx = moment(beam["x"], moment=2)**0.5  # variance -> std dev.
    sigpx = moment(beam["px"], moment=2)**0.5
    sigy = moment(beam["y"], moment=2)**0.5
    sigpy = moment(beam["py"], moment=2)**0.5
    sigt = moment(beam["t"], moment=2)**0.5
    sigpt = moment(beam["pt"], moment=2)**0.5

    epstrms = beam.cov(ddof=0)
    emittance_x = (sigx**2 * sigpx**2 - epstrms["x"]["px"]**2)**0.5
    emittance_y = (sigy**2 * sigpy**2 - epstrms["y"]["py"]**2)**0.5
    emittance_t = (sigt**2 * sigpt**2 - epstrms["t"]["pt"]**2)**0.5

    return (
        sigx, sigy, sigt,
        emittance_x, emittance_y, emittance_t)


def read_file(file_pattern):
    for filename in glob.glob(file_pattern):
        df = pd.read_csv(filename, delimiter=r"\s+")
        if 'step' not in df.columns:
            step = int(re.findall(r"[0-9]+", filename)[0])
            df['step'] = step
        yield df


def read_time_series(file_pattern):
    """Read in all CSV files from each MPI rank (and potentially OpenMP
    thread). Concatenate into one Pandas dataframe.

    Returns
    -------
    pandas.DataFrame
    """
    return pd.concat(
        read_file(file_pattern),
        axis=0,
        ignore_index=True,
    ) #.set_index('id')


# options to run this script
parser = argparse.ArgumentParser(description='Plot the FODO benchmark.')
parser.add_argument('--save-png', action="store_true",
    help='non-interactive run: save to PNGs')
args = parser.parse_args()


# IMPACT-Z data
iz_x = pd.read_csv(
    "impactz/IMPACTZ_outputX.txt",
    delimiter=r"\s+",
    skiprows=1,
    names=['s', 'mean_x', 'sigma_x', 'mean_px', 'sigma_px', 'Twiss_alpha', 'emittance_x']
).set_index('s')
iz_y = pd.read_csv(
    "impactz/IMPACTZ_outputY.txt",
    delimiter=r"\s+",
    skiprows=1,
    names=['s', 'mean_y', 'sigma_y', 'mean_py', 'sigma_py', 'Twiss_alpha', 'emittance_y']
).set_index('s').drop(columns=['Twiss_alpha'])
iz = iz_x.join(
    iz_y
)

# initial/final beam on rank zero
beam = read_time_series("diags/beam_[0-9]*.*")
ref_particle = read_time_series("diags/ref_particle.*")
#print(beam)
#print(ref_particle)

# scaling to units
millimeter = 1.e3  # m->mm
mrad = 1.e3  # ImpactX uses "static units": momenta are normalized by the magnitude of the momentum of the reference particle p0: px/p0 (rad)
#mm_mrad = 1.e6
nm_rad = 1.e9


# select a single particle by id
#particle_42 = beam[beam["id"] == 42]
#print(particle_42)


# steps & corresponding z
steps = beam.step.unique()
steps.sort()
#print(f"steps={steps}")

s = list(map(
    lambda step: ref_particle[ref_particle["step"] == step].s.values[0],
    steps
))
x = list(map(
    lambda step: ref_particle[ref_particle["step"] == step].x.values[0],
    steps
))
z = list(map(
    lambda step: ref_particle[ref_particle["step"] == step].z.values[0],
    steps
))
#print(f"z={z}")


# beam transversal size & emittance over steps
moments = list(map(
    lambda step: (step, get_moments(beam[beam["step"] == step])),
    steps
))
#print(moments)
sigx = list(map(lambda step_val: step_val[1][0] * millimeter, moments))
sigy = list(map(lambda step_val: step_val[1][1] * millimeter, moments))
emittance_x = list(map(lambda step_val: step_val[1][3] * nm_rad, moments))
emittance_y = list(map(lambda step_val: step_val[1][4] * nm_rad, moments))

#print(sigx, sigy)

# print beam position over s
f = plt.figure(figsize=(7, 2))
ax1 = f.gca()
im_x = ax1.plot(z, x, label=r'$\z-x')

ax1.set_xlabel(r"$z$ [m]")
ax1.set_ylabel(r"$x$ [m]")
#ax1.set_xlim([0, 40])
ax1.set_ylim([-8.8, 0.5])
ax1.set_yticks([-8, -4, 0])
#ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.tight_layout()
if args.save_png:
    plt.savefig("iotalattice_position.eps")
else:
    plt.show()


# print beam transversal size over steps
f = plt.figure(figsize=(7, 2))
ax1 = f.gca()

# IMPACT-Z - ImpactX
im_sigx_iz = ax1.plot(
    iz.index,
    iz.sigma_x.multiply(millimeter),
    label=r'IMPACT-Z $\sigma_x$',
    lw='5',
    color='lightsteelblue'
)
im_sigx = ax1.plot(s, sigx, label=r'$\sigma_x$')

im_sigy_iz = ax1.plot(
    iz.index,
    iz.sigma_y.multiply(millimeter),
    label=r'IMPACT-Z $\sigma_y$',
    lw='5',
    color='peachpuff'
)
im_sigy = ax1.plot(s, sigy, label=r'$\sigma_y$')

ax1.legend(
    handles=im_sigx+im_sigy,
    loc='center right',
    bbox_to_anchor=(0.83, 0.5)
)
ax1.set_xlabel(r"distance $s$ [m]")
ax1.set_ylabel("beam size\n[mm]")
ax1.set_xlim([0, 40])
ax1.set_ylim([0, None])
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.tight_layout()
if args.save_png:
    plt.savefig("iotalattice_sigma.eps")
else:
    plt.show()
