#!/usr/bin/env python3
#
# Copyright 2022 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import argparse
import glob
import re

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy.stats import moment


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


# print beam transversal size over steps
f = plt.figure(figsize=(9, 4.8))
ax1 = f.gca()
im_sigx = ax1.plot(z, sigx, label=r'$\sigma_x$')
im_sigy = ax1.plot(z, sigy, label=r'$\sigma_y$')
ax2 = ax1.twinx()
ax2._get_lines.prop_cycler = ax1._get_lines.prop_cycler
im_emittance_x = ax2.plot(z, emittance_x, ':', label=r'$\epsilon_x$')
im_emittance_y = ax2.plot(z, emittance_y, ':', label=r'$\epsilon_y$')

ax1.legend(
    handles=im_sigx+im_sigy+im_emittance_x+im_emittance_y,
    loc='lower center'
)
ax1.set_xlabel(r"$z$ [m]")
ax1.set_ylabel(r"$\sigma_{x,y}$ [mm]")
#ax2.set_ylabel(r"$\epsilon_{x,y}$ [mm-mrad]")
ax2.set_ylabel(r"$\epsilon_{x,y}$ [nm]")
ax2.set_ylim([1.5, 2.5])
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.tight_layout()
if args.save_png:
    plt.savefig("fodo_sigma.png")
else:
    plt.show()


# beam transversal scatter plot over steps
plot_ever_nstep = 25  # this is lattice.nslice in inputs
num_plots_per_row = len(steps) // plot_ever_nstep + 1
fig, axs = plt.subplots(
    3, num_plots_per_row,
    figsize=(9, 4.8),
    sharex='row', sharey='row'
)

ncol_ax = -1
for step in steps:
    # plot initial distribution & at exit of each element
    if step != 0 and step % plot_ever_nstep != 0:
        continue
    print(step)
    ncol_ax += 1

    # x-y
    ax = axs[(0, ncol_ax)]
    beam_at_step = beam[beam["step"] == step]
    ax.scatter(
        beam_at_step.x.multiply(millimeter),
        beam_at_step.y.multiply(millimeter),
        s=0.01
    )

    ax.set_title(f"$z={z[step]}$ [m]")
    ax.set_xlabel(r"$x$ [mm]")

    # x-px
    ax = axs[(1, ncol_ax)]
    beam_at_step = beam[beam["step"] == step]
    ax.scatter(
        beam_at_step.x.multiply(millimeter),
        beam_at_step.px.multiply(mrad),
        s=0.01
    )
    ax.set_xlabel(r"$x$ [mm]")

    # y-py
    ax = axs[(2, ncol_ax)]
    beam_at_step = beam[beam["step"] == step]
    ax.scatter(
        beam_at_step.y.multiply(millimeter),
        beam_at_step.py.multiply(mrad),
        s=0.01
    )
    ax.set_xlabel(r"$y$ [mm]")

axs[(0, 0)].set_ylabel(r"$y$ [mm]")
axs[(1, 0)].set_ylabel(r"$p_x$ [mrad]")
axs[(2, 0)].set_ylabel(r"$p_y$ [mrad]")
plt.tight_layout()
if args.save_png:
    plt.savefig("fodo_scatter.png")
else:
    plt.show()
