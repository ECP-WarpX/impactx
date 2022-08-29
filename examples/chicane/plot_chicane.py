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
    sigx = moment(beam["x"], moment=2) ** 0.5  # variance -> std dev.
    sigpx = moment(beam["px"], moment=2) ** 0.5
    sigy = moment(beam["y"], moment=2) ** 0.5
    sigpy = moment(beam["py"], moment=2) ** 0.5
    sigt = moment(beam["t"], moment=2) ** 0.5
    sigpt = moment(beam["pt"], moment=2) ** 0.5

    epstrms = beam.cov(ddof=0)
    emittance_x = (sigx**2 * sigpx**2 - epstrms["x"]["px"] ** 2) ** 0.5
    emittance_y = (sigy**2 * sigpy**2 - epstrms["y"]["py"] ** 2) ** 0.5
    emittance_t = (sigt**2 * sigpt**2 - epstrms["t"]["pt"] ** 2) ** 0.5

    return (sigx, sigy, sigt, emittance_x, emittance_y, emittance_t)


def read_file(file_pattern):
    for filename in glob.glob(file_pattern):
        df = pd.read_csv(filename, delimiter=r"\s+")
        if "step" not in df.columns:
            step = int(re.findall(r"[0-9]+", filename)[0])
            df["step"] = step
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
    )  # .set_index('id')


# options to run this script
parser = argparse.ArgumentParser(description="Plot the chicane benchmark.")
parser.add_argument(
    "--save-png", action="store_true", help="non-interactive run: save to PNGs"
)
args = parser.parse_args()


# initial/final beam on rank zero
beam = read_time_series("diags/beam_[0-9]*.*")
ref_particle = read_time_series("diags/ref_particle.*")
# print(beam)
# print(ref_particle)

# scaling to units
millimeter = 1.0e3  # m->mm
# for "t": the time coordinate is scaled by c, and therefore has units of length (m) by default, so we can label the axis ct (mm)
mrad = 1.0e3  # ImpactX uses "static units": momenta are normalized by the magnitude of the momentum of the reference particle p0: px/p0 (rad)
# mm_mrad = 1.e6
nm_rad = 1.0e9


# select a single particle by id
# particle_42 = beam[beam["id"] == 42]
# print(particle_42)


# steps & corresponding z
steps = beam.step.unique()
steps.sort()
# print(f"steps={steps}")

z = list(
    map(lambda step: ref_particle[ref_particle["step"] == step].z.values[0], steps)
)
x = list(
    map(lambda step: ref_particle[ref_particle["step"] == step].x.values[0], steps)
)
# print(f"z={z}")


# beam transversal size & emittance over steps
moments = list(map(lambda step: (step, get_moments(beam[beam["step"] == step])), steps))
# print(moments)
sigx = list(map(lambda step_val: step_val[1][0] * millimeter, moments))
sigt = list(map(lambda step_val: step_val[1][2] * millimeter, moments))
emittance_x = list(map(lambda step_val: step_val[1][3] * nm_rad, moments))
emittance_t = list(map(lambda step_val: step_val[1][5] * nm_rad, moments))

# print(sigx, sigt)


# print beam transversal size over steps
f = plt.figure(figsize=(9, 4.8))

f, axs = plt.subplots(
    2, 1, figsize=(9, 4.8), sharex=True, gridspec_kw={"height_ratios": [1, 2]}
)
ax0 = axs[0]
im_xz = ax0.plot(z, x, "--", lw=3, label=r"$x$")
ax0.legend(loc="upper right")
ax0.set_ylim([0, None])
ax0.set_ylabel(r"$x$ [m]")

ax1 = axs[1]
im_sigx = ax1.plot(z, sigx, label=r"$\sigma_x$")
im_sigt = ax1.plot(z, sigt, label=r"$\sigma_t$")
ax2 = ax1.twinx()
ax2._get_lines.prop_cycler = ax1._get_lines.prop_cycler
im_emittance_x = ax2.plot(z, emittance_x, ":", label=r"$\epsilon_x$")
im_emittance_t = ax2.plot(z, emittance_t, ":", label=r"$\epsilon_t$")

ax1.legend(
    handles=im_sigx + im_sigt + im_emittance_x + im_emittance_t, loc="upper right"
)
ax1.set_xlabel(r"$z$ [m]")
ax1.set_ylabel(r"$\sigma_{x,t}$ [mm]")
# ax2.set_ylabel(r"$\epsilon_{x,y}$ [mm-mrad]")
ax2.set_ylabel(r"$\epsilon_{x,t}$ [nm]")
ax1.set_ylim([0, None])
ax2.set_ylim([0, None])
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.tight_layout()
if args.save_png:
    plt.savefig("chicane_sigma.png")
else:
    plt.show()


# beam transversal scatter plot over steps
plot_ever_nstep = 25  # this is lattice.nslice in inputs
# entry and exit dipedge are part of each bending element
# lattice.elements = sbend1 dipedge1 drift1 dipedge2 sbend2 drift2 sbend2
#                    dipedge2 drift1 dipedge1 sbend1 drift3
# lattice.nslice = 25
print_steps = [0, 26, 51, 77, 102, 128, 153, 179, 204]
num_plots_per_row = len(print_steps)
fig, axs = plt.subplots(
    4, num_plots_per_row, figsize=(9, 6.4), sharex="row", sharey="row"
)

ncol_ax = -1
for step in print_steps:
    # plot initial distribution & at exit of each element
    print(step)
    ncol_ax += 1

    # t-pt
    ax = axs[(0, ncol_ax)]
    beam_at_step = beam[beam["step"] == step]
    ax.scatter(
        beam_at_step.t.multiply(millimeter), beam_at_step.pt.multiply(mrad), s=0.01
    )
    ax.set_xlabel(r"$ct$ [mm]")
    z_unit = ""
    if ncol_ax == num_plots_per_row - 1:
        z_unit = " [m]"
    ax.set_title(f"$z={z[step]:.1f}${z_unit}")

    # x-px
    ax = axs[(1, ncol_ax)]
    beam_at_step = beam[beam["step"] == step]
    ax.scatter(
        beam_at_step.x.multiply(millimeter), beam_at_step.px.multiply(mrad), s=0.01
    )
    ax.set_xlabel(r"$x$ [mm]")

    # t-x
    ax = axs[(2, ncol_ax)]
    beam_at_step = beam[beam["step"] == step]
    ax.scatter(
        beam_at_step.t.multiply(millimeter), beam_at_step.x.multiply(millimeter), s=0.01
    )
    ax.set_xlabel(r"$ct$ [mm]")

    # t-px
    ax = axs[(3, ncol_ax)]
    beam_at_step = beam[beam["step"] == step]
    ax.scatter(
        beam_at_step.t.multiply(millimeter), beam_at_step.px.multiply(mrad), s=0.01
    )
    ax.set_xlabel(r"$ct$ [mm]")

axs[(0, 0)].set_ylabel(r"$p_t$ [mrad]")
axs[(1, 0)].set_ylabel(r"$p_x$ [mrad]")
axs[(2, 0)].set_ylabel(r"$x$ [mm]")
axs[(3, 0)].set_ylabel(r"$p_x$ [mrad]")
plt.tight_layout()
if args.save_png:
    plt.savefig("chicane_scatter.png")
else:
    plt.show()
