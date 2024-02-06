#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Ryan Sandberg, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import argparse
import glob
import re

from matplotlib import pyplot as plt
import numpy as np
import openpmd_api as io
import pandas as pd
from scipy.constants import c, e, m_e


def read_all_files(file_pattern):
    """Read in all CSV files from each MPI rank (and potentially OpenMP
    thread). Concatenate into one Pandas dataframe.
    Returns
    -------
    pandas.DataFrame
    """
    return pd.concat(
        (
            pd.read_csv(filename, delimiter=r"\s+")
            for filename in glob.glob(file_pattern)
        ),
        axis=0,
        ignore_index=True,
    ).set_index("id")


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


from enum import Enum


class TCoords(Enum):
    REF = 1
    GLOBAL = 2


def to_t(
    ref_pz, ref_pt, data_arr_s, ref_z=None, coord_type=TCoords.REF
):  # x, y, t, dpx, dpy, dpt):
    """Change to fixed t coordinates

    Parameters
    ---
    ref_pz: float, reference particle momentum in z
    ref_pt: float, reference particle pt = -gamma
    data_arr_s: Nx6 array-like structure containing fixed-s particle coordinates
    ref_z: if transforming to global coordinates
    coord_type: TCoords enum, (default is in ref coordinates) whether to get particle data relative to reference coordinate or in the global frame

    Notes
    """
    if type(data_arr_s) is pd.core.frame.DataFrame:
        coordinate_columns = [
            "position_x",
            "position_y",
            "position_t",
            "momentum_x",
            "momentum_y",
            "momentum_t",
        ]
        assert all(
            val in data_arr_s.columns for val in coordinate_columns
        ), f"data_arr_s must have columns {' '.join(coordinate_columns)}"
        x, y, t, dpx, dpy, dpt = data_arr_s[coordinate_columns].to_numpy().T
        x = data_arr_s["position_x"]
        y = data_arr_s["position_y"]
        t = data_arr_s["position_t"]
        dpx = data_arr_s["momentum_x"]
        dpy = data_arr_s["momentum_y"]
        dpt = data_arr_s["momentum_t"]

    elif type(data_arr_s) is np.ndarray:
        assert (
            data_arr_s.shape[1] == 6
        ), f"data_arr_s.shape={data_arr_s.shape} but data_arr_s must be an Nx6 array"
        x, y, t, dpx, dpy, dpt = data_arr_s.T
    else:
        raise Exception(
            f"Incompatible input type {type(data_arr_s)} for data_arr_s, must be pandas DataFrame or Nx6 array-like object"
        )
    x += ref_pz * dpx * t / (ref_pt + ref_pz * dpt)
    y += ref_pz * dpy * t / (ref_pt + ref_pz * dpt)
    pz = np.sqrt(
        -1 + (ref_pt + ref_pz * dpt) ** 2 - (ref_pz * dpx) ** 2 - (ref_pz * dpy) ** 2
    )
    t *= pz / (ref_pt + ref_pz * dpt)
    if type(data_arr_s) is pd.core.frame.DataFrame:
        data_arr_s["momentum_t"] = pz - ref_pz
        dpt = data_arr_s["momentum_t"]
    else:
        dpt[:] = pz - ref_pz
    if coord_type is TCoords.REF:
        print("applying reference normalization")
        dpt /= ref_pz
    elif coord_type is TCoords.GLOBAL:
        assert (
            ref_z is not None
        ), "Reference particle z coordinate is required to transform to global coordinates"
        print("target global coordinates")
        t += ref_z
        dpx *= ref_pz
        dpy *= ref_pz
        dpt += ref_pz
    # data_arr_t = np.column_stack([xt,yt,z,dpx,dpy,dpz])
    return  # modifies data_arr_s in place


def plot_beam_df(
    beam_at_step,
    axT,
    unit=1e6,
    unit_z=1e3,
    unit_label="$\mu$m",
    unit_z_label="mm",
    alpha=1.0,
    cmap=None,
    color="k",
    size=0.1,
    t_offset=0.0,
    label=None,
    z_ticks=None,
):
    ax = axT[0][0]
    ax.scatter(
        beam_at_step.position_x.multiply(unit),
        beam_at_step.position_y.multiply(unit),
        c=color,
        s=size,
        alpha=alpha,
        cmap=cmap,
    )
    ax.set_xlabel(r"x [%s]" % unit_label)
    ax.set_ylabel(r"y [%s]" % unit_label)
    ax.axes.ticklabel_format(axis="both", style="sci", scilimits=(-2, 2))
    ###########

    ax = axT[0][1]
    ax.scatter(
        beam_at_step.position_t.multiply(unit_z) - t_offset,
        beam_at_step.position_x.multiply(unit),
        c=color,
        s=size,
        alpha=alpha,
        cmap=cmap,
    )
    ax.set_xlabel(r"%s" % unit_z_label)
    ax.set_ylabel(r"x [%s]" % unit_label)
    ax.axes.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    if z_ticks is not None:
        ax.set_xticks(z_ticks)
    ###########

    ax = axT[0][2]
    ax.scatter(
        beam_at_step.position_t.multiply(unit_z) - t_offset,
        beam_at_step.position_y.multiply(unit),
        c=color,
        s=size,
        alpha=alpha,
        cmap=cmap,
    )
    ax.set_xlabel(r"%s" % unit_z_label)
    ax.set_ylabel(r"y [%s]" % unit_label)
    ax.axes.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    if z_ticks is not None:
        ax.set_xticks(z_ticks)
    ############
    ##########
    ax = axT[1][0]
    ax.scatter(
        beam_at_step.momentum_x,
        beam_at_step.momentum_y,
        c=color,
        s=size,
        alpha=alpha,
        cmap=cmap,
    )
    ax.set_xlabel("px")
    ax.set_ylabel("py")
    ax.axes.ticklabel_format(axis="both", style="sci", scilimits=(-2, 2))
    ##########
    ax = axT[1][1]
    ax.scatter(
        beam_at_step.momentum_t,
        #         beam_at_step.position_t.multiply(unit_z)-t_offset,
        beam_at_step.momentum_x,
        c=color,
        s=size,
        alpha=alpha,
        cmap=cmap,
    )
    ax.set_xlabel("pt")
    #     ax.set_xlabel(r'%s'%unit_z_label)
    ax.set_ylabel("px")
    ax.axes.ticklabel_format(axis="both", style="sci", scilimits=(-2, 2))
    ##########
    ax = axT[1][2]
    ax.scatter(
        beam_at_step.momentum_t,
        #         beam_at_step.position_t.multiply(unit_z)-t_offset,
        beam_at_step.momentum_y,
        c=color,
        s=size,
        alpha=alpha,
        label=label,
        cmap=cmap,
    )
    if label is not None:
        ax.legend()
    #     ax.set_xlabel(r'%s'%unit_z_label)
    ax.set_xlabel("pt")
    ax.set_ylabel("py")
    ax.axes.ticklabel_format(axis="both", style="sci", scilimits=(-2, 2))
    ############
    ############
    ##########

    ax = axT[2][0]
    ax.scatter(
        beam_at_step.position_x.multiply(unit),
        beam_at_step.momentum_x,
        c=color,
        s=size,
        alpha=alpha,
        cmap=cmap,
    )
    ax.set_xlabel(r"x [%s]" % unit_label)
    ax.set_ylabel("px")
    ax.axes.ticklabel_format(axis="both", style="sci", scilimits=(-2, 2))
    ############
    ax = axT[2][1]
    ax.scatter(
        beam_at_step.position_y.multiply(unit),
        beam_at_step.momentum_y,
        c=color,
        s=size,
        alpha=alpha,
        cmap=cmap,
    )
    ax.set_xlabel(r"y [%s]" % unit_label)
    ax.set_ylabel("py")
    ax.axes.ticklabel_format(axis="both", style="sci", scilimits=(-2, 2))

    ################
    ax = axT[2][2]
    ax.scatter(
        beam_at_step.position_t.multiply(unit_z) - t_offset,
        beam_at_step.momentum_t,
        c=color,
        s=size,
        alpha=alpha,
        cmap=cmap,
    )
    ax.set_xlabel(r"%s" % unit_z_label)
    ax.set_ylabel("pt")
    ax.axes.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    if z_ticks is not None:
        ax.set_xticks(z_ticks)
    plt.tight_layout()
    # done


# options to run this script
parser = argparse.ArgumentParser(description="Plot the ML surrogate benchmark.")
parser.add_argument(
    "--save-png", action="store_true", help="non-interactive run: save to PNGs"
)
args = parser.parse_args()

impactx_surrogate_reduced_diags = read_time_series(
    "diags/reduced_beam_characteristics.*"
)
ref_gamma = np.sqrt(1 + impactx_surrogate_reduced_diags["ref_beta_gamma"] ** 2)
beam_gamma = (
    ref_gamma
    - impactx_surrogate_reduced_diags["pt_mean"]
    * impactx_surrogate_reduced_diags["ref_beta_gamma"]
)
beam_u = np.sqrt(beam_gamma**2 - 1)
emit_x = impactx_surrogate_reduced_diags["emittance_x"]
emit_nx = emit_x * beam_u
emit_y = impactx_surrogate_reduced_diags["emittance_y"]
emit_ny = emit_y * beam_u

ix_slice = [0] + [2 + 9 * i for i in range(8)]

############# plot moments ##############
fig, axT = plt.subplots(2, 2, figsize=(10, 8))
######### emittance ##########
ax = axT[0][0]
scale = 1e6
ax.plot(
    impactx_surrogate_reduced_diags["s"][ix_slice],
    emit_nx[ix_slice] * scale,
    "bo",
    label="x",
)
ax.plot(
    impactx_surrogate_reduced_diags["s"][ix_slice],
    emit_ny[ix_slice] * scale,
    "ro",
    label="y",
)
ax.legend()
ax.set_xlabel("s [m]")
ax.set_ylabel(r"emittance (mm-mrad)")
######### energy ##########
ax = axT[0][1]
scale = m_e * c**2 / e * 1e-9
ax.plot(
    impactx_surrogate_reduced_diags["s"][ix_slice],
    beam_gamma[ix_slice] * scale,
    "go",
)
ax.set_xlabel("s [m]")
ax.set_ylabel(r"mean energy (GeV)")

######### width ##########
ax = axT[1][0]
scale = 1e6
ax.plot(
    impactx_surrogate_reduced_diags["s"][ix_slice],
    impactx_surrogate_reduced_diags["sig_x"][ix_slice] * scale,
    "bo",
    label="x",
)
ax.plot(
    impactx_surrogate_reduced_diags["s"][ix_slice],
    impactx_surrogate_reduced_diags["sig_y"][ix_slice] * scale,
    "ro",
    label="y",
)
ax.legend()
ax.set_xlabel("s [m]")
ax.set_ylabel(r"beam width ($\mu$m)")

######### divergence ##########
ax = axT[1][1]
scale = 1e3
ax.semilogy(
    impactx_surrogate_reduced_diags["s"][ix_slice],
    impactx_surrogate_reduced_diags["sig_px"][ix_slice] * scale,
    "bo",
    label="x",
)
ax.semilogy(
    impactx_surrogate_reduced_diags["s"][ix_slice],
    impactx_surrogate_reduced_diags["sig_py"][ix_slice] * scale,
    "ro",
    label="y",
)
ax.legend()
ax.set_xlabel("s [m]")
ax.set_ylabel(r"divergence (mrad)")

plt.tight_layout()

if args.save_png:
    plt.savefig("lpa_ml_surrogate_moments.png")
else:
    plt.show()


######## plot phase spaces ###########
beam_impactx_surrogate_series = io.Series(
    "diags/openPMD/monitor.bp", io.Access.read_only
)
impactx_surrogate_steps = list(beam_impactx_surrogate_series.iterations)
impactx_surrogate_ref_particle = read_time_series("diags/ref_particle.*")

millimeter = 1.0e3
micron = 1.0e6

N_stage = 9
impactx_stage_end_steps = [1] + [3 + 8 * i for i in range(N_stage)]
ise = impactx_stage_end_steps

# initial

step = 1
beam_at_step = beam_impactx_surrogate_series.iterations[step].particles["beam"].to_df()
ref_part_step = impactx_surrogate_ref_particle.loc[step]
ref_u = np.sqrt(ref_part_step["pt"] ** 2 - 1)
to_t(
    ref_u,
    ref_part_step["pt"],
    beam_at_step,
    ref_z=ref_part_step["z"],
    coord_type=TCoords.GLOBAL,
)

t_offset = impactx_surrogate_ref_particle.loc[step, "t"] * micron
fig, axT = plt.subplots(3, 3, figsize=(10, 8))
fig.suptitle(f"initially, ct={impactx_surrogate_ref_particle.at[step,'t']:.2e}")

plot_beam_df(
    beam_at_step,
    axT,
    alpha=0.6,
    color="red",
    unit_z=1e6,
    unit_z_label=r"$\xi$ [$\mu$m]",
    t_offset=t_offset,
    z_ticks=[-107.3, -106.6],
)
if args.save_png:
    plt.savefig("initial_phase_spaces.png")
else:
    plt.show()

####### final ###########


stage_i = 8
step = ise[stage_i + 1]
beam_at_step = beam_impactx_surrogate_series.iterations[step].particles["beam"].to_df()
ref_part_step = impactx_surrogate_ref_particle.loc[step]
ref_u = np.sqrt(ref_part_step["pt"] ** 2 - 1)
to_t(
    ref_u,
    ref_part_step["pt"],
    beam_at_step,
    ref_z=ref_part_step["z"],
    coord_type=TCoords.GLOBAL,
)

t_offset = impactx_surrogate_ref_particle.loc[step, "t"] * micron
fig, axT = plt.subplots(3, 3, figsize=(10, 8))
fig.suptitle(f"stage {stage_i}, ct={impactx_surrogate_ref_particle.at[step,'t']:.2e}")

plot_beam_df(
    beam_at_step,
    axT,
    alpha=0.6,
    color="red",
    unit_z=1e6,
    unit_z_label=r"$\xi$ [$\mu$m]",
    t_offset=t_offset,
    z_ticks=[-107.3, -106.6],
)
if args.save_png:
    plt.savefig(f"stage_{stage_i}_phase_spaces.png")
else:
    plt.show()
