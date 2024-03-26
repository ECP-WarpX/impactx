#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Ryan Sandberg, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import impactx_read_utilities as ix_plt
import openpmd_api as io
import transformation_utilities as coord
from matplotlib import pyplot as plt

######## plot phase spaces ###########
beam_series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
beam_steps = list(beam_series.iterations)
ref_series = ix_plt.read_time_series("diags/ref_particle.*")

millimeter = 1.0e3
micron = 1.0e6
step = 1
beam_at_step = beam_series.iterations[step].particles["beam"].to_df()
ref_part_step = ref_series.loc[step]
ref_part = coord.MyRefPart(
    ref_part_step["x"],
    ref_part_step["y"],
    ref_part_step["z"],
    ref_part_step["px"],
    ref_part_step["py"],
    ref_part_step["pz"],
    ref_part_step["pt"],
)

dx_s = beam_at_step["position_x"]
dy_s = beam_at_step["position_y"]
dt = beam_at_step["position_t"]
dpx_s = beam_at_step["momentum_x"]
dpy_s = beam_at_step["momentum_y"]
dpt = beam_at_step["momentum_t"]

dx_t, dy_t, dz, dpx_t, dpy_t, dpz = coord.to_t_from_s(
    ref_part, dx_s, dy_s, dt, dpx_s, dpy_s, dpt
)
x, y, z, px, py, pz = coord.to_global_t_from_ref_part_t(
    ref_part, dx_t, dy_t, dz, dpx_t, dpy_t, dpz
)

fig, axT = plt.subplots(1, 3, figsize=(8, 4))

ax = axT[0]
ax.hist2d(x * millimeter, y * millimeter, bins=200)
ax.set_xlabel("x (mm)")
ax.set_ylabel("y (mm)")

ax = axT[1]
ax.hist2d(px, py, bins=200)
ax.set_xlabel("px")
ax.set_ylabel("py")
ax = axT[2]
ax.hist2d(z * millimeter, pz, bins=200)
ax.set_xlabel("z (mm)")
ax.set_ylabel("pz")
ax.ticklabel_format(axis="y", style="sci", scilimits=(-1, 1))

plt.tight_layout()
plt.savefig("phase_space.png")
