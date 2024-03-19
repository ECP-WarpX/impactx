#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Ryan Sandberg
# License: BSD-3-Clause-LBNL
#

import impactx_read_utilities as ix_read
import numpy as np
import openpmd_api as io
import transformation_utilities as coord

print("Initial Beam:")

beam_impactx_surrogate_series = io.Series(
    "diags/openPMD/monitor.bp", io.Access.read_only
)
impactx_surrogate_steps = list(beam_impactx_surrogate_series.iterations)
impactx_surrogate_ref_particle = ix_read.read_time_series("diags/ref_particle.*")

num_particles = int(1e5)
millimeter = 1.0e3
micron = 1.0e6
step = 1
beam_at_step = beam_impactx_surrogate_series.iterations[step].particles["beam"].to_df()
ref_part_step = impactx_surrogate_ref_particle.loc[step]
ref_part = myref = coord.MyRefPart(
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
x1, y1, z1, px1, py1, pz1 = coord.to_lab_t_from_ref_part_t(
    ref_part, dx_t, dy_t, dz, dpx_t, dpy_t, dpz
)

sigx = x1.std()
sigy = y1.std()
sigz = z1.std()
sigpx = px1.std()
sigpy = py1.std()
sigpz = pz1.std()
print(f"sig x={sigx:.2e}, sig y={sigy:.2e}, sig z={sigz:.2e}")
print(f"sig px={sigpx:.2e}, sig py={sigpy:.2e}, sig pz={sigpz:.2e}")

atol = 0.0  # ignored
rtol = 2 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigz, sigpx, sigpy, sigpz],
    [1.46e-3, 1.46e-3, 1.0e-3, 10, 10, 200],
    rtol=rtol,
    atol=atol,
)
