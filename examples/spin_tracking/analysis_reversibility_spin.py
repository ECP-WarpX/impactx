#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#


import numpy as np
import openpmd_api as io

# initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial = series.iterations[1].particles["beam"].to_df()
beam_final = series.iterations[last_step].particles["beam"]
final = beam_final.to_df()

sxi = initial["spin_x"]
syi = initial["spin_y"]
szi = initial["spin_z"]

sxf = final["spin_x"]
syf = final["spin_y"]
szf = final["spin_z"]

dspin2 = (sxf - sxi) ** 2 + (syf - syi) ** 2 + (szf - szi) ** 2
dspin = np.sqrt(dspin2)
dspinmax = dspin.max()

print("Change in the spin:")
print("||delta s||_max", dspinmax)

is_double = np.dtype(sxi.dtype) == np.dtype(np.float64)
atol = 1.5e-9
if not is_double:
    # FIXME: the spin does not return in single precision, ||delta s||_max is up to
    # 0.84 of the maximum of 2 for unit vectors
    # https://github.com/BLAST-ImpactX/impactx/issues/1484
    atol = 1.7
print(f"  atol={atol}")

assert np.allclose(
    [dspinmax],
    [
        0.0,
    ],
    atol=atol,
)
