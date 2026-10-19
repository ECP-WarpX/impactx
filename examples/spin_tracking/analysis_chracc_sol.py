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
beam_initial = series.iterations[1].particles["beam"]
initial_sort = beam_initial.to_df().set_index("id")
beam_final = series.iterations[last_step].particles["beam"]
final_sort = beam_final.to_df().set_index("id")

sxi = initial_sort["spin_x"]
syi = initial_sort["spin_y"]
szi = initial_sort["spin_z"]

sxf = final_sort["spin_x"]
syf = final_sort["spin_y"]
szf = final_sort["spin_z"]

dspin2 = (sxf - sxi) ** 2 + (syf - syi) ** 2 + (szf - szi) ** 2
dspin = np.sqrt(dspin2)
dspinmax = dspin.max()

print("Change in the spin:")
print("||delta s||_max", dspinmax)

is_double = np.dtype(sxi.dtype) == np.dtype(np.float64)
atol = 2.1e-7 if is_double else 1.0e-5
print(f"  atol={atol}")

assert np.allclose(
    [dspinmax],
    [
        0.0,
    ],
    rtol=0.0,
    atol=atol,
)
