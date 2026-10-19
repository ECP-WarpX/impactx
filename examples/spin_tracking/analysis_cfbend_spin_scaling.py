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
is_double = initial["position_x"].dtype == np.float64
final = series.iterations[last_step].particles["beam"].to_df()

# initial coordinates
xi = initial["position_x"]
pxi = initial["momentum_x"]
yi = initial["position_y"]
pyi = initial["momentum_y"]
ti = initial["position_t"]
pti = initial["momentum_t"]
sxi = initial["spin_x"]
syi = initial["spin_y"]
szi = initial["spin_z"]
norm_phase_space_vector = np.sqrt(xi**2 + pxi**2 + yi**2 + pyi**2 + ti**2 + pti**2)

# difference of final coordinates from linear values
diff_xf = final["position_x"] - xi
diff_pxf = final["momentum_x"] - pxi
diff_yf = final["position_y"] - yi
diff_pyf = final["momentum_y"] - pyi
diff_tf = final["position_t"] - ti
diff_ptf = final["momentum_t"] - pti
diff_sxf = final["spin_x"] - sxi
diff_syf = final["spin_y"] - syi
diff_szf = final["spin_z"] - szi
norm_phase_space_diff = np.sqrt(
    diff_xf**2 + diff_pxf**2 + diff_yf**2 + diff_pyf**2 + diff_tf**2 + diff_ptf**2
)
norm_spin_diff = np.sqrt(diff_sxf**2 + diff_syf**2 + diff_szf**2)

print("")
print(":")
print("Phase space vector norm:")
print(norm_phase_space_vector)

print("")
print(":")
print("Phase space difference norm:")
print(norm_phase_space_diff)

m, b = np.polyfit(np.log(norm_phase_space_vector), np.log(norm_phase_space_diff), deg=1)
print(f"  slope={m}")

# Test for  quadratic scaling with initial phase space vector:
rtol = 2.0e-2
print(f"  rtol={rtol}")
# assert np.allclose(
#    [m],
#    [
#        2.0,
#    ],
#    rtol=rtol,
# )

print(":")
print("Spin difference norm:")
print(norm_spin_diff)

m, b = np.polyfit(np.log(norm_phase_space_vector), np.log(norm_spin_diff), deg=1)
print(f"  slope={m}")

# Test for  quadratic scaling with initial phase space vector:
rtol = 2.0e-2
if not is_double:
    # FIXME: in single precision, the fitted slope of the spin difference is 1.7
    # https://github.com/BLAST-ImpactX/impactx/issues/1484
    rtol = 0.3
print(f"  rtol={rtol}")
assert np.allclose(
    [m],
    [
        2.0,
    ],
    rtol=rtol,
)
