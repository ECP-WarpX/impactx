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

# compare number of particles
num_particles = 10000
assert num_particles == len(initial)
assert num_particles == len(final)

# numerical parameters
bucket_length = 0.2
half_bucket_length = bucket_length / 2.0

# join tables on particle ID, so we can compare the same particle initial->final
beam_joined = final.join(initial, lsuffix="_final", rsuffix="_initial")

beam_joined["ttest"] = np.fmod(
    beam_joined["position_t_initial"] + half_bucket_length, bucket_length
)
beam_joined["t_predicted"] = (
    np.fmod(beam_joined["ttest"] + bucket_length, bucket_length) - half_bucket_length
)
beam_joined["dt"] = (beam_joined["position_t_final"] - beam_joined["t_predicted"]).abs()

print(beam_joined["t_predicted"].max())
print(beam_joined["t_predicted"].min())
print(beam_joined["position_t_final"].max())
print(beam_joined["position_t_final"].min())

# particle-wise comparison of t & t_predicted:
atol = 2.0e-9 if is_double else 3.0e-7
rtol = 0.0  # large number
print()
print(f"  atol={atol} (ignored: rtol~={rtol})")

print(f"  dt_max={beam_joined['dt'].max()}")
assert np.allclose(beam_joined["dt"], 0.0, rtol=rtol, atol=atol)
