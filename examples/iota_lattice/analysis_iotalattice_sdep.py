#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
import numpy as np
import openpmd_api as io
from scipy.stats import moment


def get_moments(beam):
    """Calculate mean and std dev of functions defining the IOTA invariants
    Returns
    -------
    meanH, sigH, meanI, sigI
    """
    meanH = np.mean(beam["H"])
    sigH = moment(beam["H"], moment=2) ** 0.5
    meanI = np.mean(beam["I"])
    sigI = moment(beam["I"], moment=2) ** 0.5

    return (meanH, sigH, meanI, sigI)


# initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial = series.iterations[1].particles["beam"].to_df()
final = series.iterations[last_step].particles["beam"].to_df()

# compare number of particles
num_particles = 10000
assert num_particles == len(initial)
assert num_particles == len(final)

print("Initial Beam:")
meanH, sigH, meanI, sigI = get_moments(initial)
print(f"  meanH={meanH:e} sigH={sigH:e} meanI={meanI:e} sigI={sigI:e}")

atol = 0.0  # a big number
rtol = 1.6 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [meanH, sigH, meanI, sigI],
    [7.263202e-02, 4.454371e-02, 9.288060e-02, 8.211506e-02],
    rtol=rtol,
    atol=atol,
)


print("")
print("Final Beam:")
meanH, sigH, meanI, sigI = get_moments(final)
print(f"  meanH={meanH:e} sigH={sigH:e} meanI={meanI:e} sigI={sigI:e}")

atol = 0.0  # a big number
rtol = 1.6 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [meanH, sigH, meanI, sigI],
    [7.263202e-02, 4.454371e-02, 9.288060e-02, 8.211506e-02],
    rtol=rtol,
    atol=atol,
)

# join tables on particle ID, so we can compare the same particle initial->final
beam_joined = final.join(initial, lsuffix="_final", rsuffix="_initial")
# add new columns: dH and dI
beam_joined["dH"] = (beam_joined["H_initial"] - beam_joined["H_final"]).abs()
beam_joined["dI"] = (beam_joined["I_initial"] - beam_joined["I_final"]).abs()
# print(beam_joined)

# particle-wise comparison of H & I initial to final
Hrms = np.sqrt(sigH**2 + meanH**2)
Irms = np.sqrt(sigI**2 + meanI**2)

atol = 4.2e-3 * Hrms
rtol = 0.0  # large number
print()
print(f"  atol={atol} (ignored: rtol~={rtol})")
print(f"  dH_max={beam_joined['dH'].max()}")
assert np.allclose(beam_joined["dH"], 0.0, rtol=rtol, atol=atol)

atol = 5.7e-3 * Irms
rtol = 0.0
print()
print(f"  atol={atol} (ignored: rtol~={rtol})")
print(f"  dI_max={beam_joined['dI'].max()}")
assert np.allclose(beam_joined["dI"], 0.0, rtol=rtol, atol=atol)
