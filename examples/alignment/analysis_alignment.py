#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#


import numpy as np
import openpmd_api as io
from scipy.stats import moment, tmean


def get_moments(beam):
    """Calculate standard deviations of beam position & momenta
    and emittance values

    Returns
    -------
    meanx, meany, sigx, sigy, sigt, emittance_x, emittance_y, emittance_t
    """
    meanx = tmean(beam["position_x"])
    meany = tmean(beam["position_y"])
    sigx = moment(beam["position_x"], moment=2) ** 0.5  # variance -> std dev.
    sigpx = moment(beam["momentum_x"], moment=2) ** 0.5
    sigy = moment(beam["position_y"], moment=2) ** 0.5
    sigpy = moment(beam["momentum_y"], moment=2) ** 0.5
    sigt = moment(beam["position_t"], moment=2) ** 0.5
    sigpt = moment(beam["momentum_t"], moment=2) ** 0.5

    epstrms = beam.cov(ddof=0)
    emittance_x = (sigx**2 * sigpx**2 - epstrms["position_x"]["momentum_x"] ** 2) ** 0.5
    emittance_y = (sigy**2 * sigpy**2 - epstrms["position_y"]["momentum_y"] ** 2) ** 0.5
    emittance_t = (sigt**2 * sigpt**2 - epstrms["position_t"]["momentum_t"] ** 2) ** 0.5

    return (meanx, meany, sigx, sigy, sigt, emittance_x, emittance_y, emittance_t)


# initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial = series.iterations[1].particles["beam"].to_df()
final = series.iterations[last_step].particles["beam"].to_df()

# compare number of particles
num_particles = 100000
assert num_particles == len(initial)
assert num_particles == len(final)

print("Initial Beam:")
meanx, meany, sigx, sigy, sigt, emittance_x, emittance_y, emittance_t = get_moments(
    initial
)
print(f"  meanx={meanx:e} meany={meany:e}")
print(f"  sigx={sigx:e} sigy={sigy:e} sigt={sigt:e}")
print(
    f"  emittance_x={emittance_x:e} emittance_y={emittance_y:e} emittance_t={emittance_t:e}"
)

atol = 3.0 * num_particles ** (-0.5) * sigx
print(f" atol~={atol}")

assert np.allclose(
    [meanx, meany],
    [
        0.0,
        0.0,
    ],
    atol=atol,
)

atol = 0.0  # ignored
rtol = 2.0 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt, emittance_x, emittance_y, emittance_t],
    [
        1.160982600086e-3,
        1.160982600086e-3,
        1.0e-3,
        6.73940299e-7,
        6.73940299e-7,
        2.0e-6,
    ],
    rtol=rtol,
    atol=atol,
)


print("")
print("Final Beam:")
meanx, meany, sigx, sigy, sigt, emittance_x, emittance_y, emittance_t = get_moments(
    final
)
print(f"  meanx={meanx:e} meany={meany:e}")
print(f"  sigx={sigx:e} sigy={sigy:e} sigt={sigt:e}")
print(
    f"  emittance_x={emittance_x:e} emittance_y={emittance_y:e} emittance_t={emittance_t:e}"
)


atol = 3.0 * num_particles ** (-0.5) * sigx
print(f" atol~={atol}")

assert np.allclose(
    [meanx, meany],
    [
        1.79719761842e-4,
        3.24815908981e-4,
    ],
    atol=atol,
)

atol = 0.0  # ignored
rtol = 3.0 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt, emittance_x, emittance_y, emittance_t],
    [
        1.2372883901369e-3,
        1.3772750218080e-3,
        1.027364e-03,
        7.39388142e-7,
        7.39388142e-7,
        2.0e-6,
    ],
    rtol=rtol,
    atol=atol,
)
