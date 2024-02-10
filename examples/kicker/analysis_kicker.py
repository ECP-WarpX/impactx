#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import numpy as np
import openpmd_api as io
from scipy.stats import describe, moment


def get_moments(beam):
    """Calculate standard deviations of beam position & momenta
    and emittance values

    Returns
    -------
    sigx, sigy, sigt, emittance_x, emittance_y, emittance_t, meanpx, meanpy
    """
    sigx = moment(beam["position_x"], moment=2) ** 0.5  # variance -> std dev.
    sigpx = moment(beam["momentum_x"], moment=2) ** 0.5
    sigy = moment(beam["position_y"], moment=2) ** 0.5
    sigpy = moment(beam["momentum_y"], moment=2) ** 0.5
    sigt = moment(beam["position_t"], moment=2) ** 0.5
    sigpt = moment(beam["momentum_t"], moment=2) ** 0.5

    meanpx = describe(beam["momentum_x"]).mean
    meanpy = describe(beam["momentum_y"]).mean

    epstrms = beam.cov(ddof=0)
    emittance_x = (sigx**2 * sigpx**2 - epstrms["position_x"]["momentum_x"] ** 2) ** 0.5
    emittance_y = (sigy**2 * sigpy**2 - epstrms["position_y"]["momentum_y"] ** 2) ** 0.5
    emittance_t = (sigt**2 * sigpt**2 - epstrms["position_t"]["momentum_t"] ** 2) ** 0.5

    return (sigx, sigy, sigt, emittance_x, emittance_y, emittance_t, meanpx, meanpy)


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
sigx, sigy, sigt, emittance_x, emittance_y, emittance_t, meanpx, meanpy = get_moments(
    initial
)
print(
    f"  sigx={sigx:e} sigy={sigy:e} sigt={sigt:e} meanpx={meanpx:e} meanpy={meanpy:e}"
)
print(
    f"  emittance_x={emittance_x:e} emittance_y={emittance_y:e} emittance_t={emittance_t:e}"
)

atol = 0.0  # ignored
rtol = 5.0 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt, emittance_x, emittance_y, emittance_t],
    [
        4.017554e-03,
        4.017044e-03,
        9.977588e-04,
        1.197572e-06,
        1.210501e-06,
        2.001382e-06,
    ],
    rtol=rtol,
    atol=atol,
)

atol = rtol * emittance_x / sigx  # relative to rms beam size
assert np.allclose(
    [meanpx, meanpy],
    [
        0.0,
        0.0,
    ],
    atol=atol,
)


print("")
print("Final Beam:")
sigx, sigy, sigt, emittance_x, emittance_y, emittance_t, meanpx, meanpy = get_moments(
    final
)
print(
    f"  sigx={sigx:e} sigy={sigy:e} sigt={sigt:e} meanpx={meanpx:e} meanpy={meanpy:e}"
)
print(
    f"  emittance_x={emittance_x:e} emittance_y={emittance_y:e} emittance_t={emittance_t:e}"
)

atol = 0.0  # ignored
rtol = 5.0 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt, emittance_x, emittance_y, emittance_t],
    [
        4.017554e-03,
        4.017044e-03,
        9.977588e-04,
        1.197572e-06,
        1.210501e-06,
        2.001382e-06,
    ],
    rtol=rtol,
    atol=atol,
)


atol = rtol * emittance_x / sigx  # relative to rms beam size
print(f"  atol~={atol}")

assert np.allclose(
    [meanpx, meanpy],
    [
        2.0e-3,
        3.0e-3,
    ],
    atol=atol,
)
