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
    """Calculate standard deviations of beam position & momenta
    and emittance values

    Returns
    -------
    sigx, sigy, sigt, emittance_x, emittance_y, emittance_t
    """
    sigx = moment(beam["position_x"], moment=2) ** 0.5  # variance -> std dev.
    sigpx = moment(beam["divergence_x"], moment=2) ** 0.5
    sigy = moment(beam["position_y"], moment=2) ** 0.5
    sigpy = moment(beam["divergence_y"], moment=2) ** 0.5
    sigt = moment(beam["position_t"], moment=2) ** 0.5
    sigpt = moment(beam["momentum_t"], moment=2) ** 0.5

    epstrms = beam.cov(ddof=0)
    emittance_x = (sigx**2 * sigpx**2 - epstrms["position_x"]["momentum_x"] ** 2) ** 0.5
    emittance_y = (sigy**2 * sigpy**2 - epstrms["position_y"]["momentum_y"] ** 2) ** 0.5
    emittance_t = (sigt**2 * sigpt**2 - epstrms["position_t"]["momentum_t"] ** 2) ** 0.5

    return (sigx, sigy, sigt, emittance_x, emittance_y, emittance_t)


# initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial = series.iterations[1].particles["beam"].to_df()
final = series.iterations[last_step].particles["beam"].to_df()

# compare number of particles
num_particles = 100000
assert num_particles == len(initial)
assert num_particles == len(final)

scale = (
    (1.0 - initial.momentum_t) ** 2
    + (initial.momentum_x) ** 2
    + (initial.momentum_y) ** 2
)
xp = initial.momentum_x / np.sqrt(scale)
initial["divergence_x"] = xp
yp = initial.momentum_y / np.sqrt(scale)
initial["divergence_y"] = yp

print("Initial Beam:")
sigx, sigy, sigt, emittance_x, emittance_y, emittance_t = get_moments(initial)
print(f"  sigx={sigx:e} sigy={sigy:e} sigt={sigt:e}")
print(
    f"  emittance_x={emittance_x:e} emittance_y={emittance_y:e} emittance_t={emittance_t:e}"
)

atol = 0.0  # ignored
rtol = 3.0 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt, emittance_x, emittance_y, emittance_t],
    [
        1.288697604e-6,
        1.288697604e-6,
        1.0e-6,
        5.10997388810014764e-12,
        5.10997388810014764e-12,
        1.0e-8,
    ],
    rtol=rtol,
    atol=atol,
)


scale = (
    (1.0 - final.momentum_t) ** 2 + (final.momentum_x) ** 2 + (final.momentum_y) ** 2
)
xp = final.momentum_x / np.sqrt(scale)
final["divergence_x"] = xp
yp = final.momentum_y / np.sqrt(scale)
final["divergence_y"] = yp

print("")
print("Final Beam:")
sigx, sigy, sigt, emittance_xf, emittance_yf, emittance_tf = get_moments(final)
demittance_x = 100 * (emittance_xf - emittance_x) / emittance_x
demittance_y = 100 * (emittance_yf - emittance_y) / emittance_y
demittance_t = 100 * (emittance_tf - emittance_t) / emittance_t

print(f"  sigx={sigx:e} sigy={sigy:e} sigt={sigt:e}")
print(
    f"  emittance change x (%)={demittance_x:e} emittance change y (%)={demittance_y:e} emittance change t (%)={demittance_t:e}"
)

atol = 0.0  # ignored
rtol = 26.0 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt, demittance_x, demittance_y, emittance_t],
    [
        1.245e-6,
        1.245e-6,
        1.0e-6,
        0.94,
        0.94,
        1.0e-8,
    ],
    rtol=rtol,
    atol=atol,
)
