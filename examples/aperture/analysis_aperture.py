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
    sigpx = moment(beam["momentum_x"], moment=2) ** 0.5
    sigy = moment(beam["position_y"], moment=2) ** 0.5
    sigpy = moment(beam["momentum_y"], moment=2) ** 0.5
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

series_lost = io.Series("diags/openPMD/particles_lost.h5", io.Access.read_only)
particles_lost = series_lost.iterations[0].particles["beam"].to_df()

# compare number of particles
num_particles = 10000
assert num_particles == len(initial)
# we lost particles in apertures
assert num_particles > len(final)
assert num_particles == len(particles_lost) + len(final)

print("Initial Beam:")
sigx, sigy, sigt, emittance_x, emittance_y, emittance_t = get_moments(initial)
print(f"  sigx={sigx:e} sigy={sigy:e} sigt={sigt:e}")
print(
    f"  emittance_x={emittance_x:e} emittance_y={emittance_y:e} emittance_t={emittance_t:e}"
)

atol = 0.0  # ignored
rtol = 1.8 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt, emittance_x, emittance_y, emittance_t],
    [
        1.559531175539e-3,
        2.205510139392e-3,
        1.0e-3,
        1.0e-6,
        2.0e-6,
        1.0e-6,
    ],
    rtol=rtol,
    atol=atol,
)

# particle-wise comparison against the rectangular aperture boundary
xmax = 1.0e-3
ymax = 1.5e-3

# kept particles
dx = abs(final["position_x"]) - xmax
dy = abs(final["position_y"]) - ymax

print()
print(f"  x_max={final['position_x'].max()}")
print(f"  x_min={final['position_x'].min()}")
assert np.less_equal(dx.max(), 0.0)

print(f"  y_max={final['position_y'].max()}")
print(f"  y_min={final['position_y'].min()}")
assert np.less_equal(dy.max(), 0.0)

# lost particles
dx = abs(particles_lost["position_x"]) - xmax
dy = abs(particles_lost["position_y"]) - ymax

print()
print(f"  x_max={particles_lost['position_x'].max()}")
print(f"  x_min={particles_lost['position_x'].min()}")
assert np.greater_equal(dx.max(), 0.0)

print(f"  y_max={particles_lost['position_y'].max()}")
print(f"  y_min={particles_lost['position_y'].min()}")
assert np.greater_equal(dy.max(), 0.0)

# check that s is set correctly
lost_at_s = particles_lost["s_lost"]
drift_s = np.ones_like(lost_at_s) * 0.123
assert np.allclose(lost_at_s, drift_s)
