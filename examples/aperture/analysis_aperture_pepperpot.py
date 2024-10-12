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
num_particles = 1000000
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

# particle-wise comparison against the periodic rectangular aperture boundary
xmax = 1.5e-4
ymax = 1.0e-4
repeat_x = 1.0e-3
repeat_y = 1.0e-3


# kept particles, shifted to the fundamental domain
xshifted = abs(final["position_x"]) + xmax
yshifted = abs(final["position_y"]) + ymax
u = np.fmod(xshifted,repeat_x) - xmax
v = np.fmod(yshifted,repeat_y) - ymax

# difference from maximum aperture
dx = abs(u) - xmax
dy = abs(v) - ymax 

print()
print(f"  fundamental x_max={u.max()}")
print(f"  fundamental x_min={u.min()}")
assert np.less_equal(dx.max(), 0.0)

print(f"  fundamental y_max={v.max()}")
print(f"  fundamental y_min={v.min()}")
assert np.less_equal(dy.max(), 0.0)

# lost particles, shifted to the fundamental domain
xshifted = abs(particles_lost["position_x"]) - xmax
yshifted = abs(particles_lost["position_y"]) - ymax
u = np.fmod(xshifted,repeat_x) - xmax
v = np.fmod(yshifted,repeat_y) - ymax

# difference from maximum aperture
dx = abs(u) - xmax
dy = abs(v) - ymax 

print()
print(f"  fundamental x_max={u.max()}")
print(f"  fundamental x_min={u.min()}")
assert np.greater_equal(dx.max(), 0.0)

print(f"  fundamental y_max={v.max()}")
print(f"  fundamental y_min={v.min()}")
assert np.greater_equal(dy.max(), 0.0)
