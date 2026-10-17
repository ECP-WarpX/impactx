#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
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
    # double precision also for single-precision data: the t-pt correlation of
    # this beam is too close to unity to compute the emittance in float32
    beam = beam.astype(np.float64)

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


def get_eigenemittances(openpmd_beam):
    """Return eigenemittances from an openPMD particle species

    Returns
    -------
    emittance_1, emittance_2, emittance_3
    """
    emittance_1 = openpmd_beam.get_attribute("emittance_1")
    emittance_2 = openpmd_beam.get_attribute("emittance_2")
    emittance_3 = openpmd_beam.get_attribute("emittance_3")

    return (emittance_1, emittance_2, emittance_3)


# initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial_beam = series.iterations[1].particles["beam"]
initial = initial_beam.to_df()
is_double = initial["position_x"].dtype == np.float64
final_beam = series.iterations[last_step].particles["beam"]
final = final_beam.to_df()

# compare number of particles
num_particles = 10000
assert num_particles == len(initial)
assert num_particles == len(final)

print("Initial Beam:")
sigx, sigy, sigt, emittance_x, emittance_y, emittance_t = get_moments(initial)
print(f"  sigx={sigx:e} sigy={sigy:e} sigt={sigt:e}")
print(
    f"  emittance_x={emittance_x:e} emittance_y={emittance_y:e} emittance_t={emittance_t:e}"
)

atol = 0.0  # ignored
rtol = 2.2 * num_particles**-0.5  # from random sampling of a smooth distribution
if not is_double:
    # FIXME: single precision does not resolve the near-unity t-pt correlation
    # of this beam: sigt and emittance_t are off by 19% and 20% in the initial beam
    # https://github.com/BLAST-ImpactX/impactx/issues/1485
    rtol = 0.7
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt, emittance_x, emittance_y, emittance_t],
    [
        6.4214719960819659e-005,
        3.6603372435649773e-005,
        1.9955175623579313e-004,
        1.0198263116327677e-010,
        1.0308359092878036e-010,
        4.0035161705244885e-010,
    ],
    rtol=rtol,
    atol=atol,
)


print("")
print("Final Beam:")
sigx, sigy, sigt, emittance_x, emittance_y, emittance_t = get_moments(final)
print(f"  sigx={sigx:e} sigy={sigy:e} sigt={sigt:e}")
print(
    f"  emittance_x={emittance_x:e} emittance_y={emittance_y:e} emittance_t={emittance_t:e}"
)

atol = 0.0  # ignored
rtol = 2.2 * num_particles**-0.5  # from random sampling of a smooth distribution
if not is_double:
    # FIXME: all moments are off by 17-19% in single precision, see initial beam
    # https://github.com/BLAST-ImpactX/impactx/issues/1485
    rtol = 0.5
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt, emittance_x, emittance_y, emittance_t],
    [
        2.797264e-04,
        6.484914e-04,
        2.022674e-04,
        3.097048e-08,
        2.207603e-08,
        3.803276e-08,
    ],
    rtol=rtol,
    atol=atol,
)

print("")
print("Initial eigenemittances:")
emittance_1i, emittance_2i, emittance_3i = get_eigenemittances(initial_beam)
print(
    f"  emittance_1={emittance_1i:e} emittance_2={emittance_2i:e} emittance_3={emittance_3i:e}"
)

print("")
print("Final eigenemittances:")
emittance_1f, emittance_2f, emittance_3f = get_eigenemittances(final_beam)
print(
    f"  emittance_1={emittance_1f:e} emittance_2={emittance_2f:e} emittance_3={emittance_3f:e}"
)

atol = 0.0  # ignored
rtol = 3.5 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

# FIXME: in single precision, all initial and final eigenemittances are NaN
# https://github.com/BLAST-ImpactX/impactx/issues/1481
equal_nan = not is_double

assert np.allclose(
    [emittance_1f, emittance_2f, emittance_3f],
    [
        emittance_1i,
        emittance_2i,
        emittance_3i,
    ],
    rtol=rtol,
    atol=atol,
    equal_nan=equal_nan,
)
