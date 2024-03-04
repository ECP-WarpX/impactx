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


# openPMD data series at the beam monitors
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)

# first and last step
final_step = list(series.iterations)[-1]
first_it = series.iterations[1]
final_it = series.iterations[final_step]

# initial beam & reference particle gamma
initial = first_it.particles["beam"].to_df()
initial_gamma_ref = first_it.particles["beam"].get_attribute("gamma_ref")

# final beam & reference particle gamma
final = final_it.particles["beam"].to_df()
final_gamma_ref = final_it.particles["beam"].get_attribute("gamma_ref")

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
print(f"  gamma={initial_gamma_ref:e}")

atol = 0.0  # ignored
rtol = 1.9 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt, emittance_x, emittance_y, emittance_t, initial_gamma_ref],
    [
        5.0e-04,
        5.0e-04,
        5.0e-03,
        4.952764e-09,
        5.028325e-09,
        1.997821e-08,
        40.1389432485322889,
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
print(f"  gamma={final_gamma_ref:e}")


atol = 0.0  # ignored
rtol = 1.9 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt, emittance_x, emittance_y, emittance_t, final_gamma_ref],
    [
        5.004995e-04,
        5.005865e-04,
        3.033949e-03,
        4.067876e-09,
        4.129937e-09,
        6.432081e-05,
        48.8654787469061860,
    ],
    rtol=rtol,
    atol=atol,
)
