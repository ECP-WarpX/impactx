#!/usr/bin/env python3
#
# Copyright 2022 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import glob

import numpy as np
import pandas as pd
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


def read_all_files(file_pattern):
    """Read in all CSV files from each MPI rank (and potentially OpenMP
    thread). Concatenate into one Pandas dataframe.

    Returns
    -------
    pandas.DataFrame
    """
    return pd.concat(
        (
            pd.read_csv(filename, delimiter=r"\s+")
            for filename in glob.glob(file_pattern)
        ),
        axis=0,
        ignore_index=True,
    ).set_index("id")


# initial/final beam on rank zero
initial = read_all_files("diags/nonlinear_lens_invariants_000000.*")
final = read_all_files("diags/nonlinear_lens_invariants_final.*")

# compare number of particles
num_particles = 10000
assert num_particles == len(initial)
assert num_particles == len(final)

print("Initial Beam:")
meanH, sigH, meanI, sigI = get_moments(initial)
print(f"  meanH={meanH:e} sigH={sigH:e} meanI={meanI:e} sigI={sigI:e}")

atol = 0.0  # ignored
rtol = num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [meanH, sigH, meanI, sigI],
    [1.604948e-01, 1.757985e-01, 2.882956e-01, 3.844099e-01],
    rtol=rtol,
    atol=atol,
)


print("")
print("Final Beam:")
meanH, sigH, meanI, sigI = get_moments(final)
print(f"  meanH={meanH:e} sigH={sigH:e} meanI={meanI:e} sigI={sigI:e}")

atol = 0.0  # ignored
rtol = num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [meanH, sigH, meanI, sigI],
    [1.605336e-01, 1.756235e-01, 2.880581e-01, 3.839625e-01],
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
atol = 2.0e-3
rtol = 1.0  # large number
print()
print(f"  atol={atol} (ignored: rtol~={rtol})")

print(f"  dH_max={beam_joined['dH'].max()}")
assert np.allclose(beam_joined["dH"], 0.0, rtol=rtol, atol=atol)

atol = 3.0e-3
print(f"  atol={atol} (ignored: rtol~={rtol})")
print(f"  dI_max={beam_joined['dI'].max()}")
assert np.allclose(beam_joined["dI"], 0.0, rtol=rtol, atol=atol)
