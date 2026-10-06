#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import glob
import re

import numpy as np
import pandas as pd

from impactx import Config


def read_file(file_pattern):
    for filename in glob.glob(file_pattern):
        df = pd.read_csv(filename, delimiter=r"\s+")
        if "step" not in df.columns:
            step = int(re.findall(r"[0-9]+", filename)[0])
            df["step"] = step
        yield df


def read_time_series(file_pattern):
    """Read in all CSV files from each MPI rank (and potentially OpenMP
    thread). Concatenate into one Pandas dataframe.

    Returns
    -------
    pandas.DataFrame
    """
    return pd.concat(
        read_file(file_pattern),
        axis=0,
        ignore_index=True,
    )  # .set_index('id')


# read reduced diagnostics
rbc = read_time_series("diags/reduced_beam_characteristics.*")
is_double = Config.precision != "SINGLE"

s = rbc["s"]
sigma_x = rbc["sigma_x"]
sigma_y = rbc["sigma_y"]
sigma_t = rbc["sigma_t"]
emittance_x = rbc["emittance_x"]
emittance_y = rbc["emittance_y"]
emittance_t = rbc["emittance_t"]

sigma_xi = sigma_x.iloc[0]
sigma_yi = sigma_y.iloc[0]
sigma_ti = sigma_t.iloc[0]
emittance_xi = emittance_x.iloc[0]
emittance_yi = emittance_y.iloc[0]
emittance_ti = emittance_t.iloc[0]

length = len(s) - 1

sf = s.iloc[length]
sigma_xf = sigma_x.iloc[length]
sigma_yf = sigma_y.iloc[length]
sigma_tf = sigma_t.iloc[length]
emittance_xf = emittance_x.iloc[length]
emittance_yf = emittance_y.iloc[length]
emittance_tf = emittance_t.iloc[length]


print("Initial Beam:")
print(f"  sigx={sigma_xi:e} sigy={sigma_yi:e} sigt={sigma_ti:e}")
print(
    f"  emittance_x={emittance_xi:e} emittance_y={emittance_yi:e} emittance_t={emittance_ti:e}"
)

atol = 0.0  # ignored
rtol = 1.0e-3  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigma_xi, sigma_yi, sigma_ti, emittance_xi, emittance_yi, emittance_ti],
    [
        8.590000e-04,
        8.590000e-04,
        7.071068e-07,
        1.000000e-06,
        1.000000e-06,
        1.000000e-12,
    ],
    rtol=rtol,
    atol=atol,
)


print("")
print("Final Beam:")
print(f"  sigx={sigma_xf:e} sigy={sigma_yf:e} sigt={sigma_tf:e}")
print(
    f"  emittance_x={emittance_xf:e} emittance_y={emittance_yf:e} emittance_t={emittance_tf:e}"
)

atol = 0.0  # ignored
rtol = 1.0e-3  # from random sampling of a smooth distribution
if not is_double:
    # FIXME: in single precision, the final emittance_t is off by up to 0.7% (the
    # envelope model preserves it to 4e-6 in double precision)
    # https://github.com/BLAST-ImpactX/impactx/issues/1684
    rtol = 1.0e-2
print(f"  rtol={rtol} (ignored: atol~={atol})")

# the beam returns to within 0.1% of its initial transverse size over the FODO cell:
# these are the converged values of the envelope model, reached to 1e-5 at nslice=50
assert np.allclose(
    [sigma_xf, sigma_yf, sigma_tf, emittance_xf, emittance_yf, emittance_tf],
    [
        8.585939e-04,
        8.581282e-04,
        4.140854e-05,
        1.000000e-06,
        1.000000e-06,
        1.000000e-12,
    ],
    rtol=rtol,
    atol=atol,
)
