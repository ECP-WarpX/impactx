#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Ryan Sandberg, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import argparse

import amrex.space3d as amr

try:
    import cupy as cp

    cupy_available = True
except ImportError:
    cupy_available = False

import numpy as np
import scipy.optimize as opt
from surrogate_model_definitions import surrogate_model

from impactx import (
    Config,
    CoordSystem,
    ImpactX,
    ImpactXParIter,
    coordinate_transformation,
    distribution,
    elements,
)

try:
    import torch
except ImportError:
    print("Warning: Cannot import PyTorch. Skipping test.")
    import sys

    sys.exit(0)

import zipfile
from urllib import request

parser = argparse.ArgumentParser()
parser.add_argument(
    "--num_particles",
    "-N",
    type=int,
    default=100000,
    help="number of particles to use in beam",
)
parser.add_argument(
    "--N_stages",
    "-ns",
    type=int,
    default=15,
    choices=range(1, 16),
    help="number of LPA stages to simulate",
)
args = parser.parse_args()
if Config.have_gpu and cupy_available:
    array = cp.array
    stack = cp.stack
    sqrt = cp.sqrt
    device = torch.device("cuda")
    if Config.gpu_backend == "SYCL":
        print("Warning: SYCL GPU backend not yet implemented for Python")
else:
    array = np.array
    stack = np.stack
    sqrt = np.sqrt
    device = None
if device is not None:
    print(f"device={device}")
else:
    print("device set to default, cpu")

N_stage = args.N_stages
tune_by_x_or_y = "x"
npart = args.num_particles
ebeam_lpa_z0 = -107e-6
L_plasma = 0.28
L_transport = 0.03
L_stage_period = L_plasma + L_transport
drift_after_LPA = 43e-6
L_surrogate = abs(ebeam_lpa_z0) + L_plasma + drift_after_LPA


def download_and_unzip(url, data_dir):
    request.urlretrieve(url, data_dir)
    with zipfile.ZipFile(data_dir, "r") as zip_dataset:
        zip_dataset.extractall()


print(
    "Downloading trained models from Zenodo.org - this might take a minute...",
    flush=True,
)
data_url = "https://zenodo.org/records/10810754/files/models.zip?download=1"
download_and_unzip(data_url, "models.zip")

# It was found that the PyTorch multithreaded defaults interfere with MPI-enabled AMReX
# when initializing the models: https://github.com/AMReX-Codes/pyamrex/issues/322
# So we manually set the number of threads to serial (1).
if Config.have_mpi:
    n_threads = torch.get_num_threads()
    torch.set_num_threads(1)
model_list = [
    surrogate_model(f"models/beam_stage_{stage_i}_model.pt", device=device)
    for stage_i in range(N_stage)
]
if Config.have_mpi:
    torch.set_num_threads(n_threads)

pp_amrex = amr.ParmParse("amrex")
pp_amrex.add("the_arena_init_size", 0)
pp_amrex.add("the_device_arena_init_size", 0)

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = False
sim.diagnostics = True  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 1 GeV electron beam with an initial
# unnormalized rms emittance of 1 nm
ref_u = 1957
energy_gamma = np.sqrt(1 + ref_u**2)
energy_MeV = 0.510998950 * energy_gamma  # reference energy
bunch_charge_C = 10.0e-15  # used with space charge


#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(energy_MeV)
ref.z = ebeam_lpa_z0

pc = sim.particle_container()

distr = distribution.Gaussian(
    lambdaX=0.75e-6,
    lambdaY=0.75e-6,
    lambdaT=0.1e-6,
    lambdaPx=1.33 / energy_gamma,
    lambdaPy=1.33 / energy_gamma,
    lambdaPt=1e-8,
    muxpx=0.0,
    muypy=0.0,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

n_slice = 1


class LPASurrogateStage(elements.Programmable):
    def __init__(self, stage_i, surrogate_model, surrogate_length, stage_start):
        elements.Programmable.__init__(self)
        self.stage_i = stage_i
        self.surrogate_model = surrogate_model
        self.surrogate_length = surrogate_length
        self.stage_start = stage_start
        self.push = self.surrogate_push
        self.ds = surrogate_length

    def surrogate_push(self, pc, step):
        ref_part = pc.ref_particle()
        ref_z_i = ref_part.z
        ref_z_i_LPA = ref_z_i - self.stage_start
        ref_z_f = ref_z_i + self.surrogate_length
        ref_part_tensor = torch.tensor(
            [
                ref_part.x,
                ref_part.y,
                ref_z_i_LPA,
                ref_part.px,
                ref_part.py,
                ref_part.pz,
            ],
            device=device,
            dtype=torch.float64,
        )
        ref_beta_gamma = torch.sqrt(torch.sum(ref_part_tensor[3:] ** 2))
        ref_beta_gamma = ref_beta_gamma.to(device)

        with torch.no_grad():
            ref_part_model_final = self.surrogate_model(ref_part_tensor)
        ref_uz_f = ref_part_model_final[5]
        ref_beta_gamma_final = ref_uz_f
        ref_part_final = torch.tensor(
            [0, 0, ref_z_f, 0, 0, ref_uz_f], device=device, dtype=torch.float64
        )

        coordinate_transformation(pc, CoordSystem.t)

        for lvl in range(pc.finest_level + 1):
            for pti in ImpactXParIter(pc, level=lvl):
                soa = pti.soa()
                real_arrays = soa.get_real_data()
                x = array(real_arrays[0], copy=False)
                y = array(real_arrays[1], copy=False)
                t = array(real_arrays[2], copy=False)
                px = array(real_arrays[3], copy=False)
                py = array(real_arrays[4], copy=False)
                pt = array(real_arrays[5], copy=False)
                data_arr = torch.tensor(
                    stack([x, y, t, px, py, pt], axis=1),
                    device=device,
                    dtype=torch.float64,
                )

                data_arr[:, 0] += ref_part.x
                data_arr[:, 1] += ref_part.y
                data_arr[:, 2] += ref_z_i_LPA
                data_arr[:, 3:] *= ref_beta_gamma
                data_arr[:, 3] += ref_part.px
                data_arr[:, 4] += ref_part.py
                data_arr[:, 5] += ref_part.pz

                with torch.no_grad():
                    data_arr_post_model = self.surrogate_model(data_arr)

                #  z += stage start
                data_arr_post_model[:, 2] += self.stage_start
                # back to ref particle coordinates
                for ii in range(3):
                    data_arr_post_model[:, ii] -= ref_part_final[ii]
                    data_arr_post_model[:, 3 + ii] -= ref_part_final[3 + ii]
                    data_arr_post_model[:, 3 + ii] /= ref_beta_gamma_final

                x[:] = array(data_arr_post_model[:, 0])
                y[:] = array(data_arr_post_model[:, 1])
                t[:] = array(data_arr_post_model[:, 2])
                px[:] = array(data_arr_post_model[:, 3])
                py[:] = array(data_arr_post_model[:, 4])
                pt[:] = array(data_arr_post_model[:, 5])

        # TODO this part needs to be corrected for general geometry
        # where the initial vector might not point in z
        # and even if it does, bending elements may change the direction

        ref_part.x = ref_part_final[0]
        ref_part.y = ref_part_final[1]
        ref_part.z = ref_part_final[2]
        ref_gamma = torch.sqrt(1 + ref_beta_gamma_final**2)
        ref_part.px = ref_part_final[3]
        ref_part.py = ref_part_final[4]
        ref_part.pz = ref_part_final[5]
        ref_part.pt = -ref_gamma
        ref_part.s += self.surrogate_length
        ref_part.t += self.surrogate_length

        coordinate_transformation(pc, CoordSystem.s)
        ## Done!


L_transport = 0.03
L_lens = 0.003
L_focal = 0.5 * L_transport
L_drift = 0.5 * (L_transport - L_lens)
K = np.sqrt(2.0 / L_focal / L_lens)
Kt = 1e-11  # number chosen arbitrarily since 0 isn't allowed
dL = 0

L_drift_minus_surrogate = L_drift
L_drift_1 = L_drift - drift_after_LPA - dL

L_drift_before_2nd_stage = abs(ebeam_lpa_z0)
L_drift_2 = L_drift - L_drift_before_2nd_stage + dL


def get_lattice_element_iter(sim, j):
    assert (
        0 <= j < len(sim.lattice)
    ), f"Argument j must be a nonnegative integer satisfying 0 <= j < {len(sim.lattice)}, not {j}"
    i = 0
    lat_it = sim.lattice.__iter__()
    next(lat_it)
    while i != j:
        next(lat_it)
        i += 1
    return lat_it


def lens_eqn(k, lens_length, alpha, beta, gamma):
    return np.tan(k * lens_length) + 2 * alpha / (k * beta - gamma / k)


k_list = []


class UpdateConstF(elements.Programmable):
    def __init__(self, sim, stage_i, lattice_index, x_or_y):
        elements.Programmable.__init__(self)
        self.sim = sim
        self.stage_i = stage_i
        self.lattice_index = lattice_index
        self.x_or_y = x_or_y
        self.push = self.set_lens

    def set_lens(self, pc, step):
        # get envelope parameters
        rbc = pc.reduced_beam_characteristics()
        alpha = rbc[f"alpha_{self.x_or_y}"]
        beta = rbc[f"beta_{self.x_or_y}"]
        gamma = (1 + alpha**2) / beta
        # solve for k_new
        sol = opt.root_scalar(
            lens_eqn, bracket=[100, 300], args=(L_lens, alpha, beta, gamma)
        )
        k_new = sol.root
        # set lens
        self_it = get_lattice_element_iter(self.sim, self.lattice_index)
        following_lens = next(self_it)
        k_list.append(k_new)
        following_lens.kx = k_new
        following_lens.ky = k_new


lpa_stages = []
for i in range(N_stage):
    lpa = LPASurrogateStage(i, model_list[i], L_surrogate, L_stage_period * i)
    lpa.nslice = n_slice
    lpa.ds = L_surrogate
    lpa_stages.append(lpa)

monitor = elements.BeamMonitor("monitor")
for i in range(N_stage):
    sim.lattice.extend(
        [
            monitor,
            lpa_stages[i],
        ]
    )

    if i != N_stage - 1:
        sim.lattice.extend(
            [
                monitor,
                elements.Drift(ds=L_drift_1),
                monitor,
                UpdateConstF(
                    sim=sim, stage_i=i, lattice_index=5 + 9 * i, x_or_y=tune_by_x_or_y
                ),
                elements.ConstF(ds=L_lens, kx=K, ky=K, kt=Kt),
                monitor,
                elements.Drift(ds=L_drift_2),
            ]
        )
sim.lattice.extend([monitor])

sim.evolve()
sim.finalize()
del sim
