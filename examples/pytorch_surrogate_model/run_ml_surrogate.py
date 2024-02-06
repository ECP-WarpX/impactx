#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Ryan Sandberg, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import sys
import tarfile
from urllib import request

import numpy as np
from surrogate_model_definitions import surrogate_model

try:
    import torch
except ImportError:
    print("Warning: Cannot import PyTorch. Skipping test.")
    sys.exit(0)

from impactx import (
    ImpactX,
    ImpactXParIter,
    TransformationDirection,
    coordinate_transformation,
    distribution,
    elements,
)


def download_and_unzip(url, data_dir):
    request.urlretrieve(url, data_dir)
    with tarfile.open(data_dir) as tar_dataset:
        tar_dataset.extractall()


# load models
N_stage = 9

data_url = (
    "https://zenodo.org/records/10368972/files/ml_example_inference.tar.gz?download=1"
)
download_and_unzip(data_url, "inference_dataset")

dataset_dir = "datasets/"
model_dir = "models/"

model_list = [
    surrogate_model(
        dataset_dir + f"dataset_beam_stage_{i}.pt",
        model_dir + f"beam_stage_{i}_model.pt",
    )
    for i in range(N_stage)
]

# information specific to the WarpX simulation
# for which the neural networks are surrogates
ebeam_lpa_z0 = -107e-6
L_plasma = 0.28
L_transport = 0.03
L_stage_period = L_plasma + L_transport
drift_after_LPA = 43e-6
L_surrogate = abs(ebeam_lpa_z0) + L_plasma + drift_after_LPA

# number of slices per ds in each lattice element
ns = 1


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
        array = np.array

        ref_part = pc.ref_particle()
        ref_z_i = ref_part.z
        ref_z_i_LPA = ref_z_i - self.stage_start
        ref_z_f = ref_z_i + self.surrogate_length

        ref_part_tensor = torch.tensor(
            [ref_part.x, ref_part.y, ref_z_i_LPA, ref_part.px, ref_part.py, ref_part.pz]
        )
        ref_beta_gamma = np.sqrt(torch.sum(ref_part_tensor[3:] ** 2))

        with torch.no_grad():
            ref_part_model_final = self.surrogate_model(ref_part_tensor.float())
        ref_uz_f = ref_part_model_final[5]
        ref_beta_gamma_final = (
            ref_uz_f  # NOT np.sqrt(torch.sum(ref_part_model_final[3:]**2))
        )
        ref_part_final = torch.tensor([0, 0, ref_z_f, 0, 0, ref_uz_f])

        # transform
        coordinate_transformation(pc, TransformationDirection.to_fixed_t)

        for lvl in range(pc.finest_level + 1):
            for pti in ImpactXParIter(pc, level=lvl):
                aos = pti.aos()
                aos_arr = array(aos, copy=False)

                soa = pti.soa()
                real_arrays = soa.GetRealData()
                px = array(real_arrays[0], copy=False)
                py = array(real_arrays[1], copy=False)
                pt = array(real_arrays[2], copy=False)
                data_arr = (
                    torch.tensor(
                        np.vstack(
                            [aos_arr["x"], aos_arr["y"], aos_arr["z"], real_arrays[:3]]
                        )
                    )
                    .float()
                    .T
                )

                data_arr[:, 0] += ref_part.x
                data_arr[:, 1] += ref_part.y
                data_arr[:, 2] += ref_z_i_LPA
                data_arr[:, 3:] *= ref_beta_gamma
                data_arr[:, 3] += ref_part.px
                data_arr[:, 4] += ref_part.py
                data_arr[:, 5] += ref_part.pz
                #     # TODO this part needs to be corrected for general geometry
                #     # where the initial vector might not point in z
                #     # and even if it does, bending elements may change the direction
                #     # i.e. do we need to make sure beam is pointing in the right direction?
                #     # assume for now it is

                with torch.no_grad():
                    data_arr_post_model = self.surrogate_model(data_arr.float())

                #  need to add stage start to z
                data_arr_post_model[:, 2] += self.stage_start

                # back to ref particle coordinates
                for ii in range(3):
                    data_arr_post_model[:, ii] -= ref_part_final[ii]
                    data_arr_post_model[:, 3 + ii] -= ref_part_final[3 + ii]
                    data_arr_post_model[:, 3 + ii] /= ref_beta_gamma_final

                aos_arr["x"] = data_arr_post_model[:, 0]
                aos_arr["y"] = data_arr_post_model[:, 1]
                aos_arr["z"] = data_arr_post_model[:, 2]
                px[:] = data_arr_post_model[:, 3]
                py[:] = data_arr_post_model[:, 4]
                pt[:] = data_arr_post_model[:, 5]

        # TODO this part needs to be corrected for general geometry
        # where the initial vector might not point in z
        # and even if it does, bending elements may change the direction

        ref_part.x = ref_part_final[0]
        ref_part.y = ref_part_final[1]
        ref_part.z = ref_part_final[2]
        ref_gamma = np.sqrt(1 + ref_beta_gamma_final**2)
        ref_part.px = ref_part_final[3]
        ref_part.py = ref_part_final[4]
        ref_part.pz = ref_part_final[5]
        ref_part.pt = -ref_gamma
        # for now, I am applying the hack of manually setting s=z=ct.
        # this will need to be revisited and evaluated more correctly
        # when the accelerator length is more consistently established
        ref_part.s += self.surrogate_length
        ref_part.t += self.surrogate_length
        # ref_part.s += pge1.ds
        # ref_part.t += pge1.ds / ref_beta

        coordinate_transformation(pc, TransformationDirection.to_fixed_s)
        ## Done!


lpa_stage_list = []
for i in range(N_stage):
    lpa = LPASurrogateStage(i, model_list[i], L_surrogate, L_stage_period * i)
    lpa.nslice = ns
    lpa.ds = L_surrogate
    lpa_stage_list.append(lpa)

#########
sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 1 GeV electron beam with an initial
# normalized rms emittance of 1 nm
ref_u = 1957
energy_gamma = np.sqrt(1 + ref_u**2)
energy_MeV = 0.510998950 * energy_gamma  # reference energy
bunch_charge_C = 10.0e-15  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(energy_MeV)
ref.z = ebeam_lpa_z0

#   particle bunch
distr = distribution.Gaussian(
    sigmaX=0.75e-6,
    sigmaY=0.75e-6,
    sigmaT=0.1e-6,
    sigmaPx=1.33 / energy_gamma,
    sigmaPy=1.33 / energy_gamma,
    sigmaPt=1e-8,
    muxpx=0.0,
    muypy=0.0,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)
pc = sim.particle_container()

L_transport = 0.03
L_lens = 0.003
L_focal = 0.5 * L_transport
L_drift = 0.5 * (L_transport - L_lens)
Kxy = np.sqrt(2.0 / L_focal / L_lens)
Kt = 1e-11

L_drift_minus_surrogate = L_drift
L_drift_1 = L_drift - drift_after_LPA

L_drift_before_2nd_stage = abs(ebeam_lpa_z0)
L_drift_2 = L_drift - L_drift_before_2nd_stage


#########

###
monitor = elements.BeamMonitor("monitor")
for i in range(N_stage):
    sim.lattice.extend(
        [
            monitor,
            lpa_stage_list[i],
        ]
    )

    if i != N_stage - 1:
        sim.lattice.extend(
            [
                monitor,
                elements.Drift(ds=L_drift_1),
                monitor,
                elements.ConstF(ds=L_lens, kx=Kxy, ky=Kxy, kt=Kt),
                monitor,
                elements.Drift(ds=L_drift_2),
            ]
        )
sim.lattice.extend([monitor])

sim.evolve()

# clean shutdown
sim.finalize()
