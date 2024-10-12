#!/usr/bin/env python3
#
# Copyright 2022-2024 The ImpactX Community
#
# Authors: Alex Bojanich, Chad Mitchell, Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import sys

import matplotlib.pyplot as plt
import numpy as np

np.set_printoptions(threshold=sys.maxsize)
from conftest import basepath

from impactx import Config, ImpactX, amr, wakeconvolution


def test_wake(save_png=True):
    sim = ImpactX()
    sim.load_inputs_file(basepath + "/examples/chicane/input_chicane_csr.in")
    sim.n_cell = [16, 24, 32]
    sim.slice_step_diagnostics = False

    sim.init_grids()
    sim.init_beam_distribution_from_inputs()
    sim.init_lattice_elements_from_inputs()
    sim.track_particles()

    sim.deposit_charge()

    element_has_csr = True
    R = 10.35  # Units [m]
    sigma_t = 1.9975134930563207e-05

    if element_has_csr:
        pc = sim.particle_container()
        x_min, y_min, t_min, x_max, y_max, t_max = pc.min_and_max_positions()

        is_unity_particle_weight = False
        GetNumberDensity = True

        num_bins = 150
        bin_min = t_min
        bin_max = t_max
        bin_size = (bin_max - bin_min) / (num_bins - 1)

        # Create charge_distribution and slopes as PODVector_real_std
        charge_distribution = amr.PODVector_real_std(num_bins + 1)
        slopes = amr.PODVector_real_std(num_bins)

        # Call deposit_charge with the correct type
        wakeconvolution.deposit_charge(
            pc,
            charge_distribution,
            bin_min,
            bin_size,
            is_unity_particle_weight,
        )

        # Call derivative_charge with the correct types
        wakeconvolution.derivative_charge(
            charge_distribution,
            slopes,
            bin_size,
            GetNumberDensity,
        )

        # Convert charge distribution to numpy array and plot it
        charge_distribution_np = charge_distribution.to_numpy()
        charge_distribution_s_values = np.linspace(
            bin_min, bin_max, len(charge_distribution_np)
        )
        plt.figure()
        plt.plot(charge_distribution_s_values, charge_distribution_np)
        plt.xlabel("Longitudinal Position s (m)")
        plt.ylabel("Charge density")
        plt.title("Charge density vs Longitudinal Position")
        if save_png:
            plt.savefig("density.png")
        else:
            plt.show()

        # Convert slopes to numpy array and plot it
        slopes_np = slopes.to_numpy()
        slopes_s_values = np.linspace(bin_min, bin_max, len(slopes_np))
        plt.figure()
        plt.plot(slopes_s_values, slopes_np)
        plt.xlabel("Longitudinal Position s (m)")
        plt.ylabel("Slopes")
        plt.title("Slopes vs Longitudinal Position")
        if save_png:
            plt.savefig("slopes.png")
        else:
            plt.show()

        # Create wake_function as PODVector_real_std with size 2 * len(slopes)
        wake_function = amr.PODVector_real_std(2 * len(slopes))
        wake_s_values = np.linspace(bin_min, bin_max, 2 * len(slopes))
        for i in range(len(wake_function)):
            if i < num_bins:
                s = i * bin_size
            else:
                s = (i - 2 * num_bins) * bin_size
            wake_s_values[i] = s
            wake_function[i] = wakeconvolution.w_l_csr(s, R, bin_size)

        # Convert wake_function to numpy array and plot it
        wake_function_np = wake_function.to_numpy()
        plt.figure()
        plt.plot(wake_s_values, wake_function_np)
        plt.xlabel("Longitudinal Position s (m)")
        plt.ylabel("Wake Function")
        plt.title("Wake Function vs Longitudinal Position")
        if save_png:
            plt.savefig("wake_function.png")
        else:
            plt.show()

        # Call convolve_fft with the correct types and capture the output
        convolved_wakefield = wakeconvolution.convolve_fft(
            slopes, wake_function, bin_size
        )

        # Convert the result to numpy array
        convolved_wakefield_np = convolved_wakefield.to_numpy()

        # Adjust the s_values to match the length of convolved_wakefield_np
        s_values = np.linspace(bin_min, bin_max, len(convolved_wakefield_np))
        normalized_s_values = s_values / sigma_t

        plt.plot(normalized_s_values, convolved_wakefield_np)
        plt.xlabel("Longitudinal Position s/sigma_s")
        plt.ylabel("Wakefield (V C/m)")
        plt.title("Convolved CSR Wakefield vs Longitudinal Position")
        if save_png:
            plt.savefig("convolved_wakefield.png")
        else:
            plt.show()
        plt.close("all")

    sim.finalize()


if __name__ == "__main__":
    # Call MPI_Init and MPI_Finalize only once:
    if Config.have_mpi:
        from mpi4py import MPI  # noqa

    amr.initialize([])
    test_wake(save_png=False)
