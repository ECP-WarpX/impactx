import sys

import matplotlib.pyplot as plt
import numpy as np

np.set_printoptions(threshold=sys.maxsize)
# import pdb; pdb.set_trace() - not needed when using breakpoint()
from conftest import basepath

from impactx import ImpactX, amr, wakeconvolution

# Check MPI is off for this test
try:
    from mpi4py import MPI

    mpi_enabled = True
except ImportError:
    mpi_enabled = False

if mpi_enabled:
    print("MPI is enabled. Skipping this script.")
    sys.exit(0)


def test_wake(save_png=True):
    """
    Run ImpactX, calculate wakefield convolution from particle data, and plot the output
    """

    # Initialize AMReX
    amr.initialize([])

    try:
        sim = ImpactX()

        sim.n_cell = [16, 24, 32]
        sim.load_inputs_file(
            basepath + "/../../examples/chicane/input_chicane_csr.in"
        )  # Specify example path here
        sim.slice_step_diagnostics = False

        print("Initializing grids")
        sim.init_grids()
        print("Initializing beam distribution from inputs")
        sim.init_beam_distribution_from_inputs()
        print("Initializing lattice elements from inputs")
        sim.init_lattice_elements_from_inputs()

        print("Starting evolution")
        sim.evolve()  # Add following script inside of evolve for CSR along lattice
        print("Evolution completed")

        # Deposit charge
        print("Depositing charge")
        sim.deposit_charge()

        # Check for CSR elements and perform necessary calculations
        element_has_csr = True
        R = 10.35  # Units [m] TODO: Read in value
        beam_charge = 1.0e-09  # Units [C] TODO: Read in value

        # Enter loop if lattice has bend element
        if element_has_csr:
            print("Inside CSR element check")
            # Measure beam size, extract the min, max of particle positions
            pc = sim.particle_container()
            x_min, y_min, t_min, x_max, y_max, t_max = pc.min_and_max_positions()
            print(f"x_min: {x_min}, y_min: {y_min}, t_min: {t_min}")

            # Set parameters for charge deposition
            is_unity_particle_weight = True
            GetNumberDensity = True

            num_bins = 100  # Set resolution
            bin_min = t_min
            bin_max = t_max
            bin_size = (bin_max - bin_min) / num_bins

            padding_factor = 1  # Keep set to 1
            pad_factor = 2  # Change this to change the zero-padding
            sigma_t = 1.9975134930563207e-05

            # Calculate original length of the convolution result
            original_length = num_bins

            # Calculate target length with padding
            target_length = pad_factor * (2 * num_bins - 1)

            # Zero-pad the corresponding particle distribution array to match new position range
            total_padding = target_length - original_length

            # Ensure padding is split equally on both sides
            if total_padding % 2 == 0:
                pad_before = total_padding // 2
                pad_after = total_padding // 2
            else:
                pad_before = total_padding // 2
                pad_after = total_padding // 2 + 1

            # Calculate new start and end points
            new_start = bin_min - pad_before * bin_size
            new_end = bin_max + pad_after * bin_size

            # Generate s_values and normalize
            s_values = np.linspace(new_start, new_end, target_length)
            s_values_bin_size = (new_end - new_start) / target_length

            # Allocate memory for the charge profile
            charge_distribution = np.zeros(target_length, dtype=np.double)
            print(len(charge_distribution))

            # Perform charge deposition
            print("Performing charge deposition")
            wakeconvolution.deposit_charge(
                pc,
                charge_distribution,
                target_length,
                new_start,
                s_values_bin_size,
                is_unity_particle_weight,
            )

            # Compute the derivative of the charge distribution
            print("Computing the derivative of the charge distribution")
            slopes = np.zeros(target_length - 1, dtype=np.double)
            print(len(slopes))
            wakeconvolution.derivative_charge(
                charge_distribution,
                slopes,
                target_length,
                s_values_bin_size,
                GetNumberDensity,
            )

            # Calculate the CSR wake function
            print("Calculating the CSR wake function")
            wake_function = np.array(
                [
                    wakeconvolution.w_l_csr(new_start + (i * s_values_bin_size), R)
                    for i in range(target_length)
                ],
                dtype=np.double,
            )
            print(len(wake_function))

            # Perform FFT convolution
            print("Performing FFT convolution")
            convoluted_wakefield = np.zeros(2 * target_length - 1, dtype=np.double)
            print("These are my slopes:", slopes)
            print("These are my wake function results:", wake_function)
            print("My bin size is:", s_values_bin_size)
            print("My padding factor is:", padding_factor)
            # breakpoint() #Use pdb to debug
            wakeconvolution.convolve_fft(
                slopes,
                wake_function,
                s_values_bin_size,
                convoluted_wakefield,
                padding_factor,
            )
            print("My convoluted result is:", convoluted_wakefield)

            lower_bound = 2 * s_values[0]
            upper_bound = 2 * s_values[-1]
            s_values_full_range = np.linspace(
                lower_bound, upper_bound, 2 * len(s_values) - 1
            )
            normalized_s_values = s_values_full_range / sigma_t

            # Plot convoluted wakefield
            print("Plotting convoluted wakefield")
            # breakpoint()
            plt.plot(
                normalized_s_values, convoluted_wakefield, label="Convoluted Wakefield"
            )
            # plt.xlim(left = -5)
            # plt.xlim(right = 5)
            plt.xlabel("Longitudinal Position s/sigma_s")
            plt.ylabel("Wakefield (V C/m)")
            plt.legend()
            plt.title("Convoluted Wakefield")
            if save_png:
                plt.savefig("convoluted_wakefield.png")
            else:
                plt.show()
            plt.close("all")

    finally:
        # Finalize simulation
        print("Finalizing simulation")
        sim.finalize()


if __name__ == "__main__":
    test_wake()
