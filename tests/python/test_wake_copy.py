import matplotlib.pyplot as plt
import numpy as np
from conftest import basepath

from impactx import ImpactX, amr, wakeconvolution

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
            basepath + "/../../examples_backup/chicane/input_chicane.in"
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
        R = 10.35  # Units [m] TODO: Read in value
        beam_charge = 1.0e-09  # Units [C] TODO: Read in value

        # Enter loop if lattice has bend element
        if element_has_csr:
            print("Inside CSR element check")
            # Measure beam size, extract the min, max of particle positions
            pc = sim.particle_container()
            x_min, y_min, t_min, x_max, y_max, t_max = pc.min_and_max_positions()
            print(f"x_min: {x_min}, y_min: {y_min}, t_min: {t_min}")
            print(f"x_max: {x_max}, y_max: {y_max}, t_max: {t_max}")

            # Set parameters for charge deposition
            is_unity_particle_weight = True
            GetNumberDensity = True

            num_bins = 100  # Set resolution
            bin_min = t_min
            bin_max = t_max
            bin_size = (bin_max - bin_min) / num_bins

            padding_factor = 1 # Make plotting script compatible
            sigma_t = 1.9975134930563207e-05

            # Allocate memory for the charge profile
            charge_distribution = np.zeros(num_bins, dtype=np.double)

            # Perform charge deposition
            print("Performing charge deposition")
            wakeconvolution.deposit_charge(
                pc,
                charge_distribution,
                num_bins,
                bin_min,
                bin_size,
                is_unity_particle_weight,
            )

            # Compute the derivative of the charge distribution
            print("Computing the derivative of the charge distribution")
            slopes = np.zeros(num_bins - 1, dtype=np.double)
            wakeconvolution.derivative_charge(
                charge_distribution,
                slopes,
                num_bins,
                bin_size,
                GetNumberDensity,
            )

            # Calculate the CSR wake function
            print("Calculating the CSR wake function")
            wake_function = np.array(
                [
                    wakeconvolution.w_l_csr(bin_min + (i * bin_size), R, beam_charge)
                    for i in range(num_bins)
                ],
                dtype=np.double,
            )

            # Perform FFT convolution
            print("Performing FFT convolution")
            convoluted_wakefield = np.zeros(2 * num_bins - 1, dtype=np.double)
            wakeconvolution.convolve_fft(
                slopes,
                wake_function[:-1],
                bin_size,
                convoluted_wakefield,
                padding_factor,
            )

            # Plot convoluted wakefield
            print("Plotting convoluted wakefield")

            lower_bound = 2 * bin_min
            upper_bound = 2 * bin_max
            s_values = np.linspace(lower_bound, upper_bound, len(convoluted_wakefield))
            normalized_s_values = s_values / sigma_t

            plt.figure()
            plt.plot(normalized_s_values, convoluted_wakefield, label="Convoluted Wakefield")
            plt.xlim(left = -5)
            plt.xlim(right = 5)
            plt.xlabel("Longitudinal Position s/sigma_s")
            plt.ylabel("Wakefield (V C/m)")
            plt.legend()
            plt.title("Convoluted Wakefield")
            if save_png:
                plt.savefig("convoluted_wakefield.png")
            else:
                plt.show()

    finally:
        # Finalize simulation
        print("Finalizing simulation")
        sim.finalize()


if __name__ == "__main__":
    test_wake()