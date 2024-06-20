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
            basepath + "/../../examples/chicane/input_chicane.in"
        )  # Specify example path here
        sim.slice_step_diagnostics = False

        sim.init_grids()
        sim.init_beam_distribution_from_inputs()
        sim.init_lattice_elements_from_inputs()

        sim.evolve() #Add following script inside of evolve for CSR along lattice

        # Deposit charge
        sim.deposit_charge()

        # Check for CSR elements and perform necessary calculations
        element_has_csr = True
        R = 10.35 #Units [m] TODO: Read in value
        beam_charge = 1.0e-09 #Units [C] TODO: Read in value

        # Enter loop if lattice has bend element
        if element_has_csr:
            # Measure beam size, extract the min, max of particle positions
            pc = sim.particle_container()
            x_min, y_min, t_min, x_max, y_max, t_max = pc.min_and_max_positions()

            # Set parameters for charge deposition
            is_unity_particle_weight = True
            GetNumberDensity = True

            num_bins = 100  # Set resolution
            bin_min = t_min
            bin_max = t_max
            bin_size = (bin_max - bin_min) / num_bins

            # Allocate memory for the charge profile
            charge_distribution = np.zeros(num_bins, dtype=np.double)

            # Perform charge deposition
            wakeconvolution.deposit_charge(
                pc,
                charge_distribution,
                num_bins,
                bin_min,
                bin_size,
                is_unity_particle_weight,
            )

            # Compute the derivative of the charge distribution
            slopes = np.zeros(num_bins - 1, dtype=np.double)
            wakeconvolution.derivative_charge(
                charge_distribution.tolist(),
                slopes.tolist(),
                num_bins,
                bin_size,
                GetNumberDensity,
            )

            # Calculate the CSR wake function
            wake_function = np.array(
                [
                    wakeconvolution.w_l_csr(bin_min + (i * bin_size), R, beam_charge)
                    for i in range(num_bins)
                ],
                dtype=np.double,
            )

            # Perform FFT convolution
            convoluted_wakefield = np.zeros(num_bins - 1, dtype=np.double)
            wakeconvolution.convolve_fft(
                slopes.tolist(),
                wake_function[: num_bins - 1].tolist(),
                bin_size,
                convoluted_wakefield.tolist(),
                padding_factor=1,
            )

            # Plot convoluted wakefield
            s_values = np.linspace(bin_min, bin_max, num_bins)
            plt.figure()
            plt.plot(s_values[:-1], convoluted_wakefield, label="Convoluted Wakefield")
            plt.xlabel("s [m]")
            plt.ylabel("Wakefield")
            plt.legend()
            plt.title("Convoluted Wakefield")
            if save_png:
                plt.savefig("convoluted_wakefield.png")
            else:
                plt.show()

    finally:
        # Finalize simulation
        sim.finalize()
        amr.finalize()


if __name__ == "__main__":
    test_wake()
