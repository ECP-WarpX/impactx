import math
import matplotlib.pyplot as plt
import numpy as np

import os
from conftest import basepath

from impactx import ImpactX, amr, wakeconvolution

def test_wake(save_png=True):
    """
    Run ImpactX, calculate wakefield convolution from particle data, and plot the output
    """
    sim = ImpactX()

    sim.n_cell = [16, 24, 32]
    sim.load_inputs_file(basepath + "/../../examples/chicane/input_chicane.in") #Specify example path here
    sim.slice_step_diagnostics = False

    sim.init_grids()
    sim.init_beam_distribution_from_inputs()
    sim.init_lattice_elements_from_inputs()

    sim.evolve()

    #Deposit charge
    sim.deposit_charge()

    #Extract charge distribution from simulation
    rho = sim.rho(lev=0)
    rs = rho.sum_unique(comp=0, local=False)
    
    gm = sim.Geom(lev=0)
    dr = gm.data().CellSize()
    dV = np.prod(dr)

    charge_distribution = dV * rs  #Charge [C]

    #Define the parameters for the W_L_CSR function
    R = 10.35  #Curvature radius [m]
    beam_charge = -1.0e-9  #Beam charge [C]

    num_bins = len(charge_distribution)
    bin_size = dr[2]  #Longitudinal direction s [m]
    padding_factor = 1 #Zero-padding multiplicity

    #Generate s range - TODO: Include python bindings for ChargeBinning.cpp functions for the binning here
    s_values = np.linspace((-num_bins // 2) * bin_size, (num_bins // 2) * bin_size, num_bins)

    #Call the C++ function to perform wake function calculation
    wake_function = np.array([wakeconvolution.w_l_csr(s, R, beam_charge) for s in s_values], dtype=np.double)

    #Convert charge distribution to numpy array
    charge_distribution = np.array(charge_distribution, dtype=np.double)
    convoluted_wakefield = np.zeros(num_bins, dtype=np.double)

    #Call the C++ function to perform FFT convolution
    wakeconvolution.convolve_fft(
        charge_distribution,
        wake_function,
        bin_size,
        convoluted_wakefield,
        padding_factor
    )

    #Plot convoluted wakefield
    plt.figure()
    plt.plot(s_values, convoluted_wakefield, label='Convoluted Wakefield')
    plt.xlabel('s [m]')
    plt.ylabel('Wakefield')
    plt.legend()
    plt.title('Convoluted Wakefield')
    if save_png:
        plt.savefig("convoluted_wakefield.png")
    else:
        plt.show()

    #Finalize simulation
    sim.finalize()

if __name__ == "__main__":
    test_wake()