import matplotlib.pyplot as plt
import numpy as np
import openpmd_api as io  # Install with python3 -m pip install openpmd-api
import pandas as pd
from scipy.signal import fftconvolve

"""Create arrays of number, charge density along binned midpoint longitudinal beam bunch position """

# Open data file
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)

# Create 1D longitudinal charge distribution array rho
arr_np = np.sum(arr_np, axis=0)
arr_np = np.sum(arr_np, axis=0)

"""Manual Binning"""

# Get beam parameters at a time the beam passed the beam monitor element
last_step = list(series.iterations)[-1]
beam_initial = series.iterations[1].particles["beam"].to_df()
beam_final = series.iterations[last_step].particles["beam"].to_df()

# Create bins for beam bunch position 't' coordinates and retrieve the bin edges
n_bins = 60

t_bins_initial, bin_edges_initial = pd.cut(
    beam_initial["position_t"], bins=n_bins, retbins=True
)
t_bins_final, bin_edges_final = pd.cut(
    beam_final["position_t"], bins=n_bins, retbins=True
)

# Calculate the midpoints of bins using bin edges
midpoints_initial = 0.5 * (bin_edges_initial[:-1] + bin_edges_initial[1:])
midpoints_final = 0.5 * (bin_edges_final[:-1] + bin_edges_final[1:])

# Assign bins column to the DataFrame
beam_initial["t_bin_initial"] = t_bins_initial
beam_final["t_bin_final"] = t_bins_final

# Find number of particles stored in each bin

# Group by 't' bins and count the number of particles stored in each bin
bin_counts_initial = beam_initial.groupby("t_bin_initial").size()
bin_counts_final = beam_final.groupby("t_bin_final").size()

sigma_t = 0.00019976129161569057  # original value read out from file

# Calculate charge per particle
beam_charge = 1e-9
num_particle = 1e6
charge_per_particle = beam_charge / num_particle

""" Find mean value of x and y for particles in each bin"""

# Calculate the mean of the x, y positions for each bin
mean_positions_initial_x = beam_initial.groupby("t_bin_initial")["position_x"].mean()
mean_positions_initial_y = beam_initial.groupby("t_bin_initial")["position_y"].mean()

# Create arrays of midpoint and corresponding particle charge

beam_initial_t = np.array(midpoints_initial)
beam_initial_num = np.array(bin_counts_initial)
beam_initial_charge = np.array(bin_counts_initial * charge_per_particle)

beam_mean_x = np.array(mean_positions_initial_x)
beam_mean_y = np.array(mean_positions_initial_y)

beam_mean_x[np.isnan(beam_mean_x)] = 0
beam_mean_y[np.isnan(beam_mean_y)] = 0

"""Zero-pad arrays to check for convergence"""

pad_factor = 8  # Padding multiplier - Starts at 1 corresponding to a 2x pad
target_length = pad_factor * (2 * n_bins - 1)

# Extend the original position range symmetrically

# Zero-pad the corresponding particle distribution array to match new position range
original_length = len(beam_initial_t)
total_padding = target_length - original_length

# Ensure padding is split equally on both sides
if total_padding % 2 == 0:
    pad_before = total_padding // 2
    pad_after = total_padding // 2
else:
    pad_before = total_padding // 2
    pad_after = total_padding // 2 + 1

# Determine the interval between existing data points
delta_t = beam_initial_t[1] - beam_initial_t[0]

# Calculate new start and end points
new_start = (
    beam_initial_t[0] - pad_before * delta_t
)  # Adds N x delta_t to the endpoints
new_end = beam_initial_t[-1] + pad_after * delta_t

# Apply the padding
beam_initial_t_padded = np.linspace(new_start, new_end, target_length)
beam_initial_num_padded = np.pad(beam_initial_num, (pad_before, pad_after), "constant")
beam_initial_charge_padded = np.pad(
    beam_initial_charge, (pad_before, pad_after), "constant"
)
beam_mean_x_padded = np.pad(beam_mean_x, (pad_before, pad_after), "constant")
beam_mean_y_padded = np.pad(beam_mean_y, (pad_before, pad_after), "constant")

"""Define Constants"""

c = 3e8  # Speed of light [m/s]
e = 1.602e-19  # Elementary charge [C]
m_e = 9.109e-31  # Electron mass [kg]
Z0 = 377  # Free space impedance [Ohm]
epsilon_0 = 8.854e-12  # Vacuum permittivity [F/m]
alpha_1 = 0.4648

# Parameters for Resistive Wall
a = 5.9e-3  # iris radius [m]
g = 7.49e-3  # gap [m]
L = 8.75e-3  # period length [m]

# Parameters for CSR
R = 10.35  # bend radius [m]

"""Define Wake Functions"""

# Define Unit Step Function

# 1 for s > 0 & 0 for s < 0


def unit_step(s):
    return np.heaviside(s, 0)


# Define Alpha Function


def alpha(s):
    return 1 - alpha_1 * np.sqrt(s) - (1 - 2 * alpha_1) * s


# Define Wake Function for Resistive Wall


def W_T_RF(s, a, g, L):
    s0 = (0.169 * a**1.79 * g**0.38) / L**1.17
    term = np.sqrt(np.abs(s) / s0) * np.exp(-np.sqrt(np.abs(s) / s0))
    return (4 * Z0 * c * s0 * unit_step(s)) / (np.pi * a**4) * term


def W_L_RF(s, a, g, L):
    s00 = g * (a / (alpha(g / L) * L)) ** 2 / 8
    return (Z0 * c * unit_step(s) * np.exp(-np.sqrt(np.abs(s) / s00))) / (np.pi * a**2)


# Define Wake Function for CSR

W_tilde = 0.444717830072860  # known analytical value
N = beam_charge / e
rc = e**2 / (4 * np.pi * epsilon_0 * m_e * c**2)
kappa = (2 * rc * m_e * c**2) / (3 ** (1 / 3) * R ** (2 / 3))


def W_L_CSR(s, R):
    return -(N * kappa * unit_step(s)) / np.abs(s) ** (1 / 3)


"""Calculate Charge & Number of Particle Slopes"""

beam_pdf_nopad = beam_initial_num / delta_t / num_particle
beam_pdf = beam_initial_num_padded / delta_t / num_particle

charge_slopes = []
for i in range(len(beam_initial_charge_padded) - 1):
    charge_slopes.append(
        (beam_initial_charge_padded[i + 1] - beam_initial_charge_padded[i])
        / (beam_initial_t_padded[i + 1] - beam_initial_t_padded[i])
    )

num_slopes = []
for i in range(len(beam_pdf) - 1):
    num_slopes.append(
        (beam_pdf[i + 1] - beam_pdf[i])
        / (beam_initial_t_padded[i + 1] - beam_initial_t_padded[i])
    )

num_slopes_nopad = []
for i in range(len(beam_pdf_nopad) - 1):
    num_slopes_nopad.append(
        (beam_pdf_nopad[i + 1] - beam_pdf_nopad[i])
        / (beam_initial_t[i + 1] - beam_initial_t[i])
    )

charge_slopes = np.array(charge_slopes)
num_slopes = np.array(num_slopes)
num_slopes_nopad = np.array(num_slopes_nopad)

# Create convolution function that takes wake functions as input


def convolve_with_fft(beam_profile, wake_func, delta_t):

    # Calculate the wake function values
    wake_values = wake_func

    # Perform the FFT-based convolution
    convolution_signal = fftconvolve(beam_profile, wake_values)

    return delta_t * convolution_signal


"""Perform FFT Convolution for CSR Wake Function"""


# Try analytical function solutions
def func_3(s):
    # return np.exp(-s**2 / (2 * sigma_t**2)) / (np.sqrt(2 * np.pi) * sigma_t)
    return -s * np.exp(-(s**2) / (2 * sigma_t**2)) / (np.sqrt(2 * np.pi) * sigma_t**3)


def func_4(s):
    return -N * kappa * unit_step(s) / np.abs(s) ** (1 / 3)


# Perform the CSR wakefield convolution
convolved_signal_func = convolve_with_fft(
    func_3(beam_initial_t_padded), func_4(beam_initial_t_padded), delta_t
)
convolved_signal_CSR = convolve_with_fft(
    num_slopes, W_L_CSR(beam_initial_t_padded, R), delta_t
)
convolved_signal_CSR_nopad = convolve_with_fft(
    num_slopes_nopad, W_L_CSR(beam_initial_t, R), delta_t
)

# Renormalize positions by bunch length
lower_bound = 2 * beam_initial_t_padded[0]
upper_bound = 2 * beam_initial_t_padded[-1]
t_range = np.linspace(lower_bound, upper_bound, 2 * len(beam_initial_t_padded) - 1)
normalized_t = t_range / sigma_t

lower_bound_nopad = 2 * beam_initial_t[0]
upper_bound_nopad = 2 * beam_initial_t[-1]
t_range_nopad = np.linspace(
    lower_bound_nopad, upper_bound_nopad, 2 * len(beam_initial_t) - 1
)
normalized_t_nopad = t_range_nopad / sigma_t

# Set x limits to observe the convergence
x_min_limit = beam_initial_t[0] / sigma_t
x_max_limit = beam_initial_t[-1] / sigma_t

# Plot convolved wake function
plt.plot(normalized_t_nopad[:-1], convolved_signal_CSR_nopad * 1e9, label="No Padding")
plt.plot(normalized_t[:-1], convolved_signal_CSR * 1e9, label="Zero Padding")
# plt.plot(normalized_t, convolved_signal_func * 1e9, label = 'Analytical Functions')
plt.xlabel("Longitudinal Position s/sigma")
plt.ylabel("Wakefield Convolution (V pC/mm)")
plt.title("CSR Wake Analytical Function Convolution vs Longitudinal Position")
plt.xlim(left=-5)
plt.xlim(right=5)
plt.grid(True)
plt.legend()
plt.show()

"""Perform FFT Convolution for RF Wake Function"""

# Define input signal for resistive wall wake function
signal_resistive_wall_x = beam_mean_x_padded * beam_pdf
signal_resistive_wall_y = beam_mean_y_padded * beam_pdf

# Perform the resistive wall wakefield convolution
convolved_signal_RF_x = convolve_with_fft(
    signal_resistive_wall_x, W_T_RF(beam_initial_t_padded, a, g, L), delta_t
)
convolved_signal_RF_y = convolve_with_fft(
    signal_resistive_wall_y, W_T_RF(beam_initial_t_padded, a, g, L), delta_t
)
convolved_signal_RF_z = convolve_with_fft(
    beam_initial_num_padded, W_L_RF(beam_initial_t_padded, a, g, L), delta_t
)

# Plot convolved wake function
plt.plot(normalized_t, convolved_signal_RF_z)
plt.xlim(left=x_min_limit)
plt.xlim(right=x_max_limit)
plt.title("Chicane Longitudinal RF Wake Convolution vs Longitudinal Position")
plt.legend()
plt.show()
