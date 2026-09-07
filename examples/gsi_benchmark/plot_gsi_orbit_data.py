import openpmd_api as io  # install with python3 -m pip install openpmd-api

# open the data series
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)

# from openpmd_viewer import OpenPMDTimeSeries
import matplotlib.pyplot as plt

# 2. Define the range of iterations to plot (e.g., steps 1000 to 5000 with a period of 1000)
steps = list(series.iterations)

plt.figure(figsize=(8, 6))

mm_scale = 1.0e3
for it in steps:
    # Fetch x and momentum_x for electrons
    beam = series.iterations[it].particles["beam"].to_df()
    x = beam["position_x"] * mm_scale
    px = beam["momentum_x"] * mm_scale

    # Example: Plot x vs. px as a scatter
    plt.scatter(x, px, c="red", s=3)

plt.xlabel("x [mm]", fontsize=12)
plt.ylabel("px/p0 [mrad]", fontsize=12)
plt.title("Test Particles: Horizontal Coordinates")
plt.show()
