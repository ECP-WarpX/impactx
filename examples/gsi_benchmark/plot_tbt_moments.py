import matplotlib.pyplot as plt
import numpy as np
import openpmd_api as io  # install with python3 -m pip install openpmd-api

# open the data series
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)

sigmax_data = []
emittancex_data = []
mm_scale = 1.0e3

for j in list(series.iterations):
    beamj = series.iterations[j].particles["beam"]
    sigma_x = beamj.get_attribute("sigma_x") * mm_scale
    emittance_x = beamj.get_attribute("emittance_x") * mm_scale * mm_scale
    sigmax_data.append(sigma_x)
    emittancex_data.append(emittance_x)

turn_data = np.arange(0, len(list(series.iterations)), 1)
plt.xlabel("turn", fontsize=12)
plt.ylabel("sigma_x [mm]", fontsize=12)
plt.title("Turn-by-Turn Beam Size")
plt.plot(turn_data, sigmax_data)
plt.show()
