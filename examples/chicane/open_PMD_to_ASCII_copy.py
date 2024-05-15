import openpmd_api as io  # install with python3 -m pip install openpmd-api

# open the data series
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)

# print all steps found
print(list(series.iterations))

# select the last time the beam passed the beam monitor element
last_step = list(series.iterations)[-1]
beam_initial = series.iterations[1].particles["beam"].to_df()
beam_final = series.iterations[last_step].particles["beam"].to_df()

beam_initial.to_csv("beam_initial.csv", sep=" ", header=True)
beam_final.to_csv("beam_final.csv", sep=" ", header=True)
