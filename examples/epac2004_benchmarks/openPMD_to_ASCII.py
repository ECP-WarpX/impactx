import openpmd_api as io  # install with python3 -m pip install openpmd-api

# access the initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial = series.iterations[1].particles["beam"].to_df()
final = series.iterations[last_step].particles["beam"].to_df()

# print the intial beam
initial.to_csv("beam_initial.csv", sep=" ", header=True)

# print the final beam
final.to_csv("beam_final.csv", sep=" ", header=True)
