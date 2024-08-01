import pandas as pd

from impactx import distribution
from impactx import elements

from trame.app import get_server
import asyncio
import subprocess
from trame.widgets import xterm

# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------

server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller

# -----------------------------------------------------------------------------
# Content
# -----------------------------------------------------------------------------

class analyzeFunctions:

    # -----------------------------------------------------------------------------
    # Functions for Beam Characteristic and Ref particle data table
    # -----------------------------------------------------------------------------

    @staticmethod
    def load_data(file_path):
        """
        Function to read provided file_path
        """
        df = pd.read_csv(file_path, sep=" ")
        return df

    @staticmethod
    def convert_to_dict(combined_data):
        """
        Function to convert data into dictionary format.
        Used to have correct dataType in Vuetify data table.
        """
        dictionary = combined_data.to_dict(orient="records")
        columns = combined_data.columns
        headers = [
            {"text": column.strip(), "value": column.strip()} for column in columns
        ]
        return dictionary, headers

    @staticmethod
    def combine_files(file1_name, file2_name):
        """
        Function to merge two files together.
        """
        file1 = analyzeFunctions.load_data(file1_name)
        file2 = analyzeFunctions.load_data(file2_name)
        return pd.merge(file1, file2, how="outer")

    @staticmethod
    def filter_headers(allHeaders, selected_headers):
        """
        Function to retrieve only retrieve
        user selected headers
        """
        filtered_headers = []
        for selectedHeader in allHeaders:
            if selectedHeader["value"] in selected_headers:
                filtered_headers.append(selectedHeader)
        return filtered_headers

    @staticmethod
    def filter_data(allData, selected_headers):
        """
        Function to retrieve only retrieve data for
        user selected headers
        """
        filtered_data = []
        for row in allData:
            filtered_row = {}
            for key, value in row.items():
                if key in selected_headers:
                    filtered_row[key] = value
            filtered_data.append(filtered_row)
        return filtered_data

    # -----------------------------------------------------------------------------
    # Function to print simulation output in terminal view
    # -----------------------------------------------------------------------------

    async def outputTerminal(simulation_function_name):
        ctrl.terminal_println(f"Running {simulation_function_name}...")
        ctrl.terminal_println(f"npart: {state.npart}\nkin_energy_MeV: {state.kin_energy_MeV}")

        # Define the command to run based on the simulation function name
        if simulation_function_name == "run_simulation":
            command = ["python", "-c", "from Analyze.plot_phase_space.phaseSpace import run_simulation; run_simulation()"]
        elif simulation_function_name == "run_optimize_triplet":
            command = ["python", "-c", "from Analyze.plot_phase_space.phaseSpace import run_optimize_triplet; run_optimize_triplet()"]
        else:
            ctrl.terminal_println(f"Unknown simulation function: {simulation_function_name}")
            return

        # Run the specified simulation function as a separate process
        process = await asyncio.create_subprocess_exec(
            *command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,  # Capture errors to the same stream as output
        )

        # Read output from the process and print it to the xterm widget
        while True:
            output = await process.stdout.readline()
            if output == b"" and await process.wait() is not None:
                break
            if output:
                ctrl.terminal_println(output.decode().strip())

        ctrl.terminal_println(f"{simulation_function_name} complete.")
