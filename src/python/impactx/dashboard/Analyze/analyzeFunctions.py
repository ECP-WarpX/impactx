import asyncio
import subprocess

import pandas as pd
from trame.app import get_server

# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------

server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller

# -----------------------------------------------------------------------------
# Content
# -----------------------------------------------------------------------------


class analyzeFunctions:
    """
    Helper functions for
    preparing the contents for the 'Analyze' page
    """

    # -----------------------------------------------------------------------------
    # Functions for Beam Characteristic and Ref particle data table
    # -----------------------------------------------------------------------------

    @staticmethod
    def load_data(file_path):
        """
        Reads data from the provided file path.
        :param file_path: The path to the file to be read.
        :return: A DataFrame containing the data from the file.
        """

        df = pd.read_csv(file_path, sep=" ")
        return df

    @staticmethod
    def convert_to_dict(combined_data):
        """
        Converts data to dictionary format
        for Vuetify data table.
        :param combined_data: The DataFrame to be converted.
        :return: A tuple containing the data as a list of dictionaries and the headers.
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
        Merges two files together.
        :param file1_name: The name of the first file.
        :param file2_name: The name of the second file.
        :return: A DataFrame containing the merged data from the two files.
        """

        file1 = analyzeFunctions.load_data(file1_name)
        file2 = analyzeFunctions.load_data(file2_name)
        return pd.merge(file1, file2, how="outer")

    @staticmethod
    def filter_headers(allHeaders, selected_headers):
        """
        Retrieves only user-selected headers.
        :param allHeaders: The list of all headers.
        :param selected_headers: The list of headers selected by the user.
        :return: A list of filtered headers based on user selection.
        """

        filtered_headers = []
        for selectedHeader in allHeaders:
            if selectedHeader["value"] in selected_headers:
                filtered_headers.append(selectedHeader)
        return filtered_headers

    @staticmethod
    def filter_data(allData, selected_headers):
        """
        Retrieves data for user-selected headers.
        :param allData: The list of all data.
        :param selected_headers: The list of headers selected by the user.
        :return: A list of filtered data based on user selection.
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
        """
        Function to print out simulation results in
        terminal view. (Not working as intended, 8/4/24)
        """

        ctrl.terminal_println(f"Running {simulation_function_name}...")
        ctrl.terminal_println(
            f"npart: {state.npart}\nkin_energy_MeV: {state.kin_energy_MeV}"
        )

        # Define the command to run based on the simulation function name
        if simulation_function_name == "run_simulation":
            command = [
                "python",
                "-c",
                "from Analyze.plot_PhaseSpaceProjections.phaseSpace import run_simulation; run_simulation()",
            ]
        elif simulation_function_name == "run_optimize_triplet":
            command = [
                "python",
                "-c",
                "from Analyze.plot_PhaseSpaceProjections.phaseSpace import run_optimize_triplet; run_optimize_triplet()",
            ]
        else:
            ctrl.terminal_println(
                f"Unknown simulation function: {simulation_function_name}"
            )
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
