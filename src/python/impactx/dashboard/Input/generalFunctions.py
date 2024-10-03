"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

import inspect
import os
import re
import subprocess
import webbrowser

from ..trame_setup import setup_server

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Code
# -----------------------------------------------------------------------------


class generalFunctions:
    @staticmethod
    def documentation(section_name):
        """
        Opens a tab to the specified section link in the documentation.
        :param section_name (str): The name of the documentation section to open.
        """
        url_dict = {
            "LatticeElements": "https://impactx.readthedocs.io/en/latest/usage/python.html#lattice-elements",
            "BeamDistributions": "https://impactx.readthedocs.io/en/latest/usage/python.html#initial-beam-distributions",
            "pythonParameters": "https://impactx.readthedocs.io/en/latest/usage/python.html#general",
            "space_charge_documentation": "https://impactx.readthedocs.io/en/latest/usage/parameters.html#space-charge",
            "CSR": "https://impactx.readthedocs.io/en/latest/usage/parameters.html#coherent-synchrotron-radiation-csr",
        }

        url = url_dict.get(section_name)
        if url is None:
            raise ValueError(f"Invalid section name: {section_name}")

        if "WSL_DISTRO_NAME" in os.environ:
            subprocess.run(["explorer.exe", url])
        else:
            webbrowser.open_new_tab(url)

    # -----------------------------------------------------------------------------
    # Validation functions
    # -----------------------------------------------------------------------------

    @staticmethod
    def determine_input_type(value):
        """
        Determines the type of the input value.
        :param value: The input value whose type needs to be determined.
        :return: A tuple containing the value converted to its determined type and the type itself.
        """

        try:
            return int(value), int
        except ValueError:
            try:
                return float(value), float
            except ValueError:
                return value, str

    @staticmethod
    def validate_against(input_value, value_type, additional_conditions=None):
        """
        Validates the input value against the desired type and additional conditions.
        :param input_value: The value to validate.
        :param value_type: The desired type ('int', 'float', 'str').
        :param conditions: A list of additional conditions to validate.
        :return: A list of error messages. An empty list if there are no errors.
        """
        errors = []
        value = None

        # value_type checking
        if value_type == "int":
            if input_value is None:
                errors.append("Must be an integer")
            else:
                try:
                    value = int(input_value)
                except ValueError:
                    errors.append("Must be an integer")
        elif value_type == "float":
            if input_value is None:
                errors.append("Must be a float")
            else:
                try:
                    value = float(input_value)
                except ValueError:
                    errors.append("Must be a float")
        elif value_type == "str":
            if input_value is None:
                errors.append("Must be a string")
            else:
                value = str(input_value)
        else:
            errors.append("Unknown type")

        # addition_conditions checking
        if errors == [] and additional_conditions:
            for condition in additional_conditions:
                if condition == "non_zero" and value == 0:
                    errors.append("Must be non-zero")
                if condition == "positive" and value <= 0:
                    errors.append("Must be positive")
                if condition == "negative" and value >= 0:
                    errors.append("Must be negative")

        return errors

    @staticmethod
    def update_simulation_validation_status():
        """
        Checks if any input fields are not provided with the correct input type.
        Updates the state to enable or disable the run simulation button.
        """

        error_details = []

        # Check for errors in distribution parameters
        for param in state.selectedDistributionParameters:
            if param["parameter_error_message"]:
                error_details.append(
                    f"{param['parameter_name']}: {param['parameter_error_message']}"
                )

        # Check for errors in lattice parameters
        for lattice in state.selectedLatticeList:
            for param in lattice["parameters"]:
                if param["parameter_error_message"]:
                    error_details.append(
                        f"Lattice {lattice['name']} - {param['parameter_name']}: {param['parameter_error_message']}"
                    )

        # Check for errors in input card
        if state.npart_validation:
            error_details.append(f"Number of Particles: {state.npart_validation}")
        if state.kin_energy_validation:
            error_details.append(f"Kinetic Energy: {state.kin_energy_validation}")
        if state.bunch_charge_C_validation:
            error_details.append(f"Bunch Charge: {state.bunch_charge_C_validation}")
        if state.charge_qe_validation:
            error_details.append(f"Ref. Particle Charge: {state.charge_qe_validation}")
        if state.mass_MeV_validation:
            error_details.append(f"Ref. Particle Mass: {state.mass_MeV}")

        if state.selectedLatticeList == []:
            error_details.append("LatticeListIsEmpty")

        state.disableRunSimulationButton = bool(error_details)

    # -----------------------------------------------------------------------------
    # Class, parameter, default value, and default type retrievals
    # -----------------------------------------------------------------------------

    @staticmethod
    def find_classes(module_name):
        """
        Returns a list of all classes in the given module.
        :param module_name: The module to inspect.
        :return: A list of tuples containing class names.
        """

        results = []
        for name in dir(module_name):
            attr = getattr(module_name, name)
            if inspect.isclass(attr):
                results.append((name, attr))
        return results

    @staticmethod
    def find_init_docstring_for_classes(classes):
        """
        Retrieves the __init__ docstring of the given classes.
        :param classes: A list of typles containing class names.
        :return: A dictionary with class names as keys and their __init__ docstrings as values.
        """

        if not isinstance(classes, (list, tuple)):
            raise TypeError("The 'classes' argument must be a list or tuple.")

        docstrings = {}
        for name, cls in classes:
            init_method = getattr(cls, "__init__", None)
            if init_method:
                docstring = cls.__init__.__doc__
                docstrings[name] = docstring
        return docstrings

    @staticmethod
    def extract_parameters(docstring):
        """
        Parses specific information from docstrings.
        Aimed to retrieve parameter names, values, and types.
        :param docstring: The docstring to parse.
        :return: A list of tuples containing parameter names, default values, and types.
        """

        parameters = []
        docstring = re.search(r"\((.*?)\)", docstring).group(
            1
        )  # Return class name and init signature
        docstring = docstring.split(",")

        for parameter in docstring:
            if parameter.startswith("self"):
                continue

            name = parameter
            default = None
            parameter_type = "Any"

            if ":" in parameter:
                split_by_semicolon = parameter.split(":", 1)
                name = split_by_semicolon[0].strip()
                type_and_default = split_by_semicolon[1].strip()
                if "=" in type_and_default:
                    split_by_equals = type_and_default.split("=", 1)
                    parameter_type = split_by_equals[0].strip()
                    default = split_by_equals[1].strip()
                    if default.startswith("'") and default.endswith("'"):
                        default = default[1:-1]
                else:
                    parameter_type = type_and_default

            parameters.append((name, default, parameter_type))

        return parameters

    @staticmethod
    def class_parameters_with_defaults(module_name):
        """
        Given a module name, outputs a dictionary of class names and their parameters.
        Keys are class names, and values are lists of parameter information (name, default value, type).
        :param module_name: The module to inspect.
        :return: A dictionary with class names as keys and parameter information as values.
        """

        classes = generalFunctions.find_classes(module_name)
        docstrings = generalFunctions.find_init_docstring_for_classes(classes)

        result = {}

        for class_name, docstring in docstrings.items():
            parameters = generalFunctions.extract_parameters(docstring)
            result[class_name] = parameters

        return result

    @staticmethod
    def select_classes(module_name):
        """
        Given a module name, outputs a list of all class names in the module.
        :param module_name: The module to inspect.
        :return: A list of class names.
        """

        return list(generalFunctions.class_parameters_with_defaults(module_name))

    @staticmethod
    def convert_to_correct_type(value, desired_type):
        """
        Converts the given value to the desired type.
        :param value: The value to convert.
        :param desired_type: The type to convert the value to ('int', 'float', 'str').
        :return: The value converted to the desired type.
        """

        if value is None:
            raise ValueError("Cannot convert to desired type")
        if desired_type == "int":
            return int(value)
        elif desired_type == "float":
            return float(value)
        elif desired_type == "str":
            return str(value)
        else:
            raise ValueError("Unknown type")
