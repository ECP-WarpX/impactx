
import webbrowser
import subprocess
import os

import inspect
import re

from trame.app import get_server
from impactx import distribution, elements

# -----------------------------------------------------------------------------
# Server setup
# -----------------------------------------------------------------------------

server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller

# -----------------------------------------------------------------------------
# Code
# -----------------------------------------------------------------------------
ANSI_RED = "\033[91m"
ANSI_RESET = "\033[0m"

class generalFunctions:
    
    def documentation(section_name):
        """
        Function that opens tab to section_name link
        """
        if section_name == "LatticeElements":
            url = "https://impactx.readthedocs.io/en/latest/usage/python.html#lattice-elements"
        elif section_name == "BeamDistributions":
            url = "https://impactx.readthedocs.io/en/latest/usage/python.html#initial-beam-distributions"
        elif section_name == "pythonParameters":
            url = "https://impactx.readthedocs.io/en/latest/usage/python.html#general"
        else:
            raise ValueError(f"Invalid section name: {section_name}")
        
        if 'WSL_DISTRO_NAME' in os.environ:
            subprocess.run(['explorer.exe', url])
        else:
            webbrowser.open_new_tab(url)

# -----------------------------------------------------------------------------
# Validation functions
# -----------------------------------------------------------------------------

    def determine_input_type(value):
        """"
        Used to help find out the value type
        """
        try:
            return int(value), int
        except ValueError:
            try:
                return float(value), float
            except ValueError:
                return value, str
        
    def validate_against(input_value, value_type):
        """
        Function which returns error message if
        input value type does not match desired
        type.
        """
        if value_type == "int":
            if input_value is None:
                return ["Must be an integer"]
            try:
                value = int(input_value)
                if value < -10000000000:
                    return ["Must be positive"]
                return []
            except ValueError:
                return ["Must be an integer"]

        elif value_type == "float":
            if input_value is None:
                return ["Must be a float"]
            try:
                value = float(input_value)
                if value < -1000000000:
                    return ["Must be positive"]
                return []
            except ValueError:
                return ["Must be a float"]

        elif value_type == "str":
            if input_value is None:
                return ["Must be a string"]
            return []

        else:
            return ["Unknown type"]
        
    def update_runSimulation_validation_checking():
        """
        Function to check if any input fields are not
        provided with the correct input type.
        Updates states as True or False given result.
        """
        error_details = []

        # Check for errors in distribution parameters
        for param in state.selectedDistributionParameters:
            if param["parameter_error_message"]:
                error_details.append(f"{param['parameter_name']}: {param['parameter_error_message']}")

        # Check for errors in lattice parameters
        for lattice in state.selectedLatticeList:
            for param in lattice['parameters']:
                if param['parameter_error_message']:
                    error_details.append(f"Lattice {lattice['name']} - {param['parameter_name']}: {param['parameter_error_message']}")

        # Check for errors in input card
        if state.npart_validation:
            error_details.append(f"Number of Particles: {state.npart_validation}")
        if state.kin_energy_validation:
            error_details.append(f"Kinetic Energy: {state.kin_energy_validation}")
        if state.bunch_charge_C_validation:
            error_details.append(f"Bunch Charge: {state.bunch_charge_C_validation}")
        if state.selectedLatticeList == []:
            error_details.append("LatticeListIsEmpty")

        # Print all collected error messages
        # if error_details:
        #     print("Errors detected in the following parameters:")
        #     for error in error_details:
        #         print(ANSI_RED + error + ANSI_RESET)

        state.disableRunSimulationButton = bool(error_details)

# -----------------------------------------------------------------------------
# Class, parameter, default value, and default type retrievals
# -----------------------------------------------------------------------------

    def findAllClasses(module_name):
        """
        Returns list of all classes in given module_name
        """
        results = []
        for name in dir(module_name):
            attr = getattr(module_name, name)
            if inspect.isclass(attr):
                results.append((name, attr))
        return results


    def findInitDocstringForClasses(classes):
        """
        Retrieves the __init__ docstring of given classes
        """
        docstrings = {}
        for name, cls in classes:
            init_method = getattr(cls, '__init__', None)
            if init_method:
                docstring = cls.__init__.__doc__
                docstrings[name] = docstring
        return docstrings

    def extractParameters(docstring):
        """
        Parses specific information from docstrings.
        Aimed to retrieve parameter names/values/types.
        """
        parameters = []
        docstring = re.search(r'\((.*?)\)', docstring).group(1)  # Return class name and init signature
        docstring = docstring.split(',')

        for parameter in docstring:
            if parameter.startswith('self'):
                continue
            
            name = parameter
            default = None
            parameter_type = 'Any' 

            if ':' in parameter:
                split_by_semicolon = parameter.split(':', 1)
                name = split_by_semicolon[0].strip()
                type_and_default = split_by_semicolon[1].strip()
                if '=' in type_and_default:
                    split_by_equals = type_and_default.split('=', 1)
                    parameter_type = split_by_equals[0].strip()
                    default = split_by_equals[1].strip()
                    if (default.startswith("'") and default.endswith("'")):
                        default = default[1:-1]
                else:
                    parameter_type = type_and_default

            parameters.append((name, default, parameter_type))

        return parameters
    
    def classAndParametersAndDefaultValueAndType(module_name):
        """
        Given module_name, outputs a dictionary.
        Keys are each class name of the module_name
        Values are dictionaries of parameter information,
        such as default name/value/type.
        """
        
        classes = generalFunctions.findAllClasses(module_name)
        docstrings = generalFunctions.findInitDocstringForClasses(classes)

        result = {}

        for class_name, docstring in docstrings.items():
            parameters = generalFunctions.extractParameters(docstring)
            result[class_name] = parameters

        return result

    def selectClasses(module_name):
        """
        Given module_name, outputs a list
        of all class names in module_name.
        """
        return list(generalFunctions.classAndParametersAndDefaultValueAndType(module_name))

    def convert_to_correct_type(value, desired_type):
        """
        Converts given value to the desired_type.
        Used for input parameters.
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