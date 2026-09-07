"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

import keyword
from typing import Union

from ... import state
from ..utils import GeneralFunctions

ALLOWED_INPUT_TYPES = {"int", "float", "str", "bool"}
INT_ERROR_MESSAGE = "Must be an integer"
FLOAT_ERROR_MESSAGE = "Must be a float"
BOOL_ERROR_MESSAGE = "Must be a boolean: True or False"
NON_ZERO_ERROR = "Must be non-zero"
POSITIVE_ERROR = "Must be positive"
NEGATIVE_ERROR = "Must be negative"
N_CELL_MULTIPLE_ERROR = "Must be a multiple of its blocking factor"
PYTHON_IDENTIFIER_ERROR = "Must be a valid Python identifier"

# Utilized in prob_relative validation
GREATER_THAN_THREE_ERROR = "Must be greater than 3"
GREATER_THAN_ONE_ERROR = "Must be greater than 1"
LESS_THAN_PREVIOUS_ERROR = "Must be less than previous value"


class DashboardValidation:
    """
    Contains all validation logic for the ImpactX dashboard inputs.
    """

    @staticmethod
    def update_error_message_on_ui(state_name: str, error_message: str) -> None:
        """
        Sets the error message for a given state input field.

        :param state_name: The name of the state field to attach the error message to.
        :param error_message: The error message to set.
        """
        validation_name = f"{state_name}_error_message"

        if not hasattr(state, validation_name):
            raise AttributeError(
                f"The error message state does not exist: {validation_name}'"
            )

        if getattr(state, validation_name) != error_message:
            setattr(state, validation_name, error_message)

    @staticmethod
    def is_blank(input_value) -> bool:
        """
        Check whether a value was left blank by the user.

        The literal ``"None"`` counts as blank as well: it is the signature
        default that :meth:`InputDefaultsHelper.extract_parameters` reports for
        an optional parameter.

        :param input_value: The value to check.
        :return: True if the value carries no user input.
        """
        return input_value is None or str(input_value).strip() in ("", "None")

    @staticmethod
    def validate(
        input_name: str,
        input_value: Union[float, int],
        category: str | None = None,
        parameter_type: str | None = None,
        optional: bool = False,
    ) -> list[str]:
        """
        Validates the input value against its default type and any additional conditions.

        :param input_name: The name of the input to validate.
        :param input_value: The value to validate.
        :param category: The category of validation (e.g., 'distribution', 'lattice').
        :param parameter_type: The explicit type to use ('int', 'float', 'str', 'bool'). If provided, overrides type lookup.
        :param optional: If True, a blank value is accepted and the parameter is left at its default.
        :return: A list of error messages. An empty list if there are no errors.
        """
        if optional and DashboardValidation.is_blank(input_value):
            return []

        input_type = DashboardValidation._get_input_type(
            input_name, category, parameter_type
        )

        if input_type not in ALLOWED_INPUT_TYPES:
            return [f"Unknown or unsupported type '{input_type}'"]

        if input_type == "str":
            if not DashboardValidation.is_valid_input_name(str(input_value)):
                return [PYTHON_IDENTIFIER_ERROR]
            return []

        if input_type == "bool":
            if str(input_value).strip() in ("True", "False", "true", "false"):
                return []
            return [BOOL_ERROR_MESSAGE]

        numeric_input = GeneralFunctions.convert_to_numeric(input_value)
        type_errors = DashboardValidation._validate_type(numeric_input, input_type)

        if type_errors:
            return type_errors

        additional_validation = DashboardValidation._validate_additional_conditions(
            input_name, numeric_input
        )
        return additional_validation

    @staticmethod
    def _get_input_type(
        input_name: str, category: str | None, parameter_type: str | None
    ) -> str:
        """
        Retrieves the expected input type for validation.

        :param input_name: Name of the parameter.
        :param category: The category of validation (e.g., 'distribution', 'lattice').
        :param parameter_type: The explicit type to use ('int', 'float', 'str'). If provided, overrides type lookup.
        :return: The resolved type as a string ('int', 'float', or 'str').
        """

        if parameter_type:
            return parameter_type

        input_type = None
        if category in {"distribution", "lattice"}:
            input_type = GeneralFunctions.get_default(category, "types")
        else:
            input_type = GeneralFunctions.get_default(input_name, "types")

        return input_type if input_type in ALLOWED_INPUT_TYPES else "str"

    @staticmethod
    def _validate_type(
        numeric_input: Union[float, int, None], value_type: str
    ) -> list[str]:
        """
        Validates a numeric input against the expected type ('int' or 'float').

        :param numeric_input: The value to validate (already converted to numeric or None).
        :param value_type: The expected type ('int' or 'float').
        :return: A list of error messages. Empty if valid.
        """

        if numeric_input is None:
            error_message = (
                INT_ERROR_MESSAGE if value_type == "int" else FLOAT_ERROR_MESSAGE
            )
            return [error_message]

        is_int = isinstance(numeric_input, int)
        is_float = isinstance(numeric_input, (int, float))

        if value_type == "int" and not is_int:
            return [INT_ERROR_MESSAGE]
        elif value_type == "float" and not is_float:
            return [FLOAT_ERROR_MESSAGE]

        return []

    @staticmethod
    def _validate_additional_conditions(
        input_name: str, value: Union[float, int]
    ) -> list[str]:
        """
        Validates additional numeric conditions (e.g., non-zero, positive, negative)
        based on rules defined for the input name.

        :param input_name: Name of the input parameter.
        :param value: Numeric value to validate.
        :return: List of error messages. Empty if all conditions are satisfied.
        """

        lookup_name = "lambda" if "lambda" in input_name else input_name
        additional_conditions = (
            GeneralFunctions.get_default(lookup_name, "validation_condition") or []
        )

        errors = []
        for condition in additional_conditions:
            if condition == "non_zero" and value == 0:
                errors.append(NON_ZERO_ERROR)
            elif condition == "positive" and value <= 0:
                errors.append(POSITIVE_ERROR)
            elif condition == "negative" and value >= 0:
                errors.append(NEGATIVE_ERROR)

        return errors

    @staticmethod
    def is_valid_input_name(user_input: str) -> bool:
        """
        Check if the user input is a valid Python name.
        """
        if user_input is None:
            return True
        return user_input.isidentifier() and not keyword.iskeyword(user_input)

    @staticmethod
    def update_n_cell_validation(direction: str) -> None:
        """
        Validates whether the 'n_cell_<direction>' value is a multiple of 'blocking_factor_<direction>'.

        :param direction: One of 'x', 'y', or 'z' indicating which directional value to check.
        """

        n_cell = GeneralFunctions.convert_to_numeric(
            getattr(state, f"n_cell_{direction}", None)
        )
        blocking_factor = GeneralFunctions.convert_to_numeric(
            getattr(state, f"blocking_factor_{direction}", None)
        )

        if blocking_factor is None or blocking_factor == 0:
            return

        if n_cell % blocking_factor != 0:
            DashboardValidation.update_error_message_on_ui(
                f"n_cell_{direction}", N_CELL_MULTIPLE_ERROR
            )
        else:
            DashboardValidation.update_error_message_on_ui(f"n_cell_{direction}", "")

    @staticmethod
    def validate_prob_relative_fields(index: int, prob_relative_value: float) -> str:
        """
        Validates the prob_relative_fields based on the index and solver type.

        :param index: Index of the modified prob_relative_field.
        :param prob_relative_value: User-provided value to validate.
        :return: Error message if invalid.
        """
        error_message = ""

        try:
            prob_relative_value = float(prob_relative_value)
            poisson_solver = state.poisson_solver

            if index == 0:
                if poisson_solver == "multigrid":
                    if prob_relative_value <= 3:
                        error_message = GREATER_THAN_THREE_ERROR
                elif poisson_solver == "fft":
                    if prob_relative_value <= 1:
                        error_message = GREATER_THAN_ONE_ERROR
            else:
                previous_value = float(state.prob_relative[index - 1])
                if prob_relative_value >= previous_value:
                    error_message = f"{LESS_THAN_PREVIOUS_ERROR} ({previous_value})"
                else:
                    if prob_relative_value <= 1:
                        error_message = GREATER_THAN_ONE_ERROR
        except ValueError:
            error_message = FLOAT_ERROR_MESSAGE

        return error_message
