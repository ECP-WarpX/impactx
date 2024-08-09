"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from trame.widgets import vuetify

from impactx import distribution

from ...trame_setup import setup_server
from ..generalFunctions import generalFunctions

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Helpful
# -----------------------------------------------------------------------------

DISTRIBUTIONS_MODULE_NAME = distribution

state.listOfDistributions = generalFunctions.select_classes(DISTRIBUTIONS_MODULE_NAME)
state.listOfDistributionsAndParametersAndDefault = (
    generalFunctions.class_parameters_with_defaults(DISTRIBUTIONS_MODULE_NAME)
)

# -----------------------------------------------------------------------------
# Default
# -----------------------------------------------------------------------------

state.selectedDistribution = "Waterbag"
state.selectedDistributionType = "Native"
state.selectedDistributionParameters = []

# -----------------------------------------------------------------------------
# Main Functions
# -----------------------------------------------------------------------------


def populate_distribution_parameters(selectedDistribution):
    """
    Populates distribution parameters based on the selected distribution.
    :param selectedDistribution (str): The name of the selected distribution
        whos parameters need to be populated.
    """

    selectedDistributionParameters = (
        state.listOfDistributionsAndParametersAndDefault.get(selectedDistribution, [])
    )

    state.selectedDistributionParameters = [
        {
            "parameter_name": parameter[0],
            "parameter_default_value": parameter[1],
            "parameter_type": parameter[2],
            "parameter_error_message": generalFunctions.validate_against(
                parameter[1], parameter[2]
            ),
        }
        for parameter in selectedDistributionParameters
    ]

    generalFunctions.update_simulation_validation_status()
    return selectedDistributionParameters


def update_distribution_parameters(
    parameterName, parameterValue, parameterErrorMessage
):
    """
    Updates the value of a distribution parameter and its error message.

    :param parameterName (str): The name of the parameter to update.
    :param parameterValue: The new value for the parameter.
    :param parameterErrorMessage: The error message related to the parameter's value.
    """

    for param in state.selectedDistributionParameters:
        if param["parameter_name"] == parameterName:
            param["parameter_default_value"] = parameterValue
            param["parameter_error_message"] = parameterErrorMessage

    generalFunctions.update_simulation_validation_status()
    state.dirty("selectedDistributionParameters")


# -----------------------------------------------------------------------------
# Write to file functions
# -----------------------------------------------------------------------------


def parameter_input_checker():
    """
    Helper function to check if user input is valid.
    :return: A dictionary with parameter names as keys and their validated values.
    """

    parameter_input = {}
    for param in state.selectedDistributionParameters:
        if param["parameter_error_message"] == []:
            parameter_input[param["parameter_name"]] = float(
                param["parameter_default_value"]
            )
        else:
            parameter_input[param["parameter_name"]] = 0.0

    return parameter_input


def distribution_parameters():
    """
    Writes user input for distribution parameters in suitable format for simulation code.
    :return: An instance of the selected distribution class, initialized with user-provided parameters.
    """

    distribution_name = state.selectedDistribution
    parameters = parameter_input_checker()

    distr = getattr(distribution, distribution_name)(**parameters)
    return distr


# -----------------------------------------------------------------------------
# Callbacks
# -----------------------------------------------------------------------------


@state.change("selectedDistribution")
def on_distribution_name_change(selectedDistribution, **kwargs):
    populate_distribution_parameters(selectedDistribution)


@state.change("selectedDistributionType")
def on_distribution_type_change(**kwargs):
    populate_distribution_parameters(state.selectedDistribution)


@ctrl.add("updateDistributionParameters")
def on_distribution_parameter_change(parameter_name, parameter_value, parameter_type):
    parameter_value, input_type = generalFunctions.determine_input_type(parameter_value)
    error_message = generalFunctions.validate_against(parameter_value, parameter_type)

    update_distribution_parameters(parameter_name, parameter_value, error_message)


# -----------------------------------------------------------------------------
# Content
# -----------------------------------------------------------------------------


class DistributionParameters:
    """
    User-Input section for beam distribution.
    """

    @staticmethod
    def card():
        """
        Creates UI content for beam distribution.
        """

        with vuetify.VCard(style="width: 340px; height: 300px"):
            with vuetify.VCardTitle("Distribution Parameters"):
                vuetify.VSpacer()
                vuetify.VIcon(
                    "mdi-information",
                    style="color: #00313C;",
                    click=lambda: generalFunctions.documentation("BeamDistributions"),
                )
            vuetify.VDivider()
            with vuetify.VCardText():
                with vuetify.VRow():
                    with vuetify.VCol(cols=8):
                        vuetify.VCombobox(
                            label="Select Distribution",
                            v_model=("selectedDistribution",),
                            items=("listOfDistributions",),
                            dense=True,
                        )
                    with vuetify.VCol(cols=4):
                        vuetify.VSelect(
                            v_model=("selectedDistributionType",),
                            label="Type",
                            items=(["Native", "Twiss"],),
                            # change=(ctrl.kin_energy_unit_change, "[$event]"),
                            dense=True,
                            disabled=True,
                        )
                with vuetify.VRow(classes="my-2"):
                    for i in range(3):
                        with vuetify.VCol(cols=4, classes="py-0"):
                            with vuetify.VRow(
                                v_for="(parameter, index) in selectedDistributionParameters"
                            ):
                                with vuetify.VCol(
                                    v_if=f"index % 3 == {i}", classes="py-1"
                                ):
                                    vuetify.VTextField(
                                        label=("parameter.parameter_name",),
                                        v_model=("parameter.parameter_default_value",),
                                        change=(
                                            ctrl.updateDistributionParameters,
                                            "[parameter.parameter_name, $event, parameter.parameter_type]",
                                        ),
                                        error_messages=(
                                            "parameter.parameter_error_message",
                                        ),
                                        type="number",
                                        dense=True,
                                    )
