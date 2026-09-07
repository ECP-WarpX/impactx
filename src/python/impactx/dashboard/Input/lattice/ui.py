"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from impactx import elements

from ... import ctrl, state, vuetify
from ...Input.components import (
    CardBase,
    CardComponents,
    InputComponents,
    NavigationComponents,
)
from .. import DashboardDefaults
from ..defaults import BEAM_MONITOR_DEFAULT_NAME
from ..defaults_helper import InputDefaultsHelper
from ..validation import DashboardValidation, errors_tracker
from ..visualization.lattice import update_lattice_statistics
from .utils import LatticeConfigurationHelper
from .variable_handler import LatticeVariableHandler

state.lattice_elements_using_variables = {}
state.is_selected_element_invalid = True

LATTICE_ELEMENTS_MODULE_NAME = elements
state.listOfLatticeElementParametersAndDefault = (
    InputDefaultsHelper.class_parameters_with_defaults(LATTICE_ELEMENTS_MODULE_NAME)
)

state.selected_lattice_list = []
state.nslice = ""


def add_lattice_element() -> dict:
    """
    Appends the currently selected lattice element and its parameters to the lattice list.
    """

    selected_lattice = state.selected_lattice
    parameters_data = state.listOfLatticeElementParametersAndDefault.get(
        selected_lattice, []
    )

    parameters = []
    for name, default_value, default_type in parameters_data:
        # extract_parameters() reports the string "None" for a parameter whose
        # signature default is None, and the value None for a parameter without
        # a default. The former is optional, the latter is required.
        is_optional = default_value == "None"

        value = "" if is_optional else default_value
        if default_type == "bool":
            # real boolean for the checkbox v-model (the string "False" is truthy in JS)
            value = str(value).strip().lower() == "true"
        if selected_lattice == "BeamMonitor" and name == "name" and not value:
            value = BEAM_MONITOR_DEFAULT_NAME

        error_message = DashboardValidation.validate(
            name,
            value,
            category="lattice",
            parameter_type=default_type,
            optional=is_optional,
        )

        parameters.append(
            {
                "parameter_name": name,
                "ui_input": value,
                "sim_input": value,
                "parameter_type": default_type,
                "parameter_is_optional": is_optional,
                "parameter_error_message": error_message,
            }
        )

    lattice_element = {
        "name": selected_lattice,
        "parameters": parameters,
    }

    state.selected_lattice_list.append(lattice_element)
    errors_tracker.update_simulation_validation_status()
    update_lattice_statistics()
    return lattice_element


# -----------------------------------------------------------------------------
# Callbacks
# -----------------------------------------------------------------------------


@state.change("selected_lattice_list")
def on_selected_lattice_list_change(selected_lattice_list, **kwargs):
    if selected_lattice_list == []:
        state.isSelectedLatticeListEmpty = "Please select a lattice element"
        errors_tracker.update_simulation_validation_status()
    else:
        state.isSelectedLatticeListEmpty = ""


@state.change("selected_lattice")
def on_lattice_element_name_change(selected_lattice, **kwargs):
    lattice_list = DashboardDefaults.LISTS["lattice_list"]
    state.is_selected_element_invalid = selected_lattice not in lattice_list


@ctrl.add("add_latticeElement")
def on_add_lattice_element_click():
    add_lattice_element()
    state.dirty("selected_lattice_list")


def process_if_variable(index, parameter_name, ui_input, parameter_type):
    """
    If the updated lattice parameter value uses or potentially uses a variable, this
    function returns the simulation value by lookup, adds the element to a dictionary
    which contains current or potential variables, and also returns true or false if it
    is a current or potential variable.

    :param index: The index of the lattice element in the lattice list config.
    :param parameter_name: The specific lattice element parameter name.
    :param ui_input: The value present on the UI end..
    :param parameter_type: The lattice element parameters type.
    """
    ui_input = ui_input.strip()
    sim_input = ui_input
    binding = None

    is_negative_input = ui_input.startswith("-") and ui_input != "-"
    var_name = ui_input[1:] if is_negative_input else ui_input

    is_variable, variable_index = LatticeVariableHandler.determine_if_existing_variable(
        var_name
    )
    is_potential_variable = DashboardValidation.is_valid_input_name(var_name)

    if is_variable:
        sim_value = state.variables[variable_index]["value"]
        if sim_value is not None:
            sim_input = -sim_value if is_negative_input else sim_value

    if is_variable or is_potential_variable:
        binding = {
            "element_reference": state.selected_lattice_list[index],
            "parameter_name": parameter_name,
            "ui_input": ui_input,
            "parameter_type": parameter_type,
        }

    return sim_input, binding


@ctrl.add("updateLatticeElementParameters")
def on_lattice_element_parameter_change(
    index, parameter_name, ui_input, parameter_type
):
    if parameter_type == "bool":
        # booleans come from a checkbox and cannot reference variables
        ui_input = bool(ui_input)
        sim_input, bounded_or_pending_variable = ui_input, None
    else:
        sim_input, bounded_or_pending_variable = process_if_variable(
            index, parameter_name, ui_input, parameter_type
        )

    key = (id(state.selected_lattice_list[index]), parameter_name)
    if bounded_or_pending_variable is not None:
        state.lattice_elements_using_variables[key] = bounded_or_pending_variable
    else:
        state.lattice_elements_using_variables.pop(key, None)

    for param in state.selected_lattice_list[index]["parameters"]:
        if param["parameter_name"] == parameter_name:
            param["ui_input"] = ui_input
            param["sim_input"] = sim_input
            param["parameter_error_message"] = DashboardValidation.validate(
                parameter_name,
                sim_input,
                category="lattice",
                parameter_type=parameter_type,
                optional=param.get("parameter_is_optional", False),
            )

    errors_tracker.update_simulation_validation_status()
    state.dirty("selected_lattice_list")
    update_lattice_statistics()


@ctrl.add("deleteLatticeElement")
def on_delete_LatticeElement_click(index):
    state.selected_lattice_list.pop(index)
    state.dirty("selected_lattice_list")
    update_lattice_statistics()


@ctrl.add("move_latticeElementIndex_up")
def on_move_latticeElementIndex_up_click(index):
    if index > 0:
        state.selected_lattice_list[index], state.selected_lattice_list[index - 1] = (
            state.selected_lattice_list[index - 1],
            state.selected_lattice_list[index],
        )
        state.dirty("selected_lattice_list")
        update_lattice_statistics()


@ctrl.add("move_latticeElementIndex_down")
def on_move_latticeElementIndex_down_click(index):
    if index < len(state.selected_lattice_list) - 1:
        state.selected_lattice_list[index], state.selected_lattice_list[index + 1] = (
            state.selected_lattice_list[index + 1],
            state.selected_lattice_list[index],
        )
        state.dirty("selected_lattice_list")
        update_lattice_statistics()


@ctrl.add("nsliceDefaultChange")
def update_default_value(parameter_name, new_value):
    data = InputDefaultsHelper.class_parameters_with_defaults(elements)

    for key, parameters in data.items():
        for i, param in enumerate(parameters):
            if param[0] == parameter_name:
                parameters[i] = (param[0], new_value, param[2])

    state.listOfLatticeElementParametersAndDefault = data


# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------


class LatticeConfiguration(CardBase):
    HEADER_NAME = "Lattice Configuration"

    def __init__(self):
        super().__init__()

    def init_settings_dialog(self):
        with vuetify.VDialog(
            v_model=("lattice_configuration_dialog_settings", False), width="500px"
        ):
            LatticeConfiguration.dialog_settings()

    def card_content(self):
        self.init_settings_dialog()
        with vuetify.VCard(**self.card_props):
            CardComponents.input_header(
                self.HEADER_NAME,
                additional_components={
                    "end": LatticeConfigurationHelper.settings,
                },
            )
            with vuetify.VCardText(**self.CARD_TEXT_OVERFLOW):
                with vuetify.VRow(**self.ROW_STYLE):
                    with vuetify.VCol(cols=2):
                        InputComponents.text_field(
                            label="Periods",
                        )
                    vuetify.VDivider(vertical=True)
                    with vuetify.VCol(cols=True):
                        InputComponents.autocomplete(
                            label="Select Accelerator Lattice",
                            v_model_name="selected_lattice",
                            items=("lattice_list",),
                            error_messages=("isSelectedLatticeListEmpty",),
                        )
                    with vuetify.VCol(cols="auto"):
                        vuetify.VBtn(
                            "ADD",
                            id="add_lattice_element",
                            color="primary",
                            click=ctrl.add_latticeElement,
                            disabled=("is_selected_element_invalid",),
                        )
                with vuetify.VRow(
                    **self.ROW_STYLE,
                    v_for="(latticeElement, index) in selected_lattice_list",
                    align="center",
                    style="flex-wrap: nowrap;",
                ):
                    with vuetify.VCol(cols="auto"):
                        LatticeConfigurationHelper.move_element_up()
                        LatticeConfigurationHelper.move_element_down()
                        LatticeConfigurationHelper.delete_element()
                    with vuetify.VCol(cols="auto"):
                        vuetify.VChip(
                            text=("latticeElement.name",),
                            style="justify-content: center",
                        )
                    with vuetify.VCol(
                        v_for="(parameter, parameterIndex) in latticeElement.parameters",
                        cols="auto",
                    ):
                        vuetify.VCheckbox(
                            v_if="parameter.parameter_type === 'bool'",
                            label=("parameter.parameter_name",),
                            v_model=("parameter.ui_input",),
                            id=("parameter.parameter_name + (index + 1)",),
                            update_modelValue=(
                                ctrl.updateLatticeElementParameters,
                                "[index, parameter.parameter_name, $event, parameter.parameter_type]",
                            ),
                            error_messages=("parameter.parameter_error_message",),
                            density="comfortable",
                            hide_details="auto",
                        )
                        vuetify.VTextField(
                            v_if="parameter.parameter_type !== 'bool'",
                            label=("parameter.parameter_name",),
                            v_model=("parameter.ui_input",),
                            id=("parameter.parameter_name + (index + 1)",),
                            update_modelValue=(
                                ctrl.updateLatticeElementParameters,
                                "[index, parameter.parameter_name, $event, parameter.parameter_type]",
                            ),
                            error_messages=("parameter.parameter_error_message",),
                            density="comfortable",
                            variant="underlined",
                            style="width: 100px;",
                        )

    @staticmethod
    def defaults_handler():
        """
        Displays the content for the 'Defaults' tab
        in the lattice configuration settings.

        Allows users to pre-determine default values for
        any parameter name. Example: user can set 'nslice' to 25
        and every element added thereafter will have the nslice value
        of 25 as default.
        """
        with vuetify.VCardText():
            with vuetify.VRow():
                with vuetify.VCol(cols=3):
                    InputComponents.text_field(
                        label="nslice",
                        v_model_name="nslice",
                        change=(
                            ctrl.nsliceDefaultChange,
                            "['nslice', $event]",
                        ),
                    )

    # -----------------------------------------------------------------------------
    # Dialogs
    # -----------------------------------------------------------------------------

    @staticmethod
    def dialog_settings():
        """
        Provides controls for lattice element configuration,
        allowing dashboard users to define parameter defaults.
        """
        dialog_name = "lattice_configuration_dialog_tab_settings"

        with NavigationComponents.create_dialog_tabs(
            dialog_name, 2, ["Variables", "Defaults"]
        ):
            with vuetify.VTabsWindow(v_model=(dialog_name, 0)):
                with vuetify.VTabsWindowItem():
                    LatticeVariableHandler.variable_handler()
                with vuetify.VTabsWindowItem():
                    LatticeConfiguration.defaults_handler()
