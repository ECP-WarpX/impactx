"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from trame.widgets import vuetify

from ...trame_setup import setup_server
from ..generalFunctions import generalFunctions
from .inputFunctions import InputFunctions

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Callbacks
# -----------------------------------------------------------------------------


@ctrl.add("on_input_change")
def validate_and_convert_to_correct_type(
    value, desired_type, state_name, validation_name
):
    validation_result = generalFunctions.validate_against(value, desired_type)
    setattr(state, validation_name, validation_result)
    generalFunctions.update_simulation_validation_status()

    if validation_result == []:
        converted_value = generalFunctions.convert_to_correct_type(value, desired_type)

        if getattr(state, state_name) != converted_value:
            setattr(state, state_name, converted_value)
            if state_name == "kin_energy":
                state.kin_energy_MeV = InputFunctions.value_of_kin_energy_MeV(
                    converted_value, state.kin_energy_unit
                )


@ctrl.add("kin_energy_unit_change")
def on_convert_kin_energy_change(new_unit):
    old_unit = state.old_kin_energy_unit
    if old_unit != new_unit and float(state.kin_energy) > 0:
        state.kin_energy = InputFunctions.update_kin_energy_on_display(
            old_unit, new_unit, state.kin_energy
        )
        state.kin_energy_unit = new_unit
        state.old_kin_energy_unit = new_unit
        state.kin_energy_MeV = InputFunctions.value_of_kin_energy_MeV(
            float(state.kin_energy), new_unit
        )


# -----------------------------------------------------------------------------
# Content
# -----------------------------------------------------------------------------


class InputParameters:
    """
    User-Input section for beam properties.
    """

    def __init__(self):
        state.particle_shape = 2
        state.npart = 1000
        state.kin_energy = 2.0e3
        state.kin_energy_MeV = state.kin_energy
        state.bunch_charge_C = 1.0e-9
        state.kin_energy_unit = "MeV"
        state.old_kin_energy_unit = "MeV"

        state.npart_validation = []
        state.kin_energy_validation = []
        state.bunch_charge_C_validation = []

    def card(self):
        """
        Creates UI content for beam properties.
        """

        with vuetify.VCard(style="width: 340px; height: 300px"):
            with vuetify.VCardTitle("Input Parameters"):
                vuetify.VSpacer()
                vuetify.VIcon(
                    "mdi-information",
                    style="color: #00313C;",
                    click=lambda: generalFunctions.documentation("pythonParameters"),
                )
            vuetify.VDivider()
            with vuetify.VCardText():
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(cols=12, classes="py-0"):
                        vuetify.VCheckbox(
                            label="Space Charge",
                            v_model=("space_charge", False),
                            dense=True,
                            classes="mt-0",
                            style="margin-bottom: 0.5px",
                        )
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(cols=12, classes="py-0"):
                        vuetify.VTextField(
                            v_model=("npart",),
                            label="Number of Particles",
                            error_messages=("npart_validation",),
                            change=(
                                ctrl.on_input_change,
                                "[$event, 'int','npart','npart_validation']",
                            ),
                            type="number",
                            dense=True,
                        )
                with vuetify.VRow(classes="my-2"):
                    with vuetify.VCol(cols=8, classes="py-0"):
                        vuetify.VTextField(
                            v_model=("kin_energy",),
                            label="Kinetic Energy",
                            error_messages=("kin_energy_validation",),
                            change=(
                                ctrl.on_input_change,
                                "[$event, 'float','kin_energy','kin_energy_validation']",
                            ),
                            type="number",
                            dense=True,
                            classes="mr-2",
                        )
                    with vuetify.VCol(cols=4, classes="py-0"):
                        vuetify.VSelect(
                            v_model=("kin_energy_unit",),
                            label="Unit",
                            items=(["meV", "eV", "keV", "MeV", "GeV", "TeV"],),
                            change=(ctrl.kin_energy_unit_change, "[$event]"),
                            dense=True,
                        )
                with vuetify.VRow(classes="my-2"):
                    with vuetify.VCol(cols=8, classes="py-0"):
                        vuetify.VTextField(
                            label="Bunch Charge",
                            v_model=("bunch_charge_C",),
                            error_messages=("bunch_charge_C_validation",),
                            change=(
                                ctrl.on_input_change,
                                "[$event, 'float','bunch_charge_C','bunch_charge_C_validation']",
                            ),
                            type="number",
                            dense=True,
                        )
                    with vuetify.VCol(cols=4, classes="py-0"):
                        vuetify.VTextField(
                            label="Unit",
                            value="C",
                            dense=True,
                            disabled=True,
                        )
