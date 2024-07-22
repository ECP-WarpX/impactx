from trame.app import get_server
from trame.widgets  import vuetify

from Input.generalFunctions import generalFunctions
from Input.inputParametersCard.inputFunctions import inputFunctions

# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------

server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller

# -----------------------------------------------------------------------------
# Callbacks
# -----------------------------------------------------------------------------

@ctrl.add("on_input_change")
def validate_and_convert_to_correct_type(value, desired_type, state_name, validation_name):
    validation_result = generalFunctions.validate_against(value, desired_type)
    setattr(state, validation_name, validation_result)
    generalFunctions.update_runSimulation_validation_checking()

    if validation_result == []:
        converted_value = generalFunctions.convert_to_correct_type(value, desired_type)
        print(f"{state_name} changed to {converted_value} (type: {type(converted_value)})")
        if getattr(state, state_name) != converted_value:
            setattr(state, state_name, converted_value)    
            if state_name == 'kin_energy':
                state.kin_energy_MeV = inputFunctions.value_of_kin_energy_MeV(converted_value, state.kin_energy_unit)
                print(f"Value of - state.kin_energy_MeV: {state.kin_energy_MeV}")
                
@ctrl.add("kin_energy_unit_change")
def on_convert_kin_energy_change(new_unit):
    old_unit = state.old_kin_energy_unit
    if old_unit != new_unit and float(state.kin_energy) > 0:
        state.kin_energy = inputFunctions.update_kin_energy_on_display(old_unit, new_unit, state.kin_energy)
        state.kin_energy_unit = new_unit
        state.old_kin_energy_unit = new_unit
        state.kin_energy_MeV = inputFunctions.value_of_kin_energy_MeV(float(state.kin_energy), new_unit)
        # print(f"Value of - state.kin_energy_MeV (back-end value): {state.kin_energy_MeV}")
        # print(f"Value of - state.kin_energy (front-end value): {state.kin_energy}")

        # print(f"Units were changed to {new_unit}")
    # print(f"old unit is {old_unit}")
    # print(f"new unit is {new_unit}")

# -----------------------------------------------------------------------------
# Content
# -----------------------------------------------------------------------------

class inputParameters:
    
    def __init__ (self):
        state.particle_shape = 1
        state.npart = 100
        state.kin_energy = 2.0e3
        state.kin_energy_MeV = state.kin_energy
        state.bunch_charge_C = 2e5
        state.kin_energy_unit = "MeV"
        state.old_kin_energy_unit = "MeV"

        state.npart_validation = []
        state.kin_energy_validation = []
        state.bunch_charge_C_validation = []

        print(f"Initial value of - state.kin_energy_MeV: {state.kin_energy_MeV}")

    def card(self):
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
                vuetify.VCombobox(
                    v_model=("particle_shape",),
                    label="Particle Shape",
                    items=([1, 2, 3],),
                    dense=True,
                )
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(cols=12, classes="py-0"):
                        vuetify.VTextField(
                            v_model=("npart",),
                            label="Number of Particles",
                            error_messages=("npart_validation",),
                            change=(ctrl.on_input_change, "[$event, 'int','npart','npart_validation']"),
                            type="number",
                            dense=True,
                        )
                with vuetify.VRow(classes="my-2"):
                    with vuetify.VCol(cols=8, classes="py-0"):
                        vuetify.VTextField(
                            v_model=("kin_energy",),
                            label="Kinetic Energy",
                            error_messages=("kin_energy_validation",),
                            change=(ctrl.on_input_change, "[$event, 'float','kin_energy','kin_energy_validation']"),
                            type="number",
                            dense=True,
                            classes="mr-2",
                        )
                    with vuetify.VCol(cols=4, classes="py-0"):
                        vuetify.VSelect(
                            v_model=("kin_energy_unit",),
                            label="Unit",
                            items=(["meV", "eV", "MeV", "GeV", "TeV"],),
                            change=(ctrl.kin_energy_unit_change, "[$event]"),
                            dense=True,
                        )
                with vuetify.VRow(classes="my-2"):
                    with vuetify.VCol(cols=8, classes="py-0"):
                        vuetify.VTextField(
                            label="Bunch Charge",
                            v_model=("bunch_charge_C",),
                            error_messages=("bunch_charge_C_validation",),
                            change=(ctrl.on_input_change, "[$event, 'float','bunch_charge_C','bunch_charge_C_validation']"),
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