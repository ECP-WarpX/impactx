from trame.app import get_server
from trame.widgets import vuetify

from Input.generalFunctions import generalFunctions
from impactx import elements

# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------

server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller

# -----------------------------------------------------------------------------
# Helpful
# -----------------------------------------------------------------------------

LATTICE_ELEMENTS_MODULE_NAME = elements

state.listOfLatticeElements = generalFunctions.selectClasses(LATTICE_ELEMENTS_MODULE_NAME)
state.listOfLatticeElementParametersAndDefault = generalFunctions.classAndParametersAndDefaultValueAndType(LATTICE_ELEMENTS_MODULE_NAME)

# -----------------------------------------------------------------------------
# Default
# -----------------------------------------------------------------------------

state.selectedLattice = None #  Selected lattice is Empty by default
state.selectedLatticeList = [] # Selected lattice list is Empty by default
state.nsliceDefaultValue = None

# -----------------------------------------------------------------------------
# Main Functions
# -----------------------------------------------------------------------------

def add_lattice_element():
    selectedLattice = state.selectedLattice
    selectedLatticeParameters = state.listOfLatticeElementParametersAndDefault.get(selectedLattice, [])

    selectedLatticeElement = {
        "name": selectedLattice,
        "parameters": [
            {"parameter_name": parameter[0], 
             "parameter_default_value": parameter[1], 
             "parameter_type": parameter[2], 
             "parameter_error_message": generalFunctions.validate_against(parameter[1], parameter[2]),
            }
            for parameter in selectedLatticeParameters
        ]
    }

    state.selectedLatticeList.append(selectedLatticeElement)
    generalFunctions.update_runSimulation_validation_checking()
    return selectedLatticeElement
 
def update_latticeElement_parameters(index, parameterName, parameterValue, parameterErrorMessage):
    """
    Updates parameter value and includes error message if user input is not valid
    """
    for param in state.selectedLatticeList[index]["parameters"]:
        if param["parameter_name"] == parameterName:
            param["parameter_default_value"] = parameterValue
            param["parameter_error_message"] = parameterErrorMessage

    generalFunctions.update_runSimulation_validation_checking()
    state.dirty("selectedLatticeList")
    save_latticeElements_to_file()

# -----------------------------------------------------------------------------
# Write to file functions
# -----------------------------------------------------------------------------

def parameter_input_checker_for_lattice(latticeElement):
    """
    Helper function to check if user input is valid, if yes, then will update with value, if not then set to None.
    """
    parameter_input = {}
    for parameter in latticeElement["parameters"]:
        if parameter["parameter_error_message"] == []:
            if parameter["parameter_type"] == "str":
                parameter_input[parameter["parameter_name"]] = f"'{parameter["parameter_default_value"]}'"
            else:
                parameter_input[parameter["parameter_name"]] = parameter["parameter_default_value"]
        else:
            parameter_input[parameter["parameter_name"]] = None

    return parameter_input

def save_latticeElements_to_file():
    """
    Writes users input for lattice element parameters into file in simulation code format
    """
    with open("output_latticeElements_parameters.txt", "w") as file:
        file.write("latticeElements = [\n")
        for latticeElement in state.selectedLatticeList:
            latticeElement_name = latticeElement["name"]
            parameters  =  parameter_input_checker_for_lattice(latticeElement)

            param_values = ", ".join(f"{value}" for value in parameters.values())
            file.write(f"    elements.{latticeElement_name}({param_values}),\n")

        file.write("]\n")   

# -----------------------------------------------------------------------------
# Callbacks
# -----------------------------------------------------------------------------
@state.change("selectedLatticeList")
def on_selectedLatticeList_change(selectedLatticeList, **kwargs):
    if selectedLatticeList == []:
        state.isSelectedLatticeListEmpty = "Please select a lattice element"
        generalFunctions.update_runSimulation_validation_checking()
    else:
        state.isSelectedLatticeListEmpty = ""

@state.change("selectedLattice")
def on_lattice_element_name_change(selectedLattice, **kwargs):
    return
    # print (f"Lattice Selection Changed to: {selectedLattice}")

@ctrl.add("add_latticeElement")
def on_add_lattice_element_click():
    selectedLattice = state.selectedLattice
    if selectedLattice:
        add_lattice_element()
        save_latticeElements_to_file()
        state.dirty("selectedLatticeList")
        # print(f"ADD button clicked, added: {selectedLattice}")
        # print(f"Current list of selected lattice elements: {state.selectedLatticeList}")

@ctrl.add("updateLatticeElementParameters")
def on_lattice_element_parameter_change(index, parameter_name, parameter_value, parameter_type):
    parameter_value, input_type = generalFunctions.determine_input_type(parameter_value)
    error_message = generalFunctions.validate_against(parameter_value, parameter_type)

    update_latticeElement_parameters(index, parameter_name, parameter_value, error_message)
    print(f"Lattice element {index}, {parameter_name} changed to {parameter_value} (type: {input_type})")

@ctrl.add("clear_latticeElements")
def on_clear_lattice_element_click():
    state.selectedLatticeList = []
    save_latticeElements_to_file()

@ctrl.add("deleteLatticeElement")
def on_delete_LatticeElement_click(index):
    state.selectedLatticeList.pop(index)
    state.dirty("selectedLatticeList")
    save_latticeElements_to_file()

@ctrl.add("move_latticeElementIndex_up")
def on_move_latticeElementIndex_up_click(index):
    if index > 0:
        state.selectedLatticeList[index], state.selectedLatticeList[index - 1] = state.selectedLatticeList [index - 1], state.selectedLatticeList [index]
        state.dirty("selectedLatticeList")
        save_latticeElements_to_file()

@ctrl.add("move_latticeElementIndex_down")
def on_move_latticeElementIndex_up_click(index):
    if index < len(state.selectedLatticeList) - 1:
        state.selectedLatticeList[index], state.selectedLatticeList[index + 1] = state.selectedLatticeList [index + 1], state.selectedLatticeList [index]
        state.dirty("selectedLatticeList")
        save_latticeElements_to_file()

@ctrl.add("nsliceDefaultChange")
def update_default_value(parameter_name, new_value):
    data = generalFunctions.classAndParametersAndDefaultValueAndType(elements)
    
    for key, parameters in data.items():
        for i, param in enumerate(parameters):
            if param[0] == parameter_name:
                parameters[i] = (param[0], new_value, param[2])
    
    state.listOfLatticeElementParametersAndDefault = data
# -----------------------------------------------------------------------------
# ContentSetup
# -----------------------------------------------------------------------------

class latticeConfiguration:

    def card():
        with vuetify.VDialog(v_model=("showDialog", False), width="1200px"):
            latticeConfiguration.dialog_lattice_elementList()
        
        with vuetify.VDialog(v_model=("showDialog_settings", False), width="300px"):
            latticeConfiguration.dialog_lattice_settings()
            
        with vuetify.VCard(style="width: 696px;"):
            with vuetify.VCardTitle("Lattice Configuration"):
                vuetify.VSpacer()
                vuetify.VIcon(
                    "mdi-information",
                    classes="ml-2",
                    click=lambda: generalFunctions.documentation("LatticeElements"),
                    style="color: #00313C;",
                )
            vuetify.VDivider()
            with vuetify.VCardText():
                with vuetify.VRow(align="center", no_gutters=True):
                    with vuetify.VCol(cols=8):
                        vuetify.VCombobox(
                            label="Select Accelerator Lattice",
                            v_model=("selectedLattice", None),
                            items=("listOfLatticeElements",),
                            error_messages=("isSelectedLatticeListEmpty",),
                            dense=True,
                            classes="mr-2 pt-6"
                        )
                    with vuetify.VCol(cols="auto"):
                        vuetify.VBtn(
                            "ADD",
                            color="primary",
                            dense=True,
                            classes="mr-2",
                            click=ctrl.add_latticeElement,
                        )
                    with vuetify.VCol(cols="auto"):
                        vuetify.VBtn(
                            "CLEAR",
                            color="secondary",
                            dense=True,
                            classes="mr-2",
                            click=ctrl.clear_latticeElements,
                        )
                    with vuetify.VCol(cols="auto"):
                        vuetify.VIcon(
                            "mdi-cog",
                            click="showDialog_settings = true",
                        )
                with vuetify.VRow():
                    with vuetify.VCol():       
                        with vuetify.VCard(style="height: 300px; width: 700px; overflow-y: auto;"):
                            with vuetify.VCardTitle("Elements", classes="text-subtitle-2 pa-3"):
                                vuetify.VSpacer()
                                vuetify.VIcon(
                                    "mdi-arrow-expand",
                                    color="primary",
                                    click="showDialog = true",
                                )
                            vuetify.VDivider()
                            with vuetify.VContainer(fluid=True):
                                with vuetify.VRow(v_for="(latticeElement, index) in selectedLatticeList", align="center", no_gutters=True, style="min-width: 1500px;"):
                                    with vuetify.VCol(cols="auto", classes="pa-2"):
                                        vuetify.VIcon(
                                            "mdi-menu-up",
                                            click=(ctrl.move_latticeElementIndex_up, "[index]"),
                                        )
                                        vuetify.VIcon(
                                            "mdi-menu-down",
                                            click=(ctrl.move_latticeElementIndex_down, "[index]"),
                                        )
                                        vuetify.VIcon(
                                            "mdi-delete",
                                            click=(ctrl.deleteLatticeElement,"[index]"),
                                        )
                                        vuetify.VChip(
                                            v_text=("latticeElement.name",),
                                            dense=True,
                                            classes="mr-2",
                                            style="justify-content: center"
                                        )
                                    with vuetify.VCol(v_for="(parameter, parameterIndex) in latticeElement.parameters", cols="auto", classes="pa-2"):
                                        vuetify.VTextField(
                                            label=("parameter.parameter_name",),
                                            v_model=("parameter.parameter_default_value",),
                                            change=(ctrl.updateLatticeElementParameters, "[index, parameter.parameter_name, $event, parameter.parameter_type]"),
                                            error_messages=("parameter.parameter_error_message",),
                                            dense=True,
                                            style="width: 100px;"
                                        )

    def dialog_lattice_elementList():
        with vuetify.VCard():
            with vuetify.VCardTitle("Elements", classes="text-subtitle-2 pa-3"):
                vuetify.VSpacer()
            vuetify.VDivider()
            with vuetify.VContainer(fluid=True):
                with vuetify.VRow(v_for="(latticeElement, index) in selectedLatticeList", align="center", no_gutters=True, style="min-width: 1500px;"):
                    with vuetify.VCol(cols="auto", classes="pa-2"):
                        vuetify.VIcon(
                            "mdi-delete",
                            click=(ctrl.deleteLatticeElement,"[index]"),
                        )
                        vuetify.VChip(
                            v_text=("latticeElement.name",),
                            dense=True,
                            classes="mr-2",
                            style="justify-content: center"
                        )
                    with vuetify.VCol(v_for="(parameter, parameterIndex) in latticeElement.parameters", cols="auto", classes="pa-2"):
                        vuetify.VTextField(
                            label=("parameter.parameter_name",),
                            v_model=("parameter.parameter_default_value",),
                            change=(ctrl.updateLatticeElementParameters, "[index, parameter.parameter_name, $event, parameter.parameter_type]"),
                            error_messages=("parameter.parameter_error_message",),
                            dense=True,
                            style="width: 100px;"
                        )

    def dialog_lattice_settings():
        with vuetify.VCard():
            with vuetify.VCardTitle("Settings", classes="text-subtitle-2 pa-3"):
                vuetify.VSpacer()
            vuetify.VDivider()
            with vuetify.VContainer(fluid=True):
                with vuetify.VRow(no_gutters=True, align="center"):
                    with vuetify.VCol(no_gutters=True, cols="auto"):
                        vuetify.VListItem(
                            "nslice",
                            classes="ma-0 pl-0 font-weight-bold"
                        )
                    with vuetify.VCol(no_gutters=True):
                        vuetify.VTextField(
                            v_model=("nsliceDefaultValue",),
                            change=(ctrl.nsliceDefaultChange, "['nslice', $event]"),
                            placeholder="Value",
                            dense=True,
                            outlined=True,
                            hide_details=True,
                            style="max-width: 75px",
                            classes="ma-0 pa-0"
                        )