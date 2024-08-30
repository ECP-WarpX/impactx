from trame.widgets import vuetify

from ...trame_setup import setup_server
from ..generalFunctions import generalFunctions
from .spaceChargeFunctions import SpaceChargeFunctions

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Default
# -----------------------------------------------------------------------------

state.dynamic_size = False
state.max_level = 0
state.n_cell = [32.0, 32.0, 32.0]
state.prob_relative = []
state.particle_shape = 2
state.poisson_solver = "fft"

state.prob_relative_fields = []
state.n_cell_x = ""
state.n_cell_y = ""
state.n_cell_z = ""

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------


def populate_prob_relative_fields(max_level):
    num_prob_relative_fields = int(max_level) + 1

    state.prob_relative_fields = [
        {
            "value": "",
            "error_message": SpaceChargeFunctions.validate_prob_relative_fields(i, ""),
        }
        for i in range(num_prob_relative_fields)
    ]
    state.prob_relative = [0.0] * num_prob_relative_fields
    print(f"Reset prob_relative: {state.prob_relative}")


# -----------------------------------------------------------------------------
# Decorators
# -----------------------------------------------------------------------------
@state.change("poisson_solver")
def on_poisson_solver_change(poisson_solver, **kwargs):
    updated_prob_relative_fields = []

    for i, field in enumerate(state.prob_relative_fields):
        prob_relative_value = state.prob_relative[i]
        error_message = SpaceChargeFunctions.validate_prob_relative_fields(
            i, prob_relative_value
        )

        updated_field = {"value": field["value"], "error_message": error_message}

        updated_prob_relative_fields.append(updated_field)

    state.prob_relative_fields = updated_prob_relative_fields
    state.dirty("prob_relative_fields")


@state.change("space_charge")
def on_space_charge_change(space_charge, **kwargs):
    state.dynamic_size = space_charge


@state.change("max_level")
def on_max_level_change(max_level, **kwargs):
    populate_prob_relative_fields(max_level)


@state.change("n_cell_x", "n_cell_y", "n_cell_z")
def on_nCell_value_change(n_cell_x, n_cell_y, n_cell_z, **kwargs):
    state.error_message_x = generalFunctions.validate_against(n_cell_x, "int")
    state.error_message_y = generalFunctions.validate_against(n_cell_y, "int")
    state.error_message_z = generalFunctions.validate_against(n_cell_z, "int")

    state.n_cell = [
        int(n_cell_x) if not state.error_message_x else 0,
        int(n_cell_y) if not state.error_message_y else 0,
        int(n_cell_z) if not state.error_message_z else 0,
    ]

    state.dirty("n_cell")


@ctrl.add("update_prob_relative")
def on_update_prob_relative_call(index, value):
    prob_relative_value, input_type = generalFunctions.determine_input_type(value)

    index = int(index)

    # Validate the updated value
    error_message = SpaceChargeFunctions.validate_prob_relative_fields(
        index, prob_relative_value
    )

    if index < len(state.prob_relative):
        state.prob_relative[index] = prob_relative_value if value else 0.0
        print(f"Updated prob_relative: {state.prob_relative}")
        state.prob_relative_fields[index]["error_message"] = error_message
        state.dirty("prob_relative_fields")


# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------

@ctrl.add("updateArray")
def updateArray(index, value):
    index = int(index)
    if index < len(state.prob_relative):
        state.prob_relative[index] = float(value) if value else 0.0
        state.level_fields[index]["value"] = str(state.prob_relative[index])
        print(f"Updated prob_relative: {state.prob_relative}")

class SpaceChargeConfiguration:
    @staticmethod
    def card():
        """
        Creates UI content for space charge configuration
        """

        with vuetify.VCard(v_show="space_charge", style="width: 340px; height: 300px"):
            with vuetify.VCardTitle("Space Charge"):
                vuetify.VSpacer()
                vuetify.VIcon(
                    "mdi-information",
                    classes="ml-2",
                    click=lambda: generalFunctions.documentation(
                        "space_charge_documentation"
                    ),
                    style="color: #00313C;",
                )
            vuetify.VDivider()
            with vuetify.VCardText():
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(cols=5, classes="py-0"):
                        vuetify.VSelect(
                            label="Poisson Solver",
                            v_model=("poisson_solver",),
                            items=(["multigrid", "fft"],),
                            dense=True,
                            hide_details=True,
                        )
                    with vuetify.VCol(cols=4, classes="py-0"):
                        vuetify.VSelect(
                            label="Particle Shape",
                            v_model=("particle_shape",),
                            items=([1, 2, 3],),
                            dense=True,
                        )
                    with vuetify.VCol(cols=3, classes="py-0"):
                        vuetify.VSelect(
                            label="Max Level",
                            v_model=("max_level",),
                            items=([0, 1, 2, 3, 4],),
                            dense=True,
                        )
                with vuetify.VCol(classes="pa-0"):
                    vuetify.VListItemSubtitle(
                        "nCell",
                        classes="font-weight-bold black--text",
                    )
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(cols=3, classes="py-0"):
                        vuetify.VTextField(
                            placeholder="x",
                            v_model=("n_cell_x",),
                            error_messages=("error_message_x",),
                            type="number",
                            dense=True,
                            style="margin-top: -5px",
                        )
                    with vuetify.VCol(cols=3, classes="py-0"):
                        vuetify.VTextField(
                            placeholder="y",
                            v_model=("n_cell_y",),
                            error_messages=("error_message_y",),
                            type="number",
                            dense=True,
                            style="margin-top: -5px",
                        )
                    with vuetify.VCol(cols=3, classes="py-0"):
                        vuetify.VTextField(
                            placeholder="z",
                            v_model=("n_cell_z",),
                            error_messages=("error_message_z",),
                            type="number",
                            dense=True,
                            style="margin-top: -5px",
                        )
                with vuetify.VCol(classes="pa-0"):
                    vuetify.VListItemSubtitle(
                        "prob_relative",
                        classes="font-weight-bold black--text",
                    )
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(
                        v_for=("(field, index) in prob_relative_fields",),
                        classes="py-0",
                    ):
                        vuetify.VTextField(
                            placeholder=("val."),
                            v_model=("field.value",),
                            input=(ctrl.update_prob_relative, "[index, $event]"),
                            error_messages=("field.error_message",),
                            type="number",
                            dense=True,
                            style="margin-top: -5px",
                        )
