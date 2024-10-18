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
state.n_cell = []
state.prob_relative = []
state.particle_shape = 2
state.poisson_solver = "fft"

state.prob_relative_fields = []
state.n_cell_x = 32
state.n_cell_y = 32
state.n_cell_z = 32

state.blocking_factor_x = 16
state.blocking_factor_y = 16
state.blocking_factor_z = 16

state.mlmg_relative_tolerance = 1.0e-7
state.mlmg_absolute_tolerance = 0.0
state.mlmg_max_iters = 100
state.mlmg_verbosity = 1

state.error_message_mlmg_relative_tolerance = ""
state.error_message_mlmg_absolute_tolerance = ""
state.error_message_mlmg_max_iters = ""
state.error_message_mlmg_verbosity = ""

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------


def populate_prob_relative_fields(max_level):
    num_prob_relative_fields = int(max_level) + 1

    if state.poisson_solver == "fft":
        state.prob_relative = [1.1] + [0.0] * (num_prob_relative_fields - 1)
    elif state.poisson_solver == "multigrid":
        state.prob_relative = [3.1] + [0.0] * (num_prob_relative_fields - 1)
    else:
        state.prob_relative = [0.0] * num_prob_relative_fields

    state.prob_relative_fields = [
        {
            "value": state.prob_relative[i],
            "error_message": SpaceChargeFunctions.validate_prob_relative_fields(
                i, state.prob_relative[i]
            ),
        }
        for i in range(num_prob_relative_fields)
    ]


def update_state_values_and_errors(category, kwargs):
    directions = ["x", "y", "z"]

    for state_name, value in kwargs.items():
        if any(state_name == f"{category}_{dir}" for dir in directions):
            direction = state_name.split("_")[-1]
            error_message_name = f"error_message_{category}_{direction}"

            # update error message for blocking_factor_{direction}
            updated_blocking_factor_error_message = generalFunctions.validate_against(
                value, "int"
            )
            updated_n_cell_error_message = SpaceChargeFunctions.validate_n_cell_field(
                direction
            )

            # update error message for n_cell_{direction}
            setattr(state, error_message_name, updated_blocking_factor_error_message)
            setattr(
                state, f"error_message_n_cell_{direction}", updated_n_cell_error_message
            )

            if not updated_blocking_factor_error_message:
                setattr(state, state_name, int(value))  # convert to int if not an error

    updated_attribute = [
        getattr(state, f"{category}_{direction}", 0) for direction in directions
    ]

    if category == "blocking_factor":
        state.blocking_factor = updated_attribute
    elif category == "n_cell":
        state.n_cell = updated_attribute


# -----------------------------------------------------------------------------
# Decorators
# -----------------------------------------------------------------------------
@state.change("poisson_solver")
def on_poisson_solver_change(poisson_solver, **kwargs):
    populate_prob_relative_fields(state.max_level)
    state.dirty("prob_relative_fields")


@state.change("space_charge")
def on_space_charge_change(space_charge, **kwargs):
    state.dynamic_size = space_charge


@state.change("max_level")
def on_max_level_change(max_level, **kwargs):
    populate_prob_relative_fields(max_level)


@state.change("blocking_factor_x", "blocking_factor_y", "blocking_factor_z")
def on_blocking_factor_change(**kwargs):
    update_state_values_and_errors("blocking_factor", kwargs)


@state.change("n_cell_x", "n_cell_y", "n_cell_z")
def on_n_cell_change(**kwargs):
    update_state_values_and_errors("n_cell", kwargs)


@ctrl.add("update_prob_relative")
def on_update_prob_relative_call(index, value):
    index = int(index)

    try:
        prob_relative_value = float(value)
        state.prob_relative[index] = prob_relative_value
    except ValueError:
        prob_relative_value = 0.0
        state.prob_relative[index] = prob_relative_value

    # Validate the updated value
    error_message = SpaceChargeFunctions.validate_prob_relative_fields(
        index, prob_relative_value
    )

    state.prob_relative_fields[index]["value"] = value
    state.prob_relative_fields[index]["error_message"] = error_message

    # Validate next index if it exists
    if index + 1 < len(state.prob_relative):
        next_value = state.prob_relative[index + 1]

        next_error_message = SpaceChargeFunctions.validate_prob_relative_fields(
            index + 1, next_value
        )
        state.prob_relative_fields[index + 1]["error_message"] = next_error_message

    state.dirty("prob_relative_fields")


# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------


class SpaceChargeConfiguration:
    @staticmethod
    def card():
        """
        Creates UI content for space charge configuration
        """

        with vuetify.VDialog(v_model=("showSpaceChargeDialog", False), width="500px"):
            SpaceChargeConfiguration.dialog_space_charge_settings()

        with vuetify.VCard(v_show="space_charge", style="width: 340px;"):
            with vuetify.VCardTitle("Space Charge"):
                vuetify.VSpacer()
                vuetify.VIcon(
                    "mdi-cog",
                    classes="ml-2",
                    v_if="poisson_solver == 'multigrid'",
                    click="showSpaceChargeDialog = true",
                    style="cursor: pointer;",
                )
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
                    for direction in ["x", "y", "z"]:
                        with vuetify.VCol(cols=4, classes="py-0"):
                            vuetify.VTextField(
                                placeholder=direction,
                                v_model=(f"n_cell_{direction}",),
                                error_messages=(f"error_message_n_cell_{direction}",),
                                type="number",
                                dense=True,
                                style="margin-top: -5px",
                            )
                with vuetify.VCol(classes="pa-0"):
                    vuetify.VListItemSubtitle(
                        "Blocking Factor",
                        classes="font-weight-bold black--text mt-1",
                    )
                with vuetify.VRow(classes="my-0"):
                    for direction in ["x", "y", "z"]:
                        with vuetify.VCol(cols=4, classes="py-0"):
                            vuetify.VTextField(
                                placeholder=direction,
                                v_model=(f"blocking_factor_{direction}",),
                                error_messages=(
                                    f"error_message_blocking_factor_{direction}",
                                ),
                                type="number",
                                dense=True,
                                style="margin-top: -5px",
                            )
                with vuetify.VCol(classes="pa-0"):
                    vuetify.VListItemSubtitle(
                        "prob_relative",
                        classes="font-weight-bold black--text mt-1",
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

    @staticmethod
    def dialog_space_charge_settings():
        """
        Creates UI content for space charge configuration
        settings.
        """
        with vuetify.VCard():
            with vuetify.VTabs(
                v_model=("space_charge_tab", "Advanced Multigrid Settings")
            ):
                vuetify.VTab("Settings")
            vuetify.VDivider()
            with vuetify.VTabsItems(v_model="space_charge_tab"):
                with vuetify.VTabItem():
                    with vuetify.VContainer(fluid=True):
                        with vuetify.VRow(
                            classes="my-2", v_if="poisson_solver == 'multigrid'"
                        ):
                            with vuetify.VCol(cols=6, classes="py-0"):
                                vuetify.VTextField(
                                    label="MLMG Relative Tolerance",
                                    v_model=("mlmg_relative_tolerance",),
                                    error_messages=(
                                        "error_message_mlmg_relative_tolerance",
                                    ),
                                    type="number",
                                    dense=True,
                                )
                            with vuetify.VCol(cols=6, classes="py-0"):
                                vuetify.VTextField(
                                    label="MLMG Absolute Tolerance",
                                    v_model=("mlmg_absolute_tolerance",),
                                    error_messages=(
                                        "error_message_mlmg_absolute_tolerance",
                                    ),
                                    type="number",
                                    dense=True,
                                )
                        with vuetify.VRow(
                            classes="my-0", v_if="poisson_solver == 'multigrid'"
                        ):
                            with vuetify.VCol(cols=6, classes="py-0"):
                                vuetify.VTextField(
                                    label="MLMG Max Iterations",
                                    v_model=("mlmg_max_iters",),
                                    error_messages=("error_message_mlmg_max_iters",),
                                    type="number",
                                    dense=True,
                                )
                            with vuetify.VCol(cols=6, classes="py-0"):
                                vuetify.VTextField(
                                    label="MLMG Verbosity",
                                    v_model=("mlmg_verbosity",),
                                    error_messages=("error_message_mlmg_verbosity",),
                                    type="number",
                                    dense=True,
                                )
