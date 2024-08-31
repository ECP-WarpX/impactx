from trame.widgets import vuetify

from ...trame_setup import setup_server

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Default
# -----------------------------------------------------------------------------

state.dynamic_size = False
state.max_level = 3
state.n_cell = [0.0, 0.0, 0.0]
state.prob_relative = []
state.particle_shape = 2
state.poisson_solver = "multigrid"

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
            "label": f"prob_relative{i+1}",
            "value": "",
            "type": "number",
            "error_message": "",
        }
        for i in range(num_prob_relative_fields)
    ]
    state.prob_relative = [0.0] * num_prob_relative_fields
    print(f"Reset prob_relative: {state.prob_relative}")

# -----------------------------------------------------------------------------
# Decorators
# -----------------------------------------------------------------------------

@state.change("space_charge")
def on_space_charge_change(space_charge, **kwargs):
    state.dynamic_size = space_charge

@state.change("max_level")
def on_max_level_change(max_level, **kwargs):
    populate_prob_relative_fields(max_level)

@state.change("n_cell_x", "n_cell_y", "n_cell_z")
def on_nCell_value_change(n_cell_x, n_cell_y, n_cell_z, **kwargs):
    # modify into an array of ints
    state.n_cell = [
        int(n_cell_x) if n_cell_x else 0.0,
        int(n_cell_y) if n_cell_y else 0.0,
        int(n_cell_z) if n_cell_z else 0.0,
    ]

@ctrl.add("update_prob_relative")
def on_update_prob_relative_call(index, value):
    index = int(index)
    if index < len(state.prob_relative):
        state.prob_relative[index] = float(value) if value else 0.0
        state.prob_relative_fields[index]["value"] = str(state.prob_relative[index])
        print(f"Updated prob_relative: {state.prob_relative}")


# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------


class SpaceChargeConfiguration:
    @staticmethod
    def card():
        """
        Creates UI content for space charge configuration
        """

        with vuetify.VCard(v_show="space_charge", style="width: 340px;"):
            vuetify.VCardTitle("Space Charge Configuration")
            vuetify.VDivider()
            with vuetify.VCardText():
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(cols=6, classes="py-0"):
                        vuetify.VCombobox(
                            label="Particle Shape",
                            v_model=("particle_shape",),
                            items=([1, 2, 3],),
                            dense=True,
                        )
                    with vuetify.VCol(cols=6, classes="py-0"):
                        vuetify.VCombobox(
                            label="Poisson Solver",
                            v_model=("poisson_solver",),
                            items=(["multigrid", "fft"],),
                            dense=True,
                            hide_details=True,
                        )
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(cols=4, classes="py-0"):
                        vuetify.VTextField(
                            label="nCell_X",
                            v_model=("n_cell_x",),
                            type="number",
                            dense=True,
                        )
                    with vuetify.VCol(cols=4, classes="py-0"):
                        vuetify.VTextField(
                            label="nCell_Y",
                            v_model=("n_cell_y",),
                            type="number",
                            dense=True,
                        )
                    with vuetify.VCol(cols=4, classes="py-0"):
                        vuetify.VTextField(
                            label="nCell_Z",
                            v_model=("n_cell_z",),
                            type="number",
                            dense=True,
                        )
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(cols=4, classes="py-0"):
                        vuetify.VSelect(
                            label="Max Level",
                            v_model=("max_level",),
                            items=([0, 1, 2, 3, 4],),
                            dense=True,
                        )
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(
                        v_for=("(field, index) in prob_relative_fields",), cols="auto", classes="py-0"
                    ):
                        vuetify.VTextField(
                            label=("field.label",),
                            v_model=("field.value",),
                            input=(ctrl.update_prob_relative, "[index, $event]"),
                            dense=True,
                            style="width: 125px;",
                            type="number",
                        )
