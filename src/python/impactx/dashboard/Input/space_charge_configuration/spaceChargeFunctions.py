from ...trame_setup import setup_server

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------


class SpaceChargeFunctions:
    @staticmethod
    def validate_prob_relative_fields(index, prob_relative_value):
        """
        This function checks specific validation requirements
        for prob_relative_fields.
        :param index: The index of the prob_relative_field modified.
        :param prob_relative_value: The numerical value entered by the user.
        :return: An error message. An empty string if there is no error.
        """
        error_message = ""

        try:
            prob_relative_value = float(prob_relative_value)
            poisson_solver = state.poisson_solver

            if index == 0:
                if poisson_solver == "multigrid":
                    if prob_relative_value < 3:
                        error_message = "Must be greater than 3."
                elif poisson_solver == "fft":
                    if prob_relative_value <= 1:
                        error_message = "Must be greater than 1."
            else:
                previous_value = float(state.prob_relative[index - 1])
                if prob_relative_value >= previous_value:
                    error_message = (
                        f"Must be less than previous value ({previous_value})."
                    )
                else:
                    if prob_relative_value <= 1:
                        error_message = "Must be greater than 1."
        except ValueError:
            error_message = "Must be a float."

        return error_message

    @staticmethod
    def validate_n_cell_field(direction):
        """
        Validates that n_cell_value is a multiple of blocking_factor_value.
        """

        n_cell_value = int(getattr(state, f"n_cell_{direction}", 0))
        blocking_factor_value = int(getattr(state, f"blocking_factor_{direction}", 0))

        if blocking_factor_value == 0:
            return f"Blocking factor for {direction} cannot be zero."

        if n_cell_value % blocking_factor_value != 0:
            return f"n_cell_value for {direction} is not a multiple of blocking_factor_value."

        return ""
