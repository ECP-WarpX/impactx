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

            if poisson_solver == "multigrid":
                if index == 0 and prob_relative_value < 3:
                    error_message = "Must be greater than 3."
                elif index > 0 and prob_relative_value <= 1:
                    error_message = "Must be greater than 1."
            elif poisson_solver == "fft":
                if prob_relative_value <= 1:
                    error_message = "Must be greater than 1."
        except ValueError:
            error_message = "Must be a float."

        return error_message
