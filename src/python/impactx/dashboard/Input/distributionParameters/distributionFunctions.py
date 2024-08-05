class DistributionFunctions:
    """
    Helper functions for
    User-Input section for beam distribution.
    """

    @staticmethod
    def class_parameters_with_defaults_twiss():
        """
        Reads twissParameters.txt to obtain class parameters with default.
        Temorary function until twiss parameters are in python bindings (8/4/24).
        """

        file_path = "/mnt/c/Users/parth/Downloads/vsCode/TrameApplication/Input/distributionParametersCard/twissParameters.txt"

        with open(file_path, "r") as file:
            content = file.read()
            local_vars = {}
            exec(content, {}, local_vars)
            return local_vars["twiss_parameters"]
