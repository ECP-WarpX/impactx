class distributionFunctions:

    def classAndParametersAndDefaultValueAndType_Twiss():
        file_path = "/mnt/c/Users/parth/Downloads/vsCode/TrameApplication/Input/distributionParametersCard/twissParameters.txt"
        
        with open(file_path, "r") as file:
            content = file.read()
            local_vars = {}
            exec(content, {}, local_vars)
            return local_vars['twiss_parameters']

