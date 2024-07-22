
from trame.app import get_server
import re

# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------

server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller

# -----------------------------------------------------------------------------
# Content
# -----------------------------------------------------------------------------

import re

def retrieve_state_content():
    """
    This function builds the output file based on user inputs
    """

    with open('output_distribution_parameters.txt', 'r') as file:
        file_content = file.read()

    
    lambdaX = re.search(r'lambdaX=(\d+)', file_content).group(1)
    lambdaY = re.search(r'lambdaY=(\d+)', file_content).group(1)
    lambdaT = re.search(r'lambdaT=(\d+)', file_content).group(1)
    lambdaPx = re.search(r'lambdaPx=(\d+)', file_content).group(1)
    lambdaPy = re.search(r'lambdaPy=(\d+)', file_content).group(1)
    lambdaPt = re.search(r'lambdaPt=(\d+)', file_content).group(1)
    muxpx = re.search(r'muxpx=(\d+\.\d+)', file_content).group(1)
    muypy = re.search(r'muypy=(\d+\.\d+)', file_content).group(1)
    mutpt = re.search(r'mutpt=(\d+\.\d+)', file_content).group(1)

    
    content = f"""###############################################################################
# Particle Beam(s)
###############################################################################
beam.npart = {state.npart}
beam.units = static
beam.kin_energy = {state.kin_energy_MeV}
beam.charge = {state.charge}
beam.particle = electron
beam.distribution = {state.selectedDistribution}
beam.lambdaX = {lambdaX} 
beam.lambdaY = {lambdaY} 
beam.lambdaT = {lambdaT} 
beam.lambdaPx = {lambdaPx}
beam.lambdaPy = {lambdaPy}
beam.lambdaPt = {lambdaPt}
beam.muxpx = {muxpx}
beam.muypy = {muypy}
beam.mutpt = {mutpt}


###############################################################################
# Beamline: lattice elements and segments
###############################################################################

###############################################################################
# Algorithms
###############################################################################
algo.particle_shape = {state.particle_shape}
algo.space_charge = {state.space_charge}


###############################################################################
# Diagnostics
###############################################################################
diag.slice_step_diagnostics = {state.slice_step_diagnostics}

""" 
    return content

