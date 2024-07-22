from trame.app import get_server

# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------

server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller

# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------

conversion_factors = {
    "meV": 1.0e-9,
    "eV":  1.0e-6,
    "keV": 1.0e-3,
    "MeV": 1.0,
    "GeV": 1.0e3,
    "TeV": 1.0e6,
}

class inputFunctions:
    def value_of_kin_energy_MeV(kineticEnergyOnDisplayValue, OldUnit):
        state.kin_energy_MeV = kineticEnergyOnDisplayValue * conversion_factors[OldUnit] / conversion_factors["MeV"]
        return state.kin_energy_MeV
    
    def update_kin_energy_on_display(old_unit, new_unit, kin_energy_value):
        value_in_mev = inputFunctions.value_of_kin_energy_MeV(kin_energy_value, old_unit)
        kin_energy_value_on_display = value_in_mev * conversion_factors["MeV"] / conversion_factors[new_unit]
    
        return kin_energy_value_on_display
