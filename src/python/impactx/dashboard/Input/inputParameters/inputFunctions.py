"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from ...trame_setup import setup_server

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------

CONVERSION_FACTORS = {
    "meV": 1.0e-9,
    "eV": 1.0e-6,
    "keV": 1.0e-3,
    "MeV": 1.0,
    "GeV": 1.0e3,
    "TeV": 1.0e6,
}


class InputFunctions:
    """
    Helper functions for the
    User-Input section for beam properties.
    """

    @staticmethod
    def value_of_kin_energy_MeV(kineticEnergyOnDisplayValue, OldUnit):
        """
        Converts the kinetic energy to MeV.
        :param kineticEnergyOnDisplayValue: The kinetic energy value that is displayed in the UI.
        :param OldUnit: The previous unit of kin_energy prior to change.
        :return: The kinetic energy value converted to MeV.
        """

        state.kin_energy_MeV = (
            kineticEnergyOnDisplayValue
            * CONVERSION_FACTORS[OldUnit]
            / CONVERSION_FACTORS["MeV"]
        )
        return state.kin_energy_MeV

    @staticmethod
    def update_kin_energy_on_display(old_unit, new_unit, kin_energy_value):
        """
        Updates the kinetic energy value in the UI.
        :param old_unit: The previous unit of kin_energy prior to change.
        :param new_unit: The new unit of kin_energy changed to by user.
        :param kin_energy_value: The kinetic energy value in the current unit.
        :return: The kinetic energy value converted to the new unit.
        """

        value_in_mev = InputFunctions.value_of_kin_energy_MeV(
            kin_energy_value, old_unit
        )
        kin_energy_value_on_display = (
            value_in_mev * CONVERSION_FACTORS["MeV"] / CONVERSION_FACTORS[new_unit]
        )

        return kin_energy_value_on_display
