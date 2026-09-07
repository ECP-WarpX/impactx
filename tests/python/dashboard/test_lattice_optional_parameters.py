"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

import time

from .utils import TIMEOUT, DashboardTester

# 'Quad' parameters whose signature default is None, and are therefore optional.
# 'name' is the optional parameter that every lattice element carries.
QUAD_OPTIONAL_PARAMETERS = ("name",)

# 'Quad' parameters without a signature default, and are therefore required
QUAD_REQUIRED_PARAMETERS = ("ds", "k")

# 'Quad' parameters carrying a signature default other than None. Having a
# default does not make a parameter optional: these are pre-filled, not blank.
QUAD_DEFAULTED_PARAMETERS = {
    "dx": "0.0",
    "dy": "0.0",
    "rotation": "0.0",
    "aperture_x": "0.0",
    "aperture_y": "0.0",
    "nslice": "1",
}


def add_lattice_element(dashboard: DashboardTester, element_name: str) -> None:
    """
    Adds a lattice element to a freshly reset lattice configuration.

    The reset of the autouse fixture propagates asynchronously, so wait for the
    empty lattice before adding: the helper counts the elements it expects.

    :param dashboard: The dashboard test helper.
    :param element_name: Name of the lattice element to add.
    """
    dashboard.assert_state_list_length("selected_lattice_list", 0)
    dashboard.add_lattice_element(element_name)


def lattice_parameters(dashboard: DashboardTester, index: int) -> dict:
    """
    Returns the parameters of a lattice element, keyed by parameter name.

    :param dashboard: The dashboard test helper.
    :param index: Index of the element in the lattice configuration.
    """
    lattice_list = dashboard.get_state("selected_lattice_list")
    return {
        parameter["parameter_name"]: parameter
        for parameter in lattice_list[index]["parameters"]
    }


def assert_parameter_error(
    dashboard: DashboardTester, index: int, parameter_name: str, has_error: bool
) -> None:
    """
    Asserts a lattice parameter does or does not carry an error, with retry logic.

    :param dashboard: The dashboard test helper.
    :param index: Index of the element in the lattice configuration.
    :param parameter_name: The parameter to check.
    :param has_error: Whether the parameter is expected to carry an error.
    """
    error_message = None
    for _ in range(TIMEOUT):
        error_message = lattice_parameters(dashboard, index)[parameter_name][
            "parameter_error_message"
        ]
        if bool(error_message) == has_error:
            return
        time.sleep(1)

    raise AssertionError(
        f"'{parameter_name}' error message never became "
        f"{'non-empty' if has_error else 'empty'} (got: {error_message})"
    )


def test_optional_parameters_are_blank_and_valid(dashboard) -> None:
    """
    An optional lattice element parameter starts out blank and valid, while a
    required one is flagged until the user provides a value.

    'Quad' shows all three kinds in a single element: 'ds' and 'k' are
    required, 'name' is optional, and the remaining parameters are pre-filled
    from their signature defaults.
    """
    add_lattice_element(dashboard, "Quad")
    parameters = lattice_parameters(dashboard, 0)

    for name in QUAD_OPTIONAL_PARAMETERS:
        parameter = parameters[name]
        assert parameter["parameter_is_optional"] is True, name
        assert parameter["ui_input"] == "", name
        assert parameter["parameter_error_message"] == [], name

    for name in QUAD_REQUIRED_PARAMETERS:
        parameter = parameters[name]
        assert parameter["parameter_is_optional"] is False, name
        assert parameter["parameter_error_message"] != [], name

    for name, default_value in QUAD_DEFAULTED_PARAMETERS.items():
        parameter = parameters[name]
        assert parameter["parameter_is_optional"] is False, name
        assert parameter["ui_input"] == default_value, name
        assert parameter["parameter_error_message"] == [], name


def test_optional_parameters_are_validated_once_filled(dashboard) -> None:
    """
    An optional lattice element parameter is validated as soon as it carries a
    value, and clearing it makes the error go away again.
    """
    add_lattice_element(dashboard, "Quad")

    # not a valid Python identifier
    dashboard.set_input("name1", "2bad")
    assert_parameter_error(dashboard, 0, "name", has_error=True)

    # blank is valid for an optional parameter
    dashboard.set_input("name1", "")
    assert_parameter_error(dashboard, 0, "name", has_error=False)

    # and so is a valid value
    dashboard.set_input("name1", "q1")
    assert_parameter_error(dashboard, 0, "name", has_error=False)


def test_beam_monitor_name_still_defaults(dashboard) -> None:
    """
    'BeamMonitor' takes a required 'name', which keeps its dashboard default.
    """
    add_lattice_element(dashboard, "BeamMonitor")
    parameters = lattice_parameters(dashboard, 0)

    assert parameters["name"]["parameter_is_optional"] is False
    assert parameters["name"]["ui_input"] == "DefaultName"
    assert parameters["name"]["parameter_error_message"] == []


def test_blank_optional_name_tracks(dashboard) -> None:
    """
    A lattice whose element names are all left blank does not disable the run
    button and tracks successfully: the blank names are omitted from the
    element constructors instead of being passed on as an empty string.
    """
    BEAM_PROPERTIES = {
        "tracking_mode": "Particle Tracking",
        "space_charge": "false",
        "csr": False,
        "isr": False,
        "charge_qe": -1,
        "mass_MeV": 0.510998950,
        "npart": 1000,
        "bunch_charge_C": 1e-9,
    }

    DISTRIBUTION_PARAMETERS = {
        "distribution": "Waterbag",
        "distribution_type": "Quadratic",
        "lambdaX": 3.9984884770e-5,
        "lambdaY": 3.9984884770e-5,
        "lambdaT": 1.0e-3,
        "lambdaPx": 2.6623538760e-5,
        "lambdaPy": 2.6623538760e-5,
        "lambdaPt": 2.0e-3,
        "muxpx": -0.846574929020762,
        "muypy": 0.846574929020762,
        "mutpt": 0.0,
    }

    for param_id, value in {**BEAM_PROPERTIES, **DISTRIBUTION_PARAMETERS}.items():
        dashboard.set_input(param_id, value)

    add_lattice_element(dashboard, "Drift")
    dashboard.add_lattice_element("Quad")
    dashboard.assert_state("total_elements", 2)

    # every 'name' is deliberately left blank
    LATTICE_PARAMS = {"ds1": 0.25, "ds2": 1.0, "k2": 1.0}
    for param_id, value in LATTICE_PARAMS.items():
        dashboard.set_input(param_id, value)

    dashboard.assert_state("disableRunSimulationButton", False)

    dashboard.sb.click("#Run_route")
    dashboard.sb.click("#run_simulation_button")
    dashboard.assert_state("sim_progress_status", "Complete!")
