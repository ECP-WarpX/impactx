"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from .. import state
from ..Input.distribution.utils import DistributionFunctions
from ..Input.validation import DashboardValidation

TRACKING_MODE_COMMANDS = {
    "Particle Tracking": """\
sim.add_particles(bunch_charge_C, distr, npart)
sim.track_particles()""",
    "Envelope Tracking": """\
sim.init_envelope(ref, distr)
sim.track_envelope()""",
    "Reference Tracking": "sim.track_reference(ref)",
}

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------


def build_distribution_list():
    """
    Generates an instance of distribution inputs
    as a string for exporting purposes.
    """
    distribution_name = state.distribution
    parameters = DistributionFunctions.convert_distribution_parameters_to_valid_type()

    indentation = " " * (8 if state.distribution_type == "Twiss" else 4)
    distribution_parameters = ",\n".join(
        f"{indentation}{key}={value}" for key, value in parameters.items()
    )

    if state.distribution_type == "Twiss":
        return (
            f"distr = distribution.{distribution_name}(\n"
            f"    **twiss(\n"
            f"{distribution_parameters},\n"
            f"    )\n"
            f")"
        )
    else:
        return (
            f"distr = distribution.{distribution_name}(\n{distribution_parameters},\n)"
        )


def build_lattice_list() -> str:
    """
    Constructs the Python export string defining the lattice configuration
    from the dashboard. Assumes all lattice parameters have already been validated.

    :return: A Python-formatted string for `lattice_configuration = [...]`
    """
    lattice_elements = []

    for element in state.selected_lattice_list:
        name = element["name"]
        parameter_strings = []

        for param in element["parameters"]:
            param_name = param["parameter_name"]
            param_value = param["sim_input"]
            param_type = param["parameter_type"]

            # an optional parameter left blank stays at its default: passing an
            # empty string instead would build the wrong element
            if param.get("parameter_is_optional", False):
                if DashboardValidation.is_blank(param_value):
                    continue

            if param_type == "str":
                formatted_value = f'"{param_value}"'
            elif param_type == "bool":
                formatted_value = str(param_value).strip().capitalize()
            else:
                formatted_value = param_value
            parameter_strings.append(f"{param_name}={formatted_value}")

        element_string = f"elements.{name}(" + ", ".join(parameter_strings) + ")"
        lattice_elements.append(element_string)

    result = (
        "lattice_configuration = [\n    " + ",\n    ".join(lattice_elements) + "\n]"
    )
    return result


def build_space_charge_or_csr():
    """
    Generates simulation content for space charge
    and csr.
    """
    content = ""

    if state.space_charge != "false":
        content += f"""# Space Charge
sim.space_charge = "{state.space_charge}"
sim.dynamic_size = {state.dynamic_size}
sim.poisson_solver = '{state.poisson_solver}'
sim.particle_shape = {state.particle_shape}
sim.max_level = {state.max_level}
sim.n_cell = {state.n_cell}
sim.blocking_factor_x = [{state.blocking_factor_x}]
sim.blocking_factor_y = [{state.blocking_factor_y}]
sim.blocking_factor_z = [{state.blocking_factor_z}]
sim.prob_relative = {state.prob_relative}
"""
        if state.poisson_solver == "multigrid":
            content += f"""
# Space Charge - Multigrid-Specific Numerical Options
sim.mlmg_relative_tolerance = {state.mlmg_relative_tolerance}
sim.mlmg_absolute_tolerance = {state.mlmg_absolute_tolerance}
sim.mlmg_max_iters = {state.mlmg_max_iters}
sim.mlmg_verbosity = {state.mlmg_verbosity}
"""
    if state.csr:
        content += f"""# Coherent Synchrotron Radiation
sim.csr = {state.csr}
sim.csr_bins = {state.csr_bins}
"""
        if not state.space_charge:
            content += f"""
sim.particle_shape = {state.particle_shape}
"""
    if not content:
        content = f"""
sim.particle_shape = {state.particle_shape}
"""

    return content


def build_isr():
    """
    Generates simulation content for Incoherent Synchrotron Radiation (ISR).
    """
    if state.isr:
        return f"""# Incoherent Synchrotron Radiation
sim.isr = {state.isr}
sim.isr_order = {state.isr_order}
"""
    return ""


def build_tracking_commands() -> str:
    """
    Read the user's choice from state.tracking_mode and
    return the corresponding ImpactX sim command block.
    """
    return TRACKING_MODE_COMMANDS[state.tracking_mode]


# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------
def generate_phase_space(is_exporting: bool) -> str:
    """
    Returns the plotting section of the script as a string,
    or an empty string if the script is being exported.
    """
    if is_exporting or state.tracking_mode.lower() != "particle tracking":
        return ""
    return (
        "import matplotlib.pyplot as plt\n"
        "fig = pc.plot_phasespace()\n"
        "if fig is not None:\n"
        '    fig.savefig("phase_space_plot.png")\n'
    )


def dashboard_sim_inputs(is_exporting=False) -> str:
    """
    This function creates the template to export
    dashboard user inputs into a python script.
    """
    script = f"""
from impactx import ImpactX, distribution, elements, twiss

sim = ImpactX()

{build_space_charge_or_csr()}
sim.slice_step_diagnostics = True

sim.init_grids()

# Initialize particle beam
kin_energy_MeV = {state.kin_energy_MeV}
bunch_charge_C = {state.bunch_charge_C}
npart = {state.npart}

pc = sim.particle_container()

# Reference particle
ref = pc.ref_particle()
ref.set_charge_qe({state.charge_qe}).set_mass_MeV({state.mass_MeV}).set_kin_energy_MeV(kin_energy_MeV)

{build_distribution_list()}
{build_isr()}

{build_lattice_list()}
sim.lattice.extend(lattice_configuration)
sim.periods = {state.periods}

# Simulate
{build_tracking_commands()}

{generate_phase_space(is_exporting)}
# Clean Shutdown
sim.finalize()
"""

    return script
