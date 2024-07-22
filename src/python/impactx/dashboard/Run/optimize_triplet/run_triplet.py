# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------
from trame.app import get_server

server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller

# -----------------------------------------------------------------------------
# Trame Code
# -----------------------------------------------------------------------------
from Analyze.analyzeFunctions import analyzeFunctions # to load in distribution and lattice inputs
### state.npart
### state.particle_shape - yet to be added
### state.space_charge - yet to be added
### state.slice_step_diagnostics - yet to be added
### state.bunch_charge_C
### state.kin_energy_MeV
### distr
### line (latticeElements) - not added because doesnt work if dynamically going to be changed

# -----------------------------------------------------------------------------
# Optimize Triplet Code
# -----------------------------------------------------------------------------

#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import minimize

import amrex.space3d as amr
import impactx
from impactx import ImpactX, distribution, elements

# Call MPI_Init and MPI_Finalize only once:
if impactx.Config.have_mpi:
    from mpi4py import MPI  # noqa

verbose = False


def build_lattice(parameters: tuple, write_particles: bool) -> list:
    """
    Create the quadrupole triplet.

    Parameters
    ----------
    parameters: tuple
      quadrupole strengths k of quad 1/3 and quad 2.

    write_particles: bool
      write the particles in a beam monitor at the beginning and
      end of the simulation

    Returns
    -------
    A lattice for ImpactX: a list of impactx.elements.
    """
    q1_k, q2_k = parameters

    ns = 10  # number of slices per ds in the element

    # enforce a mirror symmetry of the triplet
    line = [
        elements.Drift(ds=2.7, nslice=ns),
        elements.Quad(ds=0.1, k=q1_k, nslice=ns),
        elements.Drift(ds=1.4, nslice=ns),
        elements.Quad(ds=0.2, k=q2_k, nslice=ns),
        elements.Drift(ds=1.4, nslice=ns),
        elements.Quad(ds=0.1, k=q1_k, nslice=ns),
        elements.Drift(ds=2.7, nslice=ns),
    ]
    # line = analyzeFunctions.read_latticeElements_file()
    if write_particles:
        monitor = elements.BeamMonitor("monitor", backend="h5")
        line = [monitor] + line + [monitor]

    return line


def run(parameters: tuple, write_particles=False, write_reduced=False) -> dict:
    """
    Run an ImpactX simulation with a new set of lattice parameters.

    Parameters
    ----------
    parameters: tuple
      quadrupole strengths k of quad 1/3 and quad 2.

    write_particles: bool
      write the particles in a beam monitor at the beginning and
      end of the simulation

    write_reduced: bool
      write the reduced diagnositcs of ImpactX to a file.

    Returns
    -------
    A dictionary with reduced diagnositcs of ImpactX, characterizing
    the beam at the end of the simulation.
    """
    pp_amrex = amr.ParmParse("amrex")
    pp_amrex.add("verbose", 0)

    sim = ImpactX()

    if verbose is False:
        sim.verbose = 0

    # set numerical parameters and IO control
    sim.particle_shape = 2  # B-spline order
    sim.space_charge = False
    sim.diagnostics = write_reduced
    sim.slice_step_diagnostics = write_reduced

    # domain decomposition & space charge mesh
    sim.init_grids()

    # load a 2 GeV electron beam with an initial
    # unnormalized rms emittance of 5 nm
    kin_energy_MeV = state.kin_energy_MeV  # reference energy
    bunch_charge_C = state.bunch_charge_C  # used with space charge
    npart = state.npart  # number of macro particles

    #   reference particle
    ref = sim.particle_container().ref_particle()
    ref.set_charge_qe(1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

    #   particle bunch
    # distr = distribution.Waterbag(
    #     lambdaX=2.0e-4,
    #     lambdaY=2.0e-4,
    #     lambdaT=3.1622776602e-5,
    #     lambdaPx=1.1180339887e-5,
    #     lambdaPy=1.1180339887e-5,
    #     lambdaPt=3.1622776602e-5,
    #     muxpx=0.894427190999916,
        # muypy=-0.894427190999916,
        # mutpt=0.0,
    # )
    distr = analyzeFunctions.read_distribution_file()
    sim.add_particles(bunch_charge_C, distr, npart)

    # design the accelerator lattice
    sim.lattice.extend(build_lattice(parameters, write_particles=write_particles))

    # run simulation
    sim.evolve()

    # in situ calculate the reduced beam characteristics
    beam = sim.particle_container()
    rbc = beam.reduced_beam_characteristics()

    # clean shutdown
    sim.finalize()

    return rbc


def objective(parameters: tuple) -> float:
    """
    A function that is evaluated by the optimizer.

    Parameters
    ----------
    parameters: tuple
      quadrupole strengths k of quad 1/3 and quad 2.

    Returns
    -------
    The L2 norm of alpha and beta of the beam at the end of the
    simulation.
    """
    if verbose:
        print(f"Run objective with parameters={parameters}...")

    rbc = run(parameters, write_particles=False, write_reduced=False)
    alpha_x, alpha_y, beta_x, beta_y = (
        rbc["alpha_x"],
        rbc["alpha_y"],
        rbc["beta_x"],
        rbc["beta_y"],
    )
    if verbose:
        print(f"alpha_x={alpha_x}, alpha_y={alpha_y}, beta_x={beta_x}, beta_y={beta_y}")
    alpha_beta_is = np.array([alpha_x, alpha_y, beta_x, beta_y])

    beta_x_goal = 0.55
    beta_y_goal = beta_x_goal
    alpha_beta_goal = np.array([0, 0, beta_x_goal, beta_y_goal])

    error = np.sum((alpha_beta_is - alpha_beta_goal) ** 2)

    if np.isnan(error):
        error = 1.0e99

    return error


# if __name__ == "__main__":

def run_optimize_triplet():
    # Initial guess for the quadrople strengths
    initial_quad_strengths = np.array([-3, 3])

    # Bounds for values to test: (min, max)
    positive = (0, None)
    negative = (None, 0)
    bounds = [negative, positive]

    # optimizer specific values
    #   https://docs.scipy.org/doc/scipy/reference/optimize.minimize-neldermead.html
    #   https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html
    options = {
        "maxiter": 1000,
    }

    # Call the optimizer
    res = minimize(
        objective,
        initial_quad_strengths,
        method="Nelder-Mead",
        # method="L-BFGS-B",
        tol=1.0e-8,
        options=options,
        bounds=bounds,
    )

    # Print the optimization result
    print("Optimal parameters for k:", res.x)
    print("L2 norm of alpha & beta at the optimum:", res.fun)

    # analytical result:
    #   k: -3.5, 2.75
    #   alpha & beta: 0, 0, 0.55, 0.55

    # final run
    rbc = run(res.x, write_particles=True, write_reduced=True)
    alpha_x, alpha_y, beta_x, beta_y = (
        rbc["alpha_x"],
        rbc["alpha_y"],
        rbc["beta_x"],
        rbc["beta_y"],
    )
    print(f"alpha_x={alpha_x} alpha_y={alpha_y}\n beta_x={beta_x}     beta_y={beta_y}")
