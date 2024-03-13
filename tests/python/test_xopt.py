#!/usr/bin/env python3
#
# Copyright 2022-2023 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import importlib

import amrex.space3d as amr
import impactx
import numpy as np
import pytest
from impactx import ImpactX, distribution, elements

# configure the test
verbose = True
gen_name = "Nelder-Mead"  # TuRBO or Nelder-Mead
max_steps = 60


def build_lattice(parameters: dict, write_particles: bool) -> list:
    """
    Create the quadrupole triplet.

    Parameters
    ----------
    parameters: dict
      quadrupole strengths k of quad 1/3 and quad 2.

    write_particles: bool
      write the particles in a beam monitor at the beginning and
      end of the simulation

    Returns
    -------
    A lattice for ImpactX: a list of impactx.elements.
    """
    q1_k, q2_k = parameters["q1_k"], parameters["q2_k"]

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

    if write_particles:
        monitor = elements.BeamMonitor("monitor", backend="h5")
        line = [monitor] + line + [monitor]

    return line


def run(parameters: dict, write_particles=False, write_reduced=False) -> dict:
    """
    Run an ImpactX simulation with a new set of lattice parameters.

    Parameters
    ----------
    parameters: dict
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
    kin_energy_MeV = 2.0e3  # reference energy
    bunch_charge_C = 100.0e-12  # used with space charge
    npart = 10000  # number of macro particles

    #   reference particle
    ref = sim.particle_container().ref_particle()
    ref.set_charge_qe(1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

    #   particle bunch
    distr = distribution.Waterbag(
        sigmaX=2.0e-4,
        sigmaY=2.0e-4,
        sigmaT=3.1622776602e-5,
        sigmaPx=1.1180339887e-5,
        sigmaPy=1.1180339887e-5,
        sigmaPt=3.1622776602e-5,
        muxpx=0.894427190999916,
        muypy=-0.894427190999916,
        mutpt=0.0,
    )
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


def objective(parameters: dict) -> dict:
    """
    A function that is evaluated by the optimizer.

    Parameters
    ----------
    parameters: dict
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
        print(
            f"  -> alpha_x={alpha_x}, alpha_y={alpha_y}, beta_x={beta_x}, beta_y={beta_y}"
        )
    alpha_beta_is = np.array([alpha_x, alpha_y, beta_x, beta_y])

    beta_x_goal = 0.55
    beta_y_goal = beta_x_goal
    alpha_beta_goal = np.array([0, 0, beta_x_goal, beta_y_goal])

    error = np.sum((alpha_beta_is - alpha_beta_goal) ** 2)

    # xopt will ignore NaN results, no cleaning needed
    rbc["error"] = error

    return rbc


@pytest.mark.skipif(
    importlib.util.find_spec("xopt") is None, reason="xopt is not available"
)
def test_xopt():
    from xopt import Xopt
    from xopt.evaluator import Evaluator
    from xopt.vocs import VOCS

    if gen_name == "TuRBO":
        from xopt.generators.bayesian import UpperConfidenceBoundGenerator as Generator

        gen_args = {"turbo_controller": "optimize"}
    elif gen_name == "Nelder-Mead":
        from xopt.generators.scipy.neldermead import NelderMeadGenerator as Generator

        gen_args = {}
    else:
        raise RuntimeError(f"Unconfigured generator named '{gen_name}'")

    # Bounds for values to test: (min, max)
    positive = [0, 10]
    negative = [-10, 0]

    # define variables and function objectives
    vocs = VOCS(
        variables={
            "q1_k": negative,
            "q2_k": positive,
        },
        objectives={"error": "MINIMIZE"},
    )

    # create Xopt evaluator, generator, and Xopt objects
    evaluator = Evaluator(function=objective)
    generator = Generator(vocs=vocs, **gen_args)
    X = Xopt(evaluator=evaluator, generator=generator, vocs=vocs)

    # Initial guess for the quadrople strengths
    initial_quad_strengths = {
        "q1_k": np.array([-3]),
        "q2_k": np.array([3]),
    }
    if gen_name == "TuRBO":
        # a few random guesses
        # X.random_evaluate(3)
        # a few somewhat educated guesses
        X.evaluate_data(initial_quad_strengths)
    elif gen_name == "Nelder-Mead":
        # a few somewhat educated guesses
        X.generator.initial_point = initial_quad_strengths

    # run optimization for 60 steps (iterations)
    for i in range(max_steps):
        X.step()

    # Print all trials
    if verbose:
        print(X.data)

        # plot
        # X.vocs.normalize_inputs(X.data).plot(*X.vocs.variable_names, kind="scatter")

    # Select the best result
    best_idx, best_error = X.vocs.select_best(X.data)
    best_run = X.data.iloc[best_idx]
    best_ks = best_run[["q1_k", "q2_k"]].to_dict(orient="index")[best_idx[0]]

    # Print the optimization result
    print("Optimal parameters for k:", best_ks)
    print("L2 norm of alpha & beta at the optimum:", best_run["error"].values[0])

    # analytical result:
    #   k: -3.5, 2.75
    #   alpha & beta: 0, 0, 0.55, 0.55

    # final run w/ detailed I/O on
    rbc = run(best_ks, write_particles=True, write_reduced=True)
    alpha_x, alpha_y, beta_x, beta_y = (
        rbc["alpha_x"],
        rbc["alpha_y"],
        rbc["beta_x"],
        rbc["beta_y"],
    )
    print(f"alpha_x={alpha_x} alpha_y={alpha_y}\n beta_x={beta_x}     beta_y={beta_y}")


if __name__ == "__main__":
    # Call MPI_Init and MPI_Finalize only once:
    if impactx.Config.have_mpi:
        from mpi4py import MPI  # noqa

    test_xopt()
