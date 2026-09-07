#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import importlib.util

import numpy as np
import pytest

import amrex.space3d as amr
from impactx import Config, ImpactX, distribution, elements

# Pre-initialize MPI (if ImpactX was built with MPI support) so that MPI stays
# owned by mpi4py and is not finalized between evaluations. Each call to run()
# below constructs and finalizes an ImpactX simulation (amrex Initialize /
# Finalize); without an external MPI owner, the first sim.finalize() would
# finalize MPI and break all following evaluations. This is safe because the
# exploration below uses the "threads" libEnsemble backend (no process forking).
if Config.have_mpi:
    from mpi4py import MPI  # noqa

# configure the test
verbose = True
n_init = 4  # number of initial, random samples before Bayesian optimization starts
max_evals = 60  # total number of ImpactX simulations to evaluate


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
        elements.Drift(name="drift1", ds=2.7, nslice=ns),
        elements.Quad(name="quad1", ds=0.1, k=q1_k, nslice=ns),
        elements.Drift(name="drift2", ds=1.4, nslice=ns),
        elements.Quad(name="quad2", ds=0.2, k=q2_k, nslice=ns),
        elements.Drift(name="drift3", ds=1.4, nslice=ns),
        elements.Quad(name="quad1", ds=0.1, k=q1_k, nslice=ns),
        elements.Drift(name="drift4", ds=2.7, nslice=ns),
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
    sim.tiny_profiler = False

    # set numerical parameters and IO control
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
    ref = sim.beam.ref
    ref.set_species("positron").set_kin_energy_MeV(kin_energy_MeV)

    #   particle bunch
    distr = distribution.Waterbag(
        lambdaX=2.0e-4,
        lambdaY=2.0e-4,
        lambdaT=3.1622776602e-5,
        lambdaPx=1.1180339887e-5,
        lambdaPy=1.1180339887e-5,
        lambdaPt=3.1622776602e-5,
        muxpx=0.894427190999916,
        muypy=-0.894427190999916,
        mutpt=0.0,
    )
    sim.add_particles(bunch_charge_C, distr, npart)

    # design the accelerator lattice
    sim.lattice.extend(build_lattice(parameters, write_particles=write_particles))

    # run simulation
    sim.track_particles()

    # in situ calculate the reduced beam characteristics
    beam = sim.beam
    rbc = beam.beam_moments()

    # clean shutdown
    sim.finalize()

    return rbc


def evaluate(input_params: dict, output_params) -> None:
    """
    A single evaluation for an ``optimas`` ``FunctionEvaluator``.

    optimas calls this function once per trial. The values of the varying
    parameters are provided in ``input_params`` and the results (objective
    and observables) have to be written back into ``output_params`` in place.

    Parameters
    ----------
    input_params: dict
      values of the varying parameters, here the quadrupole strengths k of
      quad 1/3 and quad 2.

    output_params:
      a mapping to fill with the value of the objective ``f`` (the L2 norm of
      alpha and beta of the beam at the end of the simulation) and the
      observables (alpha & beta).
    """
    parameters = {
        "q1_k": input_params["q1_k"],
        "q2_k": input_params["q2_k"],
    }
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

    # objective to minimize
    output_params["f"] = error
    # additional observables to record in the exploration history
    output_params["alpha_x"] = alpha_x
    output_params["alpha_y"] = alpha_y
    output_params["beta_x"] = beta_x
    output_params["beta_y"] = beta_y


@pytest.mark.skipif(
    importlib.util.find_spec("optimas") is None, reason="optimas is not available"
)
def test_optimas():
    from gest_api.vocs import VOCS
    from optimas.diagnostics import ExplorationDiagnostics
    from optimas.evaluators import FunctionEvaluator
    from optimas.explorations import Exploration
    from optimas.generators import AxSingleFidelityGenerator

    # Define the varying parameters (with their bounds), the objective to
    # minimize and the additional observables to record.
    vocs = VOCS(
        variables={
            "q1_k": [-6.0, 0.0],
            "q2_k": [0.0, 6.0],
        },
        objectives={"f": "MINIMIZE"},
        observables=["alpha_x", "alpha_y", "beta_x", "beta_y"],
    )

    # Create the generator: single-fidelity Bayesian optimization using Ax.
    gen = AxSingleFidelityGenerator(vocs=vocs, n_init=n_init)

    # Create the evaluator: run each trial in-process via our ImpactX function.
    ev = FunctionEvaluator(function=evaluate)

    # Create the exploration.
    #
    # We use the "threads" libEnsemble backend: it runs the evaluations and the
    # Ax generator in threads of a single process instead of forking worker
    # processes. Forking would deadlock the PyTorch/BoTorch-based Ax generator
    # and is incompatible with the repeated (re)initialization of ImpactX/AMReX
    # and MPI done here. The "threads" backend supports FunctionEvaluators.
    exp = Exploration(
        generator=gen,
        evaluator=ev,
        max_evals=max_evals,
        sim_workers=1,
        libe_comms="threads",
        exploration_dir_path="./optimize_triplet",
    )

    # run the optimization
    exp.run()

    # Analyze the exploration history and select the best result.
    diags = ExplorationDiagnostics(exp)
    if verbose:
        print(diags.history)

    best = diags.get_best_evaluation()
    best_ks = {
        "q1_k": best["q1_k"].iloc[0],
        "q2_k": best["q2_k"].iloc[0],
    }

    # Print the optimization result
    print("Optimal parameters for k:", best_ks)
    print("L2 norm of alpha & beta at the optimum:", best["f"].iloc[0])

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
    test_optimas()
