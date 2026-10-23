#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
import pytest

from impactx import Config, ImpactX, distribution, elements

# a carbon-12, fully stripped: in the charge state 6+
carbon_mass_MeV = 12.0 * 931.494  # approximately, the same in every charge state
carbon_kin_energy_MeV = 12.0 * 1.0e3
carbon_charge_qe = 6.0


def _set_carbon(sim, charge_qe):
    """Set the reference particle to a carbon-12 ion in the charge state charge_qe

    :return: the charge/mass of the carbon ion in e / eV
    """
    ref = sim.beam.ref
    ref.set_mass_MeV(carbon_mass_MeV).set_charge_qe(charge_qe)
    ref.set_kin_energy_MeV(carbon_kin_energy_MeV)
    return charge_qe / (carbon_mass_MeV * 1.0e6)


def _kicks(model, charge_qe):
    """Momenta after one collective kick on a beam of carbon ions at rest

    The beam is the same in every charge state: the positions, the number of real
    particles (weights) and the mass and kinetic energy of the reference particle.
    """
    sim = ImpactX()

    sim.particle_shape = 2
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    if model == "CSR":
        sim.space_charge = False
        sim.csr = True
        sim.csr_bins = 64
    else:
        sim.space_charge = model
    if model in ("Gauss2p5D", "2p5D"):
        sim.space_charge_apply_longitudinal_kick = True
    if model == "2p5D":
        sim.poisson_solver = "fft"
        sim.max_level = 0
        sim.n_cell = [64, 64, 1]
        sim.blocking_factor_x = [16]
        sim.blocking_factor_y = [16]
        sim.blocking_factor_z = [1]
        sim.prob_relative = [1.1]
    if model == "3D":
        sim.poisson_solver = "fft" if Config.have_fft else "multigrid"
        sim.max_level = 0
        sim.n_cell = [16, 16, 32]
        sim.blocking_factor_x = [16]
        sim.blocking_factor_y = [16]
        sim.blocking_factor_z = [16]
        sim.prob_relative = [1.1]

    try:
        sim.init_grids()

        qm_eev = _set_carbon(sim, charge_qe)

        # the same transverse points at each longitudinal sample, all at rest
        rng = np.random.default_rng(seed=42)
        num_transverse = 32
        num_longitudinal = 512
        x_m = rng.normal(0.0, 4.0e-5, num_transverse)
        y_m = rng.normal(0.0, 3.0e-5, num_transverse)
        g = rng.gamma(2.0, 1.0, num_longitudinal)
        t_j = (g - 2.0) / np.sqrt(2.0) * 1.0e-3  # skewed, rms 1 mm

        x = np.tile(x_m, num_longitudinal)
        y = np.tile(y_m, num_longitudinal)
        t = np.repeat(t_j, num_transverse)
        zeros = np.zeros_like(x)
        w = np.full_like(x, 1.0e9 / x.size)  # 1e9 real particles
        sim.beam.add_n_particles(x, y, t, zeros, zeros, zeros, qm_eev, w=w)

        # one slice of a thick element: the Strang split applies one collective kick at
        # the slice midpoint, between two transport halves that do not depend on the
        # charge, so every charge state is kicked at the same positions; the drift keeps
        # the kicked px, py and pt, the bend (given by its radius) keeps pt
        if model == "CSR":
            sim.lattice.append(elements.Sbend(name="bend", ds=0.5, rc=10.0, nslice=1))
        else:
            sim.lattice.append(elements.Drift(name="drift", ds=1.0, nslice=1))
        sim.track_particles()

        # the same particle order in every run: the order of the tiles depends on the
        # number of OpenMP threads, and the DataFrame index is the position in a tile
        df = sim.beam.to_df(local=True).reset_index(drop=True).sort_values("idcpu")
        if model == "CSR":
            momenta = ["momentum_t"]
        elif model in ("Gauss3D", "3D"):
            # kicked at fixed t: back at fixed s, pt = -(gamma - gamma_0) / (beta_0 gamma_0)
            # also contains the squared transverse kicks, which scale with q^4
            momenta = ["momentum_x", "momentum_y"]
        else:
            momenta = ["momentum_x", "momentum_y", "momentum_t"]
        return {p: df[p].to_numpy().astype(np.float64) for p in momenta}
    finally:
        sim.finalize()


@pytest.mark.parametrize("charge_qe", [carbon_charge_qe, -1.0, -carbon_charge_qe])
@pytest.mark.parametrize("model", ["Gauss3D", "3D", "Gauss2p5D", "2p5D", "CSR"])
def test_collective_kick_scales_with_charge_squared(model, charge_qe):
    """
    Collective kicks scale with the square of the ion charge, for either sign.

    For the same number of real particles, the self-fields of a beam of ions
    with the charge q scale with q, and so does the force of these fields on each
    ion: the kicks of space charge and CSR scale with q^2.
    """
    if model in ("2p5D", "CSR") and not Config.have_fft:
        pytest.skip("ImpactX built without FFT support; skipping FFT-only model")

    kicks_1 = _kicks(model, 1.0)
    kicks_q = _kicks(model, charge_qe)

    # round-off only; in single precision up to ~5e-5, wrong powers or signs of q give > 0.8
    rtol = 1.0e-3 if Config.precision == "SINGLE" else 1.0e-12
    for p, kick_1 in kicks_1.items():
        kick_q = kicks_q[p]
        kick_max = np.max(np.abs(kick_q))
        deviation = np.max(np.abs(kick_q - charge_qe**2 * kick_1))
        print(
            f"{model}, q={charge_qe:+} e, {p}: max deviation / max kick = "
            f"{deviation / kick_max:.3e}"
        )

        # the beam got a kick ...
        assert kick_max > 0.0
        # ... that scales with q^2
        assert deviation <= rtol * kick_max


@pytest.mark.parametrize("charge_qe", [carbon_charge_qe, -carbon_charge_qe])
@pytest.mark.parametrize("method", ["add_n_particles", "add_particles"])
def test_bunch_charge_weights(method, charge_qe):
    """
    A bunch charge sets the weights to the number of real particles.

    Each real particle carries the charge of the reference particle, so the
    beam charge is the bunch charge, with the sign of the reference particle.
    """
    sim = ImpactX()
    sim.diagnostics = False
    bunch_charge_C = 1.0e-9
    npart = 1000

    try:
        sim.init_grids()

        qm_eev = _set_carbon(sim, charge_qe)

        if method == "add_n_particles":
            rng = np.random.default_rng(seed=42)
            x, y, t = (rng.normal(0.0, 1.0e-3, npart) for _ in range(3))
            zeros = np.zeros_like(x)
            sim.beam.add_n_particles(
                x, y, t, zeros, zeros, zeros, qm_eev, bunch_charge=bunch_charge_C
            )
        else:
            distr = distribution.Waterbag(
                lambdaX=1.0e-3,
                lambdaY=1.0e-3,
                lambdaT=1.0e-3,
                lambdaPx=1.0e-4,
                lambdaPy=1.0e-4,
                lambdaPt=1.0e-4,
            )
            sim.add_particles(bunch_charge_C, distr, npart)

        w = sim.beam.to_df(local=True)["weighting"].to_numpy()
        charge_C = sim.beam.beam_moments()["charge_C"]
    finally:
        sim.finalize()

    # round-off: of a few operations per weight, and of the sum over all weights
    eps = np.finfo(np.float32 if Config.precision == "SINGLE" else np.float64).eps
    num_real = bunch_charge_C / (abs(charge_qe) * 1.602176634e-19)
    assert w.size == npart
    assert np.allclose(w, num_real / npart, rtol=4.0 * eps, atol=0.0)
    assert np.isclose(
        charge_C, np.sign(charge_qe) * bunch_charge_C, rtol=npart * eps, atol=0.0
    )


def test_bunch_charge_needs_ref_charge():
    """A nonzero bunch charge needs the charge of the reference particle to set the weights."""
    sim = ImpactX()
    sim.diagnostics = False

    try:
        sim.init_grids()

        # the reference particle is not set up yet: its charge is zero
        zeros = np.zeros(10)
        with pytest.raises(RuntimeError, match="reference particle charge is zero"):
            sim.beam.add_n_particles(
                zeros, zeros, zeros, zeros, zeros, zeros, 1.0, bunch_charge=1.0e-9
            )
        assert sim.beam.total_number_of_particles() == 0

        # the weights do not need the reference particle charge
        sim.beam.add_n_particles(zeros, zeros, zeros, zeros, zeros, zeros, 1.0, w=zeros)
        assert sim.beam.total_number_of_particles() == 10

        # a bunch charge of zero adds particles of zero weight, e.g., test particles
        sim.beam.add_n_particles(
            zeros, zeros, zeros, zeros, zeros, zeros, 1.0, bunch_charge=0.0
        )
        assert sim.beam.total_number_of_particles() == 20
        w = sim.beam.to_df(local=True)["weighting"].to_numpy()
        assert np.all(w == 0.0)
    finally:
        sim.finalize()


def test_bunch_charge_negative():
    """The bunch charge is a magnitude: the reference particle carries the sign."""
    sim = ImpactX()
    sim.diagnostics = False

    try:
        sim.init_grids()

        qm_eev = _set_carbon(sim, -carbon_charge_qe)

        zeros = np.zeros(10)
        with pytest.raises(RuntimeError, match="must not be negative"):
            sim.beam.add_n_particles(
                zeros, zeros, zeros, zeros, zeros, zeros, qm_eev, bunch_charge=-1.0e-9
            )
        assert sim.beam.total_number_of_particles() == 0
    finally:
        sim.finalize()


@pytest.mark.parametrize("bunch_charge_C", [0.0, 1.0e-9])
def test_bunch_charge_no_particles(bunch_charge_C):
    """Adding no particles with a bunch charge adds nothing.

    MPI ranks that get no particles of a distribution add an empty beam with a
    bunch charge of zero.
    """
    sim = ImpactX()
    sim.diagnostics = False

    try:
        sim.init_grids()

        qm_eev = _set_carbon(sim, carbon_charge_qe)

        empty = np.zeros(0)
        sim.beam.add_n_particles(
            empty,
            empty,
            empty,
            empty,
            empty,
            empty,
            qm_eev,
            bunch_charge=bunch_charge_C,
        )
        assert sim.beam.total_number_of_particles() == 0
    finally:
        sim.finalize()
